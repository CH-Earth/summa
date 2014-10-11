module coupled_em_module
! numerical recipes data types
USE nrtype
! access the number of snow and soil layers
USE data_struc,only:&
                    nSnow,        & ! number of snow layers  
                    nSoil,        & ! number of soil layers  
                    nLayers         ! total number of layers
! physical constants
USE multiconst,only:&
                    Tfreeze,      & ! temperature at freezing              (K)
                    LH_fus,       & ! latent heat of fusion                (J kg-1)
                    iden_ice,     & ! intrinsic density of ice             (kg m-3)
                    iden_water      ! intrinsic density of liquid water    (kg m-3)
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
 USE data_struc,only:type_data,attr_data,forc_data,mpar_data                        ! data structures
 USE data_struc,only:bvar_data,mvar_data,indx_data                                  ! data structures
 USE var_lookup,only:iLookDECISIONS                                                 ! named variables for elements of the decision structure
 USE data_struc,only:ix_soil,ix_snow                                                ! named variables for snow and soil
 USE var_lookup,only:iLookTYPE,iLookATTR,iLookFORCE,iLookPARAM,iLookMVAR,iLookINDEX ! named variables for structure elements
 ! subroutines
 USE getDiagVar_module,only:getDiagVar      ! compute diagnostic variables that are treated as constant over a substep
 USE vegNrgFlux_module,only:wettedFrac      ! compute wetted fraction of the canopy (used in sw radiation fluxes)
 USE canopySnow_module,only:canopySnow      ! compute interception and unloading of snow from the vegetation canopy
 USE vegSWavRad_module,only:vegSWavRad      ! compute canopy sw radiation fluxes
 USE systemSolv_module,only:systemSolv      ! solve the system of thermodynamic and hydrology equations for a given substep
 USE volicePack_module,only:volicePack      ! compute change in volumetric ice content in the snowpack
 ! look-up values for the numerical method
 USE mDecisions_module,only:         &
  iterative,                         &      ! iterative
  nonIterative,                      &      ! non-iterative
  iterSurfEnergyBal                         ! iterate only on the surface energy balance
 ! look-up values for the maximum interception capacity
 USE mDecisions_module,only:         &
                       stickySnow,   &      ! maximum interception capacity an increasing function of temerature
                       lightSnow            ! maximum interception capacity an inverse function of new snow density
 implicit none
 ! define output
 real(dp),intent(inout)               :: dt_init                ! used to initialize the size of the sub-step
 integer(i4b),intent(out)             :: err                    ! error code
 character(*),intent(out)             :: message                ! error message
 ! control the length of the sub-step
 real(dp),pointer                     :: minstep                ! minimum time step (seconds)
 real(dp),pointer                     :: maxstep                ! maximum time step (seconds)
 real(dp)                             :: dt                     ! length of time step (seconds)
 real(dp)                             :: dt_sub                 ! length of the sub-step (seconds)
 real(dp)                             :: dt_done                ! length of time step completed (seconds)
 integer(i4b)                         :: nsub                   ! number of sub-steps
 integer(i4b)                         :: niter                  ! number of iterations
 integer(i4b),parameter               :: n_inc=5                ! minimum number of iterations to increase time step
 integer(i4b),parameter               :: n_dec=15               ! maximum number of iterations to decrease time step
 real(dp),parameter                   :: F_inc = 1.25_dp        ! factor used to increase time step
 real(dp),parameter                   :: F_dec = 0.90_dp        ! factor used to decrease time step
 integer(i4b)                         :: maxiter                ! maxiumum number of iterations
 integer(i4b)                         :: iSnow                  ! index for snow layers
 ! check SWE
 real(dp)                             :: oldSWE                 ! SWE at the start of the substep
 real(dp)                             :: newSWE                 ! SWE at the end of the substep
 real(dp)                             :: delSWE                 ! change in SWE over the subtep
 real(dp)                             :: effRainfall            ! effective rainfall (kg m-2 s-1)
 real(dp)                             :: effSnowfall            ! effective snowfall (kg m-2 s-1)
 real(dp)                             :: sublimation            ! sublimation of ice from the snowpack (kg m-2 s-1)
 real(dp)                             :: snwDrainage            ! drainage of liquid water from the snowpack (m s-1)
 real(dp)                             :: sfcMeltPond            ! surface melt pond (kg m-2)
 real(dp)                             :: massBalance            ! mass balance error (kg m-2)
 ! define other local variables
 character(len=256)                   :: cmessage               ! error message
 logical(lgt)                         :: computeVegFlux         ! flag to indicate if we are computing fluxes over vegetation (.false. means veg is buried with snow)
 integer(i4b)                         :: nLayersRoots           ! number of soil layers that contain roots
 real(dp)                             :: exposedVAI             ! exposed vegetation area index
 real(dp)                             :: dt_wght                ! weight applied to each sub-step, to compute time step average
 real(dp)                             :: dCanopyWetFraction_dWat ! derivative in wetted fraction w.r.t. canopy total water (kg-1 m2)
 real(dp)                             :: dCanopyWetFraction_dT   ! derivative in wetted fraction w.r.t. canopy temperature (K-1)
 real(dp),parameter                   :: varNotUsed1=-9999._dp  ! variables used to calculate derivatives (not needed here)
 real(dp),parameter                   :: varNotUsed2=-9999._dp  ! variables used to calculate derivatives (not needed here)
 integer(i4b)                         :: iLayer                 ! index of model layers
 ! ----------------------------------------------------------------------------------------------------------------------------------------------
 ! ----------------------------------------------------------------------------------------------------------------------------------------------
 ! ** local pointers to increment fluxes
 ! ----------------------------------------------------------------------------------------------------------------------------------------------
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
 ! ----------------------------------------------------------------------------------------------------------------------------------------------
 ! initialize error control
 err=0; message="coupled_em/"

 ! identify the number of snow and soil layers
 nSnow = count(indx_data%var(iLookINDEX%layerType)%dat==ix_snow)
 nSoil = count(indx_data%var(iLookINDEX%layerType)%dat==ix_soil)

 ! compute the total number of snow and soil layers
 nLayers = nSnow + nSoil

 ! assign pointers to algorithmic control parameters
 minstep => mpar_data%var(iLookPARAM%minstep)  ! minimum time step (s)
 maxstep => mpar_data%var(iLookPARAM%maxstep)  ! maximum time step (s)
 !print*, 'minstep, maxstep = ', minstep, maxstep

 ! define maximum number of iterations
 maxiter = mpar_data%var(iLookPARAM%maxiter)

 ! get the length of the time step (seconds)
 dt = data_step

 ! compute the number of layers with roots
 nLayersRoots = count(mvar_data%var(iLookMVAR%iLayerHeight)%dat(nSnow:nLayers-1) < mpar_data%var(iLookPARAM%rootingDepth)-verySmall)
 if(nLayersRoots == 0)then; err=20; message=trim(message)//'no roots within the soil profile'; return; endif

 ! define the foliage nitrogen factor
 mvar_data%var(iLookMVAR%scalarFoliageNitrogenFactor)%dat(1) = 1._dp  ! foliage nitrogen concentration (1.0 = saturated)

 ! initialize the length of the sub-step
 dt_sub  = min(dt_init,min(dt,maxstep))
 dt_done = 0._dp

 ! initialize the number of sub-steps
 nsub=0

 ! loop through sub-steps
 do  ! continuous do statement with exit clause (alternative to "while")

  ! increment the number of sub-steps
  nsub = nsub+1

  ! recompute snow depth and SWE
  if(nSnow > 0)then
   mvar_data%var(iLookMVAR%scalarSnowDepth)%dat(1) = sum(mvar_data%var(iLookMVAR%mLayerDepth)%dat(1:nSnow))
   mvar_data%var(iLookMVAR%scalarSWE)%dat(1)       = sum( (mvar_data%var(iLookMVAR%mLayerVolFracLiq)%dat(1:nSnow)*iden_water + &
                                                           mvar_data%var(iLookMVAR%mLayerVolFracIce)%dat(1:nSnow)*iden_ice) &
                                                           * mvar_data%var(iLookMVAR%mLayerDepth)%dat(1:nSnow) )
  endif

  ! save SWE (mass balance check)
  oldSWE = mvar_data%var(iLookMVAR%scalarSWE)%dat(1)
  !print*, 'oldSWE = ', oldSWE

  ! compute the temperature of the root zone: used in vegetation phenology
  mvar_data%var(iLookMVAR%scalarRootZoneTemp)%dat(1) = sum(mvar_data%var(iLookMVAR%mLayerTemp)%dat(nSnow+1:nSnow+nLayersRoots)) / real(nLayersRoots, kind(dp))

  ! (1) compute variables that are treated as constant over a model substep...
  ! --------------------------------------------------------------------------
  ! (phenology, thermal properties, snow albedo)
  call getDiagVar(&
                  dt_sub,                       & ! intent(in): time step (s)
                  computeVegFlux,               & ! intent(out): flag to indicate if we are computing fluxes over vegetation (.false. means veg is buried with snow)
                  err,cmessage)                   ! intent(out): error control
  if(err/=0)then; err=20; message=trim(message)//trim(cmessage); return; endif

  ! compute exposed vegetation area index (depends on phenology)
  exposedVAI = mvar_data%var(iLookMVAR%scalarExposedLAI)%dat(1) + mvar_data%var(iLookMVAR%scalarExposedSAI)%dat(1)

  ! compute maximum canopy liquid water (kg m-2)
  mvar_data%var(iLookMVAR%scalarCanopyLiqMax)%dat(1) = mpar_data%var(iLookPARAM%refInterceptCapRain)*exposedVAI
  print*, 'mvar_data%var(iLookMVAR%scalarCanopyLiqMax)%dat(1) = ', mvar_data%var(iLookMVAR%scalarCanopyLiqMax)%dat(1)

  ! compute maximum canopy ice content (kg m-2)
  ! NOTE 1: this is used to compute the snow fraction on the canopy, as used in *BOTH* the radiation AND canopy sublimation routines
  ! NOTE 2: this is a different variable than the max ice used in the throughfall (snow interception) calculations
  select case(model_decisions(iLookDECISIONS%snowIncept)%iDecision)
   case(lightSnow);  mvar_data%var(iLookMVAR%scalarCanopyIceMax)%dat(1) = exposedVAI*mpar_data%var(iLookPARAM%refInterceptCapSnow)       ! use maximum per unit leaf area storage capacity for snow (kg m-2)
   case(stickySnow); mvar_data%var(iLookMVAR%scalarCanopyIceMax)%dat(1) = exposedVAI*mpar_data%var(iLookPARAM%refInterceptCapSnow)*4._dp ! use maximum per unit leaf area storage capacity for snow (kg m-2)
   case default; message=trim(message)//'unable to identify option for maximum branch interception capacity'; err=20; return
  end select ! identifying option for maximum branch interception capacity

  ! compute wetted fraction of the canopy
  ! NOTE: assume that the wetted fraction is constant over the substep for the radiation calculations
  if(computeVegFlux)then

   ! compute wetted fraction of the canopy
   call wettedFrac(&
                   ! input
                   .false.,                                                      & ! flag to denote if derivatives are required
                   .false.,                                                      & ! flag to denote if derivatives are calculated numerically
                   (mvar_data%var(iLookMVAR%scalarCanopyTemp)%dat(1) < Tfreeze), & ! flag to denote if the canopy is frozen
                   varNotUsed1,                                                  & ! derivative in canopy liquid w.r.t. canopy temperature (kg m-2 K-1)
                   varNotUsed2,                                                  & ! fraction of liquid water on the canopy
                   mvar_data%var(iLookMVAR%scalarCanopyLiq)%dat(1),              & ! canopy liquid water (kg m-2)
                   mvar_data%var(iLookMVAR%scalarCanopyIce)%dat(1),              & ! canopy ice (kg m-2)
                   mvar_data%var(iLookMVAR%scalarCanopyLiqMax)%dat(1),           & ! maximum canopy liquid water (kg m-2)
                   mvar_data%var(iLookMVAR%scalarCanopyLiqMax)%dat(1),           & ! maximum canopy ice content (kg m-2)
                   ! output
                   mvar_data%var(iLookMVAR%scalarCanopyWetFraction)%dat(1),      & ! canopy wetted fraction (-)
                   dCanopyWetFraction_dWat,                                      & ! derivative in wetted fraction w.r.t. canopy liquid water content (kg-1 m2)
                   dCanopyWetFraction_dT,                                        & ! derivative in wetted fraction w.r.t. canopy liquid water content (kg-1 m2)
                   err,cmessage)
   if(err/=0)then; message=trim(message)//trim(cmessage); return; endif
  else
   mvar_data%var(iLookMVAR%scalarCanopyWetFraction)%dat(1) = 0._dp
   dCanopyWetFraction_dWat                                 = 0._dp
   dCanopyWetFraction_dT                                   = 0._dp
  endif


  ! (2) compute canopy throughfall and unloading...
  ! -----------------------------------------------
  ! NOTE 1: this needs to be done before solving the energy and liquid water equations, to account for the heat advected with precipitation (and throughfall/unloading)
  ! NOTE 2: the unloading flux is computed using canopy drip (scalarCanopyLiqDrainage) from the previous time step
  call canopySnow(&
                  dt_sub,                       & ! intent(in): time step (seconds)
                  computeVegFlux,               & ! intent(in): flag to indicate if we are computing fluxes over vegetation (.false. means veg is buried with snow)
                  err,cmessage                  ) ! intent(out): error control
  if(err/=0)then; err=20; message=trim(message)//trim(cmessage); return; endif

  ! initialize drainage and throughfall
  ! NOTE 1: this needs to be done before solving the energy and liquid water equations, to account for the heat advected with precipitation
  ! NOTE 2: this initialization needs to be done AFTER the call to canopySnow, since canopySnow uses canopy drip drom the previous time step
  if(.not.computeVegFlux)then
   mvar_data%var(iLookMVAR%scalarThroughfallRain)%dat(1)   = mvar_data%var(iLookMVAR%scalarRainfall)%dat(1)
   mvar_data%var(iLookMVAR%scalarCanopyLiqDrainage)%dat(1) = 0._dp
  else
   mvar_data%var(iLookMVAR%scalarThroughfallRain)%dat(1)   = 0._dp
   mvar_data%var(iLookMVAR%scalarCanopyLiqDrainage)%dat(1) = 0._dp
  endif


  ! (3) compute canopy sw radiation fluxes...
  ! -----------------------------------------
  call vegSWavRad(&
                  dt_sub,                       & ! intent(in): time step (s) -- only used in Noah-MP radiation, to compute albedo
                  computeVegFlux,               & ! intent(in): logical flag to compute vegetation fluxes (.false. if veg buried by snow)
                  err,cmessage)                   ! intent(out): error control
  if(err/=0)then; err=20; message=trim(message)//trim(cmessage); return; endif


  ! recompute snow depth and SWE
  if(nSnow > 0)then
   mvar_data%var(iLookMVAR%scalarSnowDepth)%dat(1) = sum(mvar_data%var(iLookMVAR%mLayerDepth)%dat(1:nSnow))
   mvar_data%var(iLookMVAR%scalarSWE)%dat(1)       = sum( (mvar_data%var(iLookMVAR%mLayerVolFracLiq)%dat(1:nSnow)*iden_water + &
                                                           mvar_data%var(iLookMVAR%mLayerVolFracIce)%dat(1:nSnow)*iden_ice) &
                                                           * mvar_data%var(iLookMVAR%mLayerDepth)%dat(1:nSnow) )
  endif

  !print*, 'before picardSolv: SWE = ', mvar_data%var(iLookMVAR%scalarSWE)%dat(1)

  ! (4) solve model equations...
  ! ----------------------------
  do  ! (multiple attempts for non-convergence)
   ! get the new solution
   call systemSolv(&
                   ! input: model control
                   dt_sub,                                 & ! intent(in): length of the model sub-step
                   maxiter,                                & ! intent(in): maximum number of iterations
                   (nsub==1),                              & ! intent(in): logical flag to denote the first substep
                   computeVegFlux,                         & ! intent(in): logical flag to compute fluxes within the vegetation canopy
                   ! input/output: data structures
                   type_data,                              & ! intent(in):    type of vegetation and soil
                   attr_data,                              & ! intent(in):    spatial attributes
                   forc_data,                              & ! intent(in):    model forcing data
                   mpar_data,                              & ! intent(in):    model parameters
                   indx_data,                              & ! intent(in):    index data
                   mvar_data,                              & ! intent(inout): model variables for a local HRU
                   bvar_data,                              & ! intent(in):    model variables for the local basin
                   model_decisions,                        & ! intent(in):    model decisions
                   ! output: model control
                   niter,                                  & ! intent(out): number of iterations
                   err,cmessage)                             ! intent(out): error code and error message
   !print*, '*************************************************'
   !print*, '*************************************************'
   !print*, '*************************************************'
   !pause
   ! check for fatal errors
   if(err>0)then; err=20; message=trim(message)//trim(cmessage); return; endif 
   if(err==0) exit  ! exit do loop if all is a-ok
   ! if warnings (non-covergence, and so forth) then reduce the time step and try again
   dt_sub = dt_sub*0.5_dp
   if(dt_sub < minstep)then
    message=trim(message)//'dt_sub is below the minimum time step'
    err=20; return
   endif
  end do  ! (multiple attempts for non-convergence)



  ! (4) use Picard iteration to solve model equations...
  ! ----------------------------------------------------
  !do
  ! ! get the new solution
  ! call picardSolv(dt_sub,maxiter,(nsub==1),computeVegFlux,&  ! input
  !                 niter,err,cmessage)                        ! output
  ! ! check errors
  ! if(err > 0)then; message=trim(message)//trim(cmessage); return; endif  ! fatal error
  ! if(err==0) exit  ! exit do loop if all is a-ok
  ! ! err < 0 means non-fatal error, spreduce time step and try again
  ! dt_sub = dt_sub*0.1_dp
  ! ! check that the step size is still appropriate -- if not, use non-iterative solution
  ! if(dt_sub < minstep)then
  !  if(err/=0)then; message=trim(message)//'dt_sub is below the minimum time step'; return; endif
  !  dt_sub  = minstep
  !  ! just iterate once
  !  call picardSolv(dt_sub,1,(nsub==1),computeVegFlux,     &  ! input
  !                  niter,err,cmessage)                       ! output
  !  if(err/=0)then; message=trim(message)//trim(cmessage); return; endif
  !  exit ! exit do loop if all is a-ok
  ! endif
  !end do

  ! compute effective rainfall input to the snowpack
  if(nSnow > 0)then
   effRainfall = mvar_data%var(iLookMVAR%scalarThroughfallRain)%dat(1) + mvar_data%var(iLookMVAR%scalarCanopyLiqDrainage)%dat(1)
  else
   effRainfall = 0._dp  ! if no snow layers, water is added to the top of the soil zone
  endif

  ! save snow drainage
  if(nSnow > 0)then
   snwDrainage = mvar_data%var(iLookMVAR%iLayerLiqFluxSnow)%dat(nSnow)*iden_water  ! m s-1 -> kg m-2 s-1
  else
   snwDrainage = 0._dp   ! (no snow drainage when no snow layers)
  endif

  !print*, 'after picardSolv: SWE = ', mvar_data%var(iLookMVAR%scalarSWE)%dat(1)

  ! (5) account for changes in volumetric ice content in the snowpack...
  ! --------------------------------------------------------------------
  ! (add snowfall; combine and sub-divide layers if necessary)
  call volicePack(&
                  dt_sub,                       & ! intent(in): time step (seconds)
                  err,cmessage)                  ! intent(out): error control
  if(err/=0)then; err=55; message=trim(message)//trim(cmessage); return; endif

  ! re-compute the number of snow and soil layers
  nSnow   = count(indx_data%var(iLookINDEX%layerType)%dat==ix_snow)
  nSoil   = count(indx_data%var(iLookINDEX%layerType)%dat==ix_soil)
  nLayers = nSnow + nSoil

  ! check SWE
  sfcMeltPond = mvar_data%var(iLookMVAR%scalarSfcMeltPond)%dat(1)
  effSnowfall = mvar_data%var(iLookMVAR%scalarThroughfallSnow)%dat(1) + mvar_data%var(iLookMVAR%scalarCanopySnowUnloading)%dat(1)
  sublimation = mvar_data%var(iLookMVAR%scalarSnowSublimation)%dat(1)
  newSWE      = mvar_data%var(iLookMVAR%scalarSWE)%dat(1)
  delSWE      = newSWE - (oldSWE - sfcMeltPond)
  massBalance = delSWE - (effSnowfall + effRainfall + sublimation - snwDrainage)*dt_sub
  if(abs(massBalance) > 1.d-6)then
   print*,                  'nSnow       = ', nSnow
   write(*,'(a,1x,f20.10)') 'dt_sub      = ', dt_sub
   write(*,'(a,1x,f20.10)') 'oldSWE      = ', oldSWE
   write(*,'(a,1x,f20.10)') 'newSWE      = ', newSWE
   write(*,'(a,1x,f20.10)') 'delSWE      = ', delSWE
   write(*,'(a,1x,f20.10)') 'effRainfall = ', effRainfall*dt_sub
   write(*,'(a,1x,f20.10)') 'effSnowfall = ', effSnowfall*dt_sub
   write(*,'(a,1x,f20.10)') 'sublimation = ', sublimation*dt_sub
   write(*,'(a,1x,f20.10)') 'snwDrainage = ', snwDrainage*dt_sub
   write(*,'(a,1x,f20.10)') 'sfcMeltPond = ', sfcMeltPond
   write(*,'(a,1x,f20.10)') 'massBalance = ', massBalance
   message=trim(message)//'SWE does not balance'
   err=20; return
  endif

  ! compute the total number of snow and soil layers
  nLayers = nSnow + nSoil

  ! increment model fluxes
  dt_wght = dt_sub/dt ! define weight applied to each sub-step
  call increment_fluxes(dt_wght,nsub==1)

  ! increment the time step increment
  dt_done = dt_done + dt_sub
  print*, '***** ', dt_done, dt_sub, niter

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
  
 end do  ! (sub-step loop)
 !pause 'completed time step'

 !print*, 'mvar_data%var(iLookMVAR%averageCanopyLiqDrainage)%dat(1) = ', mvar_data%var(iLookMVAR%averageCanopyLiqDrainage)%dat(1)

 ! save the surface temperature (just to make things easier to visualize)
 mvar_data%var(iLookMVAR%scalarSurfaceTemp)%dat(1) = mvar_data%var(iLookMVAR%mLayerTemp)%dat(1)

 iLayer = nSnow+1
 !print*, 'nsub, mLayerTemp(iLayer), mLayerVolFracIce(iLayer) = ', nsub, mLayerTemp(iLayer), mLayerVolFracIce(iLayer)
 print*, 'nsub = ', nsub
 if(nsub>2000)then
  message=trim(message)//'number of sub-steps > 2000'
  err=20; return
 endif


 ! ********************************************************************************************************************************
 ! ********************************************************************************************************************************
 ! ********************************************************************************************************************************

 contains

  ! calculate timestep-average fluxes
  subroutine increment_fluxes(dt_wght,initialize)
  real(dp),intent(in)     :: dt_wght    ! weight assigned to sub-step
  logical(lgt),intent(in) :: initialize ! flag to initialize fluxes

  ! set up pointers
  
  ! assign pointers to timestep-average model fluxes
  if(initialize)then
   averageThroughfallSnow     => mvar_data%var(iLookMVAR%averageThroughfallSnow)%dat(1)        ! snow that reaches the ground without ever touching the canopy (kg m-2 s-1)
   averageThroughfallRain     => mvar_data%var(iLookMVAR%averageThroughfallRain)%dat(1)        ! rain that reaches the ground without ever touching the canopy (kg m-2 s-1)
   averageCanopySnowUnloading => mvar_data%var(iLookMVAR%averageCanopySnowUnloading)%dat(1)    ! unloading of snow from the vegetion canopy (kg m-2 s-1)
   averageCanopyLiqDrainage   => mvar_data%var(iLookMVAR%averageCanopyLiqDrainage)%dat(1)      ! drainage of liquid water from the vegetation canopy (kg m-2 s-1)
   averageCanopyMeltFreeze    => mvar_data%var(iLookMVAR%averageCanopyMeltFreeze)%dat(1)       ! melt/freeze of water stored in the canopy (kg m-2 s-1)
   averageSnowSublimation     => mvar_data%var(iLookMVAR%averageSnowSublimation)%dat(1)        ! snow sublimation/frost - below canopy or non-vegetated (kg m-2 s-1)
   averageGroundEvaporation   => mvar_data%var(iLookMVAR%averageGroundEvaporation)%dat(1)      ! ground evaporation/condensation - below canopy or non-vegetated (kg m-2 s-1)
   averageRainPlusMelt        => mvar_data%var(iLookMVAR%averageRainPlusMelt)%dat(1)           ! rain plus melt, as input to soil before calculating surface runoff (m s-1)
   averageSurfaceRunoff       => mvar_data%var(iLookMVAR%averageSurfaceRunoff)%dat(1)          ! surface runoff (m s-1)
   averageSoilInflux          => mvar_data%var(iLookMVAR%averageSoilInflux)%dat(1)             ! influx of water at the top of the soil profile (m s-1)
   averageSoilBaseflow        => mvar_data%var(iLookMVAR%averageSoilBaseflow)%dat(1)           ! total baseflow from throughout the soil profile (m s-1)
   averageSoilDrainage        => mvar_data%var(iLookMVAR%averageSoilDrainage)%dat(1)           ! drainage from the bottom of the soil profile (m s-1)
   averageAquiferRecharge     => mvar_data%var(iLookMVAR%averageAquiferRecharge)%dat(1)        ! recharge to the aquifer (m s-1)
   averageAquiferBaseflow     => mvar_data%var(iLookMVAR%averageAquiferBaseflow)%dat(1)        ! baseflow from the aquifer (m s-1)
   averageAquiferTranspire    => mvar_data%var(iLookMVAR%averageAquiferTranspire)%dat(1)       ! transpiration from the aquifer (m s-1)
   averageColumnOutflow       => mvar_data%var(iLookMVAR%averageColumnOutflow)%dat             ! outflow from each layer in the soil profile (m3 s-1)
  endif

  ! initialize average fluxes
  if(initialize)then
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
  endif

  ! assign pointers to the model flux variables
  ! NOTE: need to do this every sub-step becaause the model structures are re-defined
  scalarThroughfallSnow     => mvar_data%var(iLookMVAR%scalarThroughfallSnow)%dat(1)           ! snow that reaches the ground without ever touching the canopy (kg m-2 s-1)
  scalarThroughfallRain     => mvar_data%var(iLookMVAR%scalarThroughfallRain)%dat(1)           ! rain that reaches the ground without ever touching the canopy (kg m-2 s-1)
  scalarCanopySnowUnloading => mvar_data%var(iLookMVAR%scalarCanopySnowUnloading)%dat(1)       ! unloading of snow from the vegetion canopy (kg m-2 s-1)
  scalarCanopyLiqDrainage   => mvar_data%var(iLookMVAR%scalarCanopyLiqDrainage)%dat(1)         ! drainage of liquid water from the vegetation canopy (kg m-2 s-1)
  scalarCanopyMeltFreeze    => mvar_data%var(iLookMVAR%scalarCanopyMeltFreeze)%dat(1)          ! melt/freeze of water stored in the canopy (kg m-2 s-1)
  scalarSnowSublimation     => mvar_data%var(iLookMVAR%scalarSnowSublimation)%dat(1)           ! snow sublimation/frost - below canopy or non-vegetated (kg m-2 s-1)
  scalarGroundEvaporation   => mvar_data%var(iLookMVAR%scalarGroundEvaporation)%dat(1)         ! ground evaporation/condensation - below canopy or non-vegetated (kg m-2 s-1)
  scalarRainPlusMelt        => mvar_data%var(iLookMVAR%scalarRainPlusMelt)%dat(1)              ! rain plus melt, as input to soil before calculating surface runoff (m s-1)
  scalarSurfaceRunoff       => mvar_data%var(iLookMVAR%scalarSurfaceRunoff)%dat(1)             ! surface runoff (m s-1)
  scalarSoilInflux          => mvar_data%var(iLookMVAR%scalarSoilInflux)%dat(1)                ! influx of water at the top of the soil profile (m s-1)
  scalarSoilBaseflow        => mvar_data%var(iLookMVAR%scalarSoilBaseflow)%dat(1)              ! total baseflow from throughout the soil profile (m s-1)
  scalarSoilDrainage        => mvar_data%var(iLookMVAR%scalarSoilDrainage)%dat(1)              ! drainage from the bottom of the soil profile (m s-1)
  scalarAquiferRecharge     => mvar_data%var(iLookMVAR%scalarAquiferRecharge)%dat(1)           ! recharge to the aquifer (m s-1)
  scalarAquiferBaseflow     => mvar_data%var(iLookMVAR%scalarAquiferBaseflow)%dat(1)           ! baseflow from the aquifer (m s-1)
  scalarAquiferTranspire    => mvar_data%var(iLookMVAR%scalarAquiferTranspire)%dat(1)          ! transpiration from the aquifer (m s-1)
  mLayerColumnOutflow       => mvar_data%var(iLookMVAR%mLayerColumnOutflow)%dat                ! total outflow from each layer in a given soil column (m3 s-1) 

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

  end subroutine increment_fluxes

  ! ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------

 
 end subroutine coupled_em

end module coupled_em_module

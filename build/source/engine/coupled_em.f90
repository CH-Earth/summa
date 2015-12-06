! SUMMA - Structure for Unifying Multiple Modeling Alternatives
! Copyright (C) 2014-2015 NCAR/RAL
!
! This file is part of SUMMA
!
! For more information see: http://www.ral.ucar.edu/projects/summa
!
! This program is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! This program is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License
! along with this program.  If not, see <http://www.gnu.org/licenses/>.

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
 ! public subroutine coupled_em: run the coupled energy-mass model for one timestep
 ! ************************************************************************************************
 subroutine coupled_em(printRestart,output_fileSuffix,dt_init,err,message)
 ! data structures and named variables
 USE data_struc,only:data_step                                                      ! time step of forcing data (s)
 USE data_struc,only:model_decisions                                                ! model decision structure
 USE data_struc,only:type_data,attr_data,forc_data,mpar_data                        ! data structures
 USE data_struc,only:bvar_data,mvar_data,indx_data,time_data                        ! data structures
 USE var_lookup,only:iLookDECISIONS                                                 ! named variables for elements of the decision structure
 USE data_struc,only:ix_soil,ix_snow                                                ! named variables for snow and soil
 USE var_lookup,only:iLookTYPE,iLookATTR,iLookFORCE,iLookPARAM,iLookMVAR,iLookINDEX ! named variables for structure elements
 ! preliminary subroutines
 USE vegPhenlgy_module,only:vegPhenlgy      ! (1) compute vegetation phenology
 USE vegNrgFlux_module,only:wettedFrac      ! (2) compute wetted fraction of the canopy (used in sw radiation fluxes)
 USE snowAlbedo_module,only:snowAlbedo      ! (3) compute snow albedo
 USE vegSWavRad_module,only:vegSWavRad      ! (4) compute canopy sw radiation fluxes
 USE canopySnow_module,only:canopySnow      ! (5) compute interception and unloading of snow from the vegetation canopy
 USE volicePack_module,only:newsnwfall      ! (6) compute change in the top snow layer due to throughfall and unloading
 USE volicePack_module,only:volicePack      ! (7) merge and sub-divide snow layers, if necessary
 USE diagn_evar_module,only:diagn_evar      ! (8) compute diagnostic energy variables -- thermal conductivity and heat capacity
 ! the model solver
 USE systemSolv_module,only:systemSolv      ! solve the system of thermodynamic and hydrology equations for a given substep
 ! additional subroutines
 USE tempAdjust_module,only:tempAdjust      ! adjust snow temperature associated with new snowfall
 USE snwDensify_module,only:snwDensify      ! snow densification (compaction and cavitation)
 USE var_derive_module,only:calcHeight      ! module to calculate height at layer interfaces and layer mid-point
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
 character(*),intent(in)              :: output_fileSuffix      ! suffix for the output file (used to write re-start files)
 logical(lgt),intent(in)              :: printRestart           ! flag to print a re-start file
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
 real(dp)                             :: snwDrainage            ! drainage of liquid water from the snowpack (m s-1 -> kg m-2 s-1)
 real(dp)                             :: sfcMeltPond            ! surface melt pond (kg m-2)
 real(dp)                             :: massBalance            ! mass balance error (kg m-2)
 ! define other local variables
 character(len=256)                   :: cmessage               ! error message
 logical(lgt)                         :: computeVegFlux         ! flag to indicate if we are computing fluxes over vegetation (.false. means veg is buried with snow)
 integer(i4b)                         :: nLayersRoots           ! number of soil layers that contain roots
 real(dp)                             :: scalarCanopyWater      ! total canopy water (kg m-2)
 real(dp)                             :: canopyDepth            ! canopy depth (m)
 real(dp)                             :: exposedVAI             ! exposed vegetation area index
 real(dp)                             :: dt_wght                ! weight applied to each sub-step, to compute time step average
 real(dp)                             :: dCanopyWetFraction_dWat ! derivative in wetted fraction w.r.t. canopy total water (kg-1 m2)
 real(dp)                             :: dCanopyWetFraction_dT   ! derivative in wetted fraction w.r.t. canopy temperature (K-1)
 real(dp),parameter                   :: varNotUsed1=-9999._dp  ! variables used to calculate derivatives (not needed here)
 real(dp),parameter                   :: varNotUsed2=-9999._dp  ! variables used to calculate derivatives (not needed here)
 integer(i4b)                         :: iLayer                 ! index of model layers
 real(dp)                             :: volSub                 ! volumetric sublimation (kg m-3)
 real(dp),parameter                   :: tinyNumber=tiny(1._dp) ! a tiny number
 real(dp)                             :: dt_solv                ! progress towards dt_sub
 real(dp)                             :: dt_temp                ! temporary sub-step length
 real(dp)                             :: dt_prog                ! progress of time step (s)
 real(dp)                             :: dt_frac                ! fraction of time step (-)
 integer(i4b)                         :: nTemp                  ! number of temporary sub-steps
 integer(i4b)                         :: nTrial                 ! number of trial sub-steps
 logical(lgt)                         :: rejectedStep           ! flag to denote if the sub-step is rejected (convergence problem, etc.)
 ! balance checks
 real(dp),pointer                     :: scalarSnowfall         ! snowfall rate
 real(dp),pointer                     :: scalarRainfall         ! rainfall rate
 real(dp)                             :: scalarCanopyWatBalError ! water balance error for the vegetation canopy (kg m-2)
 real(dp)                             :: scalarSoilWatBalError  ! water balance error (kg m-2)
 real(dp)                             :: scalarTotalSoilLiq     ! total liquid water in the soil column (kg m-2)
 real(dp)                             :: scalarTotalSoilIce     ! total ice in the soil column (kg m-2)
 real(dp)                             :: balanceCanopyWater0    ! total water stored in the vegetation canopy at the start of the step (kg m-2)
 real(dp)                             :: balanceCanopyWater1    ! total water stored in the vegetation canopy at the end of the step (kg m-2)
 real(dp)                             :: balanceSoilWater0      ! total soil storage at the start of the step (kg m-2)
 real(dp)                             :: balanceSoilWater1      ! total soil storage at the end of the step (kg m-2)
 real(dp)                             :: balanceSoilInflux      ! input to the soil zone
 real(dp)                             :: balanceSoilBaseflow    ! output from the soil zone
 real(dp)                             :: balanceSoilDrainage    ! output from the soil zone
 real(dp)                             :: balanceSoilTranspiration     ! output from the soil zone
 real(dp)                             :: balanceAquifer0        ! total aquifer storage at the start of the step (kg m-2)
 real(dp)                             :: balanceAquifer1        ! total aquifer storage at the end of the step (kg m-2)
 real(dp)                             :: xCompress              ! compression in a given layer (m)
 real(dp)                             :: xFlux0,xFlux1          ! fluxes at the layer boundaries (m)
 ! ----------------------------------------------------------------------------------------------------------------------------------------------
 ! ** local pointers to model state variables
 ! ----------------------------------------------------------------------------------------------------------------------------------------------
 real(dp),pointer                     :: mLayerDepth(:)         ! depth of each soil layer (m)
 real(dp),pointer                     :: mLayerVolFracIce(:)    ! volumetric ice content in each soil layer (-)
 real(dp),pointer                     :: mLayerVolFracLiq(:)    ! volumetric liquid water content in each soil layer (-)
 real(dp),pointer                     :: scalarAquiferStorage   ! aquifer storage (m)
 real(dp),pointer                     :: scalarCanopyLiq        ! canopy liquid water content (kg m-2)
 real(dp),pointer                     :: scalarCanopyIce        ! canopy ice content (kg m-2)
 ! ----------------------------------------------------------------------------------------------------------------------------------------------
 ! ** local pointers to increment fluxes
 ! ----------------------------------------------------------------------------------------------------------------------------------------------
 ! local pointers to flux variables
 real(dp),pointer                     :: scalarThroughfallSnow       ! snow that reaches the ground without ever touching the canopy (kg m-2 s-1)
 real(dp),pointer                     :: scalarThroughfallRain       ! rain that reaches the ground without ever touching the canopy (kg m-2 s-1)
 real(dp),pointer                     :: scalarCanopySnowUnloading   ! unloading of snow from the vegetion canopy (kg m-2 s-1)
 real(dp),pointer                     :: scalarCanopyLiqDrainage     ! drainage of liquid water from the vegetation canopy (kg m-2 s-1)
 real(dp),pointer                     :: scalarCanopyMeltFreeze      ! melt/freeze of water stored in the canopy (kg m-2 s-1)
 real(dp),pointer                     :: scalarCanopyTranspiration   ! canopy transpiration (kg m-2 s-1)
 real(dp),pointer                     :: scalarCanopyEvaporation     ! canopy evaporation/condensation (kg m-2 s-1)
 real(dp),pointer                     :: scalarCanopySublimation     ! canopy sublimation/frost (kg m-2 s-1)
 real(dp),pointer                     :: scalarSnowSublimation       ! snow sublimation/frost - below canopy or non-vegetated (kg m-2 s-1)
 real(dp),pointer                     :: scalarGroundEvaporation     ! ground evaporation/condensation - below canopy or non-vegetated (kg m-2 s-1)
 real(dp),pointer                     :: scalarRainPlusMelt          ! rain plus melt, as input to soil before calculating surface runoff (m s-1)
 real(dp),pointer                     :: scalarSurfaceRunoff         ! surface runoff (m s-1)
 real(dp),pointer                     :: scalarSoilInflux            ! influx of water at the top of the soil profile (m s-1)
 real(dp),pointer                     :: scalarSoilCompress          ! change in storage associated with compression of the soil matrix (kg m-2)
 real(dp),pointer                     :: scalarSoilBaseflow          ! total baseflow from throughout the soil profile (m s-1)
 real(dp),pointer                     :: scalarSoilDrainage          ! drainage from the bottom of the soil profile (m s-1)
 real(dp),pointer                     :: scalarAquiferRecharge       ! recharge to the aquifer (m s-1)
 real(dp),pointer                     :: scalarAquiferBaseflow       ! baseflow from the aquifer (m s-1)
 real(dp),pointer                     :: scalarAquiferTranspire      ! transpiration from the aquifer (m s-1)
 real(dp),pointer                     :: mLayerColumnOutflow(:)      ! total outflow from each layer in a given soil column (m3 s-1)
 ! local pointers to timestep-average flux variables
 real(dp),pointer                     :: totalSoilCompress           ! change in storage associated with compression of the soil matrix (kg m-2)
 real(dp),pointer                     :: averageThroughfallSnow      ! snow that reaches the ground without ever touching the canopy (kg m-2 s-1)
 real(dp),pointer                     :: averageThroughfallRain      ! rain that reaches the ground without ever touching the canopy (kg m-2 s-1)
 real(dp),pointer                     :: averageCanopySnowUnloading  ! unloading of snow from the vegetion canopy (kg m-2 s-1)
 real(dp),pointer                     :: averageCanopyLiqDrainage    ! drainage of liquid water from the vegetation canopy (kg m-2 s-1)
 real(dp),pointer                     :: averageCanopyMeltFreeze     ! melt/freeze of water stored in the canopy (kg m-2 s-1)
 real(dp),pointer                     :: averageCanopyTranspiration  ! canopy transpiration (kg m-2 s-1)
 real(dp),pointer                     :: averageCanopyEvaporation    ! canopy evaporation/condensation (kg m-2 s-1)
 real(dp),pointer                     :: averageCanopySublimation    ! canopy sublimation/frost (kg m-2 s-1)
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

 ! This is the start of a data step for a local HRU

 ! count the number of snow and soil layers
 ! NOTE: need to re-compute the number of snow and soil layers at the start of each sub-step because the number of layers may change
 !         (nSnow and nSoil are shared in the data structure)
 nSnow = count(indx_data%var(iLookINDEX%layerType)%dat==ix_snow)
 nSoil = count(indx_data%var(iLookINDEX%layerType)%dat==ix_soil)

 ! compute the total number of snow and soil layers
 nLayers = nSnow + nSoil

 scalarCanopyLiq => mvar_data%var(iLookMVAR%scalarCanopyLiq)%dat(1)
 scalarCanopyIce => mvar_data%var(iLookMVAR%scalarCanopyIce)%dat(1)
 balanceCanopyWater0 = scalarCanopyLiq + scalarCanopyIce

 ! point to model state variables
 ! NOTE: need to do this at the start of each sub-step because number of layers may change
 mLayerDepth          => mvar_data%var(iLookMVAR%mLayerDepth)%dat(nSnow+1:nLayers)           ! depth of each soil layer (m)
 mLayerVolFracIce     => mvar_data%var(iLookMVAR%mLayerVolFracIce)%dat(nSnow+1:nLayers)      ! volumetric ice content in each soil layer (-)
 mLayerVolFracLiq     => mvar_data%var(iLookMVAR%mLayerVolFracLiq)%dat(nSnow+1:nLayers)      ! volumetric liquid water content in each soil layer (-)
 scalarAquiferStorage => mvar_data%var(iLookMVAR%scalarAquiferStorage)%dat(1)                ! aquifer storage (m)

 ! compute total soil moisture and ice at the *START* of the step (kg m-2)
 scalarTotalSoilLiq = sum(iden_water*mLayerVolFracLiq(1:nSoil)*mLayerDepth(1:nSoil))
 scalarTotalSoilIce = sum(iden_water*mLayerVolFracIce(1:nSoil)*mLayerDepth(1:nSoil))  ! NOTE: no expansion and hence use iden_water
 !scalarTotalSoilIce = sum(iden_ice  *mLayerVolFracIce(1:nSoil)*mLayerDepth(1:nSoil))

 ! get the total water in the soil (liquid plus ice) at the start of the time step (kg m-2)
 balanceSoilWater0 = scalarTotalSoilLiq + scalarTotalSoilIce

 ! get the total aquifer storage at the start of the time step (kg m-2)
 balanceAquifer0 = scalarAquiferStorage*iden_water

 ! print re-start file
 if(printRestart)then
  call printRestartFile(output_fileSuffix,dt_init,time_data,mvar_data,err,cmessage)
  if(err/=0)then; err=20; message=trim(message)//trim(cmessage); return; endif
  !pause
 endif

 ! assign pointers to algorithmic control parameters
 minstep => mpar_data%var(iLookPARAM%minstep)  ! minimum time step (s)
 maxstep => mpar_data%var(iLookPARAM%maxstep)  ! maximum time step (s)
 !print*, 'minstep, maxstep = ', minstep, maxstep

 ! define maximum number of iterations
 maxiter = nint(mpar_data%var(iLookPARAM%maxiter))

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

  ! print progress
  !print*, '*** new substep'
  !write(*,'(a,3(f11.4,1x))') 'dt_sub, dt_init, dt = ', dt_sub, dt_init, dt

  ! increment the number of sub-steps
  nsub = nsub+1

  ! save SWE
  oldSWE = mvar_data%var(iLookMVAR%scalarSWE)%dat(1)
  !print*, 'nSnow = ', nSnow
  !print*, 'oldSWE = ', oldSWE


  ! (1) compute phenology...
  ! ------------------------

  ! compute the temperature of the root zone: used in vegetation phenology
  mvar_data%var(iLookMVAR%scalarRootZoneTemp)%dat(1) = sum(mvar_data%var(iLookMVAR%mLayerTemp)%dat(nSnow+1:nSnow+nLayersRoots)) / real(nLayersRoots, kind(dp))

  ! compute the exposed LAI and SAI and whether veg is buried by snow
  call vegPhenlgy(&
                  ! input/output: data structures
                  model_decisions,             & ! intent(in):    model decisions
                  type_data,                   & ! intent(in):    type of vegetation and soil
                  attr_data,                   & ! intent(in):    spatial attributes
                  mpar_data,                   & ! intent(in):    model parameters
                  mvar_data,                   & ! intent(inout): model variables for a local HRU
                  ! output
                  computeVegFlux,              & ! intent(out): flag to indicate if we are computing fluxes over vegetation (.false. means veg is buried with snow)
                  canopyDepth,                 & ! intent(out): canopy depth (m)
                  exposedVAI,                  & ! intent(out): exposed vegetation area index (m2 m-2)
                  err,cmessage)                  ! intent(out): error control
  if(err/=0)then; err=20; message=trim(message)//trim(cmessage); return; endif


  ! (2) compute wetted canopy area...
  ! ---------------------------------

  ! compute maximum canopy liquid water (kg m-2)
  mvar_data%var(iLookMVAR%scalarCanopyLiqMax)%dat(1) = mpar_data%var(iLookPARAM%refInterceptCapRain)*exposedVAI

  ! compute maximum canopy ice content (kg m-2)
  ! NOTE 1: this is used to compute the snow fraction on the canopy, as used in *BOTH* the radiation AND canopy sublimation routines
  ! NOTE 2: this is a different variable than the max ice used in the throughfall (snow interception) calculations
  select case(model_decisions(iLookDECISIONS%snowIncept)%iDecision)
   case(lightSnow);  mvar_data%var(iLookMVAR%scalarCanopyIceMax)%dat(1) = exposedVAI*mpar_data%var(iLookPARAM%refInterceptCapSnow)       ! use maximum per unit leaf area storage capacity for snow (kg m-2)
   case(stickySnow); mvar_data%var(iLookMVAR%scalarCanopyIceMax)%dat(1) = exposedVAI*mpar_data%var(iLookPARAM%refInterceptCapSnow)*4._dp ! use maximum per unit leaf area storage capacity for snow (kg m-2)
   case default; message=trim(message)//'unable to identify option for maximum branch interception capacity'; err=20; return
  end select ! identifying option for maximum branch interception capacity
  !print*, 'mvar_data%var(iLookMVAR%scalarCanopyLiqMax)%dat(1) = ', mvar_data%var(iLookMVAR%scalarCanopyLiqMax)%dat(1)
  !print*, 'mvar_data%var(iLookMVAR%scalarCanopyIceMax)%dat(1) = ', mvar_data%var(iLookMVAR%scalarCanopyIceMax)%dat(1)

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

  ! vegetation is completely buried by snow (or no veg exisits at all)
  else
   mvar_data%var(iLookMVAR%scalarCanopyWetFraction)%dat(1) = 0._dp
   dCanopyWetFraction_dWat                                 = 0._dp
   dCanopyWetFraction_dT                                   = 0._dp
  endif


  ! (3) compute snow albedo...
  ! --------------------------
  ! NOTE: this should be done before the radiation calculations
  ! NOTE: uses snowfall; should really use canopy throughfall + canopy unloading
  call snowAlbedo(&
                  ! input: model control
                  dt_sub,                      & ! intent(in): model time step (s)
                  (nSnow > 0),                 & ! intent(in): logical flag to denote if snow is present
                  ! input/output: data structures
                  model_decisions,             & ! intent(in):    model decisions
                  mpar_data,                   & ! intent(in):    model parameters
                  mvar_data,                   & ! intent(inout): model variables for a local HRU
                  ! output: error control
                  err,cmessage)                  ! intent(out): error control
  if(err/=0)then; err=20; message=trim(message)//trim(cmessage); return; endif


  ! (4) compute canopy sw radiation fluxes...
  ! -----------------------------------------
  call vegSWavRad(&
                  dt_sub,                       & ! intent(in): time step (s) -- only used in Noah-MP radiation, to compute albedo
                  computeVegFlux,               & ! intent(in): logical flag to compute vegetation fluxes (.false. if veg buried by snow)
                  err,cmessage)                   ! intent(out): error control
  if(err/=0)then; err=20; message=trim(message)//trim(cmessage); return; endif

  print*, 'dt_done=', dt_done, 'summerLAI=', mpar_data%var(48), ' scalarBelowCanopySolar=', mvar_data%var(64)%dat

  ! (5) compute canopy throughfall and unloading...
  ! -----------------------------------------------
  ! NOTE 1: this needs to be done before solving the energy and liquid water equations, to account for the heat advected with precipitation (and throughfall/unloading)
  ! NOTE 2: the unloading flux is computed using canopy drip (scalarCanopyLiqDrainage) from the previous time step
  call canopySnow(&
                  ! input: model control
                  dt_sub,                      & ! intent(in): time step (seconds)
                  exposedVAI,                  & ! intent(in): exposed vegetation area index (m2 m-2)
                  computeVegFlux,              & ! intent(in): flag to denote if computing energy flux over vegetation
                  ! input/output: data structures
                  model_decisions,             & ! intent(in):    model decisions
                  forc_data,                   & ! intent(in):    model forcing data
                  mpar_data,                   & ! intent(in):    model parameters
                  mvar_data,                   & ! intent(inout): model variables for a local HRU
                  ! output: error control
                  err,cmessage)                  ! intent(out): error control
  if(err/=0)then; err=20; message=trim(message)//trim(cmessage); return; endif
  !print*, 'canopyIce = ', mvar_data%var(iLookMVAR%scalarCanopyIce)%dat(1)

  ! adjust canopy temperature to account for new snow
  call tempAdjust(&
                  ! input: derived parameters
                  canopyDepth,                 & ! intent(in): canopy depth (m)
                  ! input/output: data structures
                  mpar_data,                   & ! intent(in):    model parameters
                  mvar_data,                   & ! intent(inout): model variables for a local HRU
                  ! output: error control
                  err,cmessage)                  ! intent(out): error control
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

  ! (6) add snowfall to the snowpack...
  ! -----------------------------------

  ! add new snowfall to the snowpack
  ! NOTE: This needs to be done AFTER the call to canopySnow, since throughfall and unloading are computed in canopySnow
  call newsnwfall(&
                 ! input: model control
                 dt_sub,                                                    & ! time step (seconds)
                 (nSnow > 0),                                               & ! logical flag if snow layers exist
                 mpar_data%var(iLookPARAM%snowfrz_scale),                   & ! freeezing curve parameter for snow (K-1)
                 ! input: diagnostic scalar variables
                 mvar_data%var(iLookMVAR%scalarSnowfallTemp)%dat(1),        & ! computed temperature of fresh snow (K)
                 mvar_data%var(iLookMVAR%scalarNewSnowDensity)%dat(1),      & ! computed density of new snow (kg m-3)
                 mvar_data%var(iLookMVAR%scalarThroughfallSnow)%dat(1),     & ! throughfall of snow through the canopy (kg m-2 s-1)
                 mvar_data%var(iLookMVAR%scalarCanopySnowUnloading)%dat(1), & ! unloading of snow from the canopy (kg m-2 s-1)
                 ! input/output: state variables
                 mvar_data%var(iLookMVAR%scalarSWE)%dat(1),                 & ! SWE (kg m-2)
                 mvar_data%var(iLookMVAR%scalarSnowDepth)%dat(1),           & ! total snow depth (m)
                 mvar_data%var(iLookMVAR%mLayerTemp)%dat(1),                & ! temperature of the top layer (K)
                 mvar_data%var(iLookMVAR%mLayerDepth)%dat(1),               & ! depth of the top layer (m)
                 mvar_data%var(iLookMVAR%mLayerVolFracIce)%dat(1),          & ! volumetric fraction of ice of the top layer (-)
                 mvar_data%var(iLookMVAR%mLayerVolFracLiq)%dat(1),          & ! volumetric fraction of liquid water of the top layer (-)
                 ! output: error control
                 err,cmessage)                                                ! error control
  if(err/=0)then; err=30; message=trim(message)//trim(cmessage); return; endif

  ! re-compute snow depth and SWE
  if(nSnow > 0)then
   mvar_data%var(iLookMVAR%scalarSnowDepth)%dat(1) = sum(mvar_data%var(iLookMVAR%mLayerDepth)%dat(1:nSnow))
   mvar_data%var(iLookMVAR%scalarSWE)%dat(1)       = sum( (mvar_data%var(iLookMVAR%mLayerVolFracLiq)%dat(1:nSnow)*iden_water + &
                                                           mvar_data%var(iLookMVAR%mLayerVolFracIce)%dat(1:nSnow)*iden_ice) &
                                                           * mvar_data%var(iLookMVAR%mLayerDepth)%dat(1:nSnow) )
  endif
  !print*, 'SWE after snowfall = ',  mvar_data%var(iLookMVAR%scalarSWE)%dat(1)

  ! update coordinate variables
  call calcHeight(&
                  ! input/output: data structures
                  indx_data,   & ! intent(in): layer type
                  mvar_data,   & ! intent(inout): model variables for a local HRU
                  ! output: error control
                  err,cmessage)
  if(err/=0)then; err=20; message=trim(message)//trim(cmessage); return; endif

  ! ****************************************************************************************************
  ! *** MAIN SOLVER ************************************************************************************
  ! ****************************************************************************************************

  ! initialize dt
  ntemp   = 0       ! number of temporary sub-steps
  ntrial  = 0       ! number of trial sub-steps
  dt_solv = 0._dp   ! progress towards dt_sub
  dt_temp = dt_sub  ! temporary substep

  ! intialize variables needed for SWE mass balance check
  effRainfall = 0._dp  ! if no snow layers, water is added to the top of the soil zone
  snwDrainage = 0._dp  ! no snow drainage when no snow layers
  sublimation = 0._dp  ! no sublimation when no snow layers
  sfcMeltPond = 0._dp  ! surface melt pond

  ! initialize the rejected step
  rejectedStep=.false.  ! always try the first time

  ! ** continuous do loop to handle any non-convergence or mass balance issues that arise
  do  ! (multiple attempts for non-convergence etc.; minstep check to avoid excessive iteration)

   ! increment trial sub-steps
   ntrial = ntrial+1

   ! if step is rejected, then no need to revise layer structure etc.
   if(.not.rejectedStep)then

    ! (7) merge/sub-divide snow layers...
    ! -----------------------------------
    call volicePack(&
                    ! input/output: model data structures
                    model_decisions,             & ! intent(in):    model decisions
                    mpar_data,                   & ! intent(in):    model parameters
                    indx_data,                   & ! intent(inout): type of each layer
                    mvar_data,                   & ! intent(inout): model variables for a local HRU
                    ! output: error control
                    err,cmessage)                  ! intent(out): error control
    if(err/=0)then; err=55; message=trim(message)//trim(cmessage); return; endif

    ! recompute the number of snow and soil layers
    ! NOTE: do this here for greater visibility
    nSnow   = count(indx_data%var(iLookINDEX%layerType)%dat==ix_snow)
    nSoil   = count(indx_data%var(iLookINDEX%layerType)%dat==ix_soil)
    nLayers = nSnow+nSoil

    ! put the data in the structures
    indx_data%var(iLookINDEX%nSnow)%dat(1)   = nSnow
    indx_data%var(iLookINDEX%nSoil)%dat(1)   = nSoil
    indx_data%var(iLookINDEX%nLayers)%dat(1) = nLayers

    ! re-compute snow depth and SWE
    if(nSnow > 0)then
     mvar_data%var(iLookMVAR%scalarSnowDepth)%dat(1) = sum(mvar_data%var(iLookMVAR%mLayerDepth)%dat(1:nSnow))
     mvar_data%var(iLookMVAR%scalarSWE)%dat(1)       = sum( (mvar_data%var(iLookMVAR%mLayerVolFracLiq)%dat(1:nSnow)*iden_water + &
                                                             mvar_data%var(iLookMVAR%mLayerVolFracIce)%dat(1:nSnow)*iden_ice) &
                                                             * mvar_data%var(iLookMVAR%mLayerDepth)%dat(1:nSnow) )
    endif


    ! (7) compute diagnostic variables for each layer...
    ! --------------------------------------------------
    ! NOTE: this needs to be done AFTER volicePack, since layers may have been sub-divided and/or merged
    call diagn_evar(&
                    ! input: control variables
                    computeVegFlux,          & ! intent(in): flag to denote if computing the vegetation flux
                    canopyDepth,             & ! intent(in): canopy depth (m)
                    ! input/output: data structures
                    mpar_data,               & ! intent(in):    model parameters
                    indx_data,               & ! intent(in):    model layer indices
                    mvar_data,               & ! intent(inout): model variables for a local HRU
                    ! output: error control
                    err,cmessage)              ! intent(out): error control
    if(err/=0)then; err=55; message=trim(message)//trim(cmessage); return; endif


    ! (8) compute melt of the "snow without a layer"...
    ! -------------------------------------------------
    ! NOTE: forms a surface melt pond, which drains into the upper-most soil layer through the time step
    ! (check for the special case of "snow without a layer")
    if(nSnow==0)then
     call implctMelt(&
                     ! input/output: integrated snowpack properties
                     mvar_data%var(iLookMVAR%scalarSWE)%dat(1),               & ! intent(inout): snow water equivalent (kg m-2)
                     mvar_data%var(iLookMVAR%scalarSnowDepth)%dat(1),         & ! intent(inout): snow depth (m)
                     mvar_data%var(iLookMVAR%scalarSfcMeltPond)%dat(1),       & ! intent(inout): surface melt pond (kg m-2)
                     ! input/output: properties of the upper-most soil layer
                     mvar_data%var(iLookMVAR%mLayerTemp)%dat(nSnow+1),        & ! intent(inout): surface layer temperature (K)
                     mvar_data%var(iLookMVAR%mLayerDepth)%dat(nSnow+1),       & ! intent(inout): surface layer depth (m)
                     mvar_data%var(iLookMVAR%mLayerVolHtCapBulk)%dat(nSnow+1),& ! intent(inout): surface layer volumetric heat capacity (J m-3 K-1)
                     ! output: error control
                     err,cmessage                                             ) ! intent(out): error control
     if(err/=0)then; err=20; message=trim(message)//trim(cmessage); return; endif
    endif

   ! ** if previous step is not rejected
   endif

   ! test: recompute snow depth and SWE
   if(nSnow > 0)then
    mvar_data%var(iLookMVAR%scalarSnowDepth)%dat(1) = sum(mvar_data%var(iLookMVAR%mLayerDepth)%dat(1:nSnow))
    mvar_data%var(iLookMVAR%scalarSWE)%dat(1)       = sum( (mvar_data%var(iLookMVAR%mLayerVolFracLiq)%dat(1:nSnow)*iden_water + &
                                                            mvar_data%var(iLookMVAR%mLayerVolFracIce)%dat(1:nSnow)*iden_ice) &
                                                            * mvar_data%var(iLookMVAR%mLayerDepth)%dat(1:nSnow) )
   endif
   !write(*,'(a,1x,2(f20.5,1x),l1)') 'b4 systemSolv: testSWE, meltPond, rejectedStep = ', mvar_data%var(iLookMVAR%scalarSWE)%dat(1), mvar_data%var(iLookMVAR%scalarSfcMeltPond)%dat(1), rejectedStep

   ! print progress
   !if(dt_temp < dt_sub-tinyNumber)then
   ! write(*,'(a,1x,i4,1x,3(f15.2,1x))') 'ntrial, dt_temp, dt_solv, dt_sub = ', ntrial, dt_temp, dt_solv, dt_sub
   ! !pause
   !endif

   ! (9) solve model equations...
   ! ----------------------------

   ! get the new solution
   call systemSolv(&
                   ! input: model control
                   dt_temp,                                & ! intent(in): length of the model sub-step
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

   ! check for fatal errors
   if(err>0)then; err=20; message=trim(message)//trim(cmessage); return; endif
   !print*, 'after solv: mvar_data%var(iLookMVAR%iLayerLiqFluxSoil)%dat(0)*dt_temp*dt_temp/dt = ', mvar_data%var(iLookMVAR%iLayerLiqFluxSoil)%dat(0)*dt_temp*dt_temp/dt

   ! test: recompute snow depth and SWE
   if(nSnow > 0)then
    mvar_data%var(iLookMVAR%scalarSnowDepth)%dat(1) = sum(mvar_data%var(iLookMVAR%mLayerDepth)%dat(1:nSnow))
    mvar_data%var(iLookMVAR%scalarSWE)%dat(1)       = sum( (mvar_data%var(iLookMVAR%mLayerVolFracLiq)%dat(1:nSnow)*iden_water + &
                                                            mvar_data%var(iLookMVAR%mLayerVolFracIce)%dat(1:nSnow)*iden_ice) &
                                                            * mvar_data%var(iLookMVAR%mLayerDepth)%dat(1:nSnow) )
   endif
   !write(*,'(a,1x,2(i4,1x),10(f15.5,1x))') 'nTrial, nTemp, dt_temp, dt_solv, dt_sub, test SWE = ', &
   !                                         nTrial, nTemp, dt_temp, dt_solv, dt_sub, mvar_data%var(iLookMVAR%scalarSWE)%dat(1)

   ! if err<0 (warnings) and hence non-convergence
   if(err<0)then
    ! (adjust time step length)
    dt_temp = dt_temp*0.5_dp ! halve the sub-step
    write(*,'(a,1x,2(f13.3,1x))') trim(cmessage), dt_temp, minstep
    rejectedStep=.true.
    ! (check that time step greater than the minimum step)
    if(dt_temp < minstep)then
     message=trim(message)//'dt_temp is below the minimum time step'
     err=20; return
    endif
    !pause 'failed step'
    ! (try again)
    cycle  ! try again
   else
    rejectedStep=.false.
    !pause 'accepted step'
   endif

   ! check that err=0 at this point (it should be)
   if(err/=0)then; message=trim(message)//'expect err=0'; return; endif


   ! (10a) compute change in canopy ice content due to sublimation...
   ! --------------------------------------------------------------
   ! NOTE: keep in continuous do loop in case insufficient water on canopy for sublimation
   if(computeVegFlux)then

    ! remove mass of ice on the canopy
    mvar_data%var(iLookMVAR%scalarCanopyIce)%dat(1) = mvar_data%var(iLookMVAR%scalarCanopyIce)%dat(1) + &
                                                      mvar_data%var(iLookMVAR%scalarCanopySublimation)%dat(1)*dt_temp

    ! if removed all ice, take the remaining sublimation from water
    if(mvar_data%var(iLookMVAR%scalarCanopyIce)%dat(1) < 0._dp)then
     mvar_data%var(iLookMVAR%scalarCanopyLiq)%dat(1) = mvar_data%var(iLookMVAR%scalarCanopyLiq)%dat(1) + mvar_data%var(iLookMVAR%scalarCanopyIce)%dat(1)
     mvar_data%var(iLookMVAR%scalarCanopyIce)%dat(1) = 0._dp
    endif

    ! check that there is sufficient canopy water to support the converged sublimation rate over the time step dt_temp
    ! NOTE we conducted checks and time step adjustments in systemSolv above so we should not get here: hence fatal error
    if(mvar_data%var(iLookMVAR%scalarCanopyLiq)%dat(1) < -tinyNumber)then
     message=trim(message)//'canopy sublimation rate over time step dt_temp depletes more than the available water'
     err=20; return
    endif

   endif  ! (if computing the vegetation flux)


   ! (10b) compute change in ice content of the top snow layer due to sublimation...
   ! -----------------------------------------------------------------------------
   ! NOTE: this is done BEFORE densification
   if(nSnow > 0)then ! snow layers exist

    ! compute volumetric sublimation (-)
    volSub = dt_temp*mvar_data%var(iLookMVAR%scalarSnowSublimation)%dat(1)/mvar_data%var(iLookMVAR%mLayerDepth)%dat(1)

    ! update volumetric fraction of ice (-)
    ! NOTE: fluxes are positive downward
    mvar_data%var(iLookMVAR%mLayerVolFracIce)%dat(1) = mvar_data%var(iLookMVAR%mLayerVolFracIce)%dat(1) + volSub/iden_ice

    ! check that there is sufficient ice in the top snow layer to support the converged sublimation rate over the time step dt_temp
    ! NOTE we conducted checks and time step adjustments in systemSolv above so we should not get here: hence fatal error
    if(mvar_data%var(iLookMVAR%mLayerVolFracIce)%dat(1) < -tinyNumber)then
     message=trim(message)//'surface sublimation rate over time step dt_temp depletes more than the available water'
     err=20; return
    endif

   ! no snow
   else

    ! no snow: check that sublimation is zero
    if(abs(mvar_data%var(iLookMVAR%scalarSnowSublimation)%dat(1)) > verySmall)then
     message=trim(message)//'sublimation of snow has been computed when no snow exists'
     err=20; return
    endif

   endif  ! (if snow layers exist)
   !print*, 'ice after sublimation: ', mvar_data%var(iLookMVAR%mLayerVolFracIce)%dat(1)*iden_ice


   ! (11) account for compaction and cavitation in the snowpack...
   ! ------------------------------------------------------------
   if(nSnow>0)then
    call snwDensify(&
                    ! intent(in): variables
                    dt_temp,                                                & ! intent(in): time step (s)
                    mvar_data%var(iLookMVAR%mLayerTemp)%dat(1:nSnow),       & ! intent(in): temperature of each layer (K)
                    mvar_data%var(iLookMVAR%mLayerMeltFreeze)%dat(1:nSnow), & ! intent(in): volumetric melt in each layer (kg m-3)
                    mvar_data%var(iLookMVAR%scalarSnowSublimation)%dat(1),  & ! intent(in): sublimation from the snow surface (kg m-2 s-1)
                    ! intent(in): parameters
                    mpar_data%var(iLookPARAM%densScalGrowth),               & ! intent(in): density scaling factor for grain growth (kg-1 m3)
                    mpar_data%var(iLookPARAM%tempScalGrowth),               & ! intent(in): temperature scaling factor for grain growth (K-1)
                    mpar_data%var(iLookPARAM%grainGrowthRate),              & ! intent(in): rate of grain growth (s-1)
                    mpar_data%var(iLookPARAM%densScalOvrbdn),               & ! intent(in): density scaling factor for overburden pressure (kg-1 m3)
                    mpar_data%var(iLookPARAM%tempScalOvrbdn),               & ! intent(in): temperature scaling factor for overburden pressure (K-1)
                    mpar_data%var(iLookPARAM%base_visc),                    & ! intent(in): viscosity coefficient at T=T_frz and snow density=0 (kg m-2 s)
                    ! intent(inout): state variables
                    mvar_data%var(iLookMVAR%mLayerDepth)%dat(1:nSnow),      & ! intent(inout): depth of each layer (m)
                    mvar_data%var(iLookMVAR%mLayerVolFracLiq)%dat(1:nSnow), & ! intent(inout):  volumetric fraction of liquid water after itertations (-)
                    mvar_data%var(iLookMVAR%mLayerVolFracIce)%dat(1:nSnow), & ! intent(inout):  volumetric fraction of ice after itertations (-)
                    ! output: error control
                    err,cmessage)                     ! intent(out): error control
    if(err/=0)then; err=55; message=trim(message)//trim(cmessage); return; endif
   endif  ! if snow layers exist

   ! update coordinate variables
   call calcHeight(&
                   ! input/output: data structures
                   indx_data,   & ! intent(in): layer type
                   mvar_data,   & ! intent(inout): model variables for a local HRU
                   ! output: error control
                   err,cmessage)
   if(err/=0)then; err=20; message=trim(message)//trim(cmessage); return; endif


   ! (12) compute sub-step averages associated with the temporary steps...
   ! ---------------------------------------------------------------------

   ! keep track of the number of temporary sub-steps
   ntemp = ntemp+1

   ! increment model fluxes
   dt_wght = dt_temp/dt ! define weight applied to each sub-step
   call increment_fluxes(dt_wght,(nsub==1 .and. ntemp==1))

   !dt_prog = dt_done+dt_solv+dt_temp  ! progress in time step (s)
   !dt_frac = dt_prog/dt               ! fraction of time step completed (-)
   !write(*,'(a,1x,3(f9.3,1x),10(e20.10,1x))') 'dt_wght, dt_prog, dt_frac, totalSoilCompress, dt_prog*averageSoilInflux/dt_frac, scalarSoilCompress, scalarSoilInflux*dt_temp = ', &
   !                                            dt_wght, dt_prog, dt_frac, totalSoilCompress, dt_prog*averageSoilInflux/dt_frac, scalarSoilCompress, mvar_data%var(iLookMVAR%iLayerLiqFluxSoil)%dat(0)*dt_temp
   !pause

   ! compute effective rainfall input and snowpack drainage to/from the snowpack (kg m-2 s-1)
   if(nSnow > 0)then
    effRainfall = effRainfall + (mvar_data%var(iLookMVAR%scalarThroughfallRain)%dat(1) + mvar_data%var(iLookMVAR%scalarCanopyLiqDrainage)%dat(1) )*dt_temp/dt_sub
    snwDrainage = snwDrainage + (mvar_data%var(iLookMVAR%iLayerLiqFluxSnow)%dat(nSnow)*iden_water )*dt_temp/dt_sub ! m s-1 -> kg m-2 s-1
    sublimation = sublimation + (mvar_data%var(iLookMVAR%scalarSnowSublimation)%dat(1))*dt_temp/dt_sub
   ! compute the surface melt pond (kg m-2)
   else
    sfcMeltPond = sfcMeltPond + mvar_data%var(iLookMVAR%scalarSfcMeltPond)%dat(1)
   endif

   ! increment sub-step
   dt_solv = dt_solv + dt_temp

   ! check that we have completed the sub-step
   if(dt_solv >= dt_sub-verySmall) exit

   ! adjust length of the sub-step (make sure that we don't exceed the step)
   dt_temp = min(dt_sub - dt_solv, dt_temp)
   !print*, 'dt_temp, dt_sub = ', dt_temp, dt_sub

  end do  ! (multiple attempts for non-convergence)
  !print*, 'after do loop: dt_sub = ', dt_sub

  ! ****************************************************************************************************
  ! *** END MAIN SOLVER ********************************************************************************
  ! ****************************************************************************************************

  ! (13) check energy and mass balance...
  ! -------------------------------------

  ! recompute snow depth and SWE
  if(nSnow > 0)then
   mvar_data%var(iLookMVAR%scalarSnowDepth)%dat(1) = sum(mvar_data%var(iLookMVAR%mLayerDepth)%dat(1:nSnow))
   mvar_data%var(iLookMVAR%scalarSWE)%dat(1)       = sum( (mvar_data%var(iLookMVAR%mLayerVolFracLiq)%dat(1:nSnow)*iden_water + &
                                                           mvar_data%var(iLookMVAR%mLayerVolFracIce)%dat(1:nSnow)*iden_ice) &
                                                           * mvar_data%var(iLookMVAR%mLayerDepth)%dat(1:nSnow) )
  endif

  ! check SWE
  effSnowfall = mvar_data%var(iLookMVAR%scalarThroughfallSnow)%dat(1) + mvar_data%var(iLookMVAR%scalarCanopySnowUnloading)%dat(1)
  newSWE      = mvar_data%var(iLookMVAR%scalarSWE)%dat(1)
  delSWE      = newSWE - (oldSWE - sfcMeltPond)
  massBalance = delSWE - (effSnowfall + effRainfall + sublimation - snwDrainage)*dt_sub
  if(abs(massBalance) > 1.d-6)then
   print*,                  'nSnow       = ', nSnow
   print*,                  'nTemp       = ', nTemp
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

  ! (14) adjust length of the substep...
  ! ------------------------------------

  ! increment the time step
  dt_done = dt_done + dt_sub
  !print*, '***** ', dt_done, dt_sub, niter
  !pause ' after increment the time step'

  ! modify the length of the time step
  if(niter<n_inc) dt_sub = min(dt_temp*F_inc,maxstep)
  if(niter>n_dec) dt_sub =     dt_temp*F_dec
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
  !print*, 'in substep loop: dt_sub = ', dt_sub

 end do  ! (sub-step loop)
 !stop 'completed time step'

 ! ---
 ! (14) balance checks...
 ! ----------------------

 ! get total canopy water
 scalarCanopyLiq => mvar_data%var(iLookMVAR%scalarCanopyLiq)%dat(1)
 scalarCanopyIce => mvar_data%var(iLookMVAR%scalarCanopyIce)%dat(1)
 balanceCanopyWater1 = scalarCanopyLiq + scalarCanopyIce

 ! get snowfall and rainfall
 scalarSnowfall => mvar_data%var(iLookMVAR%scalarSnowfall)%dat(1)        ! computed snowfall rate (kg m-2 s-1)
 scalarRainfall => mvar_data%var(iLookMVAR%scalarRainfall)%dat(1)        ! computed rainfall rate (kg m-2 s-1)

 ! print progress
 !write(*,'(a,1x,f20.10)') 'balanceCanopyWater0                          = ', balanceCanopyWater0
 !write(*,'(a,1x,f20.10)') 'balanceCanopyWater1                          = ', balanceCanopyWater1
 !write(*,'(a,1x,f20.10)') '(scalarSnowfall - averageThroughfallSnow)*dt = ', (scalarSnowfall - averageThroughfallSnow)*dt
 !write(*,'(a,1x,f20.10)') '(scalarRainfall - averageThroughfallRain)*dt = ', (scalarRainfall - averageThroughfallRain)*dt
 !write(*,'(a,1x,f20.10)') 'averageCanopySnowUnloading                   = ', averageCanopySnowUnloading*dt
 !write(*,'(a,1x,f20.10)') 'averageCanopyLiqDrainage                     = ', averageCanopyLiqDrainage*dt
 !write(*,'(a,1x,f20.10)') 'averageCanopySublimation                     = ', averageCanopySublimation*dt
 !write(*,'(a,1x,f20.10)') 'averageCanopyEvaporation                     = ', averageCanopyEvaporation*dt

 ! balance checks for the canopy
 ! NOTE: need to put the balance checks in the sub-step loop so that we can re-compute if necessary
 scalarCanopyWatBalError = balanceCanopyWater1 - (balanceCanopyWater0 + (scalarSnowfall - averageThroughfallSnow)*dt + (scalarRainfall - averageThroughfallRain)*dt &
                            - averageCanopySnowUnloading*dt - averageCanopyLiqDrainage*dt + averageCanopySublimation*dt + averageCanopyEvaporation*dt)
 if(abs(scalarCanopyWatBalError) > 1.d-1)then
  print*, '** canopy water balance error:'
  write(*,'(a,1x,f20.10)') 'dt                                           = ', dt
  write(*,'(a,1x,f20.10)') 'balanceCanopyWater0                          = ', balanceCanopyWater0
  write(*,'(a,1x,f20.10)') 'balanceCanopyWater1                          = ', balanceCanopyWater1
  write(*,'(a,1x,f20.10)') '(scalarSnowfall - averageThroughfallSnow)*dt = ', (scalarSnowfall - averageThroughfallSnow)*dt
  write(*,'(a,1x,f20.10)') '(scalarRainfall - averageThroughfallRain)*dt = ', (scalarRainfall - averageThroughfallRain)*dt
  write(*,'(a,1x,f20.10)') 'averageCanopySnowUnloading                   = ', averageCanopySnowUnloading*dt
  write(*,'(a,1x,f20.10)') 'averageCanopyLiqDrainage                     = ', averageCanopyLiqDrainage*dt
  write(*,'(a,1x,f20.10)') 'averageCanopySublimation                     = ', averageCanopySublimation*dt
  write(*,'(a,1x,f20.10)') 'averageCanopyEvaporation                     = ', averageCanopyEvaporation*dt
  write(*,'(a,1x,f20.10)') 'scalarCanopyWatBalError                      = ', scalarCanopyWatBalError
  message=trim(message)//'canopy hydrology does not balance'
  err=20; return
 endif
 !pause 'canopy hydrology does balance'

 ! point to model state variables
 ! NOTE: need to do this at the end of each sub-step because number of layers may change
 mLayerDepth          => mvar_data%var(iLookMVAR%mLayerDepth)%dat(nSnow+1:nLayers)           ! depth of each soil layer (m)
 mLayerVolFracIce     => mvar_data%var(iLookMVAR%mLayerVolFracIce)%dat(nSnow+1:nLayers)      ! volumetric ice content in each soil layer (-)
 mLayerVolFracLiq     => mvar_data%var(iLookMVAR%mLayerVolFracLiq)%dat(nSnow+1:nLayers)      ! volumetric liquid water content in each soil layer (-)
 scalarAquiferStorage => mvar_data%var(iLookMVAR%scalarAquiferStorage)%dat(1)                ! aquifer storage (m)

 ! compute the liquid water and ice content at the end of the time step
 scalarTotalSoilLiq = sum(iden_water*mLayerVolFracLiq(1:nSoil)*mLayerDepth(1:nSoil))
 scalarTotalSoilIce = sum(iden_water*mLayerVolFracIce(1:nSoil)*mLayerDepth(1:nSoil))   ! NOTE: no expansion of soil, hence use iden_water
 !scalarTotalSoilIce = sum(iden_ice  *mLayerVolFracIce(1:nSoil)*mLayerDepth(1:nSoil))

 ! get the total water in the soil (liquid plus ice) at the end of the time step (kg m-2)
 balanceSoilWater1 = scalarTotalSoilLiq + scalarTotalSoilIce

 ! get the total aquifer storage at the start of the time step (kg m-2)
 balanceAquifer1 = scalarAquiferStorage*iden_water

 ! get the input and output to/from the soil zone (kg m-2)
 balanceSoilInflux        = averageSoilInflux*iden_water*dt
 balanceSoilBaseflow      = averageSoilBaseflow*iden_water*dt
 balanceSoilDrainage      = averageSoilDrainage*iden_water*dt
 balanceSoilTranspiration = averageCanopyTranspiration*dt      ! NOTE ground evaporation included in the flux at the upper boundary

 !write(*,'(a,1x,f20.10)') 'totalSoilCompress          = ', totalSoilCompress                  ! kg m-2
 !write(*,'(a,1x,f20.10)') 'averageSoilInflux          = ', averageSoilInflux*iden_water*dt    ! m s-1  -> kg m-2
 !write(*,'(a,1x,f20.10)') 'averageSoilBaseflow        = ', averageSoilBaseflow*iden_water*dt
 !write(*,'(a,1x,f20.10)') 'averageSoilDrainage        = ', averageSoilDrainage*iden_water*dt
 !write(*,'(a,1x,f20.10)') 'averageCanopyTranspiration = ', averageCanopyTranspiration*dt
 !write(*,'(a,1x,f20.10)') 'averageGroundEvaporation   = ', averageGroundEvaporation*dt

 !print*, 'sum(mLayerBaseflow) = ', sum(mvar_data%var(iLookMVAR%mLayerBaseflow)%dat)

 ! check the soil water balance
 scalarSoilWatBalError  = balanceSoilWater1 - (balanceSoilWater0 + (balanceSoilInflux + balanceSoilTranspiration - balanceSoilBaseflow - balanceSoilDrainage - totalSoilCompress) )
 if(abs(scalarSoilWatBalError) > 1.d-2)then  ! NOTE: kg m-2, so need coarse tolerance to account for precision issues
  write(*,'(a,1x,f20.10)') 'dt                        = ', dt
  write(*,'(a,1x,f20.10)') 'totalSoilCompress         = ', totalSoilCompress
  write(*,'(a,1x,f20.10)') 'balanceSoilWater0         = ', balanceSoilWater0
  write(*,'(a,1x,f20.10)') 'balanceSoilWater1         = ', balanceSoilWater1
  write(*,'(a,1x,f20.10)') 'balanceSoilInflux         = ', balanceSoilInflux
  write(*,'(a,1x,f20.10)') 'balanceSoilBaseflow       = ', balanceSoilBaseflow
  write(*,'(a,1x,f20.10)') 'balanceSoilDrainage       = ', balanceSoilDrainage
  write(*,'(a,1x,f20.10)') 'balanceSoilTranspiration  = ', balanceSoilTranspiration
  write(*,'(a,1x,f20.10)') 'scalarSoilWatBalError     = ', scalarSoilWatBalError
  ! check the water balance in each layer
  do iLayer=1,nSoil
   xCompress = mvar_data%var(iLookMVAR%mLayerCompress)%dat(iLayer)
   xFlux0    = mvar_data%var(iLookMVAR%iLayerLiqFluxSoil)%dat(iLayer-1)*dt
   xFlux1    = mvar_data%var(iLookMVAR%iLayerLiqFluxSoil)%dat(iLayer)*dt
   write(*,'(a,1x,i4,1x,10(e20.10,1x))') 'iLayer, xFlux0, xFlux1, (xFlux1 - xFlux0)/mLayerDepth(iLayer), xCompress = ', &
                                          iLayer, xFlux0, xFlux1, (xFlux1 - xFlux0)/mLayerDepth(iLayer), xCompress
  end do

  message=trim(message)//'soil hydrology does not balance'
  err=20; return
 endif
 !pause 'soil hydrology does balance'

 !print*, 'mvar_data%var(iLookMVAR%averageCanopyLiqDrainage)%dat(1) = ', mvar_data%var(iLookMVAR%averageCanopyLiqDrainage)%dat(1)

 ! save the surface temperature (just to make things easier to visualize)
 mvar_data%var(iLookMVAR%scalarSurfaceTemp)%dat(1) = mvar_data%var(iLookMVAR%mLayerTemp)%dat(1)

 iLayer = nSnow+1
 !print*, 'nsub, mLayerTemp(iLayer), mLayerVolFracIce(iLayer) = ', nsub, mLayerTemp(iLayer), mLayerVolFracIce(iLayer)
 !print*, 'nsub = ', nsub
 if(nsub>2000)then
  message=trim(message)//'number of sub-steps > 2000'
  err=20; return
 endif


 ! ********************************************************************************************************************************
 ! ********************************************************************************************************************************
 ! ********************************************************************************************************************************

 contains


  ! *********************************************************************************************************
  ! internal subroutine increment_fluxes: calculate timestep-average fluxes
  ! *********************************************************************************************************
  subroutine increment_fluxes(dt_wght,initialize)
  real(dp),intent(in)     :: dt_wght    ! weight assigned to sub-step
  logical(lgt),intent(in) :: initialize ! flag to initialize fluxes

  ! set up pointers

  ! assign pointers to timestep-average model fluxes
  if(initialize)then
   totalSoilCompress          => mvar_data%var(iLookMVAR%totalSoilCompress)%dat(1)             ! change in storage associated with compression of the soil matrix (kg m-2)
   averageThroughfallSnow     => mvar_data%var(iLookMVAR%averageThroughfallSnow)%dat(1)        ! snow that reaches the ground without ever touching the canopy (kg m-2 s-1)
   averageThroughfallRain     => mvar_data%var(iLookMVAR%averageThroughfallRain)%dat(1)        ! rain that reaches the ground without ever touching the canopy (kg m-2 s-1)
   averageCanopySnowUnloading => mvar_data%var(iLookMVAR%averageCanopySnowUnloading)%dat(1)    ! unloading of snow from the vegetion canopy (kg m-2 s-1)
   averageCanopyLiqDrainage   => mvar_data%var(iLookMVAR%averageCanopyLiqDrainage)%dat(1)      ! drainage of liquid water from the vegetation canopy (kg m-2 s-1)
   averageCanopyMeltFreeze    => mvar_data%var(iLookMVAR%averageCanopyMeltFreeze)%dat(1)       ! melt/freeze of water stored in the canopy (kg m-2 s-1)
   averageCanopyTranspiration => mvar_data%var(iLookMVAR%averageCanopyTranspiration)%dat(1)    ! canopy transpiration (kg m-2 s-1)
   averageCanopyEvaporation   => mvar_data%var(iLookMVAR%averageCanopyEvaporation)%dat(1)      ! canopy evaporation/condensation (kg m-2 s-1)
   averageCanopySublimation   => mvar_data%var(iLookMVAR%averageCanopySublimation)%dat(1)      ! canopy sublimation/frost (kg m-2 s-1)
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
   totalSoilCompress          = 0._dp  ! change in storage associated with compression of the soil matrix (kg m-2)
   averageThroughfallSnow     = 0._dp  ! snow that reaches the ground without ever touching the canopy (kg m-2 s-1)
   averageThroughfallRain     = 0._dp  ! rain that reaches the ground without ever touching the canopy (kg m-2 s-1)
   averageCanopySnowUnloading = 0._dp  ! unloading of snow from the vegetion canopy (kg m-2 s-1)
   averageCanopyLiqDrainage   = 0._dp  ! drainage of liquid water from the vegetation canopy (kg m-2 s-1)
   averageCanopyMeltFreeze    = 0._dp  ! melt/freeze of water stored in the canopy (kg m-2 s-1)
   averageCanopyTranspiration = 0._dp  ! canopy transpiration (kg m-2 s-1)
   averageCanopyEvaporation   = 0._dp  ! canopy evaporation/condensation (kg m-2 s-1)
   averageCanopySublimation   = 0._dp  ! canopy sublimation/frost (kg m-2 s-1)
   averageSnowSublimation     = 0._dp  ! snow sublimation/frost - below canopy or non-vegetated (kg m-2 s-1)
   averageGroundEvaporation   = 0._dp  ! ground evaporation/condensation - below canopy or non-vegetated (kg m-2 s-1)
   averageRainPlusMelt        = 0._dp  ! rain plus melt, as input to soil before calculating surface runoff (m s-1)
   averageSurfaceRunoff       = 0._dp  ! surface runoff (m s-1)
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
  scalarCanopyTranspiration => mvar_data%var(iLookMVAR%scalarCanopyTranspiration)%dat(1)       ! canopy transpiration (kg m-2 s-1)
  scalarCanopyEvaporation   => mvar_data%var(iLookMVAR%scalarCanopyEvaporation)%dat(1)         ! canopy evaporation/condensation (kg m-2 s-1)
  scalarCanopySublimation   => mvar_data%var(iLookMVAR%scalarCanopySublimation)%dat(1)         ! canopy sublimation/frost (kg m-2 s-1)
  scalarSnowSublimation     => mvar_data%var(iLookMVAR%scalarSnowSublimation)%dat(1)           ! snow sublimation/frost - below canopy or non-vegetated (kg m-2 s-1)
  scalarGroundEvaporation   => mvar_data%var(iLookMVAR%scalarGroundEvaporation)%dat(1)         ! ground evaporation/condensation - below canopy or non-vegetated (kg m-2 s-1)
  scalarRainPlusMelt        => mvar_data%var(iLookMVAR%scalarRainPlusMelt)%dat(1)              ! rain plus melt, as input to soil before calculating surface runoff (m s-1)
  scalarSurfaceRunoff       => mvar_data%var(iLookMVAR%scalarSurfaceRunoff)%dat(1)             ! surface runoff (m s-1)
  scalarSoilInflux          => mvar_data%var(iLookMVAR%iLayerLiqFluxSoil)%dat(0)               ! influx of water at the top of the soil profile (m s-1)
  scalarSoilDrainage        => mvar_data%var(iLookMVAR%iLayerLiqFluxSoil)%dat(nSoil)           ! drainage from the bottom of the soil profile (m s-1)
  scalarSoilCompress        => mvar_data%var(iLookMVAR%scalarSoilCompress)%dat(1)              ! change in storage associated with compression of the soil matrix (kg m-2)
  scalarSoilBaseflow        => mvar_data%var(iLookMVAR%scalarSoilBaseflow)%dat(1)              ! total baseflow from throughout the soil profile (m s-1)
  scalarAquiferRecharge     => mvar_data%var(iLookMVAR%scalarAquiferRecharge)%dat(1)           ! recharge to the aquifer (m s-1)
  scalarAquiferBaseflow     => mvar_data%var(iLookMVAR%scalarAquiferBaseflow)%dat(1)           ! baseflow from the aquifer (m s-1)
  scalarAquiferTranspire    => mvar_data%var(iLookMVAR%scalarAquiferTranspire)%dat(1)          ! transpiration from the aquifer (m s-1)
  mLayerColumnOutflow       => mvar_data%var(iLookMVAR%mLayerColumnOutflow)%dat                ! total outflow from each layer in a given soil column (m3 s-1)

  ! increment storage over the time step
  totalSoilCompress          = totalSoilCompress          + scalarSoilCompress ! NOTE mass not rate  ! change in storage associated with compression of the soil matrix (kg m-2)

  ! increment timestep-average fluxes
  averageThroughfallSnow     = averageThroughfallSnow     + scalarThroughfallSnow     *dt_wght ! snow that reaches the ground without ever touching the canopy (kg m-2 s-1)
  averageThroughfallRain     = averageThroughfallRain     + scalarThroughfallRain     *dt_wght ! rain that reaches the ground without ever touching the canopy (kg m-2 s-1)
  averageCanopySnowUnloading = averageCanopySnowUnloading + scalarCanopySnowUnloading *dt_wght ! unloading of snow from the vegetion canopy (kg m-2 s-1)
  averageCanopyLiqDrainage   = averageCanopyLiqDrainage   + scalarCanopyLiqDrainage   *dt_wght ! drainage of liquid water from the vegetation canopy (kg m-2 s-1)
  averageCanopyMeltFreeze    = averageCanopyMeltFreeze    + scalarCanopyMeltFreeze    *dt_wght ! melt/freeze of water stored in the canopy (kg m-2 s-1)
  averageCanopyTranspiration = averageCanopyTranspiration + scalarCanopyTranspiration *dt_wght ! canopy transpiration (kg m-2 s-1)
  averageCanopyEvaporation   = averageCanopyEvaporation   + scalarCanopyEvaporation   *dt_wght ! canopy evaporation/condensation (kg m-2 s-1)
  averageCanopySublimation   = averageCanopySublimation   + scalarCanopySublimation   *dt_wght ! canopy sublimation/frost (kg m-2 s-1)
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


 end subroutine coupled_em


 ! *********************************************************************************************************
 ! private subroutine implctMelt: compute melt of the "snow without a layer"
 ! *********************************************************************************************************
 subroutine implctMelt(&
                       ! input/output: integrated snowpack properties
                       scalarSWE,         & ! intent(inout): snow water equivalent (kg m-2)
                       scalarSnowDepth,   & ! intent(inout): snow depth (m)
                       scalarSfcMeltPond, & ! intent(inout): surface melt pond (kg m-2)
                       ! input/output: properties of the upper-most soil layer
                       soilTemp,          & ! intent(inout): surface layer temperature (K)
                       soilDepth,         & ! intent(inout): surface layer depth (m)
                       soilHeatcap,       & ! intent(inout): surface layer volumetric heat capacity (J m-3 K-1)
                       ! output: error control
                       err,message        ) ! intent(out): error control
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

 ! *********************************************************************************************************
 ! private subroutine printRestartFile: print a re-start file
 ! *********************************************************************************************************
 subroutine printRestartFile(&
                             output_fileSuffix,& ! intent(in): suffix defining the model experiment
                             dt_init,          & ! intent(in): time step length (s)
                             time_data,        & ! intent(in): model time structures
                             mvar_data,        & ! intent(in): model variables for a local HRU
                             err,message)        ! intent(out): error control
 ! --------------------------------------------------------------------------------------------------------
 ! --------------------------------------------------------------------------------------------------------
 ! access the derived types to define the data structures
 USE data_struc,only:var_i                  ! data vector (i4b)
 USE data_struc,only:var_dlength            ! data vector with variable length dimension (dp)
 ! access named variables defining elements in the data structures
 USE var_lookup,only:iLookMVAR,iLookTIME    ! named variables for structure elements
 ! access file paths
 USE summaFileManager,only:OUTPUT_PATH,OUTPUT_PREFIX      ! define output file
 ! access desired modules
 USE ascii_util_module,only:file_open       ! open file
 implicit none
 ! --------------------------------------------------------------------------------------------------------
 ! input
 character(*),intent(in)         :: output_fileSuffix   ! suffix defining the model experiment
 real(dp),intent(in)             :: dt_init             ! time step length (s)
 type(var_i),intent(in)          :: time_data           ! model time structures
 type(var_dlength),intent(in)    :: mvar_data           ! model variables for a local HRU
 ! output: error control
 integer(i4b),intent(out)        :: err                 ! error code
 character(*),intent(out)        :: message             ! error message
 ! --------------------------------------------------------------------------------------------------------
 ! local variables
 integer(i4b),parameter          :: ixUnit=64           ! file unit
 logical(lgt)                    :: fileOpen            ! flag to denote if the file unit is already used
 character(len=256),parameter    :: filepref='summaRestart'  ! prefix for the restart filename
 character(len=256)              :: timeString          ! string to define the time
 character(len=256)              :: filename            ! name of the restart file
 character(len=256)              :: cmessage            ! error message of downstream routine
 integer(i4b)                    :: iLayer              ! index of the model layer
 real(dp),parameter              :: valueMissing=-999._dp  ! missing value
 ! --------------------------------------------------------------------------------------------------------
 ! initialize error control
 err=0; message='printRestartFile/'

 ! define the time string
 write(timeString,'(a,i4,3(a,i2.2))') '_',time_data%var(iLookTIME%iyyy),'-',time_data%var(iLookTIME%im),'-',time_data%var(iLookTIME%id),'-',time_data%var(iLookTIME%ih)

 ! define the file name
 filename = trim(OUTPUT_PATH)//trim(filepref)//trim(timeString)//trim(output_fileSuffix)//'.txt'
 !print*, trim(filename)
 !pause

 ! check if file unit is open already
 inquire(unit=ixUnit,opened=fileOpen)
 if(fileOpen)then; err=20; message=trim(message)//'file ixUnit is open'; return; endif

 ! open file for writing
 open(ixUnit,file=trim(filename),status="unknown",action="write",iostat=err)
 if(err/=0)then
  message=trim(message)//"OpenError['"//trim(filename)//"']"
  err=20; return
 endif

 ! write a header
 write(ixUnit,'(a)') '! This is a summa re-start file'
 write(ixUnit,'(a)') '! ---------------------------------------------------------------------------------------------------------------------'

 ! write scalar state variables
 write(ixUnit,'(a)') '<start_scalar_icond>'
 write(ixUnit,'(a25,1x,e20.10)') 'dt_init                 ', dt_init                                              ! time step length (s)
 write(ixUnit,'(a25,1x,e20.10)') 'scalarCanopyIce         ', mvar_data%var(iLookMVAR%scalarCanopyIce     )%dat(1) ! canopy ice content (kg m-2)
 write(ixUnit,'(a25,1x,e20.10)') 'scalarCanopyLiq         ', mvar_data%var(iLookMVAR%scalarCanopyLiq     )%dat(1) ! canopy liquid water content (kg m-2)
 write(ixUnit,'(a25,1x,e20.10)') 'scalarCanairTemp        ', mvar_data%var(iLookMVAR%scalarCanairTemp    )%dat(1) ! temperature of the canopy air space (K)
 write(ixUnit,'(a25,1x,e20.10)') 'scalarCanopyTemp        ', mvar_data%var(iLookMVAR%scalarCanopyTemp    )%dat(1) ! temperature of the vegetation canopy (K)
 write(ixUnit,'(a25,1x,e20.10)') 'scalarSnowAlbedo        ', mvar_data%var(iLookMVAR%scalarSnowAlbedo    )%dat(1) ! snow albedo (-)
 write(ixUnit,'(a25,1x,e20.10)') 'scalarSWE               ', mvar_data%var(iLookMVAR%scalarSWE           )%dat(1) ! snow water equivalent (kg m-2)
 write(ixUnit,'(a25,1x,e20.10)') 'scalarSnowDepth         ', mvar_data%var(iLookMVAR%scalarSnowDepth     )%dat(1) ! snow depth (m)
 write(ixUnit,'(a25,1x,e20.10)') 'scalarSfcMeltPond       ', mvar_data%var(iLookMVAR%scalarSfcMeltPond   )%dat(1) ! surface melt pond (kg m-2)
 write(ixUnit,'(a25,1x,e20.10)') 'scalarAquiferStorage    ', mvar_data%var(iLookMVAR%scalarAquiferStorage)%dat(1) ! aquifer storage (m)
 write(ixUnit,'(a)') '<end_scalar_icond>'

 ! make the file easier to read
 write(ixUnit,'(a)') '! ---------------------------------------------------------------------------------------------------------------------'
 write(ixUnit,'(a)') '! ---------------------------------------------------------------------------------------------------------------------'

 ! write layer state variables
 write(ixUnit,'(a)') '<start_layer_icond>'
 write(ixUnit,'(a)') ' layerType            iLayerHeight             mLayerDepth            mLayerTemp     mLayerVolFracIce     mLayerVolFracLiq     mLayerMatricHead'

 ! write state variables for each snow layer
 if(nSnow>0)then
  do iLayer=1,nSnow  ! loop through snow layers
   write(ixUnit,'(a10,2x,2(e22.15,2x),4(e20.10,1x))') '      snow',                                            &
                                                                     mvar_data%var(iLookMVAR%iLayerHeight    )%dat(iLayer-1), &  ! height at the top of the layer (m)
                                                                     mvar_data%var(iLookMVAR%mLayerDepth     )%dat(iLayer),   &  ! depth of each layer (m)
                                                                     mvar_data%var(iLookMVAR%mLayerTemp      )%dat(iLayer),   &  ! temperature of each layer (K)
                                                                     mvar_data%var(iLookMVAR%mLayerVolFracIce)%dat(iLayer),   &  ! volumetric ice content (-)
                                                                     mvar_data%var(iLookMVAR%mLayerVolFracLiq)%dat(iLayer),   &  ! volumetric liquid water content (-)
                                                                     valueMissing                                                ! matric head (missing for snow)
  end do  ! looping through snow layers
 endif  ! if snow layers exist

 ! write state variables for each soil layer
 do iLayer=1,nSoil  ! loop through snow layers
  write(ixUnit,'(a10,2x,2(e22.15,2x),4(e20.10,1x))') '      soil',                                                  &
                                                                    mvar_data%var(iLookMVAR%iLayerHeight    )%dat(nSnow+iLayer-1), &  ! height at the top of the layer (m)
                                                                    mvar_data%var(iLookMVAR%mLayerDepth     )%dat(nSnow+iLayer),   &  ! depth of each layer (m)
                                                                    mvar_data%var(iLookMVAR%mLayerTemp      )%dat(nSnow+iLayer),   &  ! temperature of each layer (K)
                                                                    mvar_data%var(iLookMVAR%mLayerVolFracIce)%dat(nSnow+iLayer),   &  ! volumetric ice content (-)
                                                                    mvar_data%var(iLookMVAR%mLayerVolFracLiq)%dat(nSnow+iLayer),   &  ! volumetric liquid water content (-)
                                                                    mvar_data%var(iLookMVAR%mLayerMatricHead)%dat(iLayer)             ! matric head (m)
 end do  ! looping through soil layers

 ! end definition of layer variables
 write(ixUnit,'(a)') '<end_layer_icond>'

 ! close file
 close(ixUnit)

 end subroutine printRestartFile


end module coupled_em_module

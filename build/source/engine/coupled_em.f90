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
! number of variables
integer(i4b)           :: nSnow                ! number of snow layers
integer(i4b)           :: nSoil                ! number of soil layers
integer(i4b)           :: nLayers              ! total number of layers
integer(i4b)           :: nState               ! total number of state variables
contains


 ! ************************************************************************************************
 ! public subroutine coupled_em: run the coupled energy-mass model for one timestep
 ! ************************************************************************************************
 subroutine coupled_em(&
                       ! model control
                       istep,             & ! intent(in):    index of the model time step
                       hruId,             & ! intent(in):    hruId
                       dt_init,           & ! intent(inout): used to initialize the size of the sub-step
                       computeVegFlux,    & ! intent(inout): flag to indicate if we are computing fluxes over vegetation (.false. means veg is buried with snow)
                       resumeFailSolver,  & ! flag to resume solver when it failed (not converged)
                       ! data structures (input)
                       type_data,         & ! intent(in):    local classification of soil veg etc. for each HRU
                       attr_data,         & ! intent(in):    local attributes for each HRU
                       forc_data,         & ! intent(in):    model forcing data
                       mpar_data,         & ! intent(in):    model parameters
                       bvar_data,         & ! intent(in):    basin-average variables
                       ! data structures (input-output)
                       indx_data,         & ! intent(inout): model indices
                       prog_data,         & ! intent(inout): prognostic variables for a local HRU
                       diag_data,         & ! intent(inout): diagnostic variables for a local HRU
                       flux_data,         & ! intent(inout): model fluxes for a local HRU
                       ! error control
                       err,message)         ! intent(out):   error control
 ! data types
 USE data_types,only:&
                     var_i,               & ! x%var(:)            (i4b)
                     var_d,               & ! x%var(:)            (dp)
                     var_ilength,         & ! x%var(:)%dat        (i4b)
                     var_dlength            ! x%var(:)%dat        (dp)
 ! named variables for parent structures
 USE var_lookup,only:iLookDECISIONS         ! named variables for elements of the decision structure
 USE var_lookup,only:iLookATTR              ! named variables for structure elements
 USE var_lookup,only:iLookTYPE              ! named variables for structure elements
 USE var_lookup,only:iLookPROG              ! named variables for structure elements
 USE var_lookup,only:iLookDIAG              ! named variables for structure elements
 USE var_lookup,only:iLookFLUX              ! named variables for structure elements
 USE var_lookup,only:iLookFORCE             ! named variables for structure elements
 USE var_lookup,only:iLookPARAM             ! named variables for structure elements
 USE var_lookup,only:iLookINDEX             ! named variables for structure elements
 USE globalData,only:ix_soil,ix_snow        ! named variables for snow and soil
 ! named variables for child structures
 USE var_lookup,only:childFLUX_MEAN
 ! global data
 USE globalData,only:data_step              ! time step of forcing data (s)
 USE globalData,only:model_decisions        ! model decision structure
 ! structure allocations
 USE globalData,only:averageFlux_meta       ! metadata on the timestep-average model flux structure
 USE allocspace_module,only:allocLocal      ! allocate local data structures
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
 USE indexState_module,only:indexState      ! define indices for all model state variables and layers
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
 ! model control
 integer(i4b),intent(in)              :: istep                  ! index of model time step
 integer(i4b),intent(in)              :: hruId                  ! hruId
 real(dp),intent(inout)               :: dt_init                ! used to initialize the size of the sub-step
 logical(lgt),intent(inout)           :: computeVegFlux         ! flag to indicate if we are computing fluxes over vegetation (.false. means veg is buried with snow)
 logical(lgt),intent(in)              :: resumeFailSolver       ! flag to indicate continue simulation even solver does not converge
 ! data structures (input)
 type(var_i),intent(in)               :: type_data              ! type of vegetation and soil
 type(var_d),intent(in)               :: attr_data              ! spatial attributes
 type(var_d),intent(in)               :: forc_data              ! model forcing data
 type(var_d),intent(in)               :: mpar_data              ! model parameters
 type(var_dlength),intent(in)         :: bvar_data              ! basin-average model variables
 ! data structures (input-output)
 type(var_ilength),intent(inout)      :: indx_data              ! state vector geometry
 type(var_dlength),intent(inout)      :: prog_data              ! prognostic variables for a local HRU
 type(var_dlength),intent(inout)      :: diag_data              ! diagnostic variables for a local HRU 
 type(var_dlength),intent(inout)      :: flux_data              ! model fluxes for a local HRU
 ! error control
 integer(i4b),intent(out)             :: err                    ! error code
 character(*),intent(out)             :: message                ! error message
 ! control the length of the sub-step
 real(dp)                             :: minstep                ! minimum time step (seconds)
 real(dp)                             :: maxstep                ! maximum time step (seconds)
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
 logical(lgt)                         :: resumeSubStepSolver    ! a flag to resume the simulation by using the last iteration solution when solver is not converged
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
 logical(lgt)                         :: computeVegFluxOld      ! flag to indicate if we are computing fluxes over vegetation on the previous sub step
 logical(lgt)                         :: modifiedLayers         ! flag to denote that snow layers were modified
 logical(lgt)                         :: modifiedVegState       ! flag to denote that vegetation states were modified
 type(var_dlength)                    :: flux_mean              ! timestep-average model fluxes for a local HRU
 integer(i4b)                         :: nLayersRoots           ! number of soil layers that contain roots
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
 logical(lgt)                         :: firstStep              ! flag to denote if the first time step
 logical(lgt)                         :: rejectedStep           ! flag to denote if the sub-step is rejected (convergence problem, etc.)
 logical(lgt),parameter               :: checkTimeStepping=.false.      ! flag to denote a desire to check the time stepping 
 logical(lgt),parameter               :: backwardsCompatibility=.true.  ! flag to denote a desire to ensure backwards compatibility with previous branches. 
 ! balance checks
 integer(i4b)                         :: iVar                   ! loop through model variables
 real(dp)                             :: totalSoilCompress      ! change in storage associated with compression of the soil matrix (kg m-2)
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
 ! initialize error control
 err=0; message="coupled_em/"

 ! This is the start of a data step for a local HRU

 ! define the first step
 firstStep = (istep==1)

 ! count the number of snow and soil layers
 ! NOTE: need to re-compute the number of snow and soil layers at the start of each sub-step because the number of layers may change
 !         (nSnow and nSoil are shared in the data structure)
 nSnow = count(indx_data%var(iLookINDEX%layerType)%dat==ix_snow)
 nSoil = count(indx_data%var(iLookINDEX%layerType)%dat==ix_soil)

 ! compute the total number of snow and soil layers
 nLayers = nSnow + nSoil

 ! allocate space for the local fluxes
 call allocLocal(averageFlux_meta(:)%var_info,flux_mean,nSnow,nSoil,err,cmessage)
 if(err/=0)then; err=20; message=trim(message)//trim(cmessage); return; end if

 ! initialize compression
 totalSoilCompress = 0._dp  ! change in storage associated with compression of the soil matrix (kg m-2)

 ! initialize mean fluxes
 do iVar=1,size(averageFlux_meta)
  flux_mean%var(iVar)%dat(:) = 0._dp
 end do

 ! associate local variables with information in the data structures
 associate(&
 ! state variables in the vegetation canopy
 scalarCanopyLiq      => prog_data%var(iLookPROG%scalarCanopyLiq)%dat(1)                 ,&  ! canopy liquid water (kg m-2)
 scalarCanopyIce      => prog_data%var(iLookPROG%scalarCanopyIce)%dat(1)                 ,&  ! canopy ice content (kg m-2)
 ! state variables in the soil domain
 mLayerDepth          => prog_data%var(iLookPROG%mLayerDepth)%dat(nSnow+1:nLayers)       ,&  ! depth of each soil layer (m)
 mLayerVolFracIce     => prog_data%var(iLookPROG%mLayerVolFracIce)%dat(nSnow+1:nLayers)  ,&  ! volumetric ice content in each soil layer (-)
 mLayerVolFracLiq     => prog_data%var(iLookPROG%mLayerVolFracLiq)%dat(nSnow+1:nLayers)  ,&  ! volumetric liquid water content in each soil layer (-)
 scalarAquiferStorage => prog_data%var(iLookPROG%scalarAquiferStorage)%dat(1)             &  ! aquifer storage (m)
 ) ! (association of local variables with information in the data structures

 ! compute total soil moisture and ice at the *START* of the step (kg m-2)
 scalarTotalSoilLiq = sum(iden_water*mLayerVolFracLiq(1:nSoil)*mLayerDepth(1:nSoil))
 scalarTotalSoilIce = sum(iden_water*mLayerVolFracIce(1:nSoil)*mLayerDepth(1:nSoil))  ! NOTE: no expansion and hence use iden_water

 ! compute storage of water in the canopy and the soil
 balanceCanopyWater0 = scalarCanopyLiq + scalarCanopyIce
 balanceSoilWater0   = scalarTotalSoilLiq + scalarTotalSoilIce

 ! get the total aquifer storage at the start of the time step (kg m-2)
 balanceAquifer0 = scalarAquiferStorage*iden_water

 ! end association of local variables with information in the data structures
 end associate

 ! short-cut to the algorithmic control parameters
 minstep = mpar_data%var(iLookPARAM%minstep)  ! minimum time step (s)
 maxstep = mpar_data%var(iLookPARAM%maxstep)  ! maximum time step (s)
 !print*, 'minstep, maxstep = ', minstep, maxstep

 ! define maximum number of iterations
 maxiter = nint(mpar_data%var(iLookPARAM%maxiter))

 ! get the length of the time step (seconds)
 dt = data_step

 ! compute the number of layers with roots
 nLayersRoots = count(prog_data%var(iLookPROG%iLayerHeight)%dat(nSnow:nLayers-1) < mpar_data%var(iLookPARAM%rootingDepth)-verySmall)
 if(nLayersRoots == 0)then; err=20; message=trim(message)//'no roots within the soil profile'; return; end if

 ! define the foliage nitrogen factor
 diag_data%var(iLookDIAG%scalarFoliageNitrogenFactor)%dat(1) = 1._dp  ! foliage nitrogen concentration (1.0 = saturated)

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
  oldSWE = prog_data%var(iLookPROG%scalarSWE)%dat(1)
  !print*, 'nSnow = ', nSnow
  !print*, 'oldSWE = ', oldSWE

  ! (1) compute phenology...
  ! ------------------------

  ! compute the temperature of the root zone: used in vegetation phenology
  diag_data%var(iLookDIAG%scalarRootZoneTemp)%dat(1) = sum(prog_data%var(iLookPROG%mLayerTemp)%dat(nSnow+1:nSnow+nLayersRoots)) / real(nLayersRoots, kind(dp))

  ! remember if we compute the vegetation flux on the previous sub-step
  computeVegFluxOld = computeVegFlux  

  ! compute the exposed LAI and SAI and whether veg is buried by snow
  call vegPhenlgy(&
                  ! input/output: data structures
                  model_decisions,             & ! intent(in):    model decisions
                  type_data,                   & ! intent(in):    type of vegetation and soil
                  attr_data,                   & ! intent(in):    spatial attributes
                  mpar_data,                   & ! intent(in):    model parameters
                  prog_data,                   & ! intent(in):    model prognostic variables for a local HRU
                  diag_data,                   & ! intent(inout): model diagnostic variables for a local HRU
                  ! output
                  computeVegFlux,              & ! intent(out): flag to indicate if we are computing fluxes over vegetation (.false. means veg is buried with snow)
                  canopyDepth,                 & ! intent(out): canopy depth (m)
                  exposedVAI,                  & ! intent(out): exposed vegetation area index (m2 m-2)
                  err,cmessage)                  ! intent(out): error control
  if(err/=0)then; err=20; message=trim(message)//trim(cmessage); return; end if

  ! check
  if(computeVegFlux)then
   if(canopyDepth < epsilon(canopyDepth))then
    message=trim(message)//'canopy depth is zero when computeVegFlux flag is .true.'
    err=20; return
   endif
  endif

  ! flag the case where number of vegetation states has changed
  modifiedVegState = (computeVegFlux.neqv.computeVegFluxOld)

  ! (2) compute wetted canopy area...
  ! ---------------------------------

  ! compute maximum canopy liquid water (kg m-2)
  diag_data%var(iLookDIAG%scalarCanopyLiqMax)%dat(1) = mpar_data%var(iLookPARAM%refInterceptCapRain)*exposedVAI

  ! compute maximum canopy ice content (kg m-2)
  ! NOTE 1: this is used to compute the snow fraction on the canopy, as used in *BOTH* the radiation AND canopy sublimation routines
  ! NOTE 2: this is a different variable than the max ice used in the throughfall (snow interception) calculations
  select case(model_decisions(iLookDECISIONS%snowIncept)%iDecision)
   case(lightSnow);  diag_data%var(iLookDIAG%scalarCanopyIceMax)%dat(1) = exposedVAI*mpar_data%var(iLookPARAM%refInterceptCapSnow)       ! use maximum per unit leaf area storage capacity for snow (kg m-2)
   case(stickySnow); diag_data%var(iLookDIAG%scalarCanopyIceMax)%dat(1) = exposedVAI*mpar_data%var(iLookPARAM%refInterceptCapSnow)*4._dp ! use maximum per unit leaf area storage capacity for snow (kg m-2)
   case default; message=trim(message)//'unable to identify option for maximum branch interception capacity'; err=20; return
  end select ! identifying option for maximum branch interception capacity
  !print*, 'diag_data%var(iLookDIAG%scalarCanopyLiqMax)%dat(1) = ', diag_data%var(iLookDIAG%scalarCanopyLiqMax)%dat(1)
  !print*, 'diag_data%var(iLookDIAG%scalarCanopyIceMax)%dat(1) = ', diag_data%var(iLookDIAG%scalarCanopyIceMax)%dat(1)

  ! compute wetted fraction of the canopy
  ! NOTE: assume that the wetted fraction is constant over the substep for the radiation calculations
  if(computeVegFlux)then

   ! compute wetted fraction of the canopy
   call wettedFrac(&
                   ! input
                   .false.,                                                      & ! flag to denote if derivatives are required
                   .false.,                                                      & ! flag to denote if derivatives are calculated numerically
                   (prog_data%var(iLookPROG%scalarCanopyTemp)%dat(1) < Tfreeze), & ! flag to denote if the canopy is frozen
                   varNotUsed1,                                                  & ! derivative in canopy liquid w.r.t. canopy temperature (kg m-2 K-1)
                   varNotUsed2,                                                  & ! fraction of liquid water on the canopy
                   prog_data%var(iLookPROG%scalarCanopyLiq)%dat(1),              & ! canopy liquid water (kg m-2)
                   prog_data%var(iLookPROG%scalarCanopyIce)%dat(1),              & ! canopy ice (kg m-2)
                   diag_data%var(iLookDIAG%scalarCanopyLiqMax)%dat(1),           & ! maximum canopy liquid water (kg m-2)
                   diag_data%var(iLookDIAG%scalarCanopyLiqMax)%dat(1),           & ! maximum canopy ice content (kg m-2)
                   mpar_data%var(iLookPARAM%canopyWettingFactor),                & ! maximum wetted fraction of the canopy (-)
                   mpar_data%var(iLookPARAM%canopyWettingExp),                   & ! exponent in canopy wetting function (-)
                   ! output
                   diag_data%var(iLookDIAG%scalarCanopyWetFraction)%dat(1),      & ! canopy wetted fraction (-)
                   dCanopyWetFraction_dWat,                                      & ! derivative in wetted fraction w.r.t. canopy liquid water content (kg-1 m2)
                   dCanopyWetFraction_dT,                                        & ! derivative in wetted fraction w.r.t. canopy liquid water content (kg-1 m2)
                   err,cmessage)
   if(err/=0)then; message=trim(message)//trim(cmessage); return; end if

  ! vegetation is completely buried by snow (or no veg exisits at all)
  else
   diag_data%var(iLookDIAG%scalarCanopyWetFraction)%dat(1) = 0._dp
   dCanopyWetFraction_dWat                                 = 0._dp
   dCanopyWetFraction_dT                                   = 0._dp
  end if

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
                  flux_data,                   & ! intent(in):    model flux variables
                  diag_data,                   & ! intent(inout): model diagnostic variables for a local HRU
                  prog_data,                   & ! intent(inout): model prognostic variables for a local HRU
                  ! output: error control
                  err,cmessage)                  ! intent(out): error control
  if(err/=0)then; err=20; message=trim(message)//trim(cmessage); return; end if


  ! (4) compute canopy sw radiation fluxes...
  ! -----------------------------------------
  call vegSWavRad(&
                  dt_sub,                       & ! intent(in):    time step (s) -- only used in Noah-MP radiation, to compute albedo
                  nSnow,                        & ! intent(in):    number of snow layers
                  nSoil,                        & ! intent(in):    number of soil layers
                  nLayers,                      & ! intent(in):    total number of layers
                  computeVegFlux,               & ! intent(in):    logical flag to compute vegetation fluxes (.false. if veg buried by snow)
                  type_data,                    & ! intent(in):    type of vegetation and soil
                  prog_data,                    & ! intent(inout): model prognostic variables for a local HRU
                  diag_data,                    & ! intent(inout): model diagnostic variables for a local HRU
                  flux_data,                    & ! intent(inout): model flux variables
                  err,cmessage)                   ! intent(out):   error control
  if(err/=0)then; err=20; message=trim(message)//trim(cmessage); return; end if


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
                  diag_data,                   & ! intent(in):    model diagnostic variables for a local HRU
                  prog_data,                   & ! intent(inout): model prognostic variables for a local HRU
                  flux_data,                   & ! intent(inout): model flux variables
                  ! output: error control
                  err,cmessage)                  ! intent(out): error control
  if(err/=0)then; err=20; message=trim(message)//trim(cmessage); return; end if
  !print*, 'canopyIce = ', prog_data%var(iLookPROG%scalarCanopyIce)%dat(1)

  ! adjust canopy temperature to account for new snow
  call tempAdjust(&
                  ! input: derived parameters
                  canopyDepth,                 & ! intent(in): canopy depth (m)
                  ! input/output: data structures
                  mpar_data,                   & ! intent(in):    model parameters
                  prog_data,                   & ! intent(inout): model prognostic variables for a local HRU
                  diag_data,                   & ! intent(out):   model diagnostic variables for a local HRU
                  ! output: error control
                  err,cmessage)                  ! intent(out): error control
  if(err/=0)then; err=20; message=trim(message)//trim(cmessage); return; end if

  ! initialize drainage and throughfall
  ! NOTE 1: this needs to be done before solving the energy and liquid water equations, to account for the heat advected with precipitation
  ! NOTE 2: this initialization needs to be done AFTER the call to canopySnow, since canopySnow uses canopy drip drom the previous time step
  if(.not.computeVegFlux)then
   flux_data%var(iLookFLUX%scalarThroughfallRain)%dat(1)   = flux_data%var(iLookFLUX%scalarRainfall)%dat(1)
   flux_data%var(iLookFLUX%scalarCanopyLiqDrainage)%dat(1) = 0._dp
  else
   flux_data%var(iLookFLUX%scalarThroughfallRain)%dat(1)   = 0._dp
   flux_data%var(iLookFLUX%scalarCanopyLiqDrainage)%dat(1) = 0._dp
  end if

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
                 diag_data%var(iLookDIAG%scalarSnowfallTemp)%dat(1),        & ! computed temperature of fresh snow (K)
                 diag_data%var(iLookDIAG%scalarNewSnowDensity)%dat(1),      & ! computed density of new snow (kg m-3)
                 flux_data%var(iLookFLUX%scalarThroughfallSnow)%dat(1),     & ! throughfall of snow through the canopy (kg m-2 s-1)
                 flux_data%var(iLookFLUX%scalarCanopySnowUnloading)%dat(1), & ! unloading of snow from the canopy (kg m-2 s-1)
                 ! input/output: state variables
                 prog_data%var(iLookPROG%scalarSWE)%dat(1),                 & ! SWE (kg m-2)
                 prog_data%var(iLookPROG%scalarSnowDepth)%dat(1),           & ! total snow depth (m)
                 prog_data%var(iLookPROG%mLayerTemp)%dat(1),                & ! temperature of the top layer (K)
                 prog_data%var(iLookPROG%mLayerDepth)%dat(1),               & ! depth of the top layer (m)
                 prog_data%var(iLookPROG%mLayerVolFracIce)%dat(1),          & ! volumetric fraction of ice of the top layer (-)
                 prog_data%var(iLookPROG%mLayerVolFracLiq)%dat(1),          & ! volumetric fraction of liquid water of the top layer (-)
                 ! output: error control
                 err,cmessage)                                                ! error control
  if(err/=0)then; err=30; message=trim(message)//trim(cmessage); return; end if

  ! re-compute snow depth and SWE
  if(nSnow > 0)then
   prog_data%var(iLookPROG%scalarSnowDepth)%dat(1) = sum(  prog_data%var(iLookPROG%mLayerDepth)%dat(1:nSnow))
   prog_data%var(iLookPROG%scalarSWE)%dat(1)       = sum( (prog_data%var(iLookPROG%mLayerVolFracLiq)%dat(1:nSnow)*iden_water + &
                                                           prog_data%var(iLookPROG%mLayerVolFracIce)%dat(1:nSnow)*iden_ice) &
                                                         * prog_data%var(iLookPROG%mLayerDepth)%dat(1:nSnow) )
  end if
  !print*, 'SWE after snowfall = ',  prog_data%var(iLookPROG%scalarSWE)%dat(1)

  ! update coordinate variables
  call calcHeight(&
                  ! input/output: data structures
                  indx_data,   & ! intent(in): layer type
                  prog_data,   & ! intent(inout): model variables for a local HRU
                  ! output: error control
                  err,cmessage)
  if(err/=0)then; err=20; message=trim(message)//trim(cmessage); return; end if

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
  resumeSubStepSolver=.false. 
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
                    prog_data,                   & ! intent(inout): model prognostic variables for a local HRU
                    diag_data,                   & ! intent(inout): model diagnostic variables for a local HRU
                    flux_data,                   & ! intent(inout): model fluxes for a local HRU
                    ! output
                    modifiedLayers,              & ! intent(out): flag to denote that layers were modified
                    err,cmessage)                  ! intent(out): error control
    if(err/=0)then; err=55; message=trim(message)//trim(cmessage); return; end if

    ! recompute the number of snow and soil layers
    ! NOTE: do this here for greater visibility
    nSnow   = count(indx_data%var(iLookINDEX%layerType)%dat==ix_snow)
    nSoil   = count(indx_data%var(iLookINDEX%layerType)%dat==ix_soil)
    nLayers = nSnow+nSoil

    ! put the data in the structures
    indx_data%var(iLookINDEX%nSnow)%dat(1)   = nSnow
    indx_data%var(iLookINDEX%nSoil)%dat(1)   = nSoil
    indx_data%var(iLookINDEX%nLayers)%dat(1) = nLayers

    ! compute the indices for the model state variables
    if(firstStep .or. modifiedVegState .or. modifiedLayers)then
     call indexState(computeVegFlux,          & ! intent(in):    flag to denote if computing the vegetation flux
                     indx_data,               & ! intent(inout): indices defining model states and layers
                     err,cmessage)              ! intent(out):   error control
     if(err/=0)then; message=trim(message)//trim(cmessage); return; end if
    end if

    ! define the number of state variables
    nState = indx_data%var(iLookINDEX%nState)%dat(1)

    ! re-compute snow depth and SWE
    if(nSnow > 0)then
     prog_data%var(iLookPROG%scalarSnowDepth)%dat(1) = sum(  prog_data%var(iLookPROG%mLayerDepth)%dat(1:nSnow))
     prog_data%var(iLookPROG%scalarSWE)%dat(1)       = sum( (prog_data%var(iLookPROG%mLayerVolFracLiq)%dat(1:nSnow)*iden_water + &
                                                             prog_data%var(iLookPROG%mLayerVolFracIce)%dat(1:nSnow)*iden_ice) &
                                                           * prog_data%var(iLookPROG%mLayerDepth)%dat(1:nSnow) )
    end if


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
                    prog_data,               & ! intent(in):    model prognostic variables for a local HRU
                    diag_data,               & ! intent(inout): model diagnostic variables for a local HRU
                    ! output: error control
                    err,cmessage)              ! intent(out): error control
    if(err/=0)then; err=55; message=trim(message)//trim(cmessage); return; end if


    ! (8) compute melt of the "snow without a layer"...
    ! -------------------------------------------------
    ! NOTE: forms a surface melt pond, which drains into the upper-most soil layer through the time step
    ! (check for the special case of "snow without a layer")
    if(nSnow==0)then
     call implctMelt(&
                     ! input/output: integrated snowpack properties
                     prog_data%var(iLookPROG%scalarSWE)%dat(1),               & ! intent(inout): snow water equivalent (kg m-2)
                     prog_data%var(iLookPROG%scalarSnowDepth)%dat(1),         & ! intent(inout): snow depth (m)
                     prog_data%var(iLookPROG%scalarSfcMeltPond)%dat(1),       & ! intent(inout): surface melt pond (kg m-2)
                     ! input/output: properties of the upper-most soil layer
                     prog_data%var(iLookPROG%mLayerTemp)%dat(nSnow+1),        & ! intent(inout): surface layer temperature (K)
                     prog_data%var(iLookPROG%mLayerDepth)%dat(nSnow+1),       & ! intent(inout): surface layer depth (m)
                     diag_data%var(iLookDIAG%mLayerVolHtCapBulk)%dat(nSnow+1),& ! intent(inout): surface layer volumetric heat capacity (J m-3 K-1)
                     ! output: error control
                     err,cmessage                                             ) ! intent(out): error control
     if(err/=0)then; err=20; message=trim(message)//trim(cmessage); return; end if
    end if

   ! ** if previous step is not rejected
   end if

   ! test: recompute snow depth and SWE
   if(nSnow > 0)then
    prog_data%var(iLookPROG%scalarSnowDepth)%dat(1) = sum(  prog_data%var(iLookPROG%mLayerDepth)%dat(1:nSnow))
    prog_data%var(iLookPROG%scalarSWE)%dat(1)       = sum( (prog_data%var(iLookPROG%mLayerVolFracLiq)%dat(1:nSnow)*iden_water + &
                                                            prog_data%var(iLookPROG%mLayerVolFracIce)%dat(1:nSnow)*iden_ice) &
                                                          * prog_data%var(iLookPROG%mLayerDepth)%dat(1:nSnow) )
   end if
   !write(*,'(a,1x,2(f20.5,1x),l1)') 'b4 systemSolv: testSWE, meltPond, rejectedStep = ', prog_data%var(iLookPROG%scalarSWE)%dat(1), prog_data%var(iLookPROG%scalarSfcMeltPond)%dat(1), rejectedStep

   ! print progress
   !if(dt_temp < dt_sub-tinyNumber)then
   ! write(*,'(a,1x,i4,1x,3(f15.2,1x))') 'ntrial, dt_temp, dt_solv, dt_sub = ', ntrial, dt_temp, dt_solv, dt_sub
   ! !pause
   !end if

   ! (9) solve model equations...
   ! ----------------------------

   ! get the new solution
   call systemSolv(&
                   ! input: model control
                   nSnow,                                  & ! intent(in): number of snow layers
                   nSoil,                                  & ! intent(in): number of soil layers
                   nLayers,                                & ! intent(in): total number of layers
                   nState,                                 & ! intent(in): total number of layers
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
                   prog_data,                              & ! intent(inout): model prognostic variables for a local HRU
                   diag_data,                              & ! intent(inout): model diagnostic variables for a local HRU
                   flux_data,                              & ! intent(inout): model fluxes for a local HRU
                   bvar_data,                              & ! intent(in):    model variables for the local basin
                   model_decisions,                        & ! intent(in):    model decisions
                   ! output: model control
                   niter,                                  & ! intent(out): number of iterations
                   resumeSubStepSolver,                    & ! intent(in):  resume the solver even it is not 
                   err,cmessage)                             ! intent(out): error code and error message

   ! check for fatal errors
   if(err>0)then; err=20; message=trim(message)//trim(cmessage); return; end if

   ! update first step
   firstStep=.false.

   ! if err<0 (warnings) and hence non-convergence
   if(err<0)then
    ! (adjust time step length)
    dt_temp = dt_temp*0.5_dp ! halve the sub-step
    write(*,'(a,1x,2(f13.3,1x),A,I0)') trim(cmessage), dt_temp, minstep,' at HRU ',hruId
    rejectedStep=.true.
    ! (check that time step greater than the minimum step)
    if(dt_temp < minstep)then
     if (resumeFailSolver) then 
      resumeSubStepSolver = .true.
      write(*,'(a,i0,a)') 'Solver fails convergence with minimum time step at HRU ',hruId,' but proceeds intentionally.'
     else
      message=trim(message)//'dt_temp is below the minimum time step'
      err=20; return
     end if 
    endif
    !pause 'failed step'
    ! (try again)
    cycle  ! try again
   else
    rejectedStep=.false.
    !pause 'accepted step'
   end if

   ! check that err=0 at this point (it should be)
   if(err/=0)then; message=trim(message)//'expect err=0'; return; end if


   ! (10a) compute change in canopy ice content due to sublimation...
   ! --------------------------------------------------------------
   ! NOTE: keep in continuous do loop in case insufficient water on canopy for sublimation
   if(computeVegFlux)then

    ! remove mass of ice on the canopy
    prog_data%var(iLookPROG%scalarCanopyIce)%dat(1) = prog_data%var(iLookPROG%scalarCanopyIce)%dat(1) + &
                                                      flux_data%var(iLookFLUX%scalarCanopySublimation)%dat(1)*dt_temp

    ! if removed all ice, take the remaining sublimation from water
    if(prog_data%var(iLookPROG%scalarCanopyIce)%dat(1) < 0._dp)then
     prog_data%var(iLookPROG%scalarCanopyLiq)%dat(1) = prog_data%var(iLookPROG%scalarCanopyLiq)%dat(1) + prog_data%var(iLookPROG%scalarCanopyIce)%dat(1)
     prog_data%var(iLookPROG%scalarCanopyIce)%dat(1) = 0._dp
    end if

    ! check that there is sufficient canopy water to support the converged sublimation rate over the time step dt_temp
    ! NOTE we conducted checks and time step adjustments in systemSolv above so we should not get here: hence fatal error
    if(prog_data%var(iLookPROG%scalarCanopyLiq)%dat(1) < -tinyNumber)then
     message=trim(message)//'canopy sublimation rate over time step dt_temp depletes more than the available water'
     err=20; return
    end if

   end if  ! (if computing the vegetation flux)


   ! (10b) compute change in ice content of the top snow layer due to sublimation...
   ! -----------------------------------------------------------------------------
   ! NOTE: this is done BEFORE densification
   if(nSnow > 0)then ! snow layers exist

    ! compute volumetric sublimation (-)
    volSub = dt_temp*flux_data%var(iLookFLUX%scalarSnowSublimation)%dat(1)/prog_data%var(iLookPROG%mLayerDepth)%dat(1)

    ! update volumetric fraction of ice (-)
    ! NOTE: fluxes are positive downward
    prog_data%var(iLookPROG%mLayerVolFracIce)%dat(1) = prog_data%var(iLookPROG%mLayerVolFracIce)%dat(1) + volSub/iden_ice

    ! check that there is sufficient ice in the top snow layer to support the converged sublimation rate over the time step dt_temp
    ! NOTE we conducted checks and time step adjustments in systemSolv above so we should not get here: hence fatal error
    if(prog_data%var(iLookPROG%mLayerVolFracIce)%dat(1) < -tinyNumber)then
     message=trim(message)//'surface sublimation rate over time step dt_temp depletes more than the available water'
     err=20; return
    end if

   ! no snow
   else

    ! no snow: check that sublimation is zero
    if(abs(flux_data%var(iLookFLUX%scalarSnowSublimation)%dat(1)) > verySmall)then
     message=trim(message)//'sublimation of snow has been computed when no snow exists'
     err=20; return
    end if

   end if  ! (if snow layers exist)
   !print*, 'ice after sublimation: ', prog_data%var(iLookPROG%mLayerVolFracIce)%dat(1)*iden_ice


   ! (11) account for compaction and cavitation in the snowpack...
   ! ------------------------------------------------------------
   if(nSnow>0)then
    call snwDensify(&
                    ! intent(in): variables
                    dt_temp,                                                & ! intent(in): time step (s)
                    indx_data%var(iLookINDEX%nSnow)%dat(1),                 & ! intent(in): number of snow layers
                    prog_data%var(iLookPROG%mLayerTemp)%dat(1:nSnow),       & ! intent(in): temperature of each layer (K)
                    diag_data%var(iLookDIAG%mLayerMeltFreeze)%dat(1:nSnow), & ! intent(in): volumetric melt in each layer (kg m-3)
                    flux_data%var(iLookFLUX%scalarSnowSublimation)%dat(1),  & ! intent(in): sublimation from the snow surface (kg m-2 s-1)
                    ! intent(in): parameters
                    mpar_data%var(iLookPARAM%densScalGrowth),               & ! intent(in): density scaling factor for grain growth (kg-1 m3)
                    mpar_data%var(iLookPARAM%tempScalGrowth),               & ! intent(in): temperature scaling factor for grain growth (K-1)
                    mpar_data%var(iLookPARAM%grainGrowthRate),              & ! intent(in): rate of grain growth (s-1)
                    mpar_data%var(iLookPARAM%densScalOvrbdn),               & ! intent(in): density scaling factor for overburden pressure (kg-1 m3)
                    mpar_data%var(iLookPARAM%tempScalOvrbdn),               & ! intent(in): temperature scaling factor for overburden pressure (K-1)
                    mpar_data%var(iLookPARAM%baseViscosity),                 & ! intent(in): viscosity coefficient at T=T_frz and snow density=0 (kg m-2 s)
                    ! intent(inout): state variables
                    prog_data%var(iLookPROG%mLayerDepth)%dat(1:nSnow),      & ! intent(inout): depth of each layer (m)
                    prog_data%var(iLookPROG%mLayerVolFracLiq)%dat(1:nSnow), & ! intent(inout):  volumetric fraction of liquid water after itertations (-)
                    prog_data%var(iLookPROG%mLayerVolFracIce)%dat(1:nSnow), & ! intent(inout):  volumetric fraction of ice after itertations (-)
                    ! output: error control
                    err,cmessage)                     ! intent(out): error control
    if(err/=0)then; err=55; message=trim(message)//trim(cmessage); return; end if
   end if  ! if snow layers exist

   ! update coordinate variables
   call calcHeight(&
                   ! input/output: data structures
                   indx_data,   & ! intent(in): layer type
                   prog_data,   & ! intent(inout): model variables for a local HRU
                   ! output: error control
                   err,cmessage)
   if(err/=0)then; err=20; message=trim(message)//trim(cmessage); return; end if

   ! recompute snow depth and SWE
   if(nSnow > 0)then
    prog_data%var(iLookPROG%scalarSnowDepth)%dat(1) = sum(  prog_data%var(iLookPROG%mLayerDepth)%dat(1:nSnow))
    prog_data%var(iLookPROG%scalarSWE)%dat(1)       = sum( (prog_data%var(iLookPROG%mLayerVolFracLiq)%dat(1:nSnow)*iden_water + &
                                                            prog_data%var(iLookPROG%mLayerVolFracIce)%dat(1:nSnow)*iden_ice) &
                                                          * prog_data%var(iLookPROG%mLayerDepth)%dat(1:nSnow) )
   end if


   ! (12) compute sub-step averages associated with the temporary steps...
   ! ---------------------------------------------------------------------

   ! keep track of the number of temporary sub-steps
   ntemp = ntemp+1

   ! increment model fluxes
   dt_wght = dt_temp/dt ! define weight applied to each sub-step

   ! increment fluxes
   do iVar=1,size(averageFlux_meta)
    flux_mean%var(iVar)%dat(:) = flux_mean%var(iVar)%dat(:) + flux_data%var(averageFlux_meta(iVar)%ixParent)%dat(:)*dt_wght 
   end do

   ! increment change in storage associated with compression of the soil matrix (kg m-2)
   totalSoilCompress = totalSoilCompress + diag_data%var(iLookDIAG%scalarSoilCompress)%dat(1)

   ! check time stepping
   if(checkTimestepping)then
    dt_prog = dt_done+dt_solv+dt_temp  ! progress in time step (s)
    dt_frac = dt_prog/dt               ! fraction of time step completed (-)
    write(*,'(a,1x,3(f9.3,1x),10(e20.10,1x))') 'dt_wght, dt_prog, dt_frac = ', &
                                                dt_wght, dt_prog, dt_frac
   end if

   ! compute effective rainfall input and snowpack drainage to/from the snowpack (kg m-2 s-1)
   if(nSnow > 0)then
    effRainfall = effRainfall + (flux_data%var(iLookFLUX%scalarThroughfallRain)%dat(1) + flux_data%var(iLookFLUX%scalarCanopyLiqDrainage)%dat(1) )*dt_temp/dt_sub
    snwDrainage = snwDrainage + (flux_data%var(iLookFLUX%iLayerLiqFluxSnow)%dat(nSnow)*iden_water )*dt_temp/dt_sub ! m s-1 -> kg m-2 s-1
    sublimation = sublimation + (flux_data%var(iLookFLUX%scalarSnowSublimation)%dat(1))*dt_temp/dt_sub
   ! compute the surface melt pond (kg m-2)
   else
    sfcMeltPond = sfcMeltPond + prog_data%var(iLookPROG%scalarSfcMeltPond)%dat(1)
   end if

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
   prog_data%var(iLookPROG%scalarSnowDepth)%dat(1) = sum(  prog_data%var(iLookPROG%mLayerDepth)%dat(1:nSnow))
   prog_data%var(iLookPROG%scalarSWE)%dat(1)       = sum( (prog_data%var(iLookPROG%mLayerVolFracLiq)%dat(1:nSnow)*iden_water + &
                                                           prog_data%var(iLookPROG%mLayerVolFracIce)%dat(1:nSnow)*iden_ice) &
                                                         * prog_data%var(iLookPROG%mLayerDepth)%dat(1:nSnow) )
  end if

  ! check SWE
  effSnowfall = flux_data%var(iLookFLUX%scalarThroughfallSnow)%dat(1) + flux_data%var(iLookFLUX%scalarCanopySnowUnloading)%dat(1)
  newSWE      = prog_data%var(iLookPROG%scalarSWE)%dat(1)
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
  end if

  ! (14) adjust length of the substep...
  ! ------------------------------------

  ! increment the time step
  dt_done = dt_done + dt_sub
  !print*, '***** ', dt_done, dt_sub, niter
  !pause ' after increment the time step'

  ! modify the length of the time step
  if(niter<n_inc) dt_sub = min(dt_temp*F_inc,maxstep)
  if(niter>n_dec) dt_sub =     dt_temp*F_dec
  if(dt_sub < minstep)then; message=trim(message)//'dt_sub is below the minimum time step'; return; end if

  ! save the time step to initialize the subsequent step
  if(dt_done<dt .or. nsub==1) dt_init = dt_sub
  if(dt_init < 0.00001_dp .and. nsub > 10000) then
   write(message,'(a,f13.10,a,f9.2,a,i0,a)')trim(message)//"dt < 0.00001 and nsub > 10000 [dt=",dt_init,"; dt_done=",&
         dt_done,"; nsub=",nsub,"]"
   err=20; return
  end if

  ! exit do-loop if finished
  if(dt_done>=dt)exit

  ! make sure that we don't exceed the step
  dt_sub = min(dt-dt_done, dt_sub)
  !print*, 'in substep loop: dt_sub = ', dt_sub

 end do  ! (sub-step loop)
 !pause 'completed time step'

 ! ---
 ! (14) balance checks...
 ! ----------------------

 ! associate local variables with information in the data structures
 associate(&
 ! model forcing
 scalarSnowfall             => flux_mean%var(childFLUX_MEAN(iLookFLUX%scalarSnowfall)           )%dat(1)     ,&  ! computed snowfall rate (kg m-2 s-1)
 scalarRainfall             => flux_mean%var(childFLUX_MEAN(iLookFLUX%scalarRainfall)           )%dat(1)     ,&  ! computed rainfall rate (kg m-2 s-1)
 ! canopy fluxes
 averageThroughfallSnow     => flux_mean%var(childFLUX_MEAN(iLookFLUX%scalarThroughfallSnow)    )%dat(1)     ,&  ! snow that reaches the ground without ever touching the canopy (kg m-2 s-1)
 averageThroughfallRain     => flux_mean%var(childFLUX_MEAN(iLookFLUX%scalarThroughfallRain)    )%dat(1)     ,&  ! rain that reaches the ground without ever touching the canopy (kg m-2 s-1)
 averageCanopySnowUnloading => flux_mean%var(childFLUX_MEAN(iLookFLUX%scalarCanopySnowUnloading))%dat(1)     ,&  ! unloading of snow from the vegetion canopy (kg m-2 s-1)
 averageCanopyLiqDrainage   => flux_mean%var(childFLUX_MEAN(iLookFLUX%scalarCanopyLiqDrainage)  )%dat(1)     ,&  ! drainage of liquid water from the vegetation canopy (kg m-2 s-1)
 averageCanopySublimation   => flux_mean%var(childFLUX_MEAN(iLookFLUX%scalarCanopySublimation)  )%dat(1)     ,&  ! canopy sublimation/frost (kg m-2 s-1)
 averageCanopyEvaporation   => flux_mean%var(childFLUX_MEAN(iLookFLUX%scalarCanopyEvaporation)  )%dat(1)     ,&  ! canopy evaporation/condensation (kg m-2 s-1)
 ! soil fluxes
 averageSoilInflux          => flux_mean%var(childFLUX_MEAN(iLookFLUX%iLayerLiqFluxSoil)        )%dat(0)     ,&  ! influx of water at the top of the soil profile (m s-1)
 averageSoilDrainage        => flux_mean%var(childFLUX_MEAN(iLookFLUX%iLayerLiqFluxSoil)        )%dat(nSoil) ,&  ! drainage from the bottom of the soil profile (m s-1)
 averageSoilBaseflow        => flux_mean%var(childFLUX_MEAN(iLookFLUX%scalarSoilBaseflow)       )%dat(1)     ,&  ! total baseflow from throughout the soil profile (m s-1)
 averageCanopyTranspiration => flux_mean%var(childFLUX_MEAN(iLookFLUX%scalarCanopyTranspiration))%dat(1)     ,&  ! canopy transpiration (kg m-2 s-1)
 ! state variables in the vegetation canopy
 scalarCanopyLiq            => prog_data%var(iLookPROG%scalarCanopyLiq)%dat(1)                               ,&  ! canopy liquid water (kg m-2)
 scalarCanopyIce            => prog_data%var(iLookPROG%scalarCanopyIce)%dat(1)                               ,&  ! canopy ice content (kg m-2)
 ! state variables in the soil domain
 mLayerDepth                => prog_data%var(iLookPROG%mLayerDepth)%dat(nSnow+1:nLayers)                     ,&  ! depth of each soil layer (m)
 mLayerVolFracIce           => prog_data%var(iLookPROG%mLayerVolFracIce)%dat(nSnow+1:nLayers)                ,&  ! volumetric ice content in each soil layer (-)
 mLayerVolFracLiq           => prog_data%var(iLookPROG%mLayerVolFracLiq)%dat(nSnow+1:nLayers)                ,&  ! volumetric liquid water content in each soil layer (-)
 scalarAquiferStorage       => prog_data%var(iLookPROG%scalarAquiferStorage)%dat(1)                           &  ! aquifer storage (m)
 ) ! (association of local variables with information in the data structures

 ! canopy water balance
 balanceCanopyWater1 = scalarCanopyLiq + scalarCanopyIce

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
 end if

 ! compute the liquid water and ice content at the end of the time step
 scalarTotalSoilLiq = sum(iden_water*mLayerVolFracLiq(1:nSoil)*mLayerDepth(1:nSoil))
 scalarTotalSoilIce = sum(iden_water*mLayerVolFracIce(1:nSoil)*mLayerDepth(1:nSoil))   ! NOTE: no expansion of soil, hence use iden_water

 ! get the total water in the soil (liquid plus ice) at the end of the time step (kg m-2)
 balanceSoilWater1 = scalarTotalSoilLiq + scalarTotalSoilIce

 ! get the total aquifer storage at the start of the time step (kg m-2)
 balanceAquifer1 = scalarAquiferStorage*iden_water

 ! get the input and output to/from the soil zone (kg m-2)
 balanceSoilInflux        = averageSoilInflux*iden_water*dt
 balanceSoilBaseflow      = averageSoilBaseflow*iden_water*dt
 balanceSoilDrainage      = averageSoilDrainage*iden_water*dt
 balanceSoilTranspiration = averageCanopyTranspiration*dt      ! NOTE ground evaporation included in the flux at the upper boundary

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
   xCompress = diag_data%var(iLookDIAG%mLayerCompress)%dat(iLayer)
   xFlux0    = flux_mean%var(iLookFLUX%iLayerLiqFluxSoil)%dat(iLayer-1)*dt
   xFlux1    = flux_mean%var(iLookFLUX%iLayerLiqFluxSoil)%dat(iLayer)*dt
   write(*,'(a,1x,i4,1x,10(e20.10,1x))') 'iLayer, xFlux0, xFlux1, (xFlux1 - xFlux0)/mLayerDepth(iLayer), xCompress = ', &
                                          iLayer, xFlux0, xFlux1, (xFlux1 - xFlux0)/mLayerDepth(iLayer), xCompress
  end do
  ! error control
  message=trim(message)//'soil hydrology does not balance'
  err=20; return
 end if

 ! end association of local variables with information in the data structures
 end associate

 ! save the surface temperature (just to make things easier to visualize)
 prog_data%var(iLookPROG%scalarSurfaceTemp)%dat(1) = prog_data%var(iLookPROG%mLayerTemp)%dat(1)

 ! overwrite flux data with the timestep-average value
 if(.not.backwardsCompatibility)then
  do iVar=1,size(flux_mean%var)
   flux_data%var(averageFlux_meta(iVar)%ixParent)%dat = flux_mean%var(iVar)%dat
  end do
 end if

 iLayer = nSnow+1
 !print*, 'nsub, mLayerTemp(iLayer), mLayerVolFracIce(iLayer) = ', nsub, mLayerTemp(iLayer), mLayerVolFracIce(iLayer)
 !print*, 'nsub = ', nsub
 if(nsub>2000)then
  message=trim(message)//'number of sub-steps > 2000'
  err=20; return
 end if

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
   end if
   ! update depth
   scalarSnowDepth = scalarSWE/snwDensity
   ! update temperature of the top soil layer (K)
   soilTemp =  soilTemp - (LH_fus*scalarSfcMeltPond/soilDepth)/soilHeatcap
  else  ! melt is zero if the temperature of the top soil layer is less than Tfreeze
   scalarSfcMeltPond = 0._dp  ! kg m-2
  end if ! (if the temperature of the top soil layer is greater than Tfreeze)
 else  ! melt is zero if the "snow without a layer" does not exist
  scalarSfcMeltPond = 0._dp  ! kg m-2
 end if ! (if the "snow without a layer" exists)

 end subroutine implctMelt

end module coupled_em_module

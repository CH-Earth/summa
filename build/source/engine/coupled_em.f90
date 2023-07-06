! SUMMA - Structure for Unifying Multiple Modeling Alternatives
! Copyright (C) 2014-2020 NCAR/RAL; University of Saskatchewan; University of Washington
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
                    LH_sub,       & ! latent heat of sublimation           (J kg-1)
                    iden_ice,     & ! intrinsic density of ice             (kg m-3)
                    iden_water      ! intrinsic density of liquid water    (kg m-3)

! data types
USE data_types,only:&
                    var_i,               & ! x%var(:)                (i4b)
                    var_d,               & ! x%var(:)                (rkind)
                    var_ilength,         & ! x%var(:)%dat            (i4b)
                    var_dlength,         & ! x%var(:)%dat            (rkind)
                    zLookup                ! x%z(:)%var(:)%lookup(:) (rkind)

! named variables for parent structures
USE var_lookup,only:iLookDECISIONS         ! named variables for elements of the decision structure
USE var_lookup,only:iLookPROG              ! named variables for structure elements
USE var_lookup,only:iLookDIAG              ! named variables for structure elements
USE var_lookup,only:iLookFLUX              ! named variables for structure elements
USE var_lookup,only:iLookPARAM             ! named variables for structure elements
USE var_lookup,only:iLookINDEX             ! named variables for structure elements
USE globalData,only:iname_snow             ! named variables for snow
USE globalData,only:iname_soil             ! named variables for soil

! named variables for child structures
USE var_lookup,only:childFLUX_MEAN

! metadata
USE globalData,only:flux_meta              ! metadata on the model fluxes
USE globalData,only:indx_meta              ! metadata on the model index variables
USE globalData,only:diag_meta              ! metadata on the model diagnostic variables
USE globalData,only:prog_meta              ! metadata on the model prognostic variables
USE globalData,only:averageFlux_meta       ! metadata on the timestep-average model flux structure

! global data
USE globalData,only:data_step              ! time step of forcing data (s)
USE globalData,only:model_decisions        ! model decision structure
USE globalData,only:globalPrintFlag        ! the global print flag

! look-up values for the maximum interception capacity
USE mDecisions_module,only:         &
                      stickySnow,   &      ! maximum interception capacity an increasing function of temerature
                      lightSnow            ! maximum interception capacity an inverse function of new snow density

! look-up values for the groundwater parameterization
USE mDecisions_module,only:         &
                      qbaseTopmodel,&      ! TOPMODEL-ish baseflow parameterization
                      bigBucket    ,&      ! a big bucket (lumped aquifer model)
                      noExplicit           ! no explicit groundwater parameterization

! look-up values for the spatial representation of groundwater
USE mDecisions_module,only:         &
                      localColumn  ,&      ! separate groundwater representation in each local soil column
                      singleBasin          ! single groundwater store over the entire basin

! look-up values for the numerical method
USE mDecisions_module,only:         &
                      numrec       ,&      ! home-grown backward Euler solution using free versions of Numerical recipes
                      kinsol       ,&      ! SUNDIALS backward Euler solution using Kinsol
                      ida                  ! SUNDIALS solution using IDA

! privacy
implicit none
private
public::coupled_em
! algorithmic parameters
real(rkind),parameter     :: valueMissing=-9999._rkind  ! missing value, used when diagnostic or state variables are undefined
real(rkind),parameter     :: verySmall=1.e-6_rkind   ! used as an additive constant to check if substantial difference among real numbers
real(rkind),parameter     :: mpe=1.e-6_rkind         ! prevents overflow error if division by zero
real(rkind),parameter     :: dx=1.e-6_rkind          ! finite difference increment
contains


! ************************************************************************************************
! public subroutine coupled_em: run the coupled energy-mass model for one timestep
! ************************************************************************************************
subroutine coupled_em(&
                      ! model control
                      hruId,             & ! intent(in):    hruId
                      dt_init,           & ! intent(inout): used to initialize the size of the sub-step
                      dt_init_factor,    & ! intent(in):    Used to adjust the length of the timestep in the event of a failure
                      computeVegFlux,    & ! intent(inout): flag to indicate if we are computing fluxes over vegetation (.false. means veg is buried with snow)
                      fracJulDay,        & ! intent(in):    fractional julian days since the start of year
                      yearLength,        & ! intent(in):    number of days in the current year
                      ! data structures (input)
                      type_data,         & ! intent(in):    local classification of soil veg etc. for each HRU
                      attr_data,         & ! intent(in):    local attributes for each HRU
                      forc_data,         & ! intent(in):    model forcing data
                      mpar_data,         & ! intent(in):    model parameters
                      bvar_data,         & ! intent(in):    basin-average variables
                      lookup_data,       & ! intent(in):    lookup tables
                      ! data structures (input-output)
                      indx_data,         & ! intent(inout): model indices
                      prog_data,         & ! intent(inout): prognostic variables for a local HRU
                      diag_data,         & ! intent(inout): diagnostic variables for a local HRU
                      flux_data,         & ! intent(inout): model fluxes for a local HRU
                      ! error control
                      err,message)         ! intent(out):   error control
  ! structure allocations
  USE allocspace_module,only:allocLocal      ! allocate local data structures
  USE allocspace_module,only:resizeData      ! clone a data structure
  ! simulation of fluxes and residuals given a trial state vector
  USE soil_utils_module,only:liquidHead                ! compute the liquid water matric potential
  ! preliminary subroutines
  USE vegPhenlgy_module,only:vegPhenlgy      ! compute vegetation phenology
  USE vegNrgFlux_module,only:wettedFrac      ! compute wetted fraction of the canopy (used in sw radiation fluxes)
  USE snowAlbedo_module,only:snowAlbedo      ! compute snow albedo
  USE vegSWavRad_module,only:vegSWavRad      ! compute canopy sw radiation fluxes
  USE canopySnow_module,only:canopySnow      ! compute interception and unloading of snow from the vegetation canopy
  USE volicePack_module,only:newsnwfall      ! compute change in the top snow layer due to throughfall and unloading
  USE volicePack_module,only:volicePack      ! merge and sub-divide snow layers, if necessary
  USE diagn_evar_module,only:diagn_evar      ! compute diagnostic energy variables -- thermal conductivity and heat capacity
  ! the model solver
  USE indexState_module,only:indexState      ! define indices for all model state variables and layers
  USE opSplittin_module,only:opSplittin      ! solve the system of thermodynamic and hydrology equations for a given substep
  USE time_utils_module,only:elapsedSec      ! calculate the elapsed time
  ! additional subroutines
  USE tempAdjust_module,only:tempAdjust      ! adjust snow temperature associated with new snowfall
  USE var_derive_module,only:calcHeight      ! module to calculate height at layer interfaces and layer mid-point
  USE computSnowDepth_module,only:computSnowDepth

  implicit none
  ! model control
#ifdef ACTORS_ACTIVE
  integer(4),intent(in)                :: hruId                  ! hruId
#else
  integer(8),intent(in)                :: hruId                  ! hruId
#endif
  real(rkind),intent(inout)            :: dt_init                ! used to initialize the size of the sub-step
  integer(i4b),intent(in)              :: dt_init_factor         ! Used to adjust the length of the timestep in the event of a failure
  logical(lgt),intent(inout)           :: computeVegFlux         ! flag to indicate if we are computing fluxes over vegetation (.false. means veg is buried with snow)
  ! data structures (input)
  type(var_i),intent(in)               :: type_data              ! type of vegetation and soil
  type(var_d),intent(in)               :: attr_data              ! spatial attributes
  type(var_d),intent(in)               :: forc_data              ! model forcing data
  type(var_dlength),intent(in)         :: mpar_data              ! model parameters
  type(var_dlength),intent(in)         :: bvar_data              ! basin-average model variables
  type(zLookup),intent(in)             :: lookup_data            ! lookup tables
  ! data structures (input-output)
  type(var_ilength),intent(inout)      :: indx_data              ! state vector geometry
  type(var_dlength),intent(inout)      :: prog_data              ! prognostic variables for a local HRU
  type(var_dlength),intent(inout)      :: diag_data              ! diagnostic variables for a local HRU
  type(var_dlength),intent(inout)      :: flux_data              ! model fluxes for a local HRU
  real(rkind),intent(in)               :: fracJulDay             ! fractional julian days since the start of year
  integer(i4b),intent(in)              :: yearLength             ! number of days in the current year
  ! error control
  integer(i4b),intent(out)             :: err                    ! error code
  character(*),intent(out)             :: message                ! error message
  ! =====================================================================================================================================================
  ! =====================================================================================================================================================
  ! local variables
  character(len=256)                   :: cmessage               ! error message
  integer(i4b)                         :: nSnow                  ! number of snow layers
  integer(i4b)                         :: nSoil                  ! number of soil layers
  integer(i4b)                         :: nLayers                ! total number of layers
  integer(i4b)                         :: nState                 ! total number of state variables
  real(rkind)                          :: dtSave                 ! length of last input model whole sub-step (seconds)
  real(rkind)                          :: dt_sub                 ! length of model sub-step (seconds)
  real(rkind)                          :: dt_wght                ! weight applied to model sub-step (dt_sub/data_step)
  real(rkind)                          :: dt_solv                ! seconds in the data step that have been completed
  real(rkind)                          :: dtMultiplier           ! time step multiplier (-) based on what happenned in "opSplittin"
  real(rkind)                          :: minstep,maxstep        ! minimum and maximum time step length (seconds)
  real(rkind)                          :: maxstep_op             ! maximum time step length (seconds) to run opSplittin over
  real(rkind)                          :: whole_step             ! step the surface pond drainage and sublimation calculated over
  integer(i4b)                         :: nsub                   ! number of substeps
  logical(lgt)                         :: computeVegFluxOld      ! flag to indicate if we are computing fluxes over vegetation on the previous sub step
  logical(lgt)                         :: includeAquifer         ! flag to denote that an aquifer is included
  logical(lgt)                         :: modifiedLayers         ! flag to denote that snow layers were modified
  logical(lgt)                         :: modifiedVegState       ! flag to denote that vegetation states were modified
  integer(i4b)                         :: nLayersRoots           ! number of soil layers that contain roots
  real(rkind)                          :: exposedVAI             ! exposed vegetation area index
  real(rkind)                          :: dCanopyWetFraction_dWat ! derivative in wetted fraction w.r.t. canopy total water (kg-1 m2)
  real(rkind)                          :: dCanopyWetFraction_dT   ! derivative in wetted fraction w.r.t. canopy temperature (K-1)
  real(rkind),parameter                :: varNotUsed1=-9999._rkind ! variables used to calculate derivatives (not needed here)
  real(rkind),parameter                :: varNotUsed2=-9999._rkind ! variables used to calculate derivatives (not needed here)
  integer(i4b)                         :: iSnow                  ! index of snow layers
  integer(i4b)                         :: iLayer                 ! index of model layers
  real(rkind)                          :: massLiquid             ! mass liquid water (kg m-2)
  real(rkind)                          :: superflousSub          ! superflous sublimation (kg m-2 s-1)
  real(rkind)                          :: superflousNrg          ! superflous energy that cannot be used for sublimation (W m-2 [J m-2 s-1])
  integer(i4b)                         :: ixSolution             ! solution method used by opSplittin
  logical(lgt)                         :: firstSubStep           ! flag to denote if the first time step
  logical(lgt)                         :: stepFailure            ! flag to denote the need to reduce length of the coupled step and try again
  logical(lgt)                         :: tooMuchMelt            ! flag to denote that there was too much melt in a given time step
  logical(lgt)                         :: tooMuchSublim          ! flag to denote that there was too much sublimation in a given time step
  logical(lgt)                         :: doLayerMerge           ! flag to denote the need to merge snow layers
  logical(lgt)                         :: pauseFlag              ! flag to pause execution
  logical(lgt),parameter               :: backwardsCompatibility=.true.  ! flag to denote a desire to ensure backwards compatibility with previous branches
  logical(lgt)                         :: checkMassBalance       ! flag to check the mass balance
  type(var_ilength)                    :: indx_temp              ! temporary model index variables saved only on outer loop
  type(var_ilength)                    :: indx_temp0             ! temporary model index variables saved every time
  type(var_dlength)                    :: prog_temp              ! temporary model prognostic variables
  type(var_dlength)                    :: diag_temp              ! temporary model diagnostic variables
  real(rkind),allocatable              :: mLayerVolFracIceInit(:)! initial vector for volumetric fraction of ice (-)
  ! check SWE
  real(rkind)                          :: oldSWE                 ! SWE at the start of the substep
  real(rkind)                          :: newSWE                 ! SWE at the end of the substep
  real(rkind)                          :: delSWE                 ! change in SWE over the subtep
  real(rkind)                          :: innerEffRainfall       ! inner step average effective rainfall into snow (kg m-2 s-1)
  real(rkind)                          :: effRainfall            ! timestep-average effective rainfall into snow (kg m-2 s-1)
  real(rkind)                          :: effSnowfall            ! effective snowfall (kg m-2 s-1)
  real(rkind)                          :: sfcMeltPond            ! surface melt pond (kg m-2)
  real(rkind)                          :: massBalance            ! mass balance error (kg m-2)
  ! energy fluxes
  integer(i4b)                         :: iSoil                  ! index of soil layers
  type(var_dlength)                    :: flux_mean              ! timestep-average model fluxes for a local HRU
  type(var_dlength)                    :: flux_inner             ! inner step average model fluxes for a local HRU
  real(rkind),allocatable              :: meanSoilCompress(:)    ! timestep-average soil compression by layer
  real(rkind),allocatable              :: innerSoilCompress(:)   ! inner step average soil compression by layer
  ! sublimation sums over substep and means over data_step
  real(rkind)                          :: sumCanopySublimation   ! sum of sublimation from the vegetation canopy (kg m-2 s-1) over substep
  real(rkind)                          :: sumSnowSublimation     ! sum of sublimation from the snow surface (kg m-2 s-1) over substep
  real(rkind)                          :: sumLatHeatCanopyEvap   ! sum of latent heat flux for evaporation from the canopy to the canopy air space (W m-2) over substep
  real(rkind)                          :: sumSenHeatCanopy       ! sum of sensible heat flux from the canopy to the canopy air space (W m-2) over substep
  real(rkind)                          :: meanCanopySublimation  ! timestep-average sublimation from the vegetation canopy (kg m-2 s-1)
  real(rkind)                          :: meanLatHeatCanopyEvap  ! timestep-average latent heat flux for evaporation from the canopy to the canopy air space (W m-2)
  real(rkind)                          :: meanSenHeatCanopy      ! timestep-average sensible heat flux from the canopy to the canopy air space (W m-2)
  ! balance checks
  integer(i4b)                         :: iVar                   ! loop through model variables
  real(rkind)                          :: balanceSoilCompress    ! total soil compression (kg m-2)
  real(rkind)                          :: scalarCanopyWatBalError! water balance error for the vegetation canopy (kg m-2)
  real(rkind)                          :: scalarSoilWatBalError  ! water balance error (kg m-2)
  real(rkind)                          :: scalarInitCanopyLiq    ! initial liquid water on the vegetation canopy (kg m-2)
  real(rkind)                          :: scalarInitCanopyIce    ! initial ice          on the vegetation canopy (kg m-2)
  real(rkind)                          :: balanceCanopyWater0    ! total water stored in the vegetation canopy at the start of the step (kg m-2)
  real(rkind)                          :: balanceCanopyWater1    ! total water stored in the vegetation canopy at the end of the step (kg m-2)
  real(rkind)                          :: balanceSoilWater0      ! total soil storage at the start of the step (kg m-2)
  real(rkind)                          :: balanceSoilWater1      ! total soil storage at the end of the step (kg m-2)
  real(rkind)                          :: balanceSoilInflux      ! input to the soil zone
  real(rkind)                          :: balanceSoilBaseflow    ! output from the soil zone
  real(rkind)                          :: balanceSoilDrainage    ! output from the soil zone
  real(rkind)                          :: balanceSoilET          ! output from the soil zone
  real(rkind)                          :: balanceAquifer0        ! total aquifer storage at the start of the step (kg m-2)
  real(rkind)                          :: balanceAquifer1        ! total aquifer storage at the end of the step (kg m-2)
  ! test balance checks
  logical(lgt), parameter              :: printBalance=.false.   ! flag to print the balance checks
  real(rkind), allocatable             :: liqSnowInit(:)         ! volumetric liquid water conetnt of snow at the start of the time step
  real(rkind), allocatable             :: liqSoilInit(:)         ! soil moisture at the start of the time step
  ! timing information
  real(rkind)                          :: startTime              ! start time (used to compute wall clock time)
  real(rkind)                          :: endTime                ! end time (used to compute wall clock time)
  ! outer loop control
  logical(lgt)                         :: firstInnerStep         ! flag to denote if the first time step in maxstep subStep
  logical(lgt)                         :: lastInnerStep          ! flag to denote if the last time step in maxstep subStep
  logical(lgt)                         :: do_outer               ! flag to denote if doing the outer steps surrounding the call to opSplittin
  real(rkind)                          :: dt_solvInner           ! seconds in the maxstep subStep that have been completed

  ! ----------------------------------------------------------------------------------------------------------------------------------------------
  ! initialize error control
  err=0; message="coupled_em/"

  ! This is the start of a data step for a local HRU

  ! get the start time
  call cpu_time(startTime)

  ! check that the decision is supported
  if(model_decisions(iLookDECISIONS%groundwatr)%iDecision==bigBucket .and. &
      model_decisions(iLookDECISIONS%spatial_gw)%iDecision/=localColumn)then
    message=trim(message)//'expect "spatial_gw" decision to equal localColumn when "groundwatr" decision is bigBucket'
    err=20; return
  endif

  ! check if the aquifer is included
  includeAquifer = (model_decisions(iLookDECISIONS%groundwatr)%iDecision==bigBucket)

  ! initialize the numerix tracking variables
  indx_data%var(iLookINDEX%numberFluxCalc       )%dat(1) = 0  ! number of flux calculations                     (-)
  indx_data%var(iLookINDEX%numberStateSplit     )%dat(1) = 0  ! number of state splitting solutions             (-)
  indx_data%var(iLookINDEX%numberDomainSplitNrg )%dat(1) = 0  ! number of domain splitting solutions for energy (-)
  indx_data%var(iLookINDEX%numberDomainSplitMass)%dat(1) = 0  ! number of domain splitting solutions for mass   (-)
  indx_data%var(iLookINDEX%numberScalarSolutions)%dat(1) = 0  ! number of scalar solutions                      (-)

  ! link canopy depth to the information in the data structure
  canopy: associate(canopyDepth => diag_data%var(iLookDIAG%scalarCanopyDepth)%dat(1) )  ! intent(out): [dp] canopy depth (m)

    ! start by NOT pausing
    pauseFlag=.false.

    ! start by assuming that the step is successful
    stepFailure  = .false.
    doLayerMerge = .false.

    ! initialize flags to modify the veg layers or modify snow layers
    modifiedLayers    = .false.    ! flag to denote that snow layers were modified
    modifiedVegState  = .false.    ! flag to denote that vegetation states were modified

    ! define the first step and first and last inner steps
    firstSubStep = .true.
    firstInnerStep = .true.
    lastInnerStep = .false.

    ! count the number of snow and soil layers
    ! NOTE: need to re-compute the number of snow and soil layers at the start of each sub-step because the number of layers may change
    !         (nSnow and nSoil are shared in the data structure)
    nSnow = count(indx_data%var(iLookINDEX%layerType)%dat==iname_snow)
    nSoil = count(indx_data%var(iLookINDEX%layerType)%dat==iname_soil)

    ! compute the total number of snow and soil layers
    nLayers = nSnow + nSoil

    ! create temporary data structures for prognostic variables
    call resizeData(prog_meta(:),prog_data,prog_temp,err=err,message=cmessage)
    if(err/=0)then; err=20; message=trim(message)//trim(cmessage); return; endif

    ! create temporary data structures for diagnostic variables
    call resizeData(diag_meta(:),diag_data,diag_temp,err=err,message=cmessage)
    if(err/=0)then; err=20; message=trim(message)//trim(cmessage); return; endif

    ! create temporary data structures for index variables
    call resizeData(indx_meta(:),indx_data,indx_temp,err=err,message=cmessage)
    call resizeData(indx_meta(:),indx_data,indx_temp0,err=err,message=cmessage)
    if(err/=0)then; err=20; message=trim(message)//trim(cmessage); return; endif

    ! allocate space for the local fluxes
    call allocLocal(averageFlux_meta(:)%var_info,flux_mean,nSnow,nSoil,err,cmessage)
    if(err/=0)then; err=20; message=trim(message)//trim(cmessage); return; end if
    call allocLocal(averageFlux_meta(:)%var_info,flux_inner,nSnow,nSoil,err,cmessage)
    if(err/=0)then; err=20; message=trim(message)//trim(cmessage); return; end if

    ! initialize surface melt pond
    sfcMeltPond       = 0._rkind  ! change in storage associated with the surface melt pond (kg m-2)

    ! initialize fluxes to average over data_step (averaged over substep in varSubStep)
    do iVar=1,size(averageFlux_meta)
      flux_mean%var(iVar)%dat(:) = 0._rkind
    end do
    meanCanopySublimation = 0._rkind ! mean canopy sublimation
    meanLatHeatCanopyEvap = 0._rkind ! mean latent heat flux for evaporation from the canopy
    meanSenHeatCanopy     = 0._rkind ! mean sensible heat flux from the canopy
    effRainfall           = 0._rkind ! mean total effective rainfall over snow

    ! Need mean soil compression for balance checks but it is not in flux structure so handle differently 
    !  This will be a problem if nSoil changes (currently not possible)-- then might need to not keep the average
    allocate(meanSoilCompress(nSoil))
    allocate(innerSoilCompress(nSoil))
    meanSoilCompress = 0._rkind ! mean total soil compression

    ! associate local variables with information in the data structures
    associate(&
    ! state variables in the vegetation canopy
    scalarCanopyLiq      => prog_data%var(iLookPROG%scalarCanopyLiq)%dat(1)                 ,&  ! canopy liquid water (kg m-2)
    scalarCanopyIce      => prog_data%var(iLookPROG%scalarCanopyIce)%dat(1)                 ,&  ! canopy ice content (kg m-2)
    ! state variables in the soil domain
    mLayerDepth          => prog_data%var(iLookPROG%mLayerDepth)%dat(nSnow+1:nLayers)       ,&  ! depth of each soil layer (m)
    mLayerVolFracIce     => prog_data%var(iLookPROG%mLayerVolFracIce)%dat(nSnow+1:nLayers)  ,&  ! volumetric ice content in each soil layer (-)
    mLayerVolFracLiq     => prog_data%var(iLookPROG%mLayerVolFracLiq)%dat(nSnow+1:nLayers)  ,&  ! volumetric liquid water content in each soil layer (-)
    scalarAquiferStorage => prog_data%var(iLookPROG%scalarAquiferStorage)%dat(1)            ,&  ! aquifer storage (m)
    scalarTotalSoilIce   => diag_data%var(iLookDIAG%scalarTotalSoilIce)%dat(1)              ,&  ! total ice in the soil column (kg m-2)
    scalarTotalSoilLiq   => diag_data%var(iLookDIAG%scalarTotalSoilLiq)%dat(1)              ,&  ! total liquid water in the soil column (kg m-2)
    scalarTotalSoilWat   => diag_data%var(iLookDIAG%scalarTotalSoilWat)%dat(1)               &  ! total water in the soil column (kg m-2)
    ) ! (association of local variables with information in the data structures

      ! save the liquid water and ice on the vegetation canopy
      scalarInitCanopyLiq = scalarCanopyLiq    ! initial liquid water on the vegetation canopy (kg m-2)
      scalarInitCanopyIce = scalarCanopyIce    ! initial ice          on the vegetation canopy (kg m-2)

      ! compute total soil moisture and ice at the *START* of the step (kg m-2)
      scalarTotalSoilLiq = sum(iden_water*mLayerVolFracLiq(1:nSoil)*mLayerDepth(1:nSoil))
      scalarTotalSoilIce = sum(iden_water*mLayerVolFracIce(1:nSoil)*mLayerDepth(1:nSoil))  ! NOTE: no expansion and hence use iden_water

      ! compute storage of water in the canopy and the soil
      balanceCanopyWater0 = scalarCanopyLiq + scalarCanopyIce
      balanceSoilWater0   = scalarTotalSoilLiq + scalarTotalSoilIce

      ! get the total aquifer storage at the start of the time step (kg m-2)
      balanceAquifer0 = scalarAquiferStorage*iden_water

      ! save liquid water content
      if(printBalance)then
        allocate(liqSnowInit(nSnow), liqSoilInit(nSoil), stat=err)
        if(err/=0)then
          message=trim(message)//'unable to allocate space for the initial vectors'
          err=20; return
        endif
        if(nSnow>0) liqSnowInit = prog_data%var(iLookPROG%mLayerVolFracLiq)%dat(1:nSnow)
                    liqSoilInit =                         mLayerVolFracLiq
      endif

    ! end association of local variables with information in the data structures
    end associate

    ! short-cut to the algorithmic control parameters
    ! NOTE - temporary assignment of minstep to foce something reasonable
    ! changing the maxstep parameter will make the outer and inner loop computations here in coupled_em happen more frequently
    ! changing the be_steps parameter will make the inner loop computations in opSplittin happen more frequently (e.g. be_steps = 32.0 give BE32)
    minstep = 10._rkind  ! mpar_data%var(iLookPARAM%minstep)%dat(1)  ! minimum time step (s)
    maxstep = mpar_data%var(iLookPARAM%maxstep)%dat(1)  ! maximum time step (s)
    maxstep_op = mpar_data%var(iLookPARAM%maxstep)%dat(1)/NINT(mpar_data%var(iLookPARAM%be_steps)%dat(1))  ! maximum time step (s) to run opSplittin over

    ! compute the number of layers with roots
    nLayersRoots = count(prog_data%var(iLookPROG%iLayerHeight)%dat(nSnow:nLayers-1) < mpar_data%var(iLookPARAM%rootingDepth)%dat(1)-verySmall)
    if(nLayersRoots == 0)then
      message=trim(message)//'no roots within the soil profile'
      err=20; return
    end if

    ! define the foliage nitrogen factor
    diag_data%var(iLookDIAG%scalarFoliageNitrogenFactor)%dat(1) = 1._rkind  ! foliage nitrogen concentration (1.0 = saturated)

    ! save SWE
    oldSWE = prog_data%var(iLookPROG%scalarSWE)%dat(1)

    ! *** compute phenology...
    ! ------------------------

    ! compute the temperature of the root zone: used in vegetation phenology
    diag_data%var(iLookDIAG%scalarRootZoneTemp)%dat(1) = sum(prog_data%var(iLookPROG%mLayerTemp)%dat(nSnow+1:nSnow+nLayersRoots)) / real(nLayersRoots, kind(rkind))

    ! remember if we compute the vegetation flux on the previous sub-step
    computeVegFluxOld = computeVegFlux

    ! compute the exposed LAI and SAI and whether veg is buried by snow
    call vegPhenlgy(&
                    ! model control
                    fracJulDay,                  & ! intent(in):    fractional julian days since the start of year
                    yearLength,                  & ! intent(in):    number of days in the current year
                    ! input/output: data structures
                    model_decisions,             & ! intent(in):    model decisions
                    type_data,                   & ! intent(in):    type of vegetation and soil
                    attr_data,                   & ! intent(in):    spatial attributes
                    mpar_data,                   & ! intent(in):    model parameters
                    prog_data,                   & ! intent(inout): model prognostic variables for a local HRU
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

    ! *** compute wetted canopy area...
    ! ---------------------------------

    ! compute maximum canopy liquid water (kg m-2)
    diag_data%var(iLookDIAG%scalarCanopyLiqMax)%dat(1) = mpar_data%var(iLookPARAM%refInterceptCapRain)%dat(1)*exposedVAI

    ! compute maximum canopy ice content (kg m-2)
    ! NOTE 1: this is used to compute the snow fraction on the canopy, as used in *BOTH* the radiation AND canopy sublimation routines
    ! NOTE 2: this is a different variable than the max ice used in the throughfall (snow interception) calculations
    ! NOTE 3: use maximum per unit leaf area storage capacity for snow (kg m-2)
    select case(model_decisions(iLookDECISIONS%snowIncept)%iDecision)
      case(lightSnow);  diag_data%var(iLookDIAG%scalarCanopyIceMax)%dat(1) = exposedVAI*mpar_data%var(iLookPARAM%refInterceptCapSnow)%dat(1)
      case(stickySnow); diag_data%var(iLookDIAG%scalarCanopyIceMax)%dat(1) = exposedVAI*mpar_data%var(iLookPARAM%refInterceptCapSnow)%dat(1)*4._rkind
      case default; message=trim(message)//'unable to identify option for maximum branch interception capacity'; err=20; return
    end select ! identifying option for maximum branch interception capacity

    ! compute wetted fraction of the canopy
    ! NOTE: assume that the wetted fraction is constant over the substep for the radiation calculations
    if(computeVegFlux)then

      ! compute wetted fraction of the canopy
      call wettedFrac(&
                      ! input
                      .false.,                                                      & ! flag to denote if derivatives are required
                      (prog_data%var(iLookPROG%scalarCanopyTemp)%dat(1) < Tfreeze), & ! flag to denote if the canopy is frozen
                      varNotUsed1,                                                  & ! derivative in canopy liquid w.r.t. canopy temperature (kg m-2 K-1)
                      varNotUsed2,                                                  & ! fraction of liquid water on the canopy
                      prog_data%var(iLookPROG%scalarCanopyLiq)%dat(1),              & ! canopy liquid water (kg m-2)
                      prog_data%var(iLookPROG%scalarCanopyIce)%dat(1),              & ! canopy ice (kg m-2)
                      diag_data%var(iLookDIAG%scalarCanopyLiqMax)%dat(1),           & ! maximum canopy liquid water (kg m-2)
                      diag_data%var(iLookDIAG%scalarCanopyIceMax)%dat(1),           & ! maximum canopy ice content (kg m-2)
                      mpar_data%var(iLookPARAM%canopyWettingFactor)%dat(1),         & ! maximum wetted fraction of the canopy (-)
                      mpar_data%var(iLookPARAM%canopyWettingExp)%dat(1),            & ! exponent in canopy wetting function (-)
                      ! output
                      diag_data%var(iLookDIAG%scalarCanopyWetFraction)%dat(1),      & ! canopy wetted fraction (-)
                      dCanopyWetFraction_dWat,                                      & ! derivative in wetted fraction w.r.t. canopy liquid water content (kg-1 m2)
                      dCanopyWetFraction_dT,                                        & ! derivative in wetted fraction w.r.t. canopy liquid water content (kg-1 m2)
                      err,cmessage)
      if(err/=0)then; message=trim(message)//trim(cmessage); return; end if

    ! vegetation is completely buried by snow (or no veg exists at all)
    else
      diag_data%var(iLookDIAG%scalarCanopyWetFraction)%dat(1) = 0._rkind
      dCanopyWetFraction_dWat                                 = 0._rkind
      dCanopyWetFraction_dT                                   = 0._rkind
    end if

    ! *** compute snow albedo...
    ! --------------------------
    ! NOTE: this should be done before the radiation calculations
    ! NOTE: uses snowfall; should really use canopy throughfall + canopy unloading
    call snowAlbedo(&
                    ! input: model control
                    data_step,                   & ! intent(in): model time step (s)
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


    ! *** compute canopy sw radiation fluxes...
    ! -----------------------------------------
    call vegSWavRad(&
                    data_step,                    & ! intent(in):    time step (s) -- only used in Noah-MP radiation, to compute albedo
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


    ! *** compute canopy throughfall and unloading...
    ! -----------------------------------------------
    ! NOTE 1: this needs to be done before solving the energy and liquid water equations, to account for the heat advected with precipitation (and throughfall/unloading)
    ! NOTE 2: the unloading flux is computed using canopy drip (scalarCanopyLiqDrainage) from the previous time step
    ! this changes canopy ice
    call canopySnow(&
                    ! input: model control
                    data_step,                   & ! intent(in): time step (seconds)
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

    ! adjust canopy temperature to account for new snow
    if(computeVegFlux)then ! logical flag to compute vegetation fluxes (.false. if veg buried by snow)
      call tempAdjust(&
                      ! input: derived parameters
                      canopyDepth,                 & ! intent(in):    canopy depth (m)
                      ! input/output: data structures
                      mpar_data,                   & ! intent(in):    model parameters
                      prog_data,                   & ! intent(inout): model prognostic variables for a local HRU
                      diag_data,                   & ! intent(inout): model diagnostic variables for a local HRU
                      ! output: error control
                      err,cmessage)                  ! intent(out): error control
                      if(err/=0)then; err=20; message=trim(message)//trim(cmessage); return; end if
    endif ! if computing fluxes over vegetation

    ! initialize drainage and throughfall
    ! NOTE 1: this needs to be done before solving the energy and liquid water equations, to account for the heat advected with precipitation
    ! NOTE 2: this initialization needs to be done AFTER the call to canopySnow, since canopySnow uses canopy drip drom the previous time step
    if(.not.computeVegFlux)then
      flux_data%var(iLookFLUX%scalarThroughfallRain)%dat(1) = flux_data%var(iLookFLUX%scalarRainfall)%dat(1)
    else
      flux_data%var(iLookFLUX%scalarThroughfallRain)%dat(1) = 0._rkind
    end if
    flux_data%var(iLookFLUX%scalarCanopyLiqDrainage)%dat(1) = 0._rkind

    ! ****************************************************************************************************
    ! *** MAIN SOLVER ************************************************************************************
    ! ****************************************************************************************************

    ! initialize the length of the sub-step and counters
    whole_step = maxstep
    dt_solv       = 0._rkind   ! length of time step that has been completed (s)
    dt_solvInner  = 0._rkind   ! length of time step that has been completed (s) in whole_step subStep
    dt_init = min(data_step,whole_step,maxstep_op) / dt_init_factor  ! initial substep length (s)
    dt_sub = dt_init
    dtSave  = whole_step       ! length of whole substep

    ! initialize the number of sub-steps
    nsub=0

    ! loop through sub-steps
    substeps: do  ! continuous do statement with exit clause (alternative to "while")

      dt_sub = min(data_step,whole_step,maxstep_op,dt_sub) ! adjust for possible whole_step changes

      ! print progress
      if(globalPrintFlag)then
        write(*,'(a,1x,4(f13.5,1x))') ' start of step: dt_init, dt_sub, dt_solv, data_step: ', dt_init, dt_sub, dt_solv, data_step
        print*, 'stepFailure = ', stepFailure
        print*, 'before resizeData: nSnow, nSoil = ', nSnow, nSoil
      endif

      ! increment the number of sub-steps
      nsub = nsub+1

      ! resize the "indx_data" structure
      ! NOTE: this is necessary because the length of index variables depends on a given split
      !        --> the resize here is overwritten later (in indexSplit)
      !        --> admittedly ugly, and retained for now
      if(stepFailure)then ! resize temp to current data, later in code current data is set to lastInnerStep data
        call resizeData(indx_meta(:),indx_temp,indx_data,err=err,message=cmessage)
        if(err/=0)then; err=20; message=trim(message)//trim(cmessage); return; endif
      else ! resize current data to temp0, temp0 is saved for next run
        call resizeData(indx_meta(:),indx_data,indx_temp0,err=err,message=cmessage)
        if(err/=0)then; err=20; message=trim(message)//trim(cmessage); return; endif
        do iVar=1,size(indx_data%var)
          indx_temp0%var(iVar)%dat(:) = indx_data%var(iVar)%dat(:)
        end do
      endif

      ! check if on outer loop, always do outer if after failed step and on then on reduced whole_step
      do_outer = .false.
      if(stepFailure) firstInnerStep = .true.
      if(firstInnerStep) do_outer = .true.

      if(do_outer)then

        if(.not.stepFailure)then
          call resizeData(indx_meta(:),indx_data,indx_temp,err=err,message=cmessage)
          if(err/=0)then; err=20; message=trim(message)//trim(cmessage); return; endif
        endif

        ! save/recover copies of index variables, temp saved on lastInnerStep, failed starts at lastInnerStep
        do iVar=1,size(indx_data%var)
          select case(stepFailure)
            case(.false.); indx_temp%var(iVar)%dat(:) = indx_data%var(iVar)%dat(:)
            case(.true.);  indx_data%var(iVar)%dat(:) = indx_temp%var(iVar)%dat(:)
          end select
        end do  ! looping through variables

        ! save/recover copies of prognostic variables
        do iVar=1,size(prog_data%var)
          select case(stepFailure)
            case(.false.); prog_temp%var(iVar)%dat(:) = prog_data%var(iVar)%dat(:)
            case(.true.);  prog_data%var(iVar)%dat(:) = prog_temp%var(iVar)%dat(:)
          end select
        end do  ! looping through variables

        ! save/recover copies of diagnostic variables
        do iVar=1,size(diag_data%var)
          select case(stepFailure)
            case(.false.); diag_temp%var(iVar)%dat(:) = diag_data%var(iVar)%dat(:)
            case(.true.);  diag_data%var(iVar)%dat(:) = diag_temp%var(iVar)%dat(:)
          end select
        end do  ! looping through variables

        ! re-assign dimension lengths
        nSnow   = count(indx_data%var(iLookINDEX%layerType)%dat==iname_snow)
        nSoil   = count(indx_data%var(iLookINDEX%layerType)%dat==iname_soil)
        nLayers = nSnow+nSoil

        ! *** merge/sub-divide snow layers...
        ! -----------------------------------
        call volicePack(&
                        ! input/output: model data structures
                        doLayerMerge,               & ! intent(in):    flag to force merge of snow layers
                        model_decisions,            & ! intent(in):    model decisions
                        mpar_data,                  & ! intent(in):    model parameters
                        indx_data,                  & ! intent(inout): type of each layer
                        prog_data,                  & ! intent(inout): model prognostic variables for a local HRU
                        diag_data,                  & ! intent(inout): model diagnostic variables for a local HRU
                        flux_data,                  & ! intent(inout): model fluxes for a local HRU
                        ! output
                        modifiedLayers,             & ! intent(out): flag to denote that layers were modified
                        err,cmessage)                ! intent(out): error control
        if(err/=0)then; err=55; message=trim(message)//trim(cmessage); return; end if

        ! save the number of snow and soil layers
        nSnow   = indx_data%var(iLookINDEX%nSnow)%dat(1)
        nSoil   = indx_data%var(iLookINDEX%nSoil)%dat(1)
        nLayers = indx_data%var(iLookINDEX%nLayers)%dat(1)

        ! compute the indices for the model state variables
        if(firstSubStep .or. modifiedVegState .or. modifiedLayers)then
          call indexState(computeVegFlux,         & ! intent(in):    flag to denote if computing the vegetation flux
                          includeAquifer,         & ! intent(in):    flag to denote if included the aquifer
                          nSnow,nSoil,nLayers,    & ! intent(in):    number of snow and soil layers, and total number of layers
                          indx_data,              & ! intent(inout): indices defining model states and layers
                          err,cmessage)             ! intent(out):   error control
          if(err/=0)then; message=trim(message)//trim(cmessage); return; end if

        end if

        ! recreate the temporary data structures
        ! NOTE: resizeData(meta, old, new, ..)
        if(modifiedVegState .or. modifiedLayers)then

          ! create temporary data structures for prognostic variables
          call resizeData(prog_meta(:),prog_data,prog_temp,copy=.true.,err=err,message=cmessage)
          if(err/=0)then; err=20; message=trim(message)//trim(cmessage); return; endif

          ! create temporary data structures for diagnostic variables
          call resizeData(diag_meta(:),diag_data,diag_temp,copy=.true.,err=err,message=cmessage)
          if(err/=0)then; err=20; message=trim(message)//trim(cmessage); return; endif

          ! create temporary data structures for index variables
          call resizeData(indx_meta(:),indx_data,indx_temp,copy=.true.,err=err,message=cmessage)
          if(err/=0)then; err=20; message=trim(message)//trim(cmessage); return; endif

          do iVar=1,size(indx_data%var)
            select case(stepFailure)
              case(.false.); indx_temp%var(iVar)%dat(:) = indx_data%var(iVar)%dat(:)
              case(.true.);  indx_data%var(iVar)%dat(:) = indx_temp%var(iVar)%dat(:)
            end select
          end do  ! looping through variables

        endif  ! if modified the states

        ! define the number of state variables
        nState = indx_data%var(iLookINDEX%nState)%dat(1)

        ! *** compute diagnostic variables for each layer...
        ! --------------------------------------------------
        ! NOTE: this needs to be done AFTER volicePack, since layers may have been sub-divided and/or merged
        call diagn_evar(&
                        ! input: control variables
                        computeVegFlux,         & ! intent(in): flag to denote if computing the vegetation flux
                        diag_data%var(iLookDIAG%scalarCanopyDepth)%dat(1),            & ! intent(in): canopy depth (m)
                        ! input/output: data structures
                        mpar_data,              & ! intent(in):    model parameters
                        indx_data,              & ! intent(in):    model layer indices
                        prog_data,              & ! intent(in):    model prognostic variables for a local HRU
                        diag_data,              & ! intent(inout): model diagnostic variables for a local HRU
                        ! output: error control
                        err,cmessage)             ! intent(out): error control
        if(err/=0)then; err=55; message=trim(message)//trim(cmessage); return; end if

        ! *** compute melt of the "snow without a layer"...
        ! -------------------------------------------------
        ! NOTE: forms a surface melt pond, which drains into the upper-most soil layer through the time step
        ! (check for the special case of "snow without a layer")
        ! this pond melts evenly over entire time of maxstep until it gets recomputed because based on SWE when computed
        if(nSnow==0) then
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
                          err,cmessage                                        ) ! intent(out): error control
          if(err/=0)then; err=20; message=trim(message)//trim(cmessage); return; end if
        endif

        ! save volumetric ice content at the start of the step
        ! NOTE: used for volumetric loss due to melt-freeze
        allocate(mLayerVolFracIceInit(nLayers)); mLayerVolFracIceInit = prog_data%var(iLookPROG%mLayerVolFracIce)%dat

        ! make sure have consistent state variables to start, later done in updateVars
        ! associate local variables with information in the data structures
        init: associate(&
          ! depth-varying soil parameters
          vGn_m                   => diag_data%var(iLookDIAG%scalarVGn_m)%dat               ,& ! intent(in):    [dp(:)]  van Genutchen "m" parameter (-)
          vGn_n                   => mpar_data%var(iLookPARAM%vGn_n)%dat                    ,& ! intent(in):    [dp(:)]  van Genutchen "n" parameter (-)
          vGn_alpha               => mpar_data%var(iLookPARAM%vGn_alpha)%dat                ,& ! intent(in):    [dp(:)]  van Genutchen "alpha" parameter (m-1)
          theta_sat               => mpar_data%var(iLookPARAM%theta_sat)%dat                ,& ! intent(in):    [dp(:)]  soil porosity (-)
          theta_res               => mpar_data%var(iLookPARAM%theta_res)%dat                ,& ! intent(in):    [dp(:)]  soil residual volumetric water content (-)
          ! state variables in the vegetation canopy
          scalarCanopyIce         => prog_data%var(iLookPROG%scalarCanopyIce)%dat(1)        ,& ! intent(in):    [dp]     mass of ice on the vegetation canopy (kg m-2)
          scalarCanopyLiq         => prog_data%var(iLookPROG%scalarCanopyLiq)%dat(1)        ,& ! intent(in):    [dp]     mass of liquid water on the vegetation canopy (kg m-2)
          scalarCanopyWat         => prog_data%var(iLookPROG%scalarCanopyWat)%dat(1)        ,& ! intent(out):   [dp]     mass of total water on the vegetation canopy (kg m-2)
          ! state variables in the snow and soil domains
          mLayerVolFracIce        => prog_data%var(iLookPROG%mLayerVolFracIce)%dat          ,& ! intent(in):    [dp(:)]  volumetric fraction of ice (-)
          mLayerVolFracLiq        => prog_data%var(iLookPROG%mLayerVolFracLiq)%dat          ,& ! intent(in):    [dp(:)]  volumetric fraction of liquid water (-)
          mLayerVolFracWat        => prog_data%var(iLookPROG%mLayerVolFracWat)%dat          ,& ! intent(out):   [dp(:)]  volumetric fraction of total water (-)
          mLayerMatricHead        => prog_data%var(iLookPROG%mLayerMatricHead)%dat          ,& ! intent(in):    [dp(:)]  matric head (m)
          mLayerMatricHeadLiq     => diag_data%var(iLookDIAG%mLayerMatricHeadLiq)%dat        & ! intent(out):   [dp(:)]  matric potential of liquid water (m)
          ) ! associations to variables in data structures

          ! compute the total water content in the vegetation canopy
          scalarCanopyWat = scalarCanopyLiq + scalarCanopyIce  ! kg m-2

          ! compute the total water content in snow and soil
          ! NOTE: no ice expansion allowed for soil
          if(nSnow>0)&
          mLayerVolFracWat(      1:nSnow  ) = mLayerVolFracLiq(      1:nSnow  ) + mLayerVolFracIce(      1:nSnow  )*(iden_ice/iden_water)
          mLayerVolFracWat(nSnow+1:nLayers) = mLayerVolFracLiq(nSnow+1:nLayers) + mLayerVolFracIce(nSnow+1:nLayers)

          ! compute the liquid water matric potential (m)
          ! NOTE: include ice content as part of the solid porosity - major effect of ice is to reduce the pore size; ensure that effSat=1 at saturation
          ! (from Zhao et al., J. Hydrol., 1997: Numerical analysis of simultaneous heat and mass transfer...)
          do iSoil=1,nSoil
            call liquidHead(mLayerMatricHead(iSoil),mLayerVolFracLiq(nSnow+iSoil),mLayerVolFracIce(nSnow+iSoil),  & ! input:  state variables
                      vGn_alpha(iSoil),vGn_n(iSoil),theta_sat(iSoil),theta_res(iSoil),vGn_m(iSoil),         & ! input:  parameters
                      matricHeadLiq=mLayerMatricHeadLiq(iSoil),                                             & ! output: liquid water matric potential (m)
                      err=err,message=cmessage)                                                               ! output: error control
            if(err/=0)then; message=trim(message)//trim(cmessage); return; endif
          end do  ! looping through soil layers (computing liquid water matric potential)

        end associate init

        ! correct increments (if need to redo inner step) and reset increment
        dt_solv = dt_solv - dt_solvInner
        dt_solvInner = 0._rkind
        lastInnerStep = .false.

        ! initialize sublimation sums to average over whole_step
        sumCanopySublimation = 0._rkind
        sumSnowSublimation   = 0._rkind
        sumLatHeatCanopyEvap = 0._rkind
        sumSenHeatCanopy     = 0._rkind
        ! initialize fluxes to average over whole_step (averaged over substep in varSubStep)
        do iVar=1,size(averageFlux_meta)
          flux_inner%var(iVar)%dat(:) = 0._rkind
        end do
        innerEffRainfall  = 0._rkind ! mean total effective rainfall over snow
        innerSoilCompress = 0._rkind ! mean total soil compression

      endif ! (do_outer loop)

      ! *** solve model equations...
      ! ----------------------------
      ! save input step
      dtSave = whole_step

      ! get the new solution
      call opSplittin(&
                      ! input: model control
                      nSnow,                                  & ! intent(in):    number of snow layers
                      nSoil,                                  & ! intent(in):    number of soil layers
                      nLayers,                                & ! intent(in):    total number of layers
                      nState,                                 & ! intent(in):    total number of layers
                      dt_sub,                                 & ! intent(in):    length of the model sub-step
                      whole_step,                             & ! intent(in):    length of whole step for surface drainage and average flux
                      (dt_solv<whole_step),                   & ! intent(in):    logical flag to denote the first loop of the whole_step in a data_step
                      firstInnerStep,                         & ! intent(in):    flag to denote if the first time step in maxstep subStep
                      computeVegFlux,                         & ! intent(in):    logical flag to compute fluxes within the vegetation canopy
                      ! input/output: data structures
                      type_data,                              & ! intent(in):    type of vegetation and soil
                      attr_data,                              & ! intent(in):    spatial attributes
                      forc_data,                              & ! intent(in):    model forcing data
                      mpar_data,                              & ! intent(in):    model parameters
                      indx_data,                              & ! intent(inout): index data
                      prog_data,                              & ! intent(inout): model prognostic variables for a local HRU
                      diag_data,                              & ! intent(inout): model diagnostic variables for a local HRU
                      flux_data,                              & ! intent(inout): model fluxes for a local HRU
                      bvar_data,                              & ! intent(in):    model variables for the local basin
                      lookup_data,                            & ! intent(in):    lookup tables
                      model_decisions,                        & ! intent(in):    model decisions
                      ! output: model control
                      dtMultiplier,                           & ! intent(out):   substep multiplier (-)
                      tooMuchMelt,                            & ! intent(out):   flag to denote that ice is insufficient to support melt
                      stepFailure,                            & ! intent(out):   flag to denote that the coupled step failed
                      ixSolution,                             & ! intent(out):   solution method used in this iteration
                      err,cmessage)                             ! intent(out):   error code and error message
      ! check for all errors (error recovery within opSplittin)
      if(err/=0)then; err=20; message=trim(message)//trim(cmessage); return; end if

      ! process the flag for too much melt
      if(tooMuchMelt)then
        stepFailure  = .true.
        doLayerMerge = .true.
      else
        doLayerMerge = .false.
      endif

      ! handle special case of the step failure
      ! NOTE: need to revert back to the previous state vector that we were happy with and reduce the time step
      ! TODO: ask isn't this what the actors program does without the code block below
      if(stepFailure)then
        ! halve whole_step, for more frequent outer loop updates
        whole_step = dtSave/2._rkind
        ! check that the step is not tiny
        if(whole_step < minstep)then
          print*,ixSolution
          print*, 'dtSave, dt_sub', dtSave, whole_step
          message=trim(message)//'length of the coupled step is below the minimum step length'
          err=20; return
        endif
        ! try again, restart step
        deallocate(mLayerVolFracIceInit)
        cycle substeps
      endif

      ! increment sublimation sums
      sumCanopySublimation = sumCanopySublimation + dt_sub*flux_data%var(iLookFLUX%scalarCanopySublimation)%dat(1)
      sumSnowSublimation   = sumSnowSublimation   + dt_sub*flux_data%var(iLookFLUX%scalarSnowSublimation)%dat(1)
      sumLatHeatCanopyEvap = sumLatHeatCanopyEvap + dt_sub*flux_data%var(iLookFLUX%scalarLatHeatCanopyEvap)%dat(1)
      sumSenHeatCanopy     = sumSenHeatCanopy     + dt_sub*flux_data%var(iLookFLUX%scalarSenHeatCanopy)%dat(1)

      ! update first step and first and last inner steps
      firstSubStep = .false.
      firstInnerStep = .false.
      if(dt_solvInner + dt_sub >= whole_step) lastInnerStep = .true.
      if(dt_solv + dt_sub >= data_step-verySmall) lastInnerStep = .true.

      ! check if on outer loop
      do_outer = .false.
      if(lastInnerStep) do_outer = .true.

      if(do_outer)then

        ! ***  remove ice due to sublimation and freeze calculations...
        ! NOTE: In the future this should be moved into the solver, makes a big difference
        ! --------------------------------------------------------------
        sublime: associate(&
          mLayerMeltFreeze        => diag_data%var(iLookDIAG%mLayerMeltFreeze)%dat,           & ! melt-freeze in each snow and soil layer (kg m-3)
          scalarCanopyLiq         => prog_data%var(iLookPROG%scalarCanopyLiq)%dat(1),         & ! liquid water stored on the vegetation canopy (kg m-2)
          scalarCanopyIce         => prog_data%var(iLookPROG%scalarCanopyIce)%dat(1),         & ! ice          stored on the vegetation canopy (kg m-2)
          scalarCanopyWat         => prog_data%var(iLookPROG%scalarCanopyWat)%dat(1),         & ! canopy ice content (kg m-2)
          mLayerVolFracIce        => prog_data%var(iLookPROG%mLayerVolFracIce)%dat,           & ! volumetric fraction of ice in the snow+soil domain (-)
          mLayerVolFracLiq        => prog_data%var(iLookPROG%mLayerVolFracLiq)%dat,           & ! volumetric fraction of liquid water in the snow+soil domain (-)
          mLayerVolFracWat        => prog_data%var(iLookPROG%mLayerVolFracWat)%dat,           & ! volumetric fraction of total water (-)
          mLayerDepth             => prog_data%var(iLookPROG%mLayerDepth)%dat                 & ! depth of each snow+soil layer (m)
          ) ! associations to variables in data structures

          ! compute the melt in each snow and soil layer
          if(nSnow>0)&
          mLayerMeltFreeze(1:nSnow) = -( mLayerVolFracIce(1:nSnow) - mLayerVolFracIceInit(1:nSnow) ) * iden_ice
          mLayerMeltFreeze(nSnow+1:nLayers) = -( mLayerVolFracIce(nSnow+1:nLayers) - mLayerVolFracIceInit(nSnow+1:nLayers) )*iden_water
          deallocate(mLayerVolFracIceInit)

          ! * compute change in canopy ice content due to sublimation...
          ! ------------------------------------------------------------
          if(computeVegFlux)then

            ! remove mass of ice on the canopy
            scalarCanopyIce = scalarCanopyIce + sumCanopySublimation

            ! if removed all ice, take the remaining sublimation from water
            if(scalarCanopyIce < 0._rkind)then
              scalarCanopyLiq = scalarCanopyLiq + scalarCanopyIce
              scalarCanopyIce = 0._rkind
            endif

            ! modify fluxes and mean fluxes if there is insufficient canopy water to support the converged sublimation rate over the whole time step
            if(scalarCanopyLiq < 0._rkind)then
              ! --> superfluous sublimation flux
              superflousSub = -scalarCanopyLiq/whole_step  ! kg m-2 s-1
              superflousNrg = superflousSub*LH_sub     ! W m-2 (J m-2 s-1)
              ! --> update fluxes and states
              sumCanopySublimation = sumCanopySublimation + superflousSub*whole_step
              sumLatHeatCanopyEvap = sumLatHeatCanopyEvap + superflousNrg*whole_step
              sumSenHeatCanopy     = sumSenHeatCanopy     - superflousNrg*whole_step
              scalarCanopyLiq      = 0._rkind
            endif

            ! update water
            scalarCanopyWat = scalarCanopyLiq + scalarCanopyIce

          end if  ! (if computing the vegetation flux)

          call computSnowDepth(&
                    whole_step,                               & ! intent(in)
                    nSnow,                                    & ! intent(in)
                    sumSnowSublimation/whole_step,            & ! intent(in)
                    mLayerVolFracLiq,                         & ! intent(inout)
                    mLayerVolFracIce,                         & ! intent(inout)
                    prog_data%var(iLookPROG%mLayerTemp)%dat,  & ! intent(in)
                    mLayerMeltFreeze,                         & ! intent(in)
                    mpar_data,                                & ! intent(in)
                    ! output
                    tooMuchSublim,                            & ! intent(out): flag to denote that there was too much sublimation in a given time step
                    mLayerDepth,                              & ! intent(inout)
                    ! error control
                    err,message)                                     ! intent(out):   error control
          if(err/=0)then; err=55; return; end if

          ! process the flag for too much sublimation
          if(tooMuchSublim)then
            stepFailure  = .true.
            doLayerMerge = .true.
          else
            doLayerMerge = .false.
          endif

          ! handle special case of the step failure
          ! NOTE: need to revert back to the previous state vector that we were happy with and reduce the time step
          if(stepFailure)then
            ! halve whole_step, for more frequent outer loop updates
            whole_step = dtSave/2._rkind
            ! check that the step is not tiny
            if(whole_step < minstep)then
              print*,ixSolution
              print*, 'dtSave, dt_sub', dtSave, whole_step
              message=trim(message)//'length of the coupled step is below the minimum step length'
              err=20; return
            endif
            ! try again, restart step (at end inner step)
            cycle substeps
          endif

          ! update coordinate variables
          call calcHeight(&
                     ! input/output: data structures
                     indx_data,   & ! intent(in): layer type
                     prog_data,   & ! intent(inout): model variables for a local HRU
                     ! output: error control
                     err,cmessage)
          if(err/=0)then; err=20; message=trim(message)//trim(cmessage); return; end if

          ! recompute snow depth, SWE, and layer water
          if(nSnow > 0)then
            prog_data%var(iLookPROG%scalarSnowDepth)%dat(1) = sum( mLayerDepth(1:nSnow) )
            prog_data%var(iLookPROG%scalarSWE)%dat(1)       = sum( (mLayerVolFracLiq(1:nSnow)*iden_water &
                                                              + mLayerVolFracIce(1:nSnow)*iden_ice) * mLayerDepth(1:nSnow) )
            mLayerVolFracWat(1:nSnow) = mLayerVolFracLiq(1:nSnow) + mLayerVolFracIce(1:nSnow)*iden_ice/iden_water
          endif

        end associate sublime

        ! increment change in storage associated with the surface melt pond (kg m-2)
        if(nSnow==0) sfcMeltPond = sfcMeltPond + prog_data%var(iLookPROG%scalarSfcMeltPond)%dat(1)

      endif ! (do_outer loop)

      ! ****************************************************************************************************
      ! *** END MAIN SOLVER ********************************************************************************
      ! ****************************************************************************************************

      ! increment mean fluxes, soil compression, and effective rainfall, reset on whole_step
      dt_wght = dt_sub/whole_step ! define weight applied to each sub-step
      do iVar=1,size(averageFlux_meta)
        flux_inner%var(iVar)%dat(:)    = flux_inner%var(iVar)%dat(:) + flux_data%var(averageFlux_meta(iVar)%ixParent)%dat(:)*dt_wght
      end do
      innerSoilCompress(:) = innerSoilCompress(:) + diag_data%var(iLookDIAG%mLayerCompress)%dat(:)*dt_wght
      if (nSnow>0) innerEffRainfall = innerEffRainfall + ( flux_data%var(iLookFLUX%scalarThroughfallRain)%dat(1) + flux_data%var(iLookFLUX%scalarCanopyLiqDrainage)%dat(1) )*dt_wght

      ! increment sub-step accepted step
      dt_solvInner = dt_solvInner + dt_sub
      dt_solv = dt_solv + dt_sub

      ! update first and last inner steps if did successful lastInnerStep, increment fluxes and flux variables over data_step
      if (lastInnerStep)then
        firstInnerStep = .true.
        lastInnerStep = .false.
        dt_solvInner = 0._rkind

        dt_wght = whole_step/data_step ! define weight applied to each sub-step
        do iVar=1,size(averageFlux_meta)
          flux_mean%var(iVar)%dat(:)    = flux_mean%var(iVar)%dat(:) + flux_inner%var(iVar)%dat(:)*dt_wght
        end do
        meanCanopySublimation = meanCanopySublimation + sumCanopySublimation/data_step
        meanLatHeatCanopyEvap = meanLatHeatCanopyEvap + sumLatHeatCanopyEvap/data_step
        meanSenHeatCanopy     = meanSenHeatCanopy     + sumSenHeatCanopy/data_step
        meanSoilCompress(:) = meanSoilCompress(:) + innerSoilCompress(:)*dt_wght
        effRainfall = effRainfall + innerEffRainfall*dt_wght
        flux_mean%var(childFLUX_MEAN(iLookFLUX%scalarCanopySublimation))%dat(1) = meanCanopySublimation
        flux_mean%var(childFLUX_MEAN(iLookFLUX%scalarLatHeatCanopyEvap))%dat(1) = meanLatHeatCanopyEvap
        flux_mean%var(childFLUX_MEAN(iLookFLUX%scalarSenHeatCanopy))%dat(1)     = meanSenHeatCanopy
      endif

      ! save the time step to initialize the subsequent step
      if(dt_solv<data_step .or. nsub==1) dt_init = dt_sub

      ! check
      if(globalPrintFlag)&
      write(*,'(a,1x,3(f18.5,1x))') 'dt_sub, dt_solv, data_step: ', dt_sub, dt_solv, data_step

      ! check that we have completed the sub-step
      if(dt_solv >= data_step-verySmall) then
        exit substeps
      endif

      ! adjust length of the sub-step (make sure that we don't exceed the step)
      dt_sub = min(data_step - dt_solv, dt_sub)

    end do  substeps ! (sub-step loop)

    ! *** add snowfall to the snowpack...
    ! -----------------------------------
    ! add new snowfall to the snowpack
    ! NOTE: This needs to be done AFTER the call to canopySnow, since throughfall and unloading are computed in canopySnow
    call newsnwfall(&
                  ! input: model control
                  data_step,                                                 & ! time step (seconds)
                  (nSnow > 0),                                               & ! logical flag if snow layers exist
                  mpar_data%var(iLookPARAM%snowfrz_scale)%dat(1),            & ! freeezing curve parameter for snow (K-1)
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

    ! re-compute snow depth, SWE, and top layer water
    if(nSnow > 0)then
      prog_data%var(iLookPROG%scalarSnowDepth)%dat(1) = sum(  prog_data%var(iLookPROG%mLayerDepth)%dat(1:nSnow))
      prog_data%var(iLookPROG%scalarSWE)%dat(1)       = sum( (prog_data%var(iLookPROG%mLayerVolFracLiq)%dat(1:nSnow)*iden_water + &
                                                              prog_data%var(iLookPROG%mLayerVolFracIce)%dat(1:nSnow)*iden_ice) &
                                                            * prog_data%var(iLookPROG%mLayerDepth)%dat(1:nSnow) )
      prog_data%var(iLookPROG%mLayerVolFracWat)%dat(1) = prog_data%var(iLookPROG%mLayerVolFracLiq)%dat(1) &
                                                        + prog_data%var(iLookPROG%mLayerVolFracIce)%dat(1)*iden_ice/iden_water
    end if

    ! re-assign dimension lengths
    nSnow   = count(indx_data%var(iLookINDEX%layerType)%dat==iname_snow)
    nSoil   = count(indx_data%var(iLookINDEX%layerType)%dat==iname_soil)
    nLayers = nSnow+nSoil

    ! update coordinate variables
    call calcHeight(&
                    ! input/output: data structures
                    indx_data,   & ! intent(in): layer type
                    prog_data,   & ! intent(inout): model variables for a local HRU
                    ! output: error control
                    err,cmessage)
    if(err/=0)then; err=20; message=trim(message)//trim(cmessage); return; end if

    ! overwrite flux_data and soil compression with the timestep-average value (returns timestep-average fluxes for scalar variables)
    do iVar=1,size(averageFlux_meta)
      flux_data%var(averageFlux_meta(iVar)%ixParent)%dat(:) = flux_mean%var(iVar)%dat(:)
    end do
    ! keep soil compression as an average like the fluxes, will not want to do this if nSoil can change
    diag_data%var(iLookDIAG%mLayerCompress)%dat(:) = meanSoilCompress
    diag_data%var(iLookDIAG%scalarSoilCompress)%dat(1) = sum(  meanSoilCompress(1:nSoil)*iden_water &
                                                             * prog_data%var(iLookPROG%mLayerDepth)%dat(nSnow+1:nLayers) )
    deallocate(innerSoilCompress)
    deallocate(meanSoilCompress)

    ! ***********************************************************************************************************************************
    ! ---
    ! *** balance checks...
    ! ---------------------

    ! save the average compression and melt pond storage in the data structures
    prog_data%var(iLookPROG%scalarSfcMeltPond)%dat(1)  = sfcMeltPond

    ! associate local variables with information in the data structures
    associate(&
      ! model decisions
      ixNumericalMethod          => model_decisions(iLookDECISIONS%num_method)%iDecision                          ,&  ! choice of numerical solver
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
      ! snow fluxes
      averageSnowSublimation     => flux_mean%var(childFLUX_MEAN(iLookFLUX%scalarSnowSublimation)    )%dat(1)     ,&  ! sublimation from the snow surface (kg m-2 s-1)
      averageSnowDrainage        => flux_mean%var(childFLUX_MEAN(iLookFLUX%scalarSnowDrainage)       )%dat(1)     ,&  ! drainage from the bottom of the snowpack (m s-1)
      ! soil fluxes
      averageSoilInflux          => flux_mean%var(childFLUX_MEAN(iLookFLUX%scalarInfiltration)       )%dat(1)     ,&  ! influx of water at the top of the soil profile (m s-1)
      averageSoilDrainage        => flux_mean%var(childFLUX_MEAN(iLookFLUX%scalarSoilDrainage)       )%dat(1)     ,&  ! drainage from the bottom of the soil profile (m s-1)
      averageSoilBaseflow        => flux_mean%var(childFLUX_MEAN(iLookFLUX%scalarSoilBaseflow)       )%dat(1)     ,&  ! total baseflow from throughout the soil profile (m s-1)
      averageSoilCompress        => diag_data%var(               iLookDIAG%scalarSoilCompress)        %dat(1)     ,&  ! soil compression (kg m-2 s-1)
      averageGroundEvaporation   => flux_mean%var(childFLUX_MEAN(iLookFLUX%scalarGroundEvaporation)  )%dat(1)     ,&  ! soil evaporation (kg m-2 s-1)
      averageCanopyTranspiration => flux_mean%var(childFLUX_MEAN(iLookFLUX%scalarCanopyTranspiration))%dat(1)     ,&  ! canopy transpiration (kg m-2 s-1)
      ! state variables in the vegetation canopy
      scalarCanopyWat            => prog_data%var(iLookPROG%scalarCanopyWat)%dat(1)                               ,&  ! canopy ice content (kg m-2)
      ! state variables in the soil domain
      mLayerDepth                => prog_data%var(iLookPROG%mLayerDepth)%dat(nSnow+1:nLayers)                     ,&  ! depth of each soil layer (m)
      mLayerVolFracIce           => prog_data%var(iLookPROG%mLayerVolFracIce)%dat(nSnow+1:nLayers)                ,&  ! volumetric ice content in each soil layer (-)
      mLayerVolFracLiq           => prog_data%var(iLookPROG%mLayerVolFracLiq)%dat(nSnow+1:nLayers)                ,&  ! volumetric liquid water content in each soil layer (-)
      scalarAquiferStorage       => prog_data%var(iLookPROG%scalarAquiferStorage)%dat(1)                          ,&  ! aquifer storage (m)
      ! error tolerance
      absConvTol_liquid          => mpar_data%var(iLookPARAM%absConvTol_liquid)%dat(1)                            ,&  ! absolute convergence tolerance for vol frac liq water (-)
      scalarTotalSoilIce         => diag_data%var(iLookDIAG%scalarTotalSoilIce)%dat(1)                            ,&  ! total ice in the soil column (kg m-2)
      scalarTotalSoilLiq         => diag_data%var(iLookDIAG%scalarTotalSoilLiq)%dat(1)                             &  ! total liquid water in the soil column (kg m-2)
      ) ! (association of local variables with information in the data structures

      ! identify the need to check the mass balance, both methods should work if tolerance coarse enough
      select case(ixNumericalMethod)
        case(ida);            checkMassBalance = .false. ! IDA balance agreement levels are controlled by set tolerances
        case(kinsol, numrec); checkMassBalance = .true.  ! KINSOL or numrec give finite difference dt_sub fluxes and were summed for an average flux
        case default; err=20; message=trim(message)//'expect num_method to be ida, kinsol, or numrec (or itertive, which is numrec)'; return
      end select

      ! -----
      ! * balance checks for the canopy...
      ! ----------------------------------

      ! if computing the vegetation flux
      if(computeVegFlux)then
        ! get the canopy water balance at the end of the time step
        balanceCanopyWater1 = scalarCanopyWat

        ! balance checks for the canopy
        ! NOTE: need to put the balance checks in the sub-step loop so that we can re-compute if necessary
        scalarCanopyWatBalError = balanceCanopyWater1 - (balanceCanopyWater0 + (scalarSnowfall - averageThroughfallSnow)*data_step + (scalarRainfall - averageThroughfallRain)*data_step &
                                  - averageCanopySnowUnloading*data_step - averageCanopyLiqDrainage*data_step + averageCanopySublimation*data_step + averageCanopyEvaporation*data_step)
        if(abs(scalarCanopyWatBalError) > absConvTol_liquid*iden_water*10._rkind .and. checkMassBalance)then
          write(*,'(a,1x,f20.10)') 'data_step                                    = ', data_step
          write(*,'(a,1x,f20.10)') 'balanceCanopyWater0                          = ', balanceCanopyWater0
          write(*,'(a,1x,f20.10)') 'balanceCanopyWater1                          = ', balanceCanopyWater1
          write(*,'(a,1x,f20.10)') 'scalarSnowfall                               = ', scalarSnowfall
          write(*,'(a,1x,f20.10)') 'scalarRainfall                               = ', scalarRainfall
          write(*,'(a,1x,f20.10)') '(scalarSnowfall - averageThroughfallSnow)    = ', (scalarSnowfall - averageThroughfallSnow)!*data_step
          write(*,'(a,1x,f20.10)') '(scalarRainfall - averageThroughfallRain)    = ', (scalarRainfall - averageThroughfallRain)!*data_step
          write(*,'(a,1x,f20.10)') 'averageCanopySnowUnloading                   = ', averageCanopySnowUnloading!*data_step
          write(*,'(a,1x,f20.10)') 'averageCanopyLiqDrainage                     = ', averageCanopyLiqDrainage!*data_step
          write(*,'(a,1x,f20.10)') 'averageCanopySublimation                     = ', averageCanopySublimation!*data_step
          write(*,'(a,1x,f20.10)') 'averageCanopyEvaporation                     = ', averageCanopyEvaporation!*data_step
          write(*,'(a,1x,f20.10)') 'scalarCanopyWatBalError                      = ', scalarCanopyWatBalError
          message=trim(message)//'canopy hydrology does not balance'
          err=20; return
        end if

      endif  ! if computing the vegetation flux

      ! -----
      ! * balance checks for SWE...
      ! ---------------------------

      ! check the individual layers
      if(printBalance .and. nSnow>0)then
        write(*,'(a,1x,10(f12.8,1x))') 'liqSnowInit       = ', liqSnowInit
        write(*,'(a,1x,10(f12.8,1x))') 'volFracLiq        = ', prog_data%var(iLookPROG%mLayerVolFracLiq)%dat(1:nSnow)
        write(*,'(a,1x,10(f12.8,1x))') 'iLayerLiqFluxSnow = ', flux_data%var(iLookFLUX%iLayerLiqFluxSnow)%dat*iden_water*data_step
        write(*,'(a,1x,10(f12.8,1x))') 'mLayerLiqFluxSnow = ', flux_data%var(iLookFLUX%mLayerLiqFluxSnow)%dat*data_step
        write(*,'(a,1x,10(f12.8,1x))') 'change volFracLiq = ', prog_data%var(iLookPROG%mLayerVolFracLiq)%dat(1:nSnow) - liqSnowInit
        deallocate(liqSnowInit, stat=err)
        if(err/=0)then
          message=trim(message)//'unable to deallocate space for the initial volumetric liquid water content of snow'
          err=20; return
        endif
      endif

      ! check SWE
      if(nSnow>0)then
        effSnowfall = averageThroughfallSnow + averageCanopySnowUnloading
        ! effRainfall is averageThroughfallRain + averageCanopyLiqDrainage only over snow
        newSWE      = prog_data%var(iLookPROG%scalarSWE)%dat(1)
        delSWE      = newSWE - (oldSWE - sfcMeltPond)
        massBalance = delSWE - (effSnowfall + effRainfall + averageSnowSublimation - averageSnowDrainage*iden_water)*data_step
        if(abs(massBalance) > absConvTol_liquid*iden_water*10._rkind .and. checkMassBalance)then
          print*,                  'nSnow       = ', nSnow
          print*,                  'nSub        = ', nSub
          write(*,'(a,1x,f20.10)') 'data_step   = ', data_step
          write(*,'(a,1x,f20.10)') 'oldSWE      = ', oldSWE
          write(*,'(a,1x,f20.10)') 'newSWE      = ', newSWE
          write(*,'(a,1x,f20.10)') 'delSWE      = ', delSWE
          write(*,'(a,1x,f20.10)') 'effRainfall = ', effRainfall*data_step
          write(*,'(a,1x,f20.10)') 'effSnowfall = ', effSnowfall*data_step
          write(*,'(a,1x,f20.10)') 'sublimation = ', averageSnowSublimation*data_step
          write(*,'(a,1x,f20.10)') 'snwDrainage = ', averageSnowDrainage*iden_water*data_step
          write(*,'(a,1x,f20.10)') 'sfcMeltPond = ', sfcMeltPond
          write(*,'(a,1x,f20.10)') 'massBalance = ', massBalance
          message=trim(message)//'SWE does not balance'
          err=20; return
        endif  ! if failed mass balance check
      endif  ! if snow layers exist

      ! -----
      ! * balance checks for soil...
      ! ----------------------------

      ! compute the liquid water and ice content at the end of the time step
      scalarTotalSoilLiq = sum(iden_water*mLayerVolFracLiq(1:nSoil)*mLayerDepth(1:nSoil))
      scalarTotalSoilIce = sum(iden_water*mLayerVolFracIce(1:nSoil)*mLayerDepth(1:nSoil))   ! NOTE: no expansion of soil, hence use iden_water

      ! get the total water in the soil (liquid plus ice) at the end of the time step (kg m-2)
      balanceSoilWater1 = scalarTotalSoilLiq + scalarTotalSoilIce

      ! get the total aquifer storage at the start of the time step (kg m-2)
      balanceAquifer1 = scalarAquiferStorage*iden_water

      ! get the input and output to/from the soil zone (kg m-2)
      balanceSoilInflux        = averageSoilInflux*iden_water*data_step
      balanceSoilBaseflow      = averageSoilBaseflow*iden_water*data_step
      balanceSoilDrainage      = averageSoilDrainage*iden_water*data_step
      balanceSoilET            = (averageCanopyTranspiration + averageGroundEvaporation)*data_step
      balanceSoilCompress      = averageSoilCompress*data_step

      ! check the individual layers
      if(printBalance)then
        write(*,'(a,1x,10(f12.8,1x))') 'liqSoilInit       = ', liqSoilInit
        write(*,'(a,1x,10(f12.8,1x))') 'volFracLiq        = ', mLayerVolFracLiq
        write(*,'(a,1x,10(f12.8,1x))') 'iLayerLiqFluxSoil = ', flux_data%var(iLookFLUX%iLayerLiqFluxSoil)%dat*iden_water*data_step
        write(*,'(a,1x,10(f12.8,1x))') 'mLayerLiqFluxSoil = ', flux_data%var(iLookFLUX%mLayerLiqFluxSoil)%dat*data_step
        write(*,'(a,1x,10(f12.8,1x))') 'change volFracLiq = ', mLayerVolFracLiq - liqSoilInit
        deallocate(liqSoilInit, stat=err)
        if(err/=0)then
        message=trim(message)//'unable to deallocate space for the initial soil moisture'
        err=20; return
        endif
      endif

      ! check the soil water balance
      scalarSoilWatBalError  = balanceSoilWater1 - (balanceSoilWater0 + (balanceSoilInflux + balanceSoilET - balanceSoilBaseflow - balanceSoilDrainage - balanceSoilCompress) )
      if(abs(scalarSoilWatBalError) > absConvTol_liquid*iden_water*10._rkind .and. checkMassBalance)then  ! NOTE: kg m-2, so need coarse tolerance to account for precision issues
        write(*,*)               'solution method           = ', ixSolution
        write(*,'(a,1x,f20.10)') 'data_step                 = ', data_step
        write(*,'(a,1x,f20.10)') 'balanceSoilCompress       = ', balanceSoilCompress
        write(*,'(a,1x,f20.10)') 'scalarTotalSoilLiq        = ', scalarTotalSoilLiq
        write(*,'(a,1x,f20.10)') 'scalarTotalSoilIce        = ', scalarTotalSoilIce
        write(*,'(a,1x,f20.10)') 'balanceSoilWater0         = ', balanceSoilWater0
        write(*,'(a,1x,f20.10)') 'balanceSoilWater1         = ', balanceSoilWater1
        write(*,'(a,1x,f20.10)') 'balanceSoilInflux         = ', balanceSoilInflux
        write(*,'(a,1x,f20.10)') 'balanceSoilBaseflow       = ', balanceSoilBaseflow
        write(*,'(a,1x,f20.10)') 'balanceSoilDrainage       = ', balanceSoilDrainage
        write(*,'(a,1x,f20.10)') 'balanceSoilET             = ', balanceSoilET
        write(*,'(a,1x,f20.10)') 'scalarSoilWatBalError     = ', scalarSoilWatBalError
        message=trim(message)//'soil hydrology does not balance'
        err=20; return
      end if

    ! end association of local variables with information in the data structures
    end associate

  ! end association to canopy depth
  end associate canopy

  ! Save the total soil water (Liquid+Ice)
  diag_data%var(iLookDIAG%scalarTotalSoilWat)%dat(1) = balanceSoilWater1
  ! save the surface temperature (just to make things easier to visualize)
  prog_data%var(iLookPROG%scalarSurfaceTemp)%dat(1) = prog_data%var(iLookPROG%mLayerTemp)%dat(1)

  ! overwrite flux data with timestep-average value for all flux_mean vars, hard-coded to not happen
  if(.not.backwardsCompatibility)then
    do iVar=1,size(flux_mean%var)
      flux_data%var(averageFlux_meta(iVar)%ixParent)%dat = flux_mean%var(iVar)%dat
    end do
  end if

  iLayer = nSnow+1
  if(nsub>50000)then
    write(message,'(a,i0)') trim(cmessage)//'number of sub-steps > 50000 for HRU ', hruId
    err=20; return
  end if

  ! get the end time
  call cpu_time(endTime)

  ! get the elapsed time
  diag_data%var(iLookDIAG%wallClockTime)%dat(1) = endTime - startTime

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
  real(rkind),intent(inout)    :: scalarSWE          ! snow water equivalent (kg m-2)
  real(rkind),intent(inout)    :: scalarSnowDepth    ! snow depth (m)
  real(rkind),intent(inout)    :: scalarSfcMeltPond  ! surface melt pond (kg m-2)
  ! input/output: properties of the upper-most soil layer
  real(rkind),intent(inout)    :: soilTemp           ! surface layer temperature (K)
  real(rkind),intent(inout)    :: soilDepth          ! surface layer depth (m)
  real(rkind),intent(inout)    :: soilHeatcap        ! surface layer volumetric heat capacity (J m-3 K-1)
  ! output: error control
  integer(i4b),intent(out)  :: err                ! error code
  character(*),intent(out)  :: message            ! error message
  ! local variables
  real(rkind)                  :: nrgRequired        ! energy required to melt all the snow (J m-2)
  real(rkind)                  :: nrgAvailable       ! energy available to melt the snow (J m-2)
  real(rkind)                  :: snwDensity         ! snow density (kg m-3)
  ! initialize error control
  err=0; message='implctMelt/'

  if(scalarSWE > 0._rkind)then
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
        scalarSWE          = 0._rkind
      else
        scalarSfcMeltPond  = nrgAvailable/LH_fus
        scalarSWE          = scalarSWE - scalarSfcMeltPond
      end if
      ! update depth
      scalarSnowDepth = scalarSWE/snwDensity
      ! update temperature of the top soil layer (K)
      soilTemp =  soilTemp - (LH_fus*scalarSfcMeltPond/soilDepth)/soilHeatcap
    else  ! melt is zero if the temperature of the top soil layer is less than Tfreeze
      scalarSfcMeltPond = 0._rkind  ! kg m-2
    end if ! (if the temperature of the top soil layer is greater than Tfreeze)
  else  ! melt is zero if the "snow without a layer" does not exist
    scalarSfcMeltPond = 0._rkind  ! kg m-2
  end if ! (if the "snow without a layer" exists)

end subroutine implctMelt

end module coupled_em_module

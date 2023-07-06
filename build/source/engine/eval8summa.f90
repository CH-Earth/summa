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

module eval8summa_module

! data types
USE nrtype

! access missing values
USE globalData,only:integerMissing  ! missing integer
USE globalData,only:realMissing     ! missing double precision number
USE globalData,only:quadMissing     ! missing quadruple precision number


! constants
USE multiconst,only:&
                    Tfreeze,      & ! temperature at freezing              (K)
                    iden_ice,     & ! intrinsic density of ice             (kg m-3)
                    iden_water      ! intrinsic density of liquid water    (kg m-3)

! provide access to the derived types to define the data structures
USE data_types,only:&
                    var_i,        & ! data vector (i4b)
                    var_d,        & ! data vector (rkind)
                    var_ilength,  & ! data vector with variable length dimension (i4b)
                    var_dlength,  & ! data vector with variable length dimension (rkind)
                    zLookup,      & ! data vector with variable length dimension (rkind)
                    model_options   ! defines the model decisions

! indices that define elements of the data structures
USE var_lookup,only:iLookDECISIONS               ! named variables for elements of the decision structure
USE var_lookup,only:iLookPARAM                   ! named variables for structure elements
USE var_lookup,only:iLookPROG                    ! named variables for structure elements
USE var_lookup,only:iLookINDEX                   ! named variables for structure elements
USE var_lookup,only:iLookDIAG                    ! named variables for structure elements
USE var_lookup,only:iLookFLUX                    ! named variables for structure elements
USE var_lookup,only:iLookDERIV                   ! named variables for structure elements

! look-up values for the choice of heat capacity computation
USE mDecisions_module,only:  &
 closedForm,                 & ! heat capacity using closed form, not using enthalpy
 enthalpyFD                    ! heat capacity using enthalpy

! look-up values for the numerical method
USE mDecisions_module,only:  &
 numrec,                     & ! home-grown backward Euler solution using free versions of Numerical recipes
 kinsol,                     & ! SUNDIALS backward Euler solution using Kinsol
 ida                           ! SUNDIALS solution using IDA

implicit none
private
public::eval8summa
#ifdef SUNDIALS_ACTIVE
public::eval8summa4kinsol
#endif

contains


! **********************************************************************************************************
! public subroutine eval8summa: compute the residual vector and the Jacobian matrix
! **********************************************************************************************************
subroutine eval8summa(&
                      ! input: model control
                      dt_cur,                  & ! intent(in):    current stepsize
                      dt,                      & ! intent(in):    entire time step for drainage pond rate
                      nSnow,                   & ! intent(in):    number of snow layers
                      nSoil,                   & ! intent(in):    number of soil layers
                      nLayers,                 & ! intent(in):    total number of layers
                      nState,                  & ! intent(in):    total number of state variables
                      insideSUN,               & ! intent(in):    flag to indicate if we are inside Sundials solver
                      firstSubStep,            & ! intent(in):    flag to indicate if we are processing the first sub-step
                      firstFluxCall,           & ! intent(inout): flag to indicate if we are processing the first flux call
                      firstSplitOper,          & ! intent(in):    flag to indicate if we are processing the first flux call in a splitting operation
                      computeVegFlux,          & ! intent(in):    flag to indicate if we need to compute fluxes over vegetation
                      scalarSolution,          & ! intent(in):    flag to indicate the scalar solution
                      ! input: state vectors
                      stateVec,                & ! intent(in):    model state vector
                      fScale,                  & ! intent(in):    characteristic scale of the function evaluations
                      sMul,                    & ! intent(inout): state vector multiplier (used in the residual calculations)
                      ! input: data structures
                      model_decisions,         & ! intent(in):    model decisions
                      lookup_data,             & ! intent(in):    lookup tables
                      type_data,               & ! intent(in):    type of vegetation and soil
                      attr_data,               & ! intent(in):    spatial attributes
                      mpar_data,               & ! intent(in):    model parameters
                      forc_data,               & ! intent(in):    model forcing data
                      bvar_data,               & ! intent(in):    average model variables for the entire basin
                      prog_data,               & ! intent(in):    model prognostic variables for a local HRU
                      ! input-output: data structures
                      indx_data,               & ! intent(inout): index data
                      diag_data,               & ! intent(inout): model diagnostic variables for a local HRU
                      flux_data,               & ! intent(inout): model fluxes for a local HRU
                      deriv_data,              & ! intent(inout): derivatives in model fluxes w.r.t. relevant state variables
                      ! input-output: baseflow
                      ixSaturation,            & ! intent(inout): index of the lowest saturated layer (NOTE: only computed on the first iteration)
                      dBaseflow_dMatric,       & ! intent(out):   derivative in baseflow w.r.t. matric head (s-1)
                      ! output: flux and residual vectors
                      feasible,                & ! intent(out):   flag to denote the feasibility of the solution
                      fluxVec,                 & ! intent(out):   flux vector
                      resSink,                 & ! intent(out):   additional (sink) terms on the RHS of the state equation
                      resVec,                  & ! intent(out):   residual vector
                      fEval,                   & ! intent(out):   function evaluation
                      err,message)               ! intent(out):   error control
  ! --------------------------------------------------------------------------------------------------------------------------------
  ! provide access to subroutines
  USE getVectorz_module, only:varExtract                ! extract variables from the state vector
  USE getVectorz_module, only:checkFeas                 ! check feasibility of state vector
  USE updateVars_module, only:updateVars                ! update prognostic variables
  USE t2enthalpy_module, only:t2enthalpy                ! compute enthalpy
  USE computFlux_module, only:soilCmpres                ! compute soil compression
  USE computFlux_module, only:computFlux                ! compute fluxes given a state vector
  USE computHeatCap_module,only:computHeatCap           ! recompute heat capacity and derivatives
  USE computHeatCap_module,only:computHeatCapAnalytic   ! recompute heat capacity and derivatives
  USE computHeatCap_module,only:computCm
  USE computHeatCap_module, only:computStatMult         ! recompute state multiplier
  USE computResid_module,only:computResid               ! compute residuals given a state vector
  USE computThermConduct_module,only:computThermConduct ! recompute thermal conductivity and derivatives
  implicit none
  ! --------------------------------------------------------------------------------------------------------------------------------
  ! --------------------------------------------------------------------------------------------------------------------------------
  ! input: model control
  real(rkind),intent(in)          :: dt_cur                 ! current stepsize
  real(rkind),intent(in)          :: dt                     ! entire time step for drainage pond rate
  integer(i4b),intent(in)         :: nSnow                  ! number of snow layers
  integer(i4b),intent(in)         :: nSoil                  ! number of soil layers
  integer(i4b),intent(in)         :: nLayers                ! total number of layers
  integer(i4b),intent(in)         :: nState                 ! total number of state variables
  logical(lgt),intent(in)         :: insideSUN              ! flag to indicate if we are inside Sundials solver
  logical(lgt),intent(in)         :: firstSubStep           ! flag to indicate if we are processing the first sub-step
  logical(lgt),intent(inout)      :: firstFluxCall          ! flag to indicate if we are processing the first flux call
  logical(lgt),intent(in)         :: firstSplitOper         ! flag to indicate if we are processing the first flux call in a splitting operation
  logical(lgt),intent(in)         :: computeVegFlux         ! flag to indicate if computing fluxes over vegetation
  logical(lgt),intent(in)         :: scalarSolution         ! flag to denote if implementing the scalar solution
  ! input: state vectors
  real(rkind),intent(in)          :: stateVec(:)            ! model state vector
  real(rkind),intent(in)          :: fScale(:)              ! characteristic scale of the function evaluations
  real(qp),intent(inout)          :: sMul(:)   ! NOTE: qp   ! state vector multiplier (used in the residual calculations)
  ! input: data structures
  type(model_options),intent(in)  :: model_decisions(:)     ! model decisions
  type(zLookup),      intent(in)  :: lookup_data            ! lookup tables
  type(var_i),        intent(in)  :: type_data              ! type of vegetation and soil
  type(var_d),        intent(in)  :: attr_data              ! spatial attributes
  type(var_dlength),  intent(in)  :: mpar_data              ! model parameters
  type(var_d),        intent(in)  :: forc_data              ! model forcing data
  type(var_dlength),  intent(in)  :: bvar_data              ! model variables for the local basin
  type(var_dlength),  intent(in)  :: prog_data              ! prognostic variables for a local HRU
  ! output: data structures
  type(var_ilength),intent(inout) :: indx_data              ! indices defining model states and layers
  type(var_dlength),intent(inout) :: diag_data              ! diagnostic variables for a local HRU
  type(var_dlength),intent(inout) :: flux_data              ! model fluxes for a local HRU
  type(var_dlength),intent(inout) :: deriv_data             ! derivatives in model fluxes w.r.t. relevant state variables
  ! input-output: baseflow
  integer(i4b),intent(inout)      :: ixSaturation           ! index of the lowest saturated layer (NOTE: only computed on the first iteration)
  real(rkind),intent(out)         :: dBaseflow_dMatric(:,:) ! derivative in baseflow w.r.t. matric head (s-1)
  ! output: flux and residual vectors
  logical(lgt),intent(out)        :: feasible               ! flag to denote the feasibility of the solution
  real(rkind),intent(out)         :: fluxVec(:)             ! flux vector
  real(rkind),intent(out)         :: resSink(:)             ! sink terms on the RHS of the flux equation
  real(qp),intent(out)            :: resVec(:) ! NOTE: qp   ! residual vector
  real(rkind),intent(out)         :: fEval                  ! function evaluation
  ! output: error control
  integer(i4b),intent(out)        :: err                    ! error code
  character(*),intent(out)        :: message                ! error message
  ! --------------------------------------------------------------------------------------------------------------------------------
  ! local variables
  ! --------------------------------------------------------------------------------------------------------------------------------
  ! state variables
  real(rkind)                        :: scalarCanairTempTrial     ! trial value for temperature of the canopy air space (K)
  real(rkind)                        :: scalarCanopyTempTrial     ! trial value for temperature of the vegetation canopy (K)
  real(rkind)                        :: scalarCanopyWatTrial      ! trial value for liquid water storage in the canopy (kg m-2)
  real(rkind),dimension(nLayers)     :: mLayerTempTrial           ! trial value for temperature of layers in the snow and soil domains (K)
  real(rkind),dimension(nLayers)     :: mLayerVolFracWatTrial     ! trial value for volumetric fraction of total water (-)
  real(rkind),dimension(nSoil)       :: mLayerMatricHeadTrial     ! trial value for total water matric potential (m)
  real(rkind),dimension(nSoil)       :: mLayerMatricHeadLiqTrial  ! trial value for liquid water matric potential (m)
  real(rkind)                        :: scalarAquiferStorageTrial ! trial value of storage of water in the aquifer (m)
  ! diagnostic variables
  real(rkind)                        :: scalarCanopyLiqTrial      ! trial value for mass of liquid water on the vegetation canopy (kg m-2)
  real(rkind)                        :: scalarCanopyIceTrial      ! trial value for mass of ice on the vegetation canopy (kg m-2)
  real(rkind),dimension(nLayers)     :: mLayerVolFracLiqTrial     ! trial value for volumetric fraction of liquid water (-)
  real(rkind),dimension(nLayers)     :: mLayerVolFracIceTrial     ! trial value for volumetric fraction of ice (-)
  ! enthalpy
  real(rkind)                        :: scalarCanopyEnthalpyTrial ! trial value for enthalpy of the vegetation canopy (J m-3)
  real(rkind),dimension(nLayers)     :: mLayerEnthalpyTrial       ! trial vector of enthalpy for snow+soil layers (J m-3)
  real(rkind)                        :: dCanEnthalpy_dTk          ! derivatives in canopy enthalpy w.r.t. temperature
  real(rkind)                        :: dCanEnthalpy_dWat         ! derivatives in canopy enthalpy w.r.t. water state
  real(rkind),dimension(nLayers)     :: dEnthalpy_dTk             ! derivatives in layer enthalpy w.r.t. temperature
  real(rkind),dimension(nLayers)     :: dEnthalpy_dWat            ! derivatives in layer enthalpy w.r.t. water state
  ! other local variables
  logical(lgt)                       :: checkLWBalance            ! flag to check longwave balance
  integer(i4b)                       :: iLayer                    ! index of model layer in the snow+soil domain
  integer(i4b)                       :: jState(1)                 ! index of model state for the scalar solution within the soil domain
  integer(i4b)                       :: ixBeg,ixEnd               ! index of indices for the soil compression routine
  real(rkind),dimension(nState)      :: rVecScaled                ! scaled residual vector
  character(LEN=256)                 :: cmessage                  ! error message of downwind routine
  real(rkind)                        :: scalarCanopyCmTrial       ! trial value of Cm for the canopy
  real(rkind),dimension(nLayers)     :: mLayerCmTrial             ! trial vector of Cm for snow+soil
  logical(lgt),parameter             :: updateCp=.false.          ! flag to indicate if we update Cp at each step
  logical(lgt),parameter             :: needCm=.false.            ! flag to indicate if the energy equation contains Cm = dH_T/dTheta_m

  ! --------------------------------------------------------------------------------------------------------------------------------
  ! association to variables in the data structures
  ! --------------------------------------------------------------------------------------------------------------------------------
  associate(&
    ! model decisions
    ixNumericalMethod       => model_decisions(iLookDECISIONS%num_method)%iDecision   ,& ! intent(in):  [i4b]   choice of numerical solver
    ixHowHeatCap            => model_decisions(iLookDECISIONS%howHeatCap)%iDecision   ,& ! intent(in):  [i4b]   heat capacity computation, with or without enthalpy
    ixRichards              => model_decisions(iLookDECISIONS%f_Richards)%iDecision   ,& ! intent(in):  [i4b]   index of the form of Richards' equation
    ! snow parameters
    snowfrz_scale           => mpar_data%var(iLookPARAM%snowfrz_scale)%dat(1)         ,& ! intent(in):  [dp]    scaling parameter for the snow freezing curve (K-1)
    ! soil parameters
    theta_sat               => mpar_data%var(iLookPARAM%theta_sat)%dat                ,& ! intent(in):  [dp(:)] soil porosity (-)
    specificStorage         => mpar_data%var(iLookPARAM%specificStorage)%dat(1)       ,& ! intent(in):  [dp]    specific storage coefficient (m-1)
    ! canopy and layer depth
    canopyDepth             => diag_data%var(iLookDIAG%scalarCanopyDepth)%dat(1)      ,& ! intent(in):  [dp   ] canopy depth (m)
    mLayerDepth             => prog_data%var(iLookPROG%mLayerDepth)%dat               ,& ! intent(in):  [dp(:)] depth of each layer in the snow-soil sub-domain (m)
    ! model state variables from the previous solution
    scalarCanairTemp        => prog_data%var(iLookPROG%scalarCanairTemp)%dat(1)       ,& ! intent(in):  [dp]    temperature of the canopy air space (K)
    scalarCanopyTemp        => prog_data%var(iLookPROG%scalarCanopyTemp)%dat(1)       ,& ! intent(in):  [dp]    temperature of the vegetation canopy (K)
    scalarCanopyWat         => prog_data%var(iLookPROG%scalarCanopyWat)%dat(1)        ,& ! intent(in):  [dp]    mass of total water on the vegetation canopy (kg m-2)
    mLayerTemp              => prog_data%var(iLookPROG%mLayerTemp)%dat                ,& ! intent(in):  [dp(:)] temperature of each snow/soil layer (K)
    mLayerVolFracWat        => prog_data%var(iLookPROG%mLayerVolFracWat)%dat          ,& ! intent(in):  [dp(:)] volumetric fraction of total water (-)
    mLayerMatricHead        => prog_data%var(iLookPROG%mLayerMatricHead)%dat          ,& ! intent(in):  [dp(:)] total water matric potential (m)
    mLayerMatricHeadLiq     => diag_data%var(iLookDIAG%mLayerMatricHeadLiq)%dat       ,& ! intent(in):  [dp(:)] liquid water matric potential (m)
    scalarAquiferStorage    => prog_data%var(iLookPROG%scalarAquiferStorage)%dat(1)   ,& ! intent(in):  [dp]    storage of water in the aquifer (m)
    ! model diagnostic variables from the previous solution
    scalarCanopyLiq         => prog_data%var(iLookPROG%scalarCanopyLiq)%dat(1)        ,& ! intent(in):  [dp(:)] mass of liquid water on the vegetation canopy (kg m-2)
    scalarCanopyIce         => prog_data%var(iLookPROG%scalarCanopyIce)%dat(1)        ,& ! intent(in):  [dp(:)] mass of ice on the vegetation canopy (kg m-2)
    scalarFracLiqVeg        => diag_data%var(iLookDIAG%scalarFracLiqVeg)%dat(1)       ,& ! intent(in):  [dp]    fraction of liquid water on vegetation (-)
    scalarSfcMeltPond       => prog_data%var(iLookPROG%scalarSfcMeltPond)%dat(1)      ,& ! intent(in):  [dp]    ponded water caused by melt of the "snow without a layer" (kg m-2)
    mLayerVolFracLiq        => prog_data%var(iLookPROG%mLayerVolFracLiq)%dat          ,& ! intent(in):  [dp(:)] volumetric fraction of liquid water (-)
    mLayerVolFracIce        => prog_data%var(iLookPROG%mLayerVolFracIce)%dat          ,& ! intent(in):  [dp(:)] volumetric fraction of ice (-)
    mLayerFracLiqSnow       => diag_data%var(iLookDIAG%mLayerFracLiqSnow)%dat         ,& ! intent(in):  [dp(:)] fraction of liquid water in each snow layer (-)
    ! enthalpy from the previous solution
    scalarCanairEnthalpy    => diag_data%var(iLookDIAG%scalarCanairEnthalpy)%dat(1)   ,& ! intent(in):  [dp]    enthalpy of the canopy air space (J m-3)
    scalarCanopyEnthalpy    => diag_data%var(iLookDIAG%scalarCanopyEnthalpy)%dat(1)   ,& ! intent(in):  [dp]    enthalpy of the vegetation canopy (J m-3)
    mLayerEnthalpy          => diag_data%var(iLookDIAG%mLayerEnthalpy)%dat            ,& ! intent(in):  [dp(:)] enthalpy of the snow+soil layers (J m-3)
    ! soil compression
    scalarSoilCompress      => diag_data%var(iLookDIAG%scalarSoilCompress)%dat(1)     ,& ! intent(in):  [dp]    total change in storage associated with compression of the soil matrix (kg m-2 s-1)
    mLayerCompress          => diag_data%var(iLookDIAG%mLayerCompress)%dat            ,& ! intent(in):  [dp(:)] change in volumetric water content due to compression of soil (s-1)
    ! derivatives
    dTheta_dTkCanopy        => deriv_data%var(iLookDERIV%dTheta_dTkCanopy)%dat(1)     ,& ! intent(in):  [dp]    derivative of volumetric liquid water content w.r.t. temperature
    dVolTot_dPsi0           => deriv_data%var(iLookDERIV%dVolTot_dPsi0)%dat           ,& ! intent(in):  [dp(:)] derivative in total water content w.r.t. total water matric potential
    dCompress_dPsi          => deriv_data%var(iLookDERIV%dCompress_dPsi)%dat          ,& ! intent(in):  [dp(:)] derivative in compressibility w.r.t. matric head (m-1)
    mLayerdTheta_dTk        => deriv_data%var(iLookDERIV%mLayerdTheta_dTk)%dat        ,& ! intent(in):  [dp(:)] derivative of volumetric liquid water content w.r.t. temperature
    dVolHtCapBulk_dPsi0     => deriv_data%var(iLookDERIV%dVolHtCapBulk_dPsi0)%dat     ,& ! intent(out): [dp(:)] derivative in bulk heat capacity w.r.t. matric potential
    dVolHtCapBulk_dTheta    => deriv_data%var(iLookDERIV%dVolHtCapBulk_dTheta)%dat    ,& ! intent(out): [dp(:)] derivative in bulk heat capacity w.r.t. volumetric water content
    dVolHtCapBulk_dCanWat   => deriv_data%var(iLookDERIV%dVolHtCapBulk_dCanWat)%dat(1),& ! intent(out): [dp]    derivative in bulk heat capacity w.r.t. volumetric water content
    dVolHtCapBulk_dTk       => deriv_data%var(iLookDERIV%dVolHtCapBulk_dTk)%dat       ,& ! intent(out): [dp(:)] derivative in bulk heat capacity w.r.t. temperature
    dVolHtCapBulk_dTkCanopy => deriv_data%var(iLookDERIV%dVolHtCapBulk_dTkCanopy)%dat(1),&!intent(out): [dp]    derivative in bulk heat capacity w.r.t. temperature
    dThermalC_dWatAbove     => deriv_data%var(iLookDERIV%dThermalC_dWatAbove)%dat     ,& ! intent(out): [dp(:)] derivative in the thermal conductivity w.r.t. water state in the layer above
    dThermalC_dWatBelow     => deriv_data%var(iLookDERIV%dThermalC_dWatBelow)%dat     ,& ! intent(out): [dp(:)] derivative in the thermal conductivity w.r.t. water state in the layer above
    dThermalC_dTempAbove    => deriv_data%var(iLookDERIV%dThermalC_dTempAbove)%dat    ,& ! intent(out): [dp(:)] derivative in the thermal conductivity w.r.t. energy state in the layer above
    dThermalC_dTempBelow    => deriv_data%var(iLookDERIV%dThermalC_dTempBelow)%dat    ,& ! intent(out): [dp(:)] derivative in the thermal conductivity w.r.t. energy state in the layer above
    ! mapping
    ixMapFull2Subset        => indx_data%var(iLookINDEX%ixMapFull2Subset)%dat         ,& ! intent(in):  [i4b(:)]mapping of full state vector to the state subset
    ixControlVolume         => indx_data%var(iLookINDEX%ixControlVolume)%dat          ,& ! intent(in):  [i4b(:)]index of control volume for different domains (veg, snow, soil)
    ! heat capacity
    heatCapVegTrial         =>  diag_data%var(iLookDIAG%scalarBulkVolHeatCapVeg)%dat(1),&! intent(out): [dp]    volumetric heat capacity of vegetation canopy
    mLayerHeatCapTrial      =>  diag_data%var(iLookDIAG%mLayerVolHtCapBulk)%dat        & ! intent(out): [dp(:)] heat capacity for snow and soil
    ) ! association to variables in the data structures
    ! --------------------------------------------------------------------------------------------------------------------------------
    ! initialize error control
    err=0; message="eval8summa/"

    ! check the feasibility of the solution always with BE numrec but not inside Sundials solver
    feasible=.true.
    if (.not.insideSUN) then
      call checkFeas(&
                    ! input
                    stateVec,                                  & ! intent(in):    model state vector (mixed units)
                    mpar_data,                                 & ! intent(in):    model parameters
                    prog_data,                                 & ! intent(in):    model prognostic variables for a local HRU
                    indx_data,                                 & ! intent(in):    indices defining model states and layers
                    ! output: feasibility
                    feasible,                                  & ! intent(inout): flag to denote the feasibility of the solution
                  ! output: error control
                    err,cmessage)                                ! intent(out):   error control

      ! early return for non-feasible solutions
      if(.not.feasible)then
        fluxVec(:) = realMissing
        resVec(:)  = quadMissing
        fEval      = realMissing
        message=trim(message)//trim(cmessage)//'non-feasible'
        return
      end if
    end if ! ( feasibility check )

    ! get the start and end indices for the soil compression calculations
    if(scalarSolution)then
      jState = pack(ixControlVolume, ixMapFull2Subset/=integerMissing)
      ixBeg  = jState(1)
      ixEnd  = jState(1)
    else
      ixBeg  = 1
      ixEnd  = nSoil
    endif

    ! initialize to state variable from the last update
    scalarCanairTempTrial     = scalarCanairTemp
    scalarCanopyTempTrial     = scalarCanopyTemp
    scalarCanopyWatTrial      = scalarCanopyWat
    scalarCanopyLiqTrial      = scalarCanopyLiq
    scalarCanopyIceTrial      = scalarCanopyIce
    mLayerTempTrial           = mLayerTemp
    mLayerVolFracWatTrial     = mLayerVolFracWat
    mLayerVolFracLiqTrial     = mLayerVolFracLiq
    mLayerVolFracIceTrial     = mLayerVolFracIce
    mLayerMatricHeadTrial     = mLayerMatricHead
    mLayerMatricHeadLiqTrial  = mLayerMatricHeadLiq
    scalarAquiferStorageTrial = scalarAquiferStorage

    ! extract variables from the model state vector
    call varExtract(&
                    ! input
                    stateVec,                 & ! intent(in):    model state vector (mixed units)
                    diag_data,                & ! intent(in):    model diagnostic variables for a local HRU
                    prog_data,                & ! intent(in):    model prognostic variables for a local HRU
                    indx_data,                & ! intent(in):    indices defining model states and layers
                    ! output: variables for the vegetation canopy
                    scalarCanairTempTrial,    & ! intent(inout): trial value of canopy air temperature (K)
                    scalarCanopyTempTrial,    & ! intent(inout): trial value of canopy temperature (K)
                    scalarCanopyWatTrial,     & ! intent(inout): trial value of canopy total water (kg m-2)
                    scalarCanopyLiqTrial,     & ! intent(inout): trial value of canopy liquid water (kg m-2)
                    ! output: variables for the snow-soil domain
                    mLayerTempTrial,          & ! intent(inout): trial vector of layer temperature (K)
                    mLayerVolFracWatTrial,    & ! intent(inout): trial vector of volumetric total water content (-)
                    mLayerVolFracLiqTrial,    & ! intent(inout): trial vector of volumetric liquid water content (-)
                    mLayerMatricHeadTrial,    & ! intent(inout): trial vector of total water matric potential (m)
                    mLayerMatricHeadLiqTrial, & ! intent(inout): trial vector of liquid water matric potential (m)
                    ! output: variables for the aquifer
                    scalarAquiferStorageTrial,& ! intent(inout): trial value of storage of water in the aquifer (m)
                    ! output: error control
                    err,cmessage)               ! intent(out):   error control
    if(err/=0)then; message=trim(message)//trim(cmessage); return; end if  ! (check for errors)

    ! update diagnostic variables and derivatives
    call updateVars(&
                    ! input
                    .false.,                                   & ! intent(in):    logical flag to adjust temperature to account for the energy used in melt+freeze
                    mpar_data,                                 & ! intent(in):    model parameters for a local HRU
                    indx_data,                                 & ! intent(in):    indices defining model states and layers
                    prog_data,                                 & ! intent(in):    model prognostic variables for a local HRU
                    diag_data,                                 & ! intent(inout): model diagnostic variables for a local HRU
                    deriv_data,                                & ! intent(inout): derivatives in model fluxes w.r.t. relevant state variables
                    ! output: variables for the vegetation canopy
                    scalarCanopyTempTrial,                     & ! intent(inout): trial value of canopy temperature (K)
                    scalarCanopyWatTrial,                      & ! intent(inout): trial value of canopy total water (kg m-2)
                    scalarCanopyLiqTrial,                      & ! intent(inout): trial value of canopy liquid water (kg m-2)
                    scalarCanopyIceTrial,                      & ! intent(inout): trial value of canopy ice content (kg m-2)
                    ! output: variables for the snow-soil domain
                    mLayerTempTrial,                           & ! intent(inout): trial vector of layer temperature (K)
                    mLayerVolFracWatTrial,                     & ! intent(inout): trial vector of volumetric total water content (-)
                    mLayerVolFracLiqTrial,                     & ! intent(inout): trial vector of volumetric liquid water content (-)
                    mLayerVolFracIceTrial,                     & ! intent(inout): trial vector of volumetric ice water content (-)
                    mLayerMatricHeadTrial,                     & ! intent(inout): trial vector of total water matric potential (m)
                    mLayerMatricHeadLiqTrial,                  & ! intent(inout): trial vector of liquid water matric potential (m)
                    ! output: error control
                    err,cmessage)                                ! intent(out):   error control
    if(err/=0)then; message=trim(message)//trim(cmessage); return; end if  ! (check for errors)

    if(updateCp)then
       ! *** compute volumetric heat capacity C_p
      if(ixHowHeatCap == enthalpyFD)then
        ! compute H_T without phase change
        call t2enthalpy(&
                         .false.,                    & ! intent(in):  ogical flag to not include phase change in enthalpy
                        ! input: data structures
                        diag_data,                   & ! intent(in):  model diagnostic variables for a local HRU
                        mpar_data,                   & ! intent(in):  parameter data structure
                        indx_data,                   & ! intent(in):  model indices
                        lookup_data,                 & ! intent(in):  lookup table data structure
                        ! input: state variables for the vegetation canopy
                        scalarCanairTempTrial,       & ! intent(in):  trial value of canopy air temperature (K)
                        scalarCanopyTempTrial,       & ! intent(in):  trial value of canopy temperature (K)
                        scalarCanopyWatTrial,        & ! intent(in):  trial value of canopy total water (kg m-2)
                        scalarCanopyIceTrial,        & ! intent(in):  trial value of canopy ice content (kg m-2)
                        ! input: variables for the snow-soil domain
                        mLayerTempTrial,             & ! intent(in):  trial vector of layer temperature (K)
                        mLayerVolFracWatTrial,       & ! intent(in):  trial vector of volumetric total water content (-)
                        mLayerMatricHeadTrial,       & ! intent(in):  trial vector of total water matric potential (m)
                        mLayerVolFracIceTrial,       & ! intent(in):  trial vector of volumetric fraction of ice (-)
                        ! input: pre-computed derivatives
                        dTheta_dTkCanopy,            & ! intent(in):  derivative in canopy volumetric liquid water content w.r.t. temperature (K-1)
                        scalarFracLiqVeg,            & ! intent(in):  fraction of canopy liquid water (-)
                        mLayerdTheta_dTk,            & ! intent(in):  derivative of volumetric liquid water content w.r.t. temperature (K-1)
                        mLayerFracLiqSnow,           & ! intent(in):  fraction of liquid water (-)
                        dVolTot_dPsi0,               & ! intent(in):  derivative in total water content w.r.t. total water matric potential (m-1)
                        ! output: enthalpy
                        scalarCanairEnthalpy,        & ! intent(out): enthalpy of the canopy air space (J m-3)
                        scalarCanopyEnthalpyTrial,   & ! intent(out): enthalpy of the vegetation canopy (J m-3)
                        mLayerEnthalpyTrial,         & ! intent(out): enthalpy of each snow+soil layer (J m-3)
                        dCanEnthalpy_dTk,            & ! intent(out): derivatives in canopy enthalpy w.r.t. temperature
                        dCanEnthalpy_dWat,           & ! intent(out): derivatives in canopy enthalpy w.r.t. water state
                        dEnthalpy_dTk,               & ! intent(out): derivatives in layer enthalpy w.r.t. temperature
                        dEnthalpy_dWat,              & ! intent(out): derivatives in layer enthalpy w.r.t. water state
                        ! output: error control
                        err,cmessage)                  ! intent(out): error control
        if(err/=0)then; message=trim(message)//trim(cmessage); return; endif


        ! *** compute volumetric heat capacity C_p = dH_T/dT
        call computHeatCap(&
                            ! input: control variables
                            nLayers,                   & ! intent(in):    number of layers (soil+snow)
                            computeVegFlux,            & ! intent(in):    flag to denote if computing the vegetation flux
                            canopyDepth,               & ! intent(in):    canopy depth (m)
                            ! input output data structures   
                            mpar_data,                 & ! intent(in):    model parameters
                            indx_data,                 & ! intent(in):    model layer indices
                            diag_data,                 & ! intent(inout): model diagnostic variables for a local HRU
                            ! input: state variables
                            scalarCanopyIceTrial,      & ! intent(in):    trial value for mass of ice on the vegetation canopy (kg m-2)
                            scalarCanopyLiqTrial,      & ! intent(in):    trial value for the liquid water on the vegetation canopy (kg m-2)
                            scalarCanopyTempTrial,     & ! intent(in):    trial value of canopy temperature (K)
                            scalarCanopyTemp,          & ! intent(in):    previous value of canopy temperature (K)
                            scalarCanopyEnthalpyTrial, & ! intent(in):    trial enthalpy of the vegetation canopy (J m-3)
                            scalarCanopyEnthalpy,      & ! intent(in):    previous enthalpy of the vegetation canopy (J m-3)
                            mLayerVolFracIceTrial,     & ! intent(in):    volumetric fraction of ice at the start of the sub-step (-)
                            mLayerVolFracLiqTrial,     & ! intent(in):    volumetric fraction of liquid water at the start of the sub-step (-)
                            mLayerTempTrial,           & ! intent(in):    trial temperature
                            mLayerTemp,                & ! intent(in):    previous temperature
                            mLayerEnthalpyTrial,       & ! intent(in):    trial enthalpy for snow and soil
                            mLayerEnthalpy,            & ! intent(in):    previous enthalpy for snow and soil
                            mLayerMatricHeadTrial,     & ! intent(in):    trial total water matric potential (m)
                            ! input: pre-computed derivatives   
                            dTheta_dTkCanopy,          & ! intent(in):    derivative in canopy volumetric liquid water content w.r.t. temperature (K-1)
                            scalarFracLiqVeg,          & ! intent(in):    fraction of canopy liquid water (-)
                            mLayerdTheta_dTk,          & ! intent(in):    derivative of volumetric liquid water content w.r.t. temperature (K-1)
                            mLayerFracLiqSnow,         & ! intent(in):    fraction of liquid water (-)
                            dVolTot_dPsi0,             & ! intent(in):    derivative in total water content w.r.t. total water matric potential (m-1)
                            dCanEnthalpy_dTk,          & ! intent(in):    derivatives in canopy enthalpy w.r.t. temperature
                            dCanEnthalpy_dWat,         & ! intent(in):    derivatives in canopy enthalpy w.r.t. water state
                            dEnthalpy_dTk,             & ! intent(in):    derivatives in layer enthalpy w.r.t. temperature
                            dEnthalpy_dWat,            & ! intent(in):    derivatives in layer enthalpy w.r.t. water state
                            ! output
                            heatCapVegTrial,           & ! intent(out):   volumetric heat capacity of vegetation canopy
                            mLayerHeatCapTrial,        & ! intent(out):   heat capacity for snow and soil
                            dVolHtCapBulk_dPsi0,       & ! intent(out):   derivative in bulk heat capacity w.r.t. matric potential
                            dVolHtCapBulk_dTheta,      & ! intent(out):   derivative in bulk heat capacity w.r.t. volumetric water content
                            dVolHtCapBulk_dCanWat,     & ! intent(out):   derivative in bulk heat capacity w.r.t. volumetric water content
                            dVolHtCapBulk_dTk,         & ! intent(out):   derivative in bulk heat capacity w.r.t. temperature
                            dVolHtCapBulk_dTkCanopy,   & ! intent(out):   derivative in bulk heat capacity w.r.t. temperature                  
                            ! output: error control
                            err,cmessage)                ! intent(out):   error control
        if(err/=0)then; message=trim(message)//trim(cmessage); return; endif
        ! update values
        mLayerEnthalpy = mLayerEnthalpyTrial
        scalarCanopyEnthalpy = scalarCanopyEnthalpyTrial
      else if(ixHowHeatCap == closedForm)then
        call computHeatCapAnalytic(&
                          ! input: control variables
                          computeVegFlux,              & ! intent(in):    flag to denote if computing the vegetation flux
                          canopyDepth,                 & ! intent(in):    canopy depth (m)
                          ! input: state variables 
                          scalarCanopyIceTrial,        & ! intent(in):    trial value for mass of ice on the vegetation canopy (kg m-2)
                          scalarCanopyLiqTrial,        & ! intent(in):    trial value for the liquid water on the vegetation canopy (kg m-2)
                          scalarCanopyTempTrial,       & ! intent(in):    trial value of canopy temperature (K)
                          mLayerVolFracIceTrial,       & ! intent(in):    volumetric fraction of ice at the start of the sub-step (-)
                          mLayerVolFracLiqTrial,       & ! intent(in):    fraction of liquid water at the start of the sub-step (-)
                          mLayerTempTrial,             & ! intent(in):    trial value of layer temperature (K)
                          mLayerMatricHeadTrial,       & ! intent(in):    trial total water matric potential (m)
                          ! input: pre-computed derivatives 
                          dTheta_dTkCanopy,            & ! intent(in):    derivative in canopy volumetric liquid water content w.r.t. temperature (K-1)
                          scalarFracLiqVeg,            & ! intent(in):    fraction of canopy liquid water (-)
                          mLayerdTheta_dTk,            & ! intent(in):    derivative of volumetric liquid water content w.r.t. temperature (K-1)
                          mLayerFracLiqSnow,           & ! intent(in):    fraction of liquid water (-)
                          dVolTot_dPsi0,               & ! intent(in):    derivative in total water content w.r.t. total water matric potential (m-1)
                          ! input output data structures
                          mpar_data,                   & ! intent(in):    model parameters
                          indx_data,                   & ! intent(in):    model layer indices
                          diag_data,                   & ! intent(inout): model diagnostic variables for a local HRU
                          ! output
                          heatCapVegTrial,             & ! intent(out):   volumetric heat capacity of vegetation canopy
                          mLayerHeatCapTrial,          & ! intent(out):   volumetric heat capacity of soil and snow
                          dVolHtCapBulk_dPsi0,         & ! intent(out):   derivative in bulk heat capacity w.r.t. matric potential
                          dVolHtCapBulk_dTheta,        & ! intent(out):   derivative in bulk heat capacity w.r.t. volumetric water content
                          dVolHtCapBulk_dCanWat,       & ! intent(out):   derivative in bulk heat capacity w.r.t. volumetric water content
                          dVolHtCapBulk_dTk,           & ! intent(out):   derivative in bulk heat capacity w.r.t. temperature
                          dVolHtCapBulk_dTkCanopy,     & ! intent(out):   derivative in bulk heat capacity w.r.t. temperature                  
                          ! output: error control
                          err,cmessage)                  ! intent(out):   error control
      endif !(choice of how compute heat capacity)

      ! compute multiplier of state vector
      call computStatMult(&
                    ! input
                    heatCapVegTrial,                  & ! intent(in):  volumetric heat capacity of vegetation canopy
                    mLayerHeatCapTrial,               & ! intent(in):  volumetric heat capacity of soil and snow
                    diag_data,                        & ! intent(in):  model diagnostic variables for a local HRU
                    indx_data,                        & ! intent(in):  indices defining model states and layers
                    ! output
                    sMul,                             & ! intent(out): multiplier for state vector (used in the residual calculations)
                    err,cmessage)                       ! intent(out): error control
      if(err/=0)then; message=trim(message)//trim(cmessage); return; endif  ! (check for errors)

      ! update thermal conductivity
      call computThermConduct(&
                          ! input: control variables
                          computeVegFlux,               & ! intent(in):    flag to denote if computing the vegetation flux
                          nLayers,                      & ! intent(in):    total number of layers
                          canopyDepth,                  & ! intent(in):    canopy depth (m)
                          ! input: state variables
                          scalarCanopyIceTrial,         & ! intent(in):    trial value for mass of ice on the vegetation canopy (kg m-2)
                          scalarCanopyLiqTrial,         & ! intent(in):    trial value of canopy liquid water (kg m-2)
                          mLayerTempTrial,              & ! intent(in):    trial temperature of layer temperature (K)
                          mLayerMatricHeadTrial,        & ! intent(in):    trial value for total water matric potential (m)
                          mLayerVolFracIceTrial,        & ! intent(in):    volumetric fraction of ice at the start of the sub-step (-)
                          mLayerVolFracLiqTrial,        & ! intent(in):    volumetric fraction of liquid water at the start of the sub-step (-)
                         ! input: pre-computed derivatives
                          mLayerdTheta_dTk,             & ! intent(in):    derivative in volumetric liquid water content w.r.t. temperature (K-1)
                          mLayerFracLiqSnow,            & ! intent(in):    fraction of liquid water (-)
                          ! input/output: data structures
                          mpar_data,                    & ! intent(in):    model parameters
                          indx_data,                    & ! intent(in):    model layer indices
                          prog_data,                    & ! intent(in):    model prognostic variables for a local HRU
                          diag_data,                    & ! intent(inout): model diagnostic variables for a local HRU
                          ! output: derivatives
                          dThermalC_dWatAbove,          & ! intent(out):   derivative in the thermal conductivity w.r.t. water state in the layer above
                          dThermalC_dWatBelow,          & ! intent(out):   derivative in the thermal conductivity w.r.t. water state in the layer above
                          dThermalC_dTempAbove,         & ! intent(out):   derivative in the thermal conductivity w.r.t. energy state in the layer above
                          dThermalC_dTempBelow,         & ! intent(out):   derivative in the thermal conductivity w.r.t. energy state in the layer above
                          ! output: error control
                          err,cmessage)                   ! intent(out): error control
      if(err/=0)then; err=55; message=trim(message)//trim(cmessage); return; end if
    else
      ! set heat capacity derivatives to 0 for constant through step
      dVolHtCapBulk_dPsi0     = 0._rkind
      dVolHtCapBulk_dTheta    = 0._rkind
      dVolHtCapBulk_dCanWat   = 0._rkind
      dVolHtCapBulk_dTk       = 0._rkind
      dVolHtCapBulk_dTkCanopy = 0._rkind
      dThermalC_dWatAbove     = 0._rkind
      dThermalC_dWatBelow     = 0._rkind
      dThermalC_dTempAbove    = 0._rkind
      dThermalC_dTempBelow    = 0._rkind  
    endif ! updateCp

    if(needCm)then
      ! compute C_m
      call computCm(&
                      ! input: control variables
                      computeVegFlux,           & ! intent(in):  flag to denote if computing the vegetation flux
                      ! input: state variables
                      scalarCanopyTempTrial,    & ! intent(in):  trial value of canopy temperature (K)
                      mLayerTempTrial,          & ! intent(in):  trial value of layer temperature (K)
                      mLayerMatricHeadTrial,    & ! intent(in):  trial value for total water matric potential (m)
                      ! input data structures
                      mpar_data,                & ! intent(in):  model parameters
                      indx_data,                & ! intent(in):  model layer indices
                      ! output
                      scalarCanopyCmTrial,      & ! intent(out): Cm for vegetation
                      mLayerCmTrial,            & ! intent(out): Cm for soil and snow
                      err,cmessage)               ! intent(out): error control
    else
      scalarCanopyCmTrial = 0._qp
      mLayerCmTrial = 0._qp
    endif ! needCm

    ! save the number of flux calls per time step
    indx_data%var(iLookINDEX%numberFluxCalc)%dat(1) = indx_data%var(iLookINDEX%numberFluxCalc)%dat(1) + 1

    ! only need to check longwave balance with numerical recipes solver 
    checkLWBalance = .false.
    if(ixNumericalMethod==numrec) checkLWBalance = .true.

    ! compute the fluxes for a given state vector
    call computFlux(&
                    ! input-output: model control
                    nSnow,                     & ! intent(in):    number of snow layers
                    nSoil,                     & ! intent(in):    number of soil layers
                    nLayers,                   & ! intent(in):    total number of layers
                    firstSubStep,              & ! intent(in):    flag to indicate if we are processing the first sub-step
                    firstFluxCall,             & ! intent(inout): flag to denote the first flux call
                    firstSplitOper,            & ! intent(in):    flag to indicate if we are processing the first flux call in a splitting operation
                    computeVegFlux,            & ! intent(in):    flag to indicate if we need to compute fluxes over vegetation
                    scalarSolution,            & ! intent(in):    flag to indicate the scalar solution
                    checkLWBalance,            & ! intent(in):    flag to check longwave balance
                    scalarSfcMeltPond/dt,      & ! intent(in):    drainage from the surface melt pond (kg m-2 s-1)
                    ! input: state variables
                    scalarCanairTempTrial,     & ! intent(in):    trial value for the temperature of the canopy air space (K)
                    scalarCanopyTempTrial,     & ! intent(in):    trial value for the temperature of the vegetation canopy (K)
                    mLayerTempTrial,           & ! intent(in):    trial value for the temperature of each snow and soil layer (K)
                    mLayerMatricHeadLiqTrial,  & ! intent(in):    trial value for the liquid water matric potential in each soil layer (m)
                    mLayerMatricHeadTrial,     & ! intent(in):    trial vector of total water matric potential (m)
                    scalarAquiferStorageTrial, & ! intent(in):    trial value of storage of water in the aquifer (m)
                    ! input: diagnostic variables defining the liquid water and ice content
                    scalarCanopyLiqTrial,      & ! intent(in):    trial value for the liquid water on the vegetation canopy (kg m-2)
                    scalarCanopyIceTrial,      & ! intent(in):    trial value for the ice on the vegetation canopy (kg m-2)
                    mLayerVolFracLiqTrial,     & ! intent(in):    trial value for the volumetric liquid water content in each snow and soil layer (-)
                    mLayerVolFracIceTrial,     & ! intent(in):    trial value for the volumetric ice in each snow and soil layer (-)
                    ! input: data structures
                    model_decisions,           & ! intent(in):    model decisions
                    type_data,                 & ! intent(in):    type of vegetation and soil
                    attr_data,                 & ! intent(in):    spatial attributes
                    mpar_data,                 & ! intent(in):    model parameters
                    forc_data,                 & ! intent(in):    model forcing data
                    bvar_data,                 & ! intent(in):    average model variables for the entire basin
                    prog_data,                 & ! intent(in):    model prognostic variables for a local HRU
                    indx_data,                 & ! intent(in):    index data
                    ! input-output: data structures
                    diag_data,                 & ! intent(inout): model diagnostic variables for a local HRU
                    flux_data,                 & ! intent(inout): model fluxes for a local HRU
                    deriv_data,                & ! intent(out):   derivatives in model fluxes w.r.t. relevant state variables
                    ! input-output: flux vector and baseflow derivatives
                    ixSaturation,              & ! intent(inout): index of the lowest saturated layer (NOTE: only computed on the first iteration)
                    dBaseflow_dMatric,         & ! intent(out):   derivative in baseflow w.r.t. matric head (s-1)
                    fluxVec,                   & ! intent(out):   flux vector (mixed units)
                    ! output: error control
                    err,cmessage)                ! intent(out):   error code and error message
    if(err/=0)then; message=trim(message)//trim(cmessage); return; end if  ! (check for errors)

    ! compute soil compressibility (-) and its derivative w.r.t. matric head (m)
    ! NOTE: we already extracted trial matrix head and volumetric liquid water as part of the flux calculations
    ! use non-prime version
    call soilCmpres(&
                    ! input:
                    dt_cur,                                 & ! intent(in):    length of the time step (seconds)
                    ixRichards,                             & ! intent(in):    choice of option for Richards' equation
                    ixBeg,ixEnd,                            & ! intent(in):    start and end indices defining desired layers
                    mLayerMatricHead(1:nSoil),              & ! intent(in):    matric head at the start of the time step (m)
                    mLayerMatricHeadTrial(1:nSoil),         & ! intent(in):    trial value of matric head (m)
                    mLayerVolFracLiqTrial(nSnow+1:nLayers), & ! intent(in):    trial value for the volumetric liquid water content in each soil layer (-)
                    mLayerVolFracIceTrial(nSnow+1:nLayers), & ! intent(in):    trial value for the volumetric ice content in each soil layer (-)
                    specificStorage,                        & ! intent(in):    specific storage coefficient (m-1)
                    theta_sat,                              & ! intent(in):    soil porosity (-)
                    ! output:
                    mLayerCompress,                         & ! intent(inout): compressibility of the soil matrix (-)
                    dCompress_dPsi,                         & ! intent(inout): derivative in compressibility w.r.t. matric head (m-1)
                    err,cmessage)                             ! intent(out):   error code and error message
    if(err/=0)then; message=trim(message)//trim(cmessage); return; end if  ! (check for errors)

    ! compute the total change in storage associated with compression of the soil matrix (kg m-2 s-1)
    scalarSoilCompress = sum(mLayerCompress(1:nSoil)*mLayerDepth(nSnow+1:nLayers))*iden_water

    ! compute the residual vector
    call computResid(&
                      ! input: model control
                      dt_cur,                    & ! intent(in):  length of the time step (seconds)
                      nSnow,                     & ! intent(in):  number of snow layers
                      nSoil,                     & ! intent(in):  number of soil layers
                      nLayers,                   & ! intent(in):  total number of layers
                      ! input: flux vectors
                      sMul,                      & ! intent(in):  state vector multiplier (used in the residual calculations)
                      fluxVec,                   & ! intent(in):  flux vector
                      ! input: state variables (already disaggregated into scalars and vectors)
                      scalarCanairTempTrial,     & ! intent(in):  trial value for the temperature of the canopy air space (K)
                      scalarCanopyTempTrial,     & ! intent(in):  trial value for the temperature of the vegetation canopy (K)
                      scalarCanopyWatTrial,      & ! intent(in):  trial value for the water on the vegetation canopy (kg m-2)
                      mLayerTempTrial,           & ! intent(in):  trial value for the temperature of each snow and soil layer (K)
                      scalarAquiferStorageTrial, & ! intent(in):  trial value of storage of water in the aquifer (m)
                      ! input: diagnostic variables defining the liquid water and ice content (function of state variables)
                      scalarCanopyIceTrial,      & ! intent(in):  trial value for the ice on the vegetation canopy (kg m-2)
                      scalarCanopyLiqTrial,      & ! intent(in):  trial value for the liq on the vegetation canopy (kg m-2)
                      mLayerVolFracIceTrial,     & ! intent(in):  trial value for the volumetric ice in each snow and soil layer (-)
                      mLayerVolFracWatTrial,     & ! intent(in):  trial value for the volumetric water in each snow and soil layer (-)
                      mLayerVolFracLiqTrial,     & ! intent(in):  trial value for the volumetric liq in each snow and soil layer (-)
                      scalarCanopyCmTrial,       & ! intent(in):  Cm of vegetation canopy (-)
                      mLayerCmTrial,             & ! intent(in):  Cm of each snow and soil layer (-)
                      ! input: data structures
                      prog_data,                 & ! intent(in):  model prognostic variables for a local HRU
                      diag_data,                 & ! intent(in):  model diagnostic variables for a local HRU
                      flux_data,                 & ! intent(in):  model fluxes for a local HRU
                      indx_data,                 & ! intent(in):  index data
                      ! output
                      resSink,                   & ! intent(out): additional (sink) terms on the RHS of the state equation
                      resVec,                    & ! intent(out): residual vector
                      err,cmessage)                ! intent(out): error control
    if(err/=0)then; message=trim(message)//trim(cmessage); return; end if  ! (check for errors)

    ! compute the function evaluation
    rVecScaled = fScale(:)*real(resVec(:), rkind)   ! scale the residual vector (NOTE: residual vector is in quadruple precision)
    fEval      = 0.5_rkind*dot_product(rVecScaled,rVecScaled)

  ! end association with the information in the data structures
  end associate

end subroutine eval8summa

#ifdef SUNDIALS_ACTIVE
! **********************************************************************************************************
! public function eval8summa4kinsol: compute the residual vector F(t,y) required for KINSOL solver
! **********************************************************************************************************
! Return values:
!    0 = success,
!    1 = recoverable error,
!   -1 = non-recoverable error
! ----------------------------------------------------------------
integer(c_int) function eval8summa4kinsol(sunvec_y, sunvec_r, user_data) &
      result(ierr) bind(C,name='eval8summa4kinsol')

  !======= Inclusions ===========
  use, intrinsic :: iso_c_binding
  use fsundials_nvector_mod
  use fnvector_serial_mod
  use type4kinsol

  !======= Declarations =========
  implicit none

  ! calling variables
  type(N_Vector)              :: sunvec_y  ! solution N_Vector    y
  type(N_Vector)              :: sunvec_r  ! residual N_Vector    F(t,y)
  type(c_ptr), value          :: user_data ! user-defined data

  ! pointers to data in SUNDIALS vectors
  type(data4kinsol), pointer  :: eqns_data ! equations data
  real(rkind), pointer        :: stateVec(:)
  real(rkind), pointer        :: rVec(:)
  logical(lgt)                :: feasible
  real(rkind)                 :: fNew      ! function values, not needed here
  !======= Internals ============

  ! get equations data from user-defined data
  call c_f_pointer(user_data, eqns_data)

  ! get data arrays from SUNDIALS vectors
  stateVec(1:eqns_data%nState)  => FN_VGetArrayPointer(sunvec_y)
  rVec(1:eqns_data%nState)  => FN_VGetArrayPointer(sunvec_r)

  ! compute the flux and the residual vector for a given state vector
  call eval8summa(&
                ! input: model control
                eqns_data%dt_cur,                  & ! intent(in):    current stepsize
                eqns_data%dt,                      & ! intent(in):    data step
                eqns_data%nSnow,                   & ! intent(in):    number of snow layers
                eqns_data%nSoil,                   & ! intent(in):    number of soil layers
                eqns_data%nLayers,                 & ! intent(in):    number of layers
                eqns_data%nState,                  & ! intent(in):    number of state variables in the current subset
                .true.,                            & ! intent(in):    inside SUNDIALS solver
                eqns_data%firstSubStep,            & ! intent(in):    flag to indicate if we are processing the first sub-step
                eqns_data%firstFluxCall,           & ! intent(inout): flag to indicate if we are processing the first flux call
                eqns_data%firstSplitOper,          & ! intent(in):    flag to indicate if we are processing the first flux call in a splitting operation
                eqns_data%computeVegFlux,          & ! intent(in):    flag to indicate if we need to compute fluxes over vegetation
                eqns_data%scalarSolution,          & ! intent(in):    flag to indicate the scalar solution
                ! input: state vectors
                stateVec,                          & ! intent(in):    model state vector
                eqns_data%fScale,                  & ! intent(in):    characteristic scale of the function evaluations
                eqns_data%sMul,                    & ! intent(inout): state vector multiplier (used in the residual calculations)
                ! input: data structures
                eqns_data%model_decisions,         & ! intent(in):    model decisions
                eqns_data%lookup_data,             & ! intent(in):    lookup data
                eqns_data%type_data,               & ! intent(in):    type of vegetation and soil
                eqns_data%attr_data,               & ! intent(in):    spatial attributes
                eqns_data%mpar_data,               & ! intent(in):    model parameters
                eqns_data%forc_data,               & ! intent(in):    model forcing data
                eqns_data%bvar_data,               & ! intent(in):    average model variables for the entire basin
                eqns_data%prog_data,               & ! intent(in):    model prognostic variables for a local HRU
                ! input-output: data structures
                eqns_data%indx_data,               & ! intent(inout): index data
                eqns_data%diag_data,               & ! intent(inout): model diagnostic variables for a local HRU
                eqns_data%flux_data,               & ! intent(inout): model fluxes for a local HRU (initial flux structure)
                eqns_data%deriv_data,              & ! intent(inout): derivatives in model fluxes w.r.t. relevant state variables
                 ! input-output: baseflow
                eqns_data%ixSaturation,            & ! intent(inout): index of the lowest saturated layer
                eqns_data%dBaseflow_dMatric,       & ! intent(out):   derivative in baseflow w.r.t. matric head (s-1)
                 ! output: flux and residual vectors
                feasible,                          & ! intent(out):   flag to denote the feasibility of the solution always true inside SUNDIALS
                eqns_data%fluxVec,                 & ! intent(out):   flux vector
                eqns_data%resSink,                 & ! intent(out):   additional (sink) terms on the RHS of the state equation
                rVec,                              & ! intent(out):   residual vector
                fNew,                              & ! intent(out):   new function evaluation
                eqns_data%err,eqns_data%message)     ! intent(out):   error control
  if(eqns_data%err > 0)then; eqns_data%message=trim(eqns_data%message); ierr=-1; return; endif
  if(eqns_data%err < 0)then; eqns_data%message=trim(eqns_data%message); ierr=1; return; endif
  
  ! return success
  ierr = 0
  return

end function eval8summa4kinsol
#endif

end module eval8summa_module

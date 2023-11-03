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

! named variables to describe the state variable type
USE globalData,only:iname_watLayer  ! named variable defining the total water state variable for snow+soil layers
USE globalData,only:iname_liqLayer  ! named variable defining the liquid  water state variable for snow+soil layers

! constants
USE multiconst,only:&
                    Tfreeze,      & ! temperature at freezing              (K)
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
public::imposeConstraints

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
    dCm_dTk                 => deriv_data%var(iLookDERIV%dCm_dTk)%dat                 ,& ! intent(out): [dp(:)] derivative in heat capacity w.r.t. temperature (J kg-1 K-2)
    dCm_dTkCanopy           => deriv_data%var(iLookDERIV%dCm_dTkCanopy)%dat(1)        ,& ! intent(out): [dp   ] derivative in heat capacity w.r.t. canopy temperature (J kg-1 K-2)
    ! mapping
    ixMapFull2Subset        => indx_data%var(iLookINDEX%ixMapFull2Subset)%dat         ,& ! intent(in):  [i4b(:)]mapping of full state vector to the state subset
    ixControlVolume         => indx_data%var(iLookINDEX%ixControlVolume)%dat          ,& ! intent(in):  [i4b(:)]index of control volume for different domains (veg, snow, soil)
    ! heat capacity
    heatCapVegTrial         => diag_data%var(iLookDIAG%scalarBulkVolHeatCapVeg)%dat(1),& ! intent(out): [dp]    volumetric heat capacity of vegetation canopy
    mLayerHeatCapTrial      => diag_data%var(iLookDIAG%mLayerVolHtCapBulk)%dat        ,& ! intent(out): [dp(:)] heat capacity for snow and soil
    ! Cm
    canopyCmTrial           => diag_data%var(iLookDIAG%scalarCanopyCm)%dat(1)         ,& ! intent(out): [dp]    Cm of the canopy
    mLayerCmTrial           => diag_data%var(iLookDIAG%mLayerCm)%dat                   & ! intent(out): [dp(:)] Cm of snow and soil
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
                      mLayerMatricHeadTrial,    & ! intent(in):  trial value for total water matric potential (-)
                      ! input data structures
                      mpar_data,                & ! intent(in):  model parameters
                      indx_data,                & ! intent(in):  model layer indices
                      ! output
                      canopyCmTrial,            & ! intent(out): Cm for vegetation (J kg K-1)
                      mLayerCmTrial,            & ! intent(out): Cm for soil and snow (J kg K-1)
                      dCm_dTk,                  & ! intent(out): derivative in Cm w.r.t. temperature (J kg K-2)
                      dCm_dTkCanopy,            & ! intent(out): derivative in Cm w.r.t. temperature (J kg K-2)
                      err,cmessage)               ! intent(out): error control
    else
      canopyCmTrial = 0._qp
      mLayerCmTrial = 0._qp
      dCm_dTk       = 0._rkind
      dCm_dTkCanopy = 0._rkind
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
                      canopyCmTrial,             & ! intent(in):  Cm of vegetation canopy (-)
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
  type(N_Vector)              :: sunvec_y    ! solution N_Vector    y
  type(N_Vector)              :: sunvec_r    ! residual N_Vector    F(t,y)
  type(c_ptr), value          :: user_data   ! user-defined data

  ! pointers to data in SUNDIALS vectors
  type(data4kinsol), pointer  :: eqns_data   ! equations data
  real(rkind), pointer        :: stateVec(:) ! solution vector
  real(rkind), pointer        :: rVec(:)     ! residual vector
  logical(lgt)                :: feasible    ! feasibility of state vector
  real(rkind)                 :: fNew        ! function values, not needed here
  integer(i4b)                :: err         ! error in imposeConstraints
  character(len=256)          :: message     ! error message of downwind routine
  !======= Internals ============

  ! get equations data from user-defined data
  call c_f_pointer(user_data, eqns_data)

  ! get data arrays from SUNDIALS vectors
  stateVec(1:eqns_data%nState)  => FN_VGetArrayPointer(sunvec_y)
  rVec(1:eqns_data%nState)  => FN_VGetArrayPointer(sunvec_r)

  ! increment the proposed iteration for simple error control if needed
  if (eqns_data%firstStateiteration) then
    eqns_data%firstStateIteration = .false.
  else
    call imposeConstraints(eqns_data%model_decisions,eqns_data%indx_data,eqns_data%prog_data,eqns_data%mpar_data,stateVec(:), &
      eqns_data%stateVecPrev, eqns_data%nState, eqns_data%nSoil, eqns_data%nSnow, message, err)
     if(err/=0)then; ierr=1; message="eval8summa4kinsol/"//trim(message); print*, message; return; end if  ! (check for errors)
  endif
  eqns_data%stateVecPrev = stateVec(:)  ! save the state vector for the next iteration
  
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

! ***************************************************************************************************************************************
! public subroutine imposeConstraints: impose solution constraints 
!   This is simple error control to reduce possible temperature increments, cross over freezing point events, and keep the state feasible
!   Imposed after the internal call of KINSOL incrementing the linesearch
! ***************************************************************************************************************************************
subroutine imposeConstraints(model_decisions,indx_data, prog_data, mpar_data, stateVec, stateVecPrev,&
    nState, nSoil, nSnow, message, err)
  ! external functions
  USE snow_utils_module,only:fracliquid                           ! compute the fraction of liquid water at a given temperature (snow)
  USE soil_utils_module,only:crit_soilT                           ! compute the critical temperature below which ice exists
  implicit none
  
  type(model_options),intent(in)           :: model_decisions(:)  ! model decisions
  type(var_ilength),intent(in)             :: indx_data           ! indices defining model states and layers
  type(var_dlength),intent(in)             :: prog_data           ! prognostic variables for a local HRU
  type(var_dlength),intent(in)             :: mpar_data           ! model parameters
  real(rkind), intent(inout)               :: stateVec(:)         ! state vector
  real(rkind), intent(in)                  :: stateVecPrev(:)     ! previous state vector
  integer(i4b),intent(in)                  :: nState              ! total number of state variables
  integer(i4b),intent(in)                  :: nSoil               ! number of soil layers
  integer(i4b),intent(in)                  :: nSnow               ! number of snow layers
  integer(i4b),intent(out)                 :: err                 ! error code
  character(len=256),intent(out)           :: message             ! error message
  ! -----------------------------------------------------------------------------------------------------
  ! temporary variables for model constraints
  real(rkind)                              :: cInc                       ! constrained temperature increment (K) -- simplified bi-section
  real(qp),dimension(nState)               :: xInc                       ! iteration increment
  real(rkind)                              :: xIncFactor                 ! scaling factor for the iteration increment (-)
  integer(i4b)                             :: iMax(1)                    ! index of maximum temperature
  real(rkind)                              :: scalarTemp                 ! temperature of an individual snow layer (K)
  real(rkind)                              :: scalarIce                  ! volumetric ice content of an individual layer (-)
  real(rkind)                              :: volFracLiq                 ! volumetric liquid water content of an individual layer (-)
  logical(lgt),dimension(nSoil)            :: crosFlag                   ! flag to denote temperature crossing from unfrozen to frozen (or vice-versa)
  logical(lgt)                             :: crosTempVeg                ! flag to denoote where temperature crosses the freezing point
  real(rkind)                              :: xPsi00                     ! matric head after applying the iteration increment (m)
  real(rkind)                              :: TcSoil                     ! critical point when soil begins to freeze (K)
  real(rkind)                              :: critDiff                   ! temperature difference from critical (K)
  real(rkind)                              :: epsT                       ! small interval above/below critical (K)
  real(rkind)                              :: zMaxTempIncrement          ! maximum temperature increment (K)
  real(rkind)                              :: zMaxMatricIncrement        ! maximum matric head increment (m)
  ! indices of model state variables
  integer(i4b)                             :: iState                     ! index of state within a specific variable type
  integer(i4b)                             :: ixNrg,ixLiq                ! index of energy and mass state variables in full state vector
  ! indices of model layers
  integer(i4b)                             :: iLayer                     ! index of model layer
  ! choice of constraints to impose
  logical(lgt)                             :: small_delTemp              ! flag to constain temperature change to be less than zMaxTempIncrement
  logical(lgt)                             :: small_delMatric            ! flag to constain matric head change to be less than zMaxMatricIncrement
  logical(lgt)                             :: detect_events              ! flag to do freezing point event detection and cross-over with epsT
  logical(lgt)                             :: water_bounds               ! flag to force water to not go above or below physical bounds
  
  ! -----------------------------------------------------------------------------------------------------
  ! association to variables in the data structures
  associate(&
    ! model decisions
    ixNumericalMethod       => model_decisions(iLookDECISIONS%num_method)%iDecision   ,& ! intent(in):  [i4b]   choice of numerical solver
    ! indices of model state variables
    ixNrgOnly               => indx_data%var(iLookINDEX%ixNrgOnly)%dat                ,& ! intent(in): [i4b(:)] list of indices in the state subset for energy states
    ixHydOnly               => indx_data%var(iLookINDEX%ixHydOnly)%dat                ,& ! intent(in): [i4b(:)] list of indices in the state subset for hydrology states
    ixMatOnly               => indx_data%var(iLookINDEX%ixMatOnly)%dat                ,& ! intent(in): [i4b(:)] list of indices in the state subset for matric head states
    ixMassOnly              => indx_data%var(iLookINDEX%ixMassOnly)%dat               ,& ! intent(in): [i4b(:)] list of indices in the state subset for canopy storage states
    ixHydType               => indx_data%var(iLookINDEX%ixHydType)%dat                ,& ! intent(in): [i4b(:)] index of the type of hydrology states in snow+soil domain
    ixStateType_subset      => indx_data%var(iLookINDEX%ixStateType_subset)%dat       ,& ! intent(in): [i4b(:)] named variables defining the states in the subset
    ! indices for specific state variables
    ixCasNrg                => indx_data%var(iLookINDEX%ixCasNrg)%dat(1)              ,& ! intent(in): [i4b]    index of canopy air space energy state variable
    ixVegNrg                => indx_data%var(iLookINDEX%ixVegNrg)%dat(1)              ,& ! intent(in): [i4b]    index of canopy energy state variable
    ixVegHyd                => indx_data%var(iLookINDEX%ixVegHyd)%dat(1)              ,& ! intent(in): [i4b]    index of canopy hydrology state variable (mass)
    ixTopNrg                => indx_data%var(iLookINDEX%ixTopNrg)%dat(1)              ,& ! intent(in): [i4b]    index of upper-most energy state in the snow-soil subdomain
    ixTopHyd                => indx_data%var(iLookINDEX%ixTopHyd)%dat(1)              ,& ! intent(in): [i4b]    index of upper-most hydrology state in the snow-soil subdomain
    ! vector of energy indices for the snow and soil domains
    ! NOTE: states not in the subset are equal to integerMissing
    ixSnowSoilNrg           => indx_data%var(iLookINDEX%ixSnowSoilNrg)%dat            ,& ! intent(in): [i4b(:)] index in the state subset for energy state variables in the snow+soil domain
    ixSnowOnlyNrg           => indx_data%var(iLookINDEX%ixSnowOnlyNrg)%dat            ,& ! intent(in): [i4b(:)] index in the state subset for energy state variables in the snow domain
    ixSoilOnlyNrg           => indx_data%var(iLookINDEX%ixSoilOnlyNrg)%dat            ,& ! intent(in): [i4b(:)] index in the state subset for energy state variables in the soil domain
    ! vector of hydrology indices for the snow and soil domains
    ! NOTE: states not in the subset are equal to integerMissing
    ixSnowSoilHyd           => indx_data%var(iLookINDEX%ixSnowSoilHyd)%dat            ,& ! intent(in): [i4b(:)] index in the state subset for hydrology state variables in the snow+soil domain
    ixSnowOnlyHyd           => indx_data%var(iLookINDEX%ixSnowOnlyHyd)%dat            ,& ! intent(in): [i4b(:)] index in the state subset for hydrology state variables in the snow domain
    ixSoilOnlyHyd           => indx_data%var(iLookINDEX%ixSoilOnlyHyd)%dat            ,& ! intent(in): [i4b(:)] index in the state subset for hydrology state variables in the soil domain
    ! number of state variables of a specific type
    nSnowSoilNrg            => indx_data%var(iLookINDEX%nSnowSoilNrg )%dat(1)         ,& ! intent(in): [i4b]    number of energy state variables in the snow+soil domain
    nSnowOnlyNrg            => indx_data%var(iLookINDEX%nSnowOnlyNrg )%dat(1)         ,& ! intent(in): [i4b]    number of energy state variables in the snow domain
    nSoilOnlyNrg            => indx_data%var(iLookINDEX%nSoilOnlyNrg )%dat(1)         ,& ! intent(in): [i4b]    number of energy state variables in the soil domain
    nSnowSoilHyd            => indx_data%var(iLookINDEX%nSnowSoilHyd )%dat(1)         ,& ! intent(in): [i4b]    number of hydrology variables in the snow+soil domain
    nSnowOnlyHyd            => indx_data%var(iLookINDEX%nSnowOnlyHyd )%dat(1)         ,& ! intent(in): [i4b]    number of hydrology variables in the snow domain
    nSoilOnlyHyd            => indx_data%var(iLookINDEX%nSoilOnlyHyd )%dat(1)         ,& ! intent(in): [i4b]    number of hydrology variables in the soil domain
  ! soil parameters
    theta_sat               => mpar_data%var(iLookPARAM%theta_sat)%dat                ,& ! intent(in): [dp(:)]  soil porosity (-)
    theta_res               => mpar_data%var(iLookPARAM%theta_res)%dat                ,& ! intent(in): [dp(:)]  residual volumetric water content (-)
    ! state variables at the start of the time step
    mLayerMatricHead        => prog_data%var(iLookPROG%mLayerMatricHead)%dat          ,& ! intent(in): [dp(:)] matric head (m)
    mLayerVolFracIce        => prog_data%var(iLookPROG%mLayerVolFracIce)%dat           & ! intent(in): [dp(:)] volumetric fraction of ice (-)
    ) ! associating variables with indices of model state variables
    ! -----------------------------------------------------------------------------------------------------

    ! initialize error control
    err=0; message='imposeConstraints/'

    ! calculate proposed increment in state vector
    xInc(1:nState) = stateVec(1:nState)*1._qp - stateVecPrev(1:nState)*1._qp
  
    ! identify which constraints to impose
    select case(ixNumericalMethod)
    case(ida); err=20; message=trim(message)//'should not be imposing constraints for IDA solver'; return
      case(kinsol)
        small_delTemp       = .true.      ! flag to constain temperature change to be less than zMaxTempIncrement
        zMaxTempIncrement   = 10._rkind   ! maximum temperature increment (K)
        small_delMatric     = .true.      ! flag to constain matric head change to be less than zMaxMatricIncrement
        zMaxMatricIncrement = 10._rkind   ! maximum matric head increment (m)
        detect_events       = .true.      ! flag to do freezing point event detection and cross-over with epsT, works best if on
        epsT                = 1.e-3_rkind ! small interval above/below critical (K), works better if larger
        water_bounds        = .true.      ! flag to force water bounds, works best if on
      case(numrec)
        small_delTemp       = .true.      ! flag to constain temperature change to be less than zMaxTempIncrement
        zMaxTempIncrement   = 10._rkind   ! maximum temperature increment (K)
        small_delMatric     = .true.      ! flag to constain matric head change to be less than zMaxMatricIncrement
        zMaxMatricIncrement = 10._rkind   ! maximum matric head increment (m)
        detect_events       = .true.      ! flag to do freezing point event detection and cross-over with epsT
        epsT                = 1.e-7_rkind ! small interval above/below critical (K)
        water_bounds        = .true.      ! flag to force water bounds
      case default; err=20; message=trim(message)//'expect num_method to be ida, kinsol, or numrec (or itertive, which is numrec)'; return
    end select
   
    ! ** limit temperature increment to zMaxTempIncrement
    ! NOTE: this can cause problems especially from a cold start when far from the solution
    if(small_delTemp)then
      if(any(abs(xInc(ixNrgOnly)) > zMaxTempIncrement))then
        iMax       = maxloc( abs(xInc(ixNrgOnly)) )                     ! index of maximum temperature increment
        xIncFactor = abs( zMaxTempIncrement/xInc(ixNrgOnly(iMax(1))) )  ! scaling factor for the iteration increment (-)
        xInc       = xIncFactor*xInc
      endif
    endif ! (small temperature change)

    ! ** limit soil water (matric head) increment to zMaxMatricIncrement if starting positive
    if(small_delMatric)then
      if(size(ixMatOnly)>0)then
        ! loop through soil layers
        do iState=1,size(ixMatOnly)
          ! define index of the hydrology state variable within the state subset
          ixLiq = ixMatOnly(iState)
          ! place constraint for matric head
          if(xInc(ixLiq) > zMaxMatricIncrement .and. stateVecPrev(ixLiq) > 0._rkind)then
            xInc(ixLiq) = zMaxMatricIncrement
          endif
        end do ! (loop through soil layers)
      endif
    endif ! (small matric head change)

    ! ** stop just above or just below the freezing point if crossing
    if(detect_events)then

      ! crossing freezing point event for vegetation
      if(ixVegNrg/=integerMissing)then
        ! initialize
        critDiff    = Tfreeze - stateVecPrev(ixVegNrg)
        crosTempVeg = .false.
        ! initially frozen (T < Tfreeze)
        if(critDiff > 0._rkind)then
          if(xInc(ixVegNrg) > critDiff)then
            crosTempVeg = .true.
            cInc        = critDiff + epsT  ! constrained temperature increment (K)
          end if
        ! initially unfrozen (T > Tfreeze)
        else
          if(xInc(ixVegNrg) < critDiff)then
            crosTempVeg = .true.
            cInc        = critDiff - epsT  ! constrained temperature increment (K)
          end if
        end if  ! switch between frozen and unfrozen
        ! scale iterations
        if(crosTempVeg)then
          xIncFactor  = cInc/xInc(ixVegNrg) ! scaling factor for the iteration increment (-)
          xInc        = xIncFactor*xInc     ! scale iteration increments
        endif
      endif  ! if the state variable for canopy temperature is included within the state subset

      ! crossing freezing point event for snow
      if(nSnowOnlyNrg > 0)then
        do iLayer=1,nSnow 
          ! check if energy state is included
          if(ixSnowOnlyNrg(iLayer)==integerMissing) cycle
          ! check temperatures, and, if necessary, scale iteration increment
          iState = ixSnowOnlyNrg(iLayer)
          if(stateVecPrev(iState) + xInc(iState) > Tfreeze)then
          ! scale iteration increment
            cInc       = 0.5_rkind*(Tfreeze - stateVecPrev(iState) ) ! constrained temperature increment (K) -- simplified bi-section
            xIncFactor = cInc/xInc(iState)                           ! scaling factor for the iteration increment (-)
            xInc       = xIncFactor*xInc
          endif  ! (if snow temperature > freezing)
        end do ! (loop through snow layers)
      endif  ! (if there are state variables for energy in the snow domain)

      ! crossing freezing point event for soil
      if(nSoilOnlyNrg>0)then
        do iLayer=1,nSoil
          ! check if energy state is included
          if(ixSoilOnlyNrg(iLayer)==integerMissing) cycle
          ! define index of the state variables within the state subset
          ixNrg = ixSoilOnlyNrg(iLayer)
          ixLiq = ixSoilOnlyHyd(iLayer)
          ! get the matric potential of total water
          if(ixLiq/=integerMissing)then
            xPsi00 = stateVecPrev(ixLiq) + xInc(ixLiq)
          else
            xPsi00 = mLayerMatricHead(iLayer)
          endif
          ! identify the critical point when soil begins to freeze (TcSoil)
          TcSoil = crit_soilT(xPsi00)
          ! get the difference from the current state and the crossing point (K)
          critDiff = TcSoil - stateVecPrev(ixNrg)
          ! initially frozen (T < TcSoil)
          if(critDiff > 0._rkind)then
            ! (check crossing above zero)
            if(xInc(ixNrg) > critDiff)then
              crosFlag(iLayer) = .true.
              xInc(ixNrg) = critDiff + epsT  ! set iteration increment to slightly above critical temperature
            endif
          ! initially unfrozen (T > TcSoil)
          else 
            ! (check crossing below zero)
            if(xInc(ixNrg) < critDiff)then
              crosFlag(iLayer) = .true.
              xInc(ixNrg) = critDiff - epsT  ! set iteration increment to slightly below critical temperature
            endif
          endif ! (switch between initially frozen and initially unfrozen)
        end do ! (loop through soil layers)
      endif ! (if there are both energy and liquid water state variables)

    endif ! (detect events)

    ! ** ensure water is within bounds
    if(water_bounds)then

      ! impose positivity for canopy liquid water
      if(ixVegHyd/=integerMissing)then
        ! check if new value of storage will be negative
        if(stateVecPrev(ixVegHyd)+xInc(ixVegHyd) < 0._rkind)then
          ! scale iteration increment
          cInc       = -0.5_rkind*stateVecPrev(ixVegHyd)  ! constrained iteration increment (K) -- simplified bi-section
          xIncFactor = cInc/xInc(ixVegHyd)                ! scaling factor for the iteration increment (-)
          xInc       = xIncFactor*xInc                    ! new iteration increment
        end if
      endif ! (if the state variable for canopy water is included within the state subset)
          
      ! impose bounds for snow water, change in total water is only due to liquid flux
      if(nSnowOnlyHyd>0)then
        ! loop through snow layers
        do iLayer=1,nSnow
          ! check if the layer is included
          if(ixSnowOnlyHyd(iLayer)==integerMissing) cycle
          if(ixSnowOnlyNrg(iLayer)/=integerMissing)then
            ! get the layer temperature (from stateVecPrev if ixSnowOnlyNrg(iLayer) is within the state vector
            scalarTemp = stateVecPrev( ixSnowOnlyNrg(iLayer) )
          else ! get the layer temperature from the last update
            scalarTemp = prog_data%var(iLookPROG%mLayerTemp)%dat(iLayer)
          endif
          ! get the volumetric fraction of liquid water and ice
          select case( ixStateType_subset( ixSnowOnlyHyd(iLayer) ) )
            case(iname_watLayer); volFracLiq = fracliquid(scalarTemp,mpar_data%var(iLookPARAM%snowfrz_scale)%dat(1)) * stateVecPrev(ixSnowOnlyHyd(iLayer))
            case(iname_liqLayer); volFracLiq = stateVecPrev(ixSnowOnlyHyd(iLayer))
            case default; err=20; message=trim(message)//'expect ixStateType_subset to be iname_watLayer or iname_liqLayer for snow hydrology'; return
          end select
          scalarIce = mLayerVolFracIce(iLayer)
          ! checking if drain more than what is available or add more than possible
          if(-xInc(ixSnowOnlyHyd(iLayer)) > volFracLiq)then 
            xInc(ixSnowOnlyHyd(iLayer)) = -0.5_rkind*volFracLiq
          elseif(xInc(ixSnowOnlyHyd(iLayer)) > 1._rkind - scalarIce - volFracLiq)then
            xInc(ixSnowOnlyHyd(iLayer)) = 0.5_rkind*(1._rkind - scalarIce - volFracLiq)
          endif
        end do ! (looping through snow layers)
      endif ! (if there are state variables for liquid water in the snow domain)
       
      ! impose bounds for soil water, change in total water is only due to liquid flux
      if(nSoilOnlyHyd>0)then
        ! loop through soil layers
        do iLayer=1,nSoil
          ! check if the layer is included
          if(ixSoilOnlyHyd(iLayer)==integerMissing) cycle
          if(ixHydType(iLayer+nSnow)==iname_watLayer .or. ixHydType(iLayer+nSnow)==iname_liqLayer)then
            ! get the volumetric fraction of liquid water
            volFracLiq = stateVecPrev(ixSoilOnlyHyd(iLayer))
            scalarIce = merge(0._rkind, mLayerVolFracIce(iLayer+nSnow), ixHydType(iLayer+nSnow)==iname_watLayer)
            ! checking if drain more than what is available or add more than possible
            if(-xInc(ixSoilOnlyHyd(iLayer)) > volFracLiq - theta_res(iLayer))then 
              xInc(ixSoilOnlyHyd(iLayer)) = -0.5_rkind*(volFracLiq - theta_res(iLayer))
            elseif(xInc(ixSoilOnlyHyd(iLayer)) > theta_sat(iLayer) - scalarIce - volFracLiq)then
              xInc(ixSoilOnlyHyd(iLayer)) = -0.5_rkind*(theta_sat(iLayer) - scalarIce - volFracLiq)
            endif
          endif ! (if the state variable is not matric head)
        end do ! (looping through soil layers)
      endif ! (if there are state variables for liquid water in the soil domain)

    endif ! (water bounds)  

    ! Update the state vector with the modified iteration increment
    stateVec(:) = stateVecPrev(:) + xInc(:)

  ! end association with variables with indices of model state variables
  end associate

end subroutine imposeConstraints


end module eval8summa_module

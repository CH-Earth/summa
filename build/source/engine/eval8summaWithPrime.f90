
module eval8summaWithPrime_module

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
                    zLookup,      & ! lookup tables
                    model_options   ! defines the model decisions

! indices that define elements of the data structures
USE var_lookup,only:iLookDECISIONS               ! named variables for elements of the decision structure
USE var_lookup,only:iLookPARAM                   ! named variables for structure elements
USE var_lookup,only:iLookPROG                    ! named variables for structure elements
USE var_lookup,only:iLookINDEX                   ! named variables for structure elements
USE var_lookup,only:iLookDIAG                    ! named variables for structure elements
USE var_lookup,only:iLookFLUX                    ! named variables for structure elements
USE var_lookup,only:iLookDERIV                   ! named variables for structure elements

! look-up values for the choice of variable in energy equations (BE residual or IDA state variable)
USE mDecisions_module,only:  &
 closedForm,                 & ! use temperature with closed form heat capacity
 enthalpyFormLU,             & ! use enthalpy with soil temperature-enthalpy lookup tables
 enthalpyForm                  ! use enthalpy with soil temperature-enthalpy analytical solution

implicit none
private
public::eval8summaWithPrime
public::eval8summa4ida


contains

! **********************************************************************************************************
! public subroutine eval8summaWithPrime: compute the residual vector
! **********************************************************************************************************
subroutine eval8summaWithPrime(&
                      ! input: model control
                      dt,                            & ! intent(in):    entire time step for drainage pond rate
                      nSnow,                         & ! intent(in):    number of snow layers
                      nSoil,                         & ! intent(in):    number of soil layers
                      nLayers,                       & ! intent(in):    total number of layers
                      nState,                        & ! intent(in):    total number of state variables
                      insideSUN,                     & ! intent(in):    flag to indicate if we are inside Sundials solver
                      firstSubStep,                  & ! intent(in):    flag to indicate if we are processing the first sub-step
                      firstFluxCall,                 & ! intent(inout): flag to indicate if we are processing the first flux call
                      firstSplitOper,                & ! intent(inout): flag to indicate if we are processing the first flux call in a splitting operation
                      computeVegFlux,                & ! intent(in):    flag to indicate if we need to compute fluxes over vegetation
                      scalarSolution,                & ! intent(in):    flag to indicate the scalar solution
                      ! input: state vectors        
                      stateVec,                      & ! intent(in):    model state vector
                      stateVecPrime,                 & ! intent(in):    derivative of model state vector
                      sMul,                          & ! intent(inout): state vector multiplier (used in the residual calculations)
                      ! input: data structures
                      model_decisions,               & ! intent(in):    model decisions
                      lookup_data,                   & ! intent(in):    lookup table data structure
                      type_data,                     & ! intent(in):    type of vegetation and soil
                      attr_data,                     & ! intent(in):    spatial attributes
                      mpar_data,                     & ! intent(in):    model parameters
                      forc_data,                     & ! intent(in):    model forcing data
                      bvar_data,                     & ! intent(in):    average model variables for the entire basin
                      prog_data,                     & ! intent(in):    model prognostic variables for a local HRU
                      ! input-output: data stuctures
                      indx_data,                     & ! intent(inout): index data
                      diag_data,                     & ! intent(inout): model diagnostic variables for a local HRU
                      flux_data,                     & ! intent(inout): model fluxes for a local HRU
                      deriv_data,                    & ! intent(inout): derivatives in model fluxes w.r.t. relevant state variables
                      ! input-output: values needed in case canopy gets buried
                      scalarCanopyEnthalpyTrial,     & ! intent(inout): trial value for enthalpy of the vegetation canopy (J m-3)
                      scalarCanopyTempTrial,         & ! intent(inout): trial value for temperature of the vegetation canopy (K), also used to start enthalpy calculations
                      scalarCanopyWatTrial,          & ! intent(inout): trial value for total water content of the vegetation canopy (kg m-2)
                      ! output: new values of variables needed in data window outside of internal IDA for rootfinding and to start enthalpy calculations
                      mLayerTempTrial,               & ! intent(inout): trial vector of layer temperature (K)
                      mLayerMatricHeadTrial,         & ! intent(out):   trial value for total water matric potential (m)
                      ! output: new prime values of variables needed in data window outside of internal IDA for Jacobian
                      scalarCanopyTempPrime,         & ! intent(out):   prime value for temperature of the vegetation canopy (K s-1)
                      scalarCanopyWatPrime,          & ! intent(out):   prime value for total water content of the vegetation canopy (kg m-2 s-1)
                      mLayerTempPrime,               & ! intent(out):   prime vector of temperature of each snow and soil layer (K s-1)
                      mLayerMatricHeadPrime,         & ! intent(out):   prime vector of matric head of each snow and soil layer (m s-1)
                      mLayerVolFracWatPrime,         & ! intent(out):   prime vector of volumetric total water content of each snow and soil layer (s-1)
                      ! input-output: baseflow    
                      ixSaturation,                  & ! intent(inout): index of the lowest saturated layer
                      dBaseflow_dMatric,             & ! intent(out):   derivative in baseflow w.r.t. matric head (s-1)
                      ! output: flux and residual vectors
                      feasible,                      & ! intent(out):   flag to denote the feasibility of the solution
                      fluxVec,                       & ! intent(out):   flux vector
                      resSink,                       & ! intent(out):   sink terms on the RHS of the flux equation
                      resVec,                        & ! intent(out):   residual vector
                      err,message)                     ! intent(out):   error control
  ! --------------------------------------------------------------------------------------------------------------------------------
  ! provide access to subroutines
  USE getVectorz_module, only:varExtract                    ! extract variables from the state vector
  USE getVectorz_module, only:checkFeas                     ! check feasibility of state vector
  USE updateVarsWithPrime_module, only:updateVarsWithPrime  ! update variables
  USE computFlux_module, only:soilCmpresPrime               ! compute soil compression
  USE computFlux_module, only:computFlux                    ! compute fluxes given a state vector
  USE computHeatCap_module,only:computHeatCapAnalytic       ! recompute closed form heat capacity (Cp) and derivatives
  USE computHeatCap_module,only:computCm                    ! compute Cm and derivatives
  USE computHeatCap_module, only:computStatMult             ! recompute state multiplier
  USE computResidWithPrime_module,only:computResidWithPrime ! compute residuals given a state vector
  USE computThermConduct_module,only:computThermConduct     ! recompute thermal conductivity and derivatives
  implicit none
  ! --------------------------------------------------------------------------------------------------------------------------------
  ! --------------------------------------------------------------------------------------------------------------------------------
  ! input: model control
  real(rkind),intent(in)          :: dt                          ! entire time step for drainage pond rate
  integer(i4b),intent(in)         :: nSnow                       ! number of snow layers
  integer(i4b),intent(in)         :: nSoil                       ! number of soil layers
  integer(i4b),intent(in)         :: nLayers                     ! total number of layers
  integer,intent(in)              :: nState                      ! total number of state variables
  logical(lgt),intent(in)         :: insideSUN                   ! flag to indicate if we are inside Sundials solver
  logical(lgt),intent(in)         :: firstSubStep                ! flag to indicate if we are processing the first sub-step
  logical(lgt),intent(inout)      :: firstFluxCall               ! flag to indicate if we are processing the first flux call
  logical(lgt),intent(inout)      :: firstSplitOper              ! flag to indicate if we are processing the first flux call in a splitting operation
  logical(lgt),intent(in)         :: computeVegFlux              ! flag to indicate if computing fluxes over vegetation
  logical(lgt),intent(in)         :: scalarSolution              ! flag to denote if implementing the scalar solution
  ! input: state vectors    
  real(rkind),intent(in)          :: stateVec(:)                 ! model state vector
  real(rkind),intent(in)          :: stateVecPrime(:)            ! model state vector
  real(qp),intent(inout)          :: sMul(:)   ! NOTE: qp        ! state vector multiplier (used in the residual calculations)
  ! input: data structures
  type(model_options),intent(in)  :: model_decisions(:)          ! model decisions
  type(zLookup),      intent(in)  :: lookup_data                 ! lookup tables
  type(var_i),        intent(in)  :: type_data                   ! type of vegetation and soil
  type(var_d),        intent(in)  :: attr_data                   ! spatial attributes
  type(var_dlength),  intent(in)  :: mpar_data                   ! model parameters
  type(var_d),        intent(in)  :: forc_data                   ! model forcing data
  type(var_dlength),  intent(in)  :: bvar_data                   ! model variables for the local basin
  type(var_dlength),  intent(in)  :: prog_data                   ! prognostic variables for a local HRU
   ! output: data structures
  type(var_ilength),intent(inout) :: indx_data                   ! indices defining model states and layers
  type(var_dlength),intent(inout) :: diag_data                   ! diagnostic variables for a local HRU
  type(var_dlength),intent(inout) :: flux_data                   ! model fluxes for a local HRU
  type(var_dlength),intent(inout) :: deriv_data                  ! derivatives in model fluxes w.r.t. relevant state variables
  ! input-output: values needed in case canopy gets buried
  real(rkind),intent(inout)       :: scalarCanopyEnthalpyTrial   ! trial value for enthalpy of the vegetation canopy (J m-3)
  real(rkind),intent(inout)       :: scalarCanopyTempTrial       ! trial value for temperature of the vegetation canopy (K)
  real(rkind),intent(inout)       :: scalarCanopyWatTrial        ! trial value for total water content of the vegetation canopy (kg m-2), also used to start enthalpy calculations
  ! output: new values of variables needed in data window outside of internal IDA for rootfinding and to start enthalpy calculations
  real(rkind),intent(inout)       :: mLayerTempTrial(:)          ! trial vector of layer temperature (K)
  real(rkind),intent(out)         :: mLayerMatricHeadTrial(:)    ! trial vector for total water matric potential (m)
  ! output: new prime values of variables needed in data window outside of internal IDA for Jacobian
  real(rkind),intent(out)         :: scalarCanopyTempPrime       ! prime value for temperature of the vegetation canopy (K s-1)
  real(rkind),intent(out)         :: scalarCanopyWatPrime        ! prime value for total water content of the vegetation canopy (kg m-2 s-1)
  real(rkind),intent(out)         :: mLayerTempPrime(:)          ! prime vector of temperature of each snow and soil layer (K s-1)
  real(rkind),intent(out)         :: mLayerMatricHeadPrime(:)    ! prime vector of matric head of each snow and soil layer (m s-1)
  real(rkind),intent(out)         :: mLayerVolFracWatPrime(:)    ! prime vector of volumetric total water content of each snow and soil layer (s-1)
  ! input-output: baseflow    
  integer(i4b),intent(inout)      :: ixSaturation                ! index of the lowest saturated layer
  real(rkind),intent(out)         :: dBaseflow_dMatric(:,:)      ! derivative in baseflow w.r.t. matric head (s-1)
  ! output: flux and residual vectors
  logical(lgt),intent(out)        :: feasible                    ! flag to denote the feasibility of the solution
  real(rkind),intent(out)         :: fluxVec(:)                  ! flux vector
  real(rkind),intent(out)         :: resSink(:)                  ! sink terms on the RHS of the flux equation
  real(qp),intent(out)            :: resVec(:) ! NOTE: qp        ! residual vector
  ! output: error control
  integer(i4b),intent(out)        :: err                         ! error code
  character(*),intent(out)        :: message                     ! error message
  ! --------------------------------------------------------------------------------------------------------------------------------
  ! local variables
  ! --------------------------------------------------------------------------------------------------------------------------------
  real(rkind)                     :: dt1                         ! residual step size
  ! state variables
  real(rkind)                     :: scalarCanairEnthalpyTrial   ! trial value for enthalpy of the canopy air space (J m-3)
  real(rkind)                     :: scalarCanairTempTrial       ! trial value for temperature of the canopy air space (K)
  real(rkind)                     :: scalarCanopyLiqTrial        ! trial value for liquid water storage in the canopy (kg m-2)
  real(rkind)                     :: scalarCanopyIceTrial        ! trial value for ice storage in the canopy (kg m-2)
  real(rkind),dimension(nSoil)    :: mLayerMatricHeadLiqTrial    ! trial value for liquid water matric potential (m)
  real(rkind),dimension(nLayers)  :: mLayerEnthalpyTrial         ! trial vector of enthalpy of each snow and soil layer (J m-3)
  real(rkind),dimension(nLayers)  :: mLayerVolFracWatTrial       ! trial vector of volumetric total water content (-)
  real(rkind),dimension(nLayers)  :: mLayerVolFracLiqTrial       ! trial vector of volumetric liquid water content (-)
  real(rkind),dimension(nLayers)  :: mLayerVolFracIceTrial       ! trial vector of volumetric fraction of ice (-)
  real(rkind)                     :: scalarAquiferStorageTrial   ! trial value for storage of water in the aquifer (m)
  ! prime state variables
  real(rkind)                     :: scalarCanairEnthalpyPrime   ! prime value for enthalpy of the canopy air space (W m-3)
  real(rkind)                     :: scalarCanairTempPrime       ! prime value for temperature of the canopy air space (K s-1)
  real(rkind)                     :: scalarCanopyEnthalpyPrime   ! prime value for enthalpy of the vegetation canopy (W m-3)
  real(rkind)                     :: scalarCanopyLiqPrime        ! prime value for liquid water storage in the canopy (kg m-2 s-1)
  real(rkind)                     :: scalarCanopyIcePrime        ! prime value for mass of ice on the vegetation canopy (kg m-2 s-1)
  real(rkind),dimension(nLayers)  :: mLayerEnthalpyPrime         ! prime vector of enthalpy of each snow and soil layer (W m-3)
  real(rkind),dimension(nLayers)  :: mLayerVolFracLiqPrime       ! prime vector of volumetric liquid water content (s-1)
  real(rkind),dimension(nLayers)  :: mLayerVolFracIcePrime       ! prime vector of volumetric fraction of ice (s-1)
  real(rkind),dimension(nSoil)    :: mLayerMatricHeadLiqPrime    ! prime vector of liquid water matric potential (m s-1)
  real(rkind)                     :: scalarAquiferStoragePrime   ! prime value of storage of water in the aquifer (m s-1)
  ! dummy state variables
  real(rkind)                     :: scalarCanairNrgTrial        ! trial value for energy of the canopy air space
  real(rkind)                     :: scalarCanopyNrgTrial        ! trial value for energy of the vegetation canopy
  real(rkind),dimension(nLayers)  :: mLayerNrgTrial              ! trial vector of energy of each snow and soil layer
  real(rkind)                     :: scalarCanairNrgPrime        ! prime value for energy of the canopy air space
  real(rkind)                     :: scalarCanopyNrgPrime        ! prime value for energy of the vegetation canopy
  real(rkind),dimension(nLayers)  :: mLayerNrgPrime              ! prime vector of energy of each snow and soil layer
  ! other local variables
  integer(i4b)                    :: iLayer                      ! index of model layer in the snow+soil domain
  integer(i4b)                    :: jState(1)                   ! index of model state for the scalar solution within the soil domain
  integer(i4b)                    :: ixBeg,ixEnd                 ! index of indices for the soil compression routine
  character(LEN=256)              :: cmessage                    ! error message of downwind routine
  logical(lgt)                    :: updateStateCp               ! flag to indicate if we update Cp at each step for LHS, set with nrgConserv choice and updateCp_closedForm flag
  logical(lgt)                    :: updateFluxCp                ! flag to indicate if we update Cp at each step for RHS, set with nrgConserv choice and updateCp_closedForm flag
  logical(lgt)                    :: needStateCm                 ! flag to indicate if the energy equation contains LHS Cm = dH_T/dTheta_m,, set with nrgConserv choice and needStateCm_closedForm flag
  logical(lgt),parameter          :: updateCp_closedForm=.true. ! nrgConserv = closedForm flag to indicate if we update Cp at each step
  logical(lgt),parameter          :: needCm_closedForm=.true.   ! nrgConserv = closedForm flag to indicate if the energy equation contains Cm = dH_T/dTheta_m

  ! --------------------------------------------------------------------------------------------------------------------------------
  ! association to variables in the data structures
  ! --------------------------------------------------------------------------------------------------------------------------------
  associate(&
    ! model decisions
    ixNrgConserv              => model_decisions(iLookDECISIONS%nrgConserv)%iDecision      ,& ! intent(in):  [i4b]    choice of variable in either energy backward Euler residual or IDA state variable
    ixRichards                => model_decisions(iLookDECISIONS%f_Richards)%iDecision      ,& ! intent(in):  [i4b]    index of the form of Richards' equation
    ! snow parameters
    snowfrz_scale             => mpar_data%var(iLookPARAM%snowfrz_scale)%dat(1)            ,& ! intent(in):  [dp]     scaling parameter for the snow freezing curve (K-1)
    ! soil parameters
    theta_sat                 => mpar_data%var(iLookPARAM%theta_sat)%dat                   ,& ! intent(in):  [dp(:)]  soil porosity (-)
    specificStorage           => mpar_data%var(iLookPARAM%specificStorage)%dat(1)          ,& ! intent(in):  [dp]     specific storage coefficient (m-1)
    ! canopy and layer depth
    canopyDepth               => diag_data%var(iLookDIAG%scalarCanopyDepth)%dat(1)         ,& ! intent(in):  [dp   ]  canopy depth (m)
    mLayerDepth               => prog_data%var(iLookPROG%mLayerDepth)%dat                  ,& ! intent(in):  [dp(:)]  depth of each layer in the snow-soil sub-domain (m)
    ! model diagnostic variables, will be updated before used
    scalarFracLiqVeg          => diag_data%var(iLookDIAG%scalarFracLiqVeg)%dat(1)          ,& ! intent(in):  [dp]     fraction of liquid water on vegetation (-)
    scalarSfcMeltPond         => prog_data%var(iLookPROG%scalarSfcMeltPond)%dat(1)         ,& ! intent(in):  [dp]     ponded water caused by melt of the "snow without a layer" (kg m-2)
    mLayerFracLiqSnow         => diag_data%var(iLookDIAG%mLayerFracLiqSnow)%dat            ,& ! intent(in):  [dp(:)]  fraction of liquid water in each snow layer (-)
    ! soil compression
    scalarSoilCompress        => diag_data%var(iLookDIAG%scalarSoilCompress)%dat(1)        ,& ! intent(in):  [dp]     total change in storage associated with compression of the soil matrix (kg m-2 s-1)
    mLayerCompress            => diag_data%var(iLookDIAG%mLayerCompress)%dat               ,& ! intent(in):  [dp(:)]  change in volumetric water content due to compression of soil (s-1)
    ! derivatives
    dTheta_dTkCanopy          => deriv_data%var(iLookDERIV%dTheta_dTkCanopy)%dat(1)        ,& ! intent(in):  [dp]     derivative of volumetric liquid water content w.r.t. temperature
    dVolTot_dPsi0             => deriv_data%var(iLookDERIV%dVolTot_dPsi0)%dat              ,& ! intent(in):  [dp(:)]  derivative in total water content w.r.t. total water matric potential
    dCompress_dPsi            => deriv_data%var(iLookDERIV%dCompress_dPsi)%dat             ,& ! intent(in):  [dp(:)]  derivative in compressibility w.r.t. matric head (m-1)
    mLayerdTheta_dTk          => deriv_data%var(iLookDERIV%mLayerdTheta_dTk)%dat           ,& ! intent(in):  [dp(:)]  derivative of volumetric liquid water content w.r.t. temperature
    dVolHtCapBulk_dPsi0       => deriv_data%var(iLookDERIV%dVolHtCapBulk_dPsi0)%dat        ,& ! intent(out): [dp(:)]  derivative in bulk heat capacity w.r.t. matric potential
    dVolHtCapBulk_dTheta      => deriv_data%var(iLookDERIV%dVolHtCapBulk_dTheta)%dat       ,& ! intent(out): [dp(:)]  derivative in bulk heat capacity w.r.t. volumetric water content
    dVolHtCapBulk_dCanWat     => deriv_data%var(iLookDERIV%dVolHtCapBulk_dCanWat)%dat(1)   ,& ! intent(out): [dp]     derivative in bulk heat capacity w.r.t. volumetric water content
    dVolHtCapBulk_dTk         => deriv_data%var(iLookDERIV%dVolHtCapBulk_dTk)%dat          ,& ! intent(out): [dp(:)]  derivative in bulk heat capacity w.r.t. temperature
    dVolHtCapBulk_dTkCanopy   => deriv_data%var(iLookDERIV%dVolHtCapBulk_dTkCanopy)%dat(1) ,& ! intent(out): [dp]     derivative in bulk heat capacity w.r.t. temperature
    dThermalC_dWatAbove       => deriv_data%var(iLookDERIV%dThermalC_dWatAbove)%dat        ,& ! intent(out): [dp(:)]  derivative in the thermal conductivity w.r.t. water state in the layer above
    dThermalC_dWatBelow       => deriv_data%var(iLookDERIV%dThermalC_dWatBelow)%dat        ,& ! intent(out): [dp(:)]  derivative in the thermal conductivity w.r.t. water state in the layer above
    dThermalC_dTempAbove      => deriv_data%var(iLookDERIV%dThermalC_dTempAbove)%dat       ,& ! intent(out): [dp(:)]  derivative in the thermal conductivity w.r.t. energy state in the layer above
    dThermalC_dTempBelow      => deriv_data%var(iLookDERIV%dThermalC_dTempBelow)%dat       ,& ! intent(out): [dp(:)]  derivative in the thermal conductivity w.r.t. energy state in the layer above
    dCm_dTk                   => deriv_data%var(iLookDERIV%dCm_dTk)%dat                    ,& ! intent(out): [dp(:)]  derivative in heat capacity w.r.t. temperature (J kg-1 K-2)
    dCm_dTkCanopy             => deriv_data%var(iLookDERIV%dCm_dTkCanopy)%dat(1)           ,& ! intent(out): [dp   ]  derivative in heat capacity w.r.t. canopy temperature (J kg-1 K-2)
    ! mapping
    ixMapFull2Subset          => indx_data%var(iLookINDEX%ixMapFull2Subset)%dat            ,& ! intent(in):  [i4b(:)] mapping of full state vector to the state subset
    ixControlVolume           => indx_data%var(iLookINDEX%ixControlVolume)%dat             ,& ! intent(in):  [i4b(:)] index of control volume for different domains (veg, snow, soil)
    ! heat capacity
    heatCapVegTrial           => diag_data%var(iLookDIAG%scalarBulkVolHeatCapVeg)%dat(1)   ,& ! intent(out): [dp]     volumetric heat capacity of vegetation canopy
    mLayerHeatCapTrial        => diag_data%var(iLookDIAG%mLayerVolHtCapBulk)%dat           ,& ! intent(out): [dp(:)]  heat capacity for snow and soil
    ! Cm
    canopyCmTrial             => diag_data%var(iLookDIAG%scalarCanopyCm)%dat(1)            ,& ! intent(out): [dp]     Cm of the canopy
    mLayerCmTrial             => diag_data%var(iLookDIAG%mLayerCm)%dat                      & ! intent(out): [dp(:)]  Cm of snow and soil
    ) ! association to variables in the data structures
    ! --------------------------------------------------------------------------------------------------------------------------------
    ! initialize error control
    err=0; message="eval8summaWithPrime/"

    ! check the feasibility of the solution only if not inside Sundials solver
    feasible=.true.
    if (.not.insideSUN) then
      call checkFeas(&
                    ! input
                    stateVec,                   & ! intent(in):    model state vector (mixed units)
                    mpar_data,                  & ! intent(in):    model parameters
                    prog_data,                  & ! intent(in):    model prognostic variables for a local HRU
                    indx_data,                  & ! intent(in):    indices defining model states and layers
                    ixNrgConserv.ne.closedForm, & ! intent(in):    flag to indicate if we are using enthalpy as state variable
                    ! output: feasibility
                    feasible,                   & ! intent(inout): flag to denote the feasibility of the solution
                  ! output: error control
                    err,cmessage)                 ! intent(out):   error control

      ! early return for non-feasible solutions
      if(.not.feasible)then
        fluxVec(:) = realMissing
        resVec(:)  = quadMissing
        message=trim(message)//trim(cmessage)//'non-feasible'
        return
      end if
    end if ! ( feasibility check )

    if(ixNrgConserv == enthalpyForm .or. ixNrgConserv == enthalpyFormLU)then ! use enthalpy as state variable, do not need state terms but do need flux term
      updateStateCp = .false.
      updateFluxCp  = .true.
      needStateCm   = .false.
    else if(ixNrgConserv == closedForm)then ! have a choice, temperature the state variable
      updateStateCp = updateCp_closedForm
      updateFluxCp  = updateCp_closedForm
      needStateCm   = needCm_closedForm
    else
      message=trim(message)//'unknown choice of variable in energy conservation backward Euler residual'
      err=1; return
    end if

    ! get the start and end indices for the soil compression calculations
    if(scalarSolution)then
      jState = pack(ixControlVolume, ixMapFull2Subset/=integerMissing)
      ixBeg  = jState(1)
      ixEnd  = jState(1)
    else
      ixBeg  = 1
      ixEnd  = nSoil
    endif

    ! Canopy layer can disappear even without splitting (snow burial), so need to take last values
    if(ixNrgConserv== enthalpyForm .or. ixNrgConserv == enthalpyFormLU)then ! use state variable as enthalpy, need to compute temperature
      scalarCanopyNrgTrial = scalarCanopyEnthalpyTrial
    else ! use state variable as temperature
      scalarCanopyNrgTrial = scalarCanopyTempTrial
    endif !(choice of how conservation of energy is implemented)

   ! Placeholder: if we decide to use splitting, we need to pass all the previous values of the state variables
    scalarCanairNrgTrial      = realMissing
    scalarCanopyLiqTrial      = realMissing
    scalarCanopyIceTrial      = realMissing
    mLayerNrgTrial            = realMissing
    mLayerVolFracWatTrial     = realMissing
    mLayerVolFracLiqTrial     = realMissing
    mLayerVolFracIceTrial     = realMissing
    mLayerMatricHeadTrial     = realMissing
    mLayerMatricHeadLiqTrial  = realMissing
    scalarAquiferStorageTrial = realMissing

    ! extract states from the state vector
    call varExtract(&
                    ! input
                    stateVec,                  & ! intent(in):    model state vector (mixed units)
                    diag_data,                 & ! intent(in):    model diagnostic variables for a local HRU
                    prog_data,                 & ! intent(in):    model prognostic variables for a local HRU
                    indx_data,                 & ! intent(in):    indices defining model states and layers
                    ! output: variables for the vegetation canopy
                    scalarCanairNrgTrial,      & ! intent(inout): trial value of energy of the canopy air space, temperature (K) or enthalpy (J m-3)
                    scalarCanopyNrgTrial,      & ! intent(inout): trial value of energy of the vegetation canopy, temperature (K) or enthalpy (J m-3)
                    scalarCanopyWatTrial,      & ! intent(inout): trial value of canopy total water (kg m-2)
                    scalarCanopyLiqTrial,      & ! intent(inout): trial value of canopy liquid water (kg m-2)
                    ! output: variables for the snow-soil domain
                    mLayerNrgTrial,            & ! intent(inout): trial vector of energy, temperature (K) or enthalpy (J m-3)
                    mLayerVolFracWatTrial,     & ! intent(inout): trial vector of volumetric total water content (-)
                    mLayerVolFracLiqTrial,     & ! intent(inout): trial vector of volumetric liquid water content (-)
                    mLayerMatricHeadTrial,     & ! intent(inout): trial vector of total water matric potential (m)
                    mLayerMatricHeadLiqTrial,  & ! intent(inout): trial vector of liquid water matric potential (m)
                    ! output: variables for the aquifer
                    scalarAquiferStorageTrial, & ! intent(inout):   trial value of storage of water in the aquifer (m)
                    ! output: error control
                    err,cmessage)               ! intent(out):   error control
    if(err/=0)then; message=trim(message)//trim(cmessage); return; end if  ! (check for errors)

    ! Placeholder: if we decide to use splitting, we need to pass all the previous values of the state variables
    scalarCanairNrgPrime      = realMissing
    scalarCanopyNrgPrime      = realMissing
    scalarCanopyWatPrime      = realMissing
    scalarCanopyLiqPrime      = realMissing
    scalarCanopyIcePrime      = realMissing
    mLayerNrgPrime            = realMissing
    mLayerVolFracWatPrime     = realMissing
    mLayerVolFracLiqPrime     = realMissing
    mLayerVolFracIcePrime     = realMissing
    mLayerMatricHeadPrime     = realMissing
    mLayerMatricHeadLiqPrime  = realMissing
    scalarAquiferStoragePrime = realMissing

    call varExtract(&
                  ! input
                  stateVecPrime,             & ! intent(in):    derivative of model state vector (mixed units)
                  diag_data,                 & ! intent(in):    model diagnostic variables for a local HRU
                  prog_data,                 & ! intent(in):    model prognostic variables for a local HRU
                  indx_data,                 & ! intent(in):    indices defining model states and layers
                  ! output: variables for the vegetation canopy
                  scalarCanairNrgPrime,      & ! intent(inout): derivative of energy of the canopy air space, temperature (K s-1) or enthalpy (W m-3)
                  scalarCanopyNrgPrime,      & ! intent(inout): derivative of energy of the vegetation canopy, temperature (K s-1) or enthalpy (W m-3)
                  scalarCanopyWatPrime,      & ! intent(inout): derivative of canopy total water (kg m-2 s-1)
                  scalarCanopyLiqPrime,      & ! intent(inout): derivative of canopy liquid water (kg m-2 s-1)
                  ! output: variables for the snow-soil domain
                  mLayerNrgPrime,            & ! intent(inout): derivative of energy of each snow and soil layer, temperature (K s-1) or enthalpy (W m-3)
                  mLayerVolFracWatPrime,     & ! intent(inout): derivative of volumetric total water content (s-1)
                  mLayerVolFracLiqPrime,     & ! intent(inout): derivative of volumetric liquid water content (s-1)
                  mLayerMatricHeadPrime,     & ! intent(inout): derivative of total water matric potential (m s-1)
                  mLayerMatricHeadLiqPrime,  & ! intent(inout): derivative of liquid water matric potential (m s-1)
                  ! output: variables for the aquifer
                  scalarAquiferStoragePrime, & ! intent(inout): derivative of storage of water in the aquifer (m s-1)
                  ! output: error control
                  err,cmessage)                ! intent(out):   error control
    if(err/=0)then; message=trim(message)//trim(cmessage); return; end if  ! (check for errors)

    if(ixNrgConserv== enthalpyForm .or. ixNrgConserv == enthalpyFormLU)then ! use state variable as enthalpy, need to compute temperature
      scalarCanairEnthalpyTrial = scalarCanairNrgTrial
      scalarCanopyEnthalpyTrial = scalarCanopyNrgTrial
      mLayerEnthalpyTrial       = mLayerNrgTrial
      scalarCanairEnthalpyPrime = scalarCanairNrgPrime
      scalarCanopyEnthalpyPrime = scalarCanopyNrgPrime
      mLayerEnthalpyPrime       = mLayerNrgPrime
      ! do not use these variables
      scalarCanairTempPrime = realMissing
      scalarCanopyTempPrime = realMissing
      mLayerTempPrime       = realMissing
    else ! use state variable as temperature
      scalarCanairTempTrial = scalarCanairNrgTrial
      scalarCanopyTempTrial = scalarCanopyNrgTrial
      mLayerTempTrial       = mLayerNrgTrial
      scalarCanairTempPrime = scalarCanairNrgPrime
      scalarCanopyTempPrime = scalarCanopyNrgPrime
      mLayerTempPrime       = mLayerNrgPrime
      ! do not use these variables
      scalarCanairEnthalpyTrial = realMissing
      scalarCanopyEnthalpyTrial = realMissing
      mLayerEnthalpyTrial       = realMissing
      scalarCanairEnthalpyPrime = realMissing
      scalarCanopyEnthalpyPrime = realMissing
      mLayerEnthalpyPrime       = realMissing     
    endif !(choice of how conservation of energy is implemented)

    ! update diagnostic variables and derivatives
    ! NOTE: if we are using enthalpy as a state variable, currently all *TempPrime, *IcePrime, and *LiqPrime are set to realMissing
    !       This possibly could cause problems (?) if we use splitting, but we are not using splitting at the moment
    call updateVarsWithPrime(&
                    ! input
                    ixNrgConserv.ne.closedForm,   & ! intent(in):    flag if need to update temperature from enthalpy
                    ixNrgConserv==enthalpyFormLU, & ! intent(in):    flag to use the lookup table for soil temperature-enthalpy
                    .true.,                       & ! intent(in):    flag if computing for Jacobian update
                    .false.,                      & ! intent(in):    flag to adjust temperature to account for the energy
                    mpar_data,                    & ! intent(in):    model parameters for a local HRU
                    indx_data,                    & ! intent(in):    indices defining model states and layers
                    prog_data,                    & ! intent(in):    model prognostic variables for a local HRU
                    diag_data,                    & ! intent(inout): model diagnostic variables for a local HRU
                    deriv_data,                   & ! intent(inout): derivatives in model fluxes w.r.t. relevant state variables
                    lookup_data,                  & ! intent(in):    lookup table data structure
                    ! input: enthalpy state variables  
                    scalarCanairEnthalpyTrial,    & ! intent(in):    trial value for enthalpy of the canopy air space (J m-3)
                    scalarCanopyEnthalpyTrial,    & ! intent(in):    trial value for enthalpy of the vegetation canopy (J m-3)
                    mLayerEnthalpyTrial,          & ! intent(in):    trial vector of enthalpy of each snow+soil layer (J m-3)                      
                    ! output: variables for the vegetation canopy
                    scalarCanairTempTrial,        & ! intent(inout): trial value of canopy air space temperature (K)
                    scalarCanopyTempTrial,        & ! intent(inout): trial value of canopy temperature (K)
                    scalarCanopyWatTrial,         & ! intent(inout): trial value of canopy total water (kg m-2)
                    scalarCanopyLiqTrial,         & ! intent(inout): trial value of canopy liquid water (kg m-2)
                    scalarCanopyIceTrial,         & ! intent(inout): trial value of canopy ice content (kg m-2)
                    scalarCanopyTempPrime,        & ! intent(inout): trial value of time derivative canopy temperature (K s-1)
                    scalarCanopyWatPrime,         & ! intent(inout): trial value of time derivative canopy total water (kg m-2 s-1)
                    scalarCanopyLiqPrime,         & ! intent(inout): trial value of time derivative canopy liquid water (kg m-2 s-1)
                    scalarCanopyIcePrime,         & ! intent(inout): trial value of time derivative canopy ice content (kg m-2 s-1)
                    ! output: variables for th snow-soil domain
                    mLayerTempTrial,              & ! intent(inout): trial vector of layer temperature (K)
                    mLayerVolFracWatTrial,        & ! intent(inout): trial vector of volumetric total water content (-)
                    mLayerVolFracLiqTrial,        & ! intent(inout): trial vector of volumetric liquid water content (-)
                    mLayerVolFracIceTrial,        & ! intent(inout): trial vector of volumetric ice water content (-)
                    mLayerMatricHeadTrial,        & ! intent(inout): trial vector of total water matric potential (m)
                    mLayerMatricHeadLiqTrial,     & ! intent(inout): trial vector of liquid water matric potential (m)
                    mLayerTempPrime,              & ! intent(inout): trial vector of time derivative layer temperature (K s-1)
                    mLayerVolFracWatPrime,        & ! intent(inout): trial vector of time derivative volumetric total water content (s-1)
                    mLayerVolFracLiqPrime,        & ! intent(inout): trial vector of time derivative volumetric liquid water content (s-1)
                    mLayerVolFracIcePrime,        & ! intent(inout): trial vector of time derivative volumetric ice water content (s-1)
                    mLayerMatricHeadPrime,        & ! intent(inout): trial vector of time derivative total water matric potential (m s-1)
                    mLayerMatricHeadLiqPrime,     & ! intent(inout): trial vector of time derivative liquid water matric potential (m s-1)
                    ! output: error control
                    err,cmessage)                   ! intent(out):   error control
    if(err/=0)then; message=trim(message)//trim(cmessage); return; end if  ! (check for errors)

    if(updateStateCp)then
      ! *** compute volumetric heat capacity C_p
      call computHeatCapAnalytic(&
                  ! input: state variables
                  canopyDepth,             & ! intent(in):    canopy depth (m)
                  scalarCanopyIceTrial,    & ! intent(in):    trial value for mass of ice on the vegetation canopy (kg m-2)
                  scalarCanopyLiqTrial,    & ! intent(in):    trial value for the liquid water on the vegetation canopy (kg m-2)
                  scalarCanopyTempTrial,   & ! intent(in):    trial value of canopy temperature (K)
                  mLayerVolFracIceTrial,   & ! intent(in):    volumetric fraction of ice at the start of the sub-step (-)
                  mLayerVolFracLiqTrial,   & ! intent(in):    fraction of liquid water at the start of the sub-step (-)
                  mLayerTempTrial,         & ! intent(in):    trial value of layer temperature (K)
                  mLayerMatricHeadTrial,   & ! intent(in):    trial total water matric potential (m)
                  ! input: pre-computed derivatives
                  dTheta_dTkCanopy,        & ! intent(in):    derivative in canopy volumetric liquid water content w.r.t. temperature (K-1)
                  scalarFracLiqVeg,        & ! intent(in):    fraction of canopy liquid water (-)
                  mLayerdTheta_dTk,        & ! intent(in):    derivative of volumetric liquid water content w.r.t. temperature (K-1)
                  mLayerFracLiqSnow,       & ! intent(in):    fraction of liquid water (-)
                  dVolTot_dPsi0,           & ! intent(in):    derivative in total water content w.r.t. total water matric potential (m-1)
                  ! input output data structures
                  mpar_data,               & ! intent(in):    model parameters
                  indx_data,               & ! intent(in):    model layer indices
                  ! output
                  heatCapVegTrial,         & ! intent(inout): volumetric heat capacity of vegetation canopy
                  mLayerHeatCapTrial,      & ! intent(inout): volumetric heat capacity of soil and snow
                  dVolHtCapBulk_dPsi0,     & ! intent(inout): derivative in bulk heat capacity w.r.t. matric potential
                  dVolHtCapBulk_dTheta,    & ! intent(inout): derivative in bulk heat capacity w.r.t. volumetric water content
                  dVolHtCapBulk_dCanWat,   & ! intent(inout): derivative in bulk heat capacity w.r.t. volumetric water content
                  dVolHtCapBulk_dTk,       & ! intent(inout): derivative in bulk heat capacity w.r.t. temperature
                  dVolHtCapBulk_dTkCanopy, & ! intent(inout): derivative in bulk heat capacity w.r.t. temperature                  
                  ! output: error control
                  err,cmessage)                  ! intent(out):  error control
      if(err/=0)then; message=trim(message)//trim(cmessage); return; endif

      ! compute multiplier of state vector
      call computStatMult(&
                    ! input
                    heatCapVegTrial,    & ! intent(in):  volumetric heat capacity of vegetation canopy
                    mLayerHeatCapTrial, & ! intent(in):  volumetric heat capacity of soil and snow
                    indx_data,          & ! intent(in):  indices defining model states and layers
                    ! output
                    sMul,               & ! intent(out): multiplier for state vector (used in the residual calculations)
                    err,cmessage)         ! intent(out): error control
      if(err/=0)then; message=trim(message)//trim(cmessage); return; endif  ! (check for errors)
    else
      ! set state heat capacity derivatives to 0 for constant through step
      dVolHtCapBulk_dPsi0     = 0._rkind
      dVolHtCapBulk_dTheta    = 0._rkind
      dVolHtCapBulk_dCanWat   = 0._rkind
      dVolHtCapBulk_dTk       = 0._rkind
      dVolHtCapBulk_dTkCanopy = 0._rkind
    endif ! updateStateCp

    if(updateFluxCp)then
      ! update thermal conductivity
      call computThermConduct(&
                          ! input: control variables
                          computeVegFlux,        & ! intent(in):    flag to denote if computing the vegetation flux
                          nLayers,               & ! intent(in):    total number of layers
                          canopyDepth,           & ! intent(in):    canopy depth (m)
                          ! input: state variables
                          scalarCanopyIceTrial,  & ! intent(in):    trial value for mass of ice on the vegetation canopy (kg m-2)
                          scalarCanopyLiqTrial,  & ! intent(in):    trial value of canopy liquid water (kg m-2)
                          mLayerTempTrial,       & ! intent(in):    trial temperature of layer temperature (K)
                          mLayerMatricHeadTrial, & ! intent(in):    trial value for total water matric potential (m)
                          mLayerVolFracIceTrial, & ! intent(in):    volumetric fraction of ice at the start of the sub-step (-)
                          mLayerVolFracLiqTrial, & ! intent(in):    volumetric fraction of liquid water at the start of the sub-step (-)
                         ! input: pre-computed derivatives
                          mLayerdTheta_dTk,      & ! intent(in):    derivative in volumetric liquid water content w.r.t. temperature (K-1)
                          mLayerFracLiqSnow,     & ! intent(in):    fraction of liquid water (-)
                          ! input/output: data structures
                          mpar_data,             & ! intent(in):    model parameters
                          indx_data,             & ! intent(in):    model layer indices
                          prog_data,             & ! intent(in):    model prognostic variables for a local HRU
                          diag_data,             & ! intent(inout): model diagnostic variables for a local HRU
                          ! output: derivative
                          dThermalC_dWatAbove,   & ! intent(out):   derivative in the thermal conductivity w.r.t. water state in the layer above
                          dThermalC_dWatBelow,   & ! intent(out):   derivative in the thermal conductivity w.r.t. water state in the layer above
                          dThermalC_dTempAbove,  & ! intent(out):   derivative in the thermal conductivity w.r.t. energy state in the layer above
                          dThermalC_dTempBelow,  & ! intent(out):   derivative in the thermal conductivity w.r.t. energy state in the layer above
                          ! output: error control
                          ! output: error control
                          err,cmessage)                   ! intent(out): error control
      if(err/=0)then; err=55; message=trim(message)//trim(cmessage); return; end if
    else
      ! set flux heat capacity derivatives to 0 for constant through step
      dThermalC_dWatAbove  = 0._rkind
      dThermalC_dWatBelow  = 0._rkind
      dThermalC_dTempAbove = 0._rkind
      dThermalC_dTempBelow = 0._rkind  
    endif ! updateFluxCp

    if(needStateCm)then
      ! compute C_m
      call computCm(&
                 ! input: state variables
                 scalarCanopyTempTrial, & ! intent(in):    trial value of canopy temperature (K)
                 mLayerTempTrial,       & ! intent(in):    trial value of layer temperature (K)
                 mLayerMatricHeadTrial, & ! intent(in):    trial value for total water matric potential (-)
                 ! input data structures
                 mpar_data,             & ! intent(in):    model parameters
                 indx_data,             & ! intent(in):    model layer indices
                 ! output
                 canopyCmTrial,         & ! intent(inout): Cm for vegetation (J kg K-1)
                 mLayerCmTrial,         & ! intent(inout): Cm for soil and snow (J kg K-1)
                 dCm_dTk,               & ! intent(inout): derivative in Cm w.r.t. temperature (J kg K-2)
                 dCm_dTkCanopy,         & ! intent(inout): derivative in Cm w.r.t. temperature (J kg K-2)
                 err,cmessage)            ! intent(inout): error control
    else
      canopyCmTrial = 0._qp
      mLayerCmTrial = 0._qp
      dCm_dTk       = 0._rkind
      dCm_dTkCanopy = 0._rkind
    endif ! needStateCm

    ! save the number of flux calls per time step
    indx_data%var(iLookINDEX%numberFluxCalc)%dat(1) = indx_data%var(iLookINDEX%numberFluxCalc)%dat(1) + 1

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
                    .false.,                   & ! intent(in):    do not check longwave balance
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
                    dBaseflow_dMatric,         & ! intent(out):   derivative in baseflow w.r.t. matric head (s-1), we will use it later in computJacobWithPrime
                    fluxVec,                   & ! intent(out):   flux vector (mixed units)
                    ! output: error control
                    err,cmessage)                ! intent(out):   error code and error message
    if(err/=0)then; message=trim(message)//trim(cmessage); return; end if  ! (check for errors)

    firstSplitOper = .false. ! after call computFlux once in dt, no longer firstSplitOper

    ! compute soil compressibility (-) and its derivative w.r.t. matric head (m)
    ! NOTE: we already extracted trial matrix head and volumetric liquid water as part of the flux calculations
    call soilCmpresPrime(&
                    ! input:
                    ixRichards,                             & ! intent(in):    choice of option for Richards' equation
                    ixBeg,ixEnd,                            & ! intent(in):    start and end indices defining desired layers
                    mLayerMatricHeadPrime(1:nSoil),         & ! intent(in):    matric head at the start of the time step (m s-1)
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
    if (insideSUN)then
      dt1 = 1._qp ! always 1 for IDA since using Prime derivatives

      call computResidWithPrime(&
                       ! input: model control
                      dt1,                        & ! intent(in):  length of the residual time step (seconds)
                      nSnow,                      & ! intent(in):  number of snow layers
                      nSoil,                      & ! intent(in):  number of soil layers
                      nLayers,                    & ! intent(in):  total number of layers
                      ixNrgConserv.ne.closedForm, & ! intent(in):  flag if enthalpy is state variable
                      ! input: flux vectors
                      sMul,                       & ! intent(in):  state vector multiplier (used in the residual calculations)
                      fluxVec,                    & ! intent(in):  flux vector
                      ! input: state variables (already disaggregated into scalars and vectors)
                      scalarCanairTempPrime,      & ! intent(in):  prime value for the temperature of the canopy air space (K s-1)
                      scalarCanopyTempPrime,      & ! intent(in):  prime value for the temperature of the vegetation canopy (K s-1)
                      scalarCanopyWatPrime,       & ! intent(in):  prime value for the water on the vegetation canopy (kg m-2 s-1)
                      mLayerTempPrime,            & ! intent(in):  prime vector of the temperature of each snow and soil layer (K s-1)
                      scalarAquiferStoragePrime,  & ! intent(in):  prime value for storage of water in the aquifer (m s-1)
                      ! input: diagnostic variables defining the liquid water and ice content (function of state variables)
                      scalarCanopyIcePrime,       & ! intent(in):  prime value for the ice on the vegetation canopy (kg m-2 s-1)
                      scalarCanopyLiqPrime,       & ! intent(in):  prime value for the liq on the vegetation canopy (kg m-2 s-1)
                      mLayerVolFracIcePrime,      & ! intent(in):  prime vector of the volumetric ice in each snow and soil layer (s-1)
                      mLayerVolFracWatPrime,      & ! intent(in):  prime vector of the volumetric water in each snow and soil layer (s-1)
                      mLayerVolFracLiqPrime,      & ! intent(in):  prime vector of the volumetric liq in each snow and soil layer (s-1)
                      ! input: enthalpy terms
                      canopyCmTrial,              & ! intent(in):  Cm of vegetation canopy (-)
                      mLayerCmTrial,              & ! intent(in):  Cm of each snow and soil layer (-)
                      scalarCanairEnthalpyPrime,  & ! intent(in):  prime value for the enthalpy of the canopy air space (W m-3)
                      scalarCanopyEnthalpyPrime,  & ! intent(in):  prime value for the of enthalpy of the vegetation canopy (W m-3)
                      mLayerEnthalpyPrime,        & ! intent(in):  prime vector of the of enthalpy of each snow and soil layer (W m-3)
                      ! input: data structures
                      prog_data,                  & ! intent(in):  model prognostic variables for a local HRU
                      diag_data,                  & ! intent(in):  model diagnostic variables for a local HRU
                      flux_data,                  & ! intent(in):  model fluxes for a local HRU
                      indx_data,                  & ! intent(in):  index data
                      ! output
                      resSink,                    & ! intent(out): additional (sink) terms on the RHS of the state equation
                      resVec,                     & ! intent(out): residual vector
                      err,cmessage)                 ! intent(out): error control
      if(err/=0)then; message=trim(message)//trim(cmessage); return; end if  ! (check for errors)

    else ! currently not using residuals outside Sundials!
      dt1 = 1._qp
    endif

  ! end association with the information in the data structures
  end associate

end subroutine eval8summaWithPrime


! **********************************************************************************************************
! public function eval8summa4ida: compute the residual vector F(t,y,y') required for IDA solver
! **********************************************************************************************************
! Return values:
!    0 = success,
!    1 = recoverable error,
!   -1 = non-recoverable error
! ----------------------------------------------------------------
integer(c_int) function eval8summa4ida(tres, sunvec_y, sunvec_yp, sunvec_r, user_data) &
      result(ierr) bind(C,name='eval8summa4ida')

  !======= Inclusions ===========
  use, intrinsic :: iso_c_binding
  use fsundials_core_mod
  use fnvector_serial_mod
  use type4ida

  !======= Declarations =========
  implicit none

  ! calling variables
  real(rkind), value          :: tres             ! current time         t
  type(N_Vector)              :: sunvec_y         ! solution N_Vector    y
  type(N_Vector)              :: sunvec_yp        ! derivative N_Vector  y'
  type(N_Vector)              :: sunvec_r         ! residual N_Vector    F(t,y,y')
  type(c_ptr), value          :: user_data        ! user-defined data

  ! pointers to data in SUNDIALS vectors
  type(data4ida), pointer     :: eqns_data        ! equations data
  real(rkind), pointer        :: stateVec(:)      ! solution vector
  real(rkind), pointer        :: stateVecPrime(:) ! derivative vector
  real(rkind), pointer        :: rVec(:)          ! residual vector
  logical(lgt)                :: feasible         ! feasibility of state vector
  !======= Internals ============

  ! get equations data from user-defined data
  call c_f_pointer(user_data, eqns_data)

  ! get data arrays from SUNDIALS vectors
  stateVec(1:eqns_data%nState)  => FN_VGetArrayPointer(sunvec_y)
  stateVecPrime(1:eqns_data%nState) => FN_VGetArrayPointer(sunvec_yp)
  rVec(1:eqns_data%nState)  => FN_VGetArrayPointer(sunvec_r)

  ! compute the flux and the residual vector for a given state vector
  call eval8summaWithPrime(&
                ! input: model control
                eqns_data%dt,                            & ! intent(in):    data step
                eqns_data%nSnow,                         & ! intent(in):    number of snow layers
                eqns_data%nSoil,                         & ! intent(in):    number of soil layers
                eqns_data%nLayers,                       & ! intent(in):    number of layers
                eqns_data%nState,                        & ! intent(in):    number of state variables in the current subset
                .true.,                                  & ! intent(in):    inside SUNDIALS solver
                eqns_data%firstSubStep,                  & ! intent(in):    flag to indicate if we are processing the first sub-step
                eqns_data%firstFluxCall,                 & ! intent(inout): flag to indicate if we are processing the first flux call
                eqns_data%firstSplitOper,                & ! intent(inout): flag to indicate if we are processing the first flux call in a splitting operation
                eqns_data%computeVegFlux,                & ! intent(in):    flag to indicate if we need to compute fluxes over vegetation
                eqns_data%scalarSolution,                & ! intent(in):    flag to indicate the scalar solution
                ! input: state vectors
                stateVec,                                & ! intent(in):    model state vector
                stateVecPrime,                           & ! intent(in):    model state vector
                eqns_data%sMul,                          & ! intent(inout): state vector multiplier (used in the residual calculations)
                ! input: data structures
                eqns_data%model_decisions,               & ! intent(in):    model decisions
                eqns_data%lookup_data,                   & ! intent(in):    lookup data
                eqns_data%type_data,                     & ! intent(in):    type of vegetation and soil
                eqns_data%attr_data,                     & ! intent(in):    spatial attributes
                eqns_data%mpar_data,                     & ! intent(in):    model parameters
                eqns_data%forc_data,                     & ! intent(in):    model forcing data
                eqns_data%bvar_data,                     & ! intent(in):    average model variables for the entire basin
                eqns_data%prog_data,                     & ! intent(in):    model prognostic variables for a local HRU
                ! input-output: data structures
                eqns_data%indx_data,                     & ! intent(inout): index data
                eqns_data%diag_data,                     & ! intent(inout): model diagnostic variables for a local HRU
                eqns_data%flux_data,                     & ! intent(inout): model fluxes for a local HRU (initial flux structure)
                eqns_data%deriv_data,                    & ! intent(inout): derivatives in model fluxes w.r.t. relevant state variables
                ! input-output: values needed in case canopy gets buried
                eqns_data%scalarCanopyEnthalpyTrial,     & ! intent(inout): trial value for enthalpy of the vegetation canopy (J m-3)
                eqns_data%scalarCanopyTempTrial,         & ! intent(inout): trial value for temperature of the vegetation canopy (K), also used to start enthalpy calculations
                eqns_data%scalarCanopyWatTrial,          & ! intent(inout): trial value for total water content of the vegetation canopy (kg m-2)
                ! output: new values of variables needed in data window outside of internal IDA for rootfinding and to start enthalpy calculations
                eqns_data%mLayerTempTrial,               & ! intent(inout): trial vector of layer temperature (K)
                eqns_data%mLayerMatricHeadTrial,         & ! intent(out):   trial value for total water matric potential (m)
                ! output: new prime values of variables needed in data window outside of internal IDA
                eqns_data%scalarCanopyTempPrime,         & ! intent(out):   prime value for temperature of the vegetation canopy (K s-1)
                eqns_data%scalarCanopyWatPrime,          & ! intent(out):   prime value for total water content of the vegetation canopy (kg m-2 s-1)
                eqns_data%mLayerTempPrime,               & ! intent(out):   prime vector of temperature of each snow and soil layer (K s-1)
                eqns_data%mLayerMatricHeadPrime,         & ! intent(out):   prime vector of matric head of each snow and soil layer (m s-1)
                eqns_data%mLayerVolFracWatPrime,         & ! intent(out):   prime vector of volumetric total water content of each snow and soil layer (s-1)
                ! input-output: baseflow
                eqns_data%ixSaturation,                  & ! intent(inout): index of the lowest saturated layer
                eqns_data%dBaseflow_dMatric,             & ! intent(out):   derivative in baseflow w.r.t. matric head (s-1)
                 ! output: flux and residual vectors
                feasible,                                & ! intent(out):   flag to denote the feasibility of the solution always true inside SUNDIALS
                eqns_data%fluxVec,                       & ! intent(out):   flux vector
                eqns_data%resSink,                       & ! intent(out):   additional (sink) terms on the RHS of the state equation
                rVec,                                    & ! intent(out):   residual vector
                eqns_data%err,eqns_data%message)           ! intent(out):   error control
  if(eqns_data%err > 0)then; eqns_data%message=trim(eqns_data%message); ierr=-1; return; endif
  if(eqns_data%err < 0)then; eqns_data%message=trim(eqns_data%message); ierr=1; return; endif

  ! save residual and return success
  eqns_data%resVec = rVec
  ierr = 0
  return

end function eval8summa4ida


end module eval8summaWithPrime_module

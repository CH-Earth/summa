
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
                    zlookup,      & ! lookup tables
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
                      lookup_data,                   & ! intent(in):    lookup data
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
                      ! input: previous values of variables needed in data window outside of internal IDA
                      scalarCanopyTempPrev,          & ! intent(in):    previous value for temperature of the vegetation canopy (K)
                      mLayerTempPrev,                & ! intent(in):    previous vector of layer temperature (K)
                      mLayerMatricHeadPrev,          & ! intent(in):    previous value for total water matric potential (m)
                      ! output: new values of variables needed in data window outside of internal IDA
                      scalarCanopyTempTrial,         & ! intent(out):   trial value for temperature of the vegetation canopy (K)
                      scalarCanopyWatTrial,          & ! intent(out):   trial value for mass of total water on the vegetation canopy (kg m-2)
                      mLayerTempTrial,               & ! intent(out):   trial vector of layer temperature (K)
                      mLayerMatricHeadTrial,         & ! intent(out):   trial value for total water matric potential (m)
                      mLayerVolFracWatTrial,         & ! intent(out):   trial vector of volumetric total water content (-)
                      ! output: new prime values of variables needed in data window outside of internal IDA
                      scalarCanairTempPrime,         & ! intent(out):   prime value for temperature of the canopy air space (K s-1)
                      scalarCanopyTempPrime,         & ! intent(out):   prime value for temperature of the vegetation canopy (K s-1)
                      scalarCanopyWatPrime,          & ! intent(out):   prime value for total water content of the vegetation canopy (kg m-2 s-1)
                      scalarCanopyIcePrime,          & ! intent(out):   prime value for mass of ice on the vegetation canopy (kg m-2 s-1)
                      mLayerTempPrime,               & ! intent(out):   prime value for temperature of each snow and soil layer (K s-1)
                      mLayerMatricHeadPrime,         & ! intent(out):   prime value for matric head of each snow and soil layer (m s-1)
                      mLayerMatricHeadLiqPrime,      & ! intent(out):   prime value for liquid water matric potential (m s-1)
                      mLayerVolFracWatPrime,         & ! intent(out):   prime value for volumetric total water content of each snow and soil layer (s-1)
                      mLayerVolFracIcePrime,         & ! intent(out):   prime value for volumetric fraction of ice (s-1)
                      ! output: enthalpy prime values    
                      scalarCanairEnthalpyPrime,     & ! intent(out):   prime value for enthalpy of the canopy air space (J m-3 s-1)
                      scalarCanopyEnthalpyPrime,     & ! intent(out):   prime value of enthalpy of the vegetation canopy (J m-3 s-1)
                      mLayerEnthalpyPrime,           & ! intent(out):   prime vector of enthalpy for snow+soil layers (J m-3 s-1)
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
  USE getVectorz_module, only:varExtract                         ! extract variables from the state vector
  USE getVectorz_module, only:checkFeas                          ! check feasibility of state vector
  USE updateVarsWithPrime_module, only:updateVarsWithPrime       ! update variables
  USE t2enthalpyAddPrime_module, only:t2enthalpyPrime            ! compute enthalpy prime and derivatives                     
  USE computFlux_module, only:soilCmpresPrime                    ! compute soil compression
  USE computFlux_module, only:computFlux                         ! compute fluxes given a state vector
  USE computHeatCapWithPrime_module,only:computHeatCapWithPrime  ! recompute heat capacity (Cp) and derivatives
  USE computHeatCap_module,only:computHeatCapAnalytic            ! recompute closed form heat capacity (Cp) and derivatives
  USE computHeatCap_module,only:computCm                         ! compute Cm and derivatives
  USE computHeatCap_module, only:computStatMult                  ! recompute state multiplier
  USE computResidWithPrime_module,only:computResidWithPrime      ! compute residuals given a state vector
  USE computThermConduct_module,only:computThermConduct          ! recompute thermal conductivity and derivatives
  implicit none
  ! --------------------------------------------------------------------------------------------------------------------------------
  ! --------------------------------------------------------------------------------------------------------------------------------
  ! input: model control
  real(rkind),intent(in)          :: dt                              ! entire time step for drainage pond rate
  integer(i4b),intent(in)         :: nSnow                           ! number of snow layers
  integer(i4b),intent(in)         :: nSoil                           ! number of soil layers
  integer(i4b),intent(in)         :: nLayers                         ! total number of layers
  integer,intent(in)              :: nState                          ! total number of state variables
  logical(lgt),intent(in)         :: insideSUN                       ! flag to indicate if we are inside Sundials solver
  logical(lgt),intent(in)         :: firstSubStep                    ! flag to indicate if we are processing the first sub-step
  logical(lgt),intent(inout)      :: firstFluxCall                   ! flag to indicate if we are processing the first flux call
  logical(lgt),intent(inout)      :: firstSplitOper                  ! flag to indicate if we are processing the first flux call in a splitting operation
  logical(lgt),intent(in)         :: computeVegFlux                  ! flag to indicate if computing fluxes over vegetation
  logical(lgt),intent(in)         :: scalarSolution                  ! flag to denote if implementing the scalar solution
  ! input: state vectors    
  real(rkind),intent(in)          :: stateVec(:)                     ! model state vector
  real(rkind),intent(in)          :: stateVecPrime(:)                ! model state vector
  real(qp),intent(inout)          :: sMul(:)   ! NOTE: qp            ! state vector multiplier (used in the residual calculations)
  ! input: data structures
  type(model_options),intent(in)  :: model_decisions(:)              ! model decisions
  type(zLookup),      intent(in)  :: lookup_data                     ! lookup tables
  type(var_i),        intent(in)  :: type_data                       ! type of vegetation and soil
  type(var_d),        intent(in)  :: attr_data                       ! spatial attributes
  type(var_dlength),  intent(in)  :: mpar_data                       ! model parameters
  type(var_d),        intent(in)  :: forc_data                       ! model forcing data
  type(var_dlength),  intent(in)  :: bvar_data                       ! model variables for the local basin
  type(var_dlength),  intent(in)  :: prog_data                       ! prognostic variables for a local HRU
  ! output: data structures
  type(var_ilength),intent(inout) :: indx_data                       ! indices defining model states and layers
  type(var_dlength),intent(inout) :: diag_data                       ! diagnostic variables for a local HRU
  type(var_dlength),intent(inout) :: flux_data                       ! model fluxes for a local HRU
  type(var_dlength),intent(inout) :: deriv_data                      ! derivatives in model fluxes w.r.t. relevant state variables
  ! input: previous values of variables needed in data window outside of internal IDA
  real(rkind),intent(in)          :: scalarCanopyTempPrev            ! previous value for temperature of the vegetation canopy (K)
  real(rkind),intent(in)          :: mLayerTempPrev(:)               ! previous vector of layer temperature (K)
  real(rkind),intent(in)          :: mLayerMatricHeadPrev(:)         ! previous value for total water matric potential (m)
  ! output: new values of variables needed in data window outside of internal IDA
  real(rkind),intent(out)         :: scalarCanopyTempTrial           ! trial value for temperature of the vegetation canopy (K)
  real(rkind),intent(out)         :: scalarCanopyWatTrial            ! trial value for liquid water storage in the canopy (kg m-2)
  real(rkind),intent(out)         :: mLayerTempTrial(:)              ! trial vector of layer temperature (K)
  real(rkind),intent(out)         :: mLayerMatricHeadTrial(:)        ! trial vector for total water matric potential (m)
  real(rkind),intent(out)         :: mLayerVolFracWatTrial(:)        ! trial vector of volumetric total water content (-)
  ! output: new prime values of variables needed in data window outside of internal IDA
  real(rkind),intent(out)         :: scalarCanairTempPrime           ! prime value for temperature of the canopy air space (K s-1)
  real(rkind),intent(out)         :: scalarCanopyTempPrime           ! prime value for temperature of the vegetation canopy (K s-1)
  real(rkind),intent(out)         :: scalarCanopyWatPrime            ! prime value for total water content of the vegetation canopy (kg m-2 s-1)
  real(rkind),intent(out)         :: scalarCanopyIcePrime            ! prime value for mass of ice on the vegetation canopy (kg m-2 s-1)
  real(rkind),intent(out)         :: mLayerTempPrime(:)              ! prime vector for temperature of each snow and soil layer (K s-1)
  real(rkind),intent(out)         :: mLayerMatricHeadPrime(:)        ! prime vector for matric head of each snow and soil layer (m s-1)
  real(rkind),intent(out)         :: mLayerMatricHeadLiqPrime(:)     ! prime vector for liquid water matric potential (m s-1)
  real(rkind),intent(out)         :: mLayerVolFracWatPrime(:)        ! prime vector for volumetric total water content of each snow and soil layer (s-1)
  real(rkind),intent(out)         :: mLayerVolFracIcePrime(:)        ! prime vector for volumetric fraction of ice (s-1)
  ! output: enthalpy prime values    
  real(rkind),intent(out)         :: scalarCanairEnthalpyPrime       ! prime value for enthalpy of the canopy air space (J m-3 s-1)
  real(rkind),intent(out)         :: scalarCanopyEnthalpyPrime       ! prime value of enthalpy of the vegetation canopy (J m-3 s-1)
  real(rkind),intent(out)         :: mLayerEnthalpyPrime(:)          ! prime vector of enthalpy for snow+soil layers (J m-3 s-1)
  ! input-output: baseflow    
  integer(i4b),intent(inout)      :: ixSaturation                    ! index of the lowest saturated layer
  real(rkind),intent(out)         :: dBaseflow_dMatric(:,:)          ! derivative in baseflow w.r.t. matric head (s-1)
  ! output: flux and residual vectors
  logical(lgt),intent(out)        :: feasible                        ! flag to denote the feasibility of the solution
  real(rkind),intent(out)         :: fluxVec(:)                      ! flux vector
  real(rkind),intent(out)         :: resSink(:)                      ! sink terms on the RHS of the flux equation
  real(qp),intent(out)            :: resVec(:) ! NOTE: qp            ! residual vector
  ! output: error control
  integer(i4b),intent(out)        :: err                             ! error code
  character(*),intent(out)        :: message                         ! error message
  ! --------------------------------------------------------------------------------------------------------------------------------
  ! local variables
  ! --------------------------------------------------------------------------------------------------------------------------------
  real(rkind)                     :: dt1                             ! residual step size
  ! state variables
  real(rkind)                     :: scalarCanairTempTrial           ! trial value for temperature of the canopy air space (K)
  real(rkind)                     :: scalarCanopyLiqTrial            ! trial value for liquid water storage in the canopy (kg m-2)
  real(rkind)                     :: scalarCanopyIceTrial            ! trial value for ice storage in the canopy (kg m-2)
  real(rkind),dimension(nSoil)    :: mLayerMatricHeadLiqTrial        ! trial value for liquid water matric potential (m)
  real(rkind),dimension(nLayers)  :: mLayerVolFracLiqTrial           ! trial vector for volumetric liquid water content (-)
  real(rkind),dimension(nLayers)  :: mLayerVolFracIceTrial           ! trial vector for volumetric ice content (-)
  real(rkind)                     :: scalarAquiferStorageTrial       ! trial value of storage of water in the aquifer (m)
  real(rkind)                     :: scalarCanopyLiqPrime            ! prime value for liquid water storage in the canopy (kg m-2 s-1)
  real(rkind),dimension(nLayers)  :: mLayerVolFracLiqPrime           ! prime vector for volumetric liquid water content (s-1)
  real(rkind)                     :: scalarAquiferStoragePrime       ! prime value of storage of water in the aquifer (m s-1)
  ! enthalpy derivatives
  real(rkind)                     :: dCanEnthalpyPrime_dTk           ! derivatives in prime canopy enthalpy w.r.t. temperature
  real(rkind)                     :: dCanEnthalpyPrime_dWat          ! derivatives in prime canopy enthalpy w.r.t. water state
  real(rkind),dimension(nLayers)  :: dEnthalpyPrime_dTk              ! derivatives in prime layer enthalpy w.r.t. temperature
  real(rkind),dimension(nLayers)  :: dEnthalpyPrime_dWat             ! derivatives in prime layer enthalpy w.r.t. water state
  real(rkind)                     :: dCanEnthalpyPrime_dTkPrime      ! derivatives in prime canopy enthalpy w.r.t. prime temperature
  real(rkind)                     :: dCanEnthalpyPrime_dWatPrime     ! derivatives in prime canopy enthalpy w.r.t. prime water state
  real(rkind),dimension(nLayers)  :: dEnthalpyPrime_dTkPrime         ! derivatives in prime layer enthalpy w.r.t. prime temperature
  real(rkind),dimension(nLayers)  :: dEnthalpyPrime_dWatPrime        ! derivatives in prime layer enthalpy w.r.t. prime water state
  ! other local variables
  integer(i4b)                    :: iLayer                          ! index of model layer in the snow+soil domain
  integer(i4b)                    :: jState(1)                       ! index of model state for the scalar solution within the soil domain
  integer(i4b)                    :: ixBeg,ixEnd                     ! index of indices for the soil compression routine
  character(LEN=256)              :: cmessage                        ! error message of downwind routine
  real(rkind)                     :: scalarCanopyCmTrial             ! trial value of Cm for the canopy
  real(rkind),dimension(nLayers)  :: mLayerCmTrial                   ! trial vector of Cm for snow+soil
  logical(lgt),parameter          :: updateCp=.true.                 ! flag to indicate if we update Cp at each step
  logical(lgt),parameter          :: needCm=.true.                   ! flag to indicate if the energy equation contains Cm = dH_T/dTheta_m

  ! --------------------------------------------------------------------------------------------------------------------------------
  ! association to variables in the data structures
  ! --------------------------------------------------------------------------------------------------------------------------------
  associate(&
    ! model decisions
    ixHowHeatCap              => model_decisions(iLookDECISIONS%howHeatCap)%iDecision         ,& ! intent(in):  [i4b]    heat capacity computation, with or without enthalpy
    ixRichards                => model_decisions(iLookDECISIONS%f_Richards)%iDecision         ,& ! intent(in):  [i4b]    index of the form of Richards' equation
    ! snow parameters
    snowfrz_scale             => mpar_data%var(iLookPARAM%snowfrz_scale)%dat(1)               ,& ! intent(in):  [dp]     scaling parameter for the snow freezing curve (K-1)
    ! soil parameters
    theta_sat                 => mpar_data%var(iLookPARAM%theta_sat)%dat                      ,& ! intent(in):  [dp(:)]  soil porosity (-)
    specificStorage           => mpar_data%var(iLookPARAM%specificStorage)%dat(1)             ,& ! intent(in):  [dp]     specific storage coefficient (m-1)
      ! canopy and layer depth
    canopyDepth               => diag_data%var(iLookDIAG%scalarCanopyDepth)%dat(1)            ,& ! intent(in):  [dp   ]  canopy depth (m)
    mLayerDepth               => prog_data%var(iLookPROG%mLayerDepth)%dat                     ,& ! intent(in):  [dp(:)]  depth of each layer in the snow-soil sub-domain (m)
    ! model diagnostic variables from a previous solution
    scalarFracLiqVeg          => diag_data%var(iLookDIAG%scalarFracLiqVeg)%dat(1)             ,& ! intent(in):  [dp]     fraction of liquid water on vegetation (-)
    scalarSfcMeltPond         => prog_data%var(iLookPROG%scalarSfcMeltPond)%dat(1)            ,& ! intent(in):  [dp]     ponded water caused by melt of the "snow without a layer" (kg m-2)
    mLayerFracLiqSnow         => diag_data%var(iLookDIAG%mLayerFracLiqSnow)%dat               ,& ! intent(in):  [dp(:)]  fraction of liquid water in each snow layer (-)
    ! soil compression
    scalarSoilCompress        => diag_data%var(iLookDIAG%scalarSoilCompress)%dat(1)           ,& ! intent(in):  [dp]     total change in storage associated with compression of the soil matrix (kg m-2 s-1)
    mLayerCompress            => diag_data%var(iLookDIAG%mLayerCompress)%dat                  ,& ! intent(in):  [dp(:)]  change in volumetric water content due to compression of soil (s-1)
    ! derivatives
    dTheta_dTkCanopy          => deriv_data%var(iLookDERIV%dTheta_dTkCanopy)%dat(1)           ,& ! intent(in):  [dp]     derivative of volumetric liquid water content w.r.t. temperature
    dVolTot_dPsi0             => deriv_data%var(iLookDERIV%dVolTot_dPsi0)%dat                 ,& ! intent(in):  [dp(:)]  derivative in total water content w.r.t. total water matric potential
    d2VolTot_dPsi02           => deriv_data%var(iLookDERIV%d2VolTot_dPsi02)%dat               ,& ! intent(in):  [dp(:)]  second derivative in total water content w.r.t. total water matric potential
    dCompress_dPsi            => deriv_data%var(iLookDERIV%dCompress_dPsi)%dat                ,& ! intent(in):  [dp(:)]  derivative in compressibility w.r.t. matric head (m-1)
    mLayerdTheta_dTk          => deriv_data%var(iLookDERIV%mLayerdTheta_dTk)%dat              ,& ! intent(in):  [dp(:)]  derivative of volumetric liquid water content w.r.t. temperature
    dVolHtCapBulk_dPsi0       => deriv_data%var(iLookDERIV%dVolHtCapBulk_dPsi0)%dat           ,& ! intent(out): [dp(:)]  derivative in bulk heat capacity w.r.t. matric potential
    dVolHtCapBulk_dTheta      => deriv_data%var(iLookDERIV%dVolHtCapBulk_dTheta)%dat          ,& ! intent(out): [dp(:)]  derivative in bulk heat capacity w.r.t. volumetric water content
    dVolHtCapBulk_dCanWat     => deriv_data%var(iLookDERIV%dVolHtCapBulk_dCanWat)%dat(1)      ,& ! intent(out): [dp]     derivative in bulk heat capacity w.r.t. volumetric water content
    dVolHtCapBulk_dTk         => deriv_data%var(iLookDERIV%dVolHtCapBulk_dTk)%dat             ,& ! intent(out): [dp(:)]  derivative in bulk heat capacity w.r.t. temperature
    dVolHtCapBulk_dTkCanopy   => deriv_data%var(iLookDERIV%dVolHtCapBulk_dTkCanopy)%dat(1)    ,& ! intent(out): [dp]     derivative in bulk heat capacity w.r.t. temperature
    dVolHtCapBulk_dPsi0Prime  => deriv_data%var(iLookDERIV%dVolHtCapBulk_dPsi0Prime)%dat      ,& ! intent(out): [dp(:)]  derivative in bulk heat capacity w.r.t. prime matric potential
    dVolHtCapBulk_dThetaPrime => deriv_data%var(iLookDERIV%dVolHtCapBulk_dThetaPrime)%dat     ,& ! intent(out): [dp(:)]  derivative in bulk heat capacity w.r.t. prime volumetric water content
    dVolHtCapBulk_dCanWatPrime=> deriv_data%var(iLookDERIV%dVolHtCapBulk_dCanWatPrime)%dat(1) ,& ! intent(out): [dp]     derivative in bulk heat capacity w.r.t. prime volumetric water content
    dVolHtCapBulk_dTkPrime    => deriv_data%var(iLookDERIV%dVolHtCapBulk_dTkPrime)%dat        ,& ! intent(out): [dp(:)]  derivative in bulk heat capacity w.r.t. prime temperature
    dVolHtCapBulk_dTkCanPrime => deriv_data%var(iLookDERIV%dVolHtCapBulk_dTkCanPrime)%dat(1)  ,& ! intent(out): [dp]     derivative in bulk heat capacity w.r.t. prime temperature       
    dThermalC_dWatAbove       => deriv_data%var(iLookDERIV%dThermalC_dWatAbove)%dat           ,& ! intent(out): [dp(:)]  derivative in the thermal conductivity w.r.t. water state in the layer above
    dThermalC_dWatBelow       => deriv_data%var(iLookDERIV%dThermalC_dWatBelow)%dat           ,& ! intent(out): [dp(:)]  derivative in the thermal conductivity w.r.t. water state in the layer above
    dThermalC_dTempAbove      => deriv_data%var(iLookDERIV%dThermalC_dTempAbove)%dat          ,& ! intent(out): [dp(:)]  derivative in the thermal conductivity w.r.t. energy state in the layer above
    dThermalC_dTempBelow      => deriv_data%var(iLookDERIV%dThermalC_dTempBelow)%dat          ,& ! intent(out): [dp(:)]  derivative in the thermal conductivity w.r.t. energy state in the layer above
    dCm_dTk                   => deriv_data%var(iLookDERIV%dCm_dTk)%dat                       ,& ! intent(out): [dp(:)]  derivative in heat capacity w.r.t. temperature (J kg-1 K-2)
    dCm_dTkCanopy             => deriv_data%var(iLookDERIV%dCm_dTkCanopy)%dat(1)              ,& ! intent(out): [dp   ]  derivative in heat capacity w.r.t. canopy temperature (J kg-1 K-2)
    ! mapping
    ixMapFull2Subset          => indx_data%var(iLookINDEX%ixMapFull2Subset)%dat               ,& ! intent(in):  [i4b(:)] mapping of full state vector to the state subset
    ixControlVolume           => indx_data%var(iLookINDEX%ixControlVolume)%dat                ,& ! intent(in):  [i4b(:)] index of control volume for different domains (veg, snow, soil)
    ! heat capacity
    heatCapVegTrial           => diag_data%var(iLookDIAG%scalarBulkVolHeatCapVeg)%dat(1)      ,& ! intent(out): [dp]     volumetric heat capacity of vegetation canopy
    mLayerHeatCapTrial        => diag_data%var(iLookDIAG%mLayerVolHtCapBulk)%dat              ,& ! intent(out): [dp(:)]  heat capacity for snow and soil
    ! Cm
    canopyCmTrial             => diag_data%var(iLookDIAG%scalarCanopyCm)%dat(1)               ,& ! intent(out): [dp]     Cm of the canopy
    mLayerCmTrial             => diag_data%var(iLookDIAG%mLayerCm)%dat                         & ! intent(out): [dp(:)]  Cm of snow and soil
    ) ! association to variables in the data structures
    ! --------------------------------------------------------------------------------------------------------------------------------
    ! initialize error control
    err=0; message="eval8summaWithPrime/"

    ! check the feasibility of the solution only if not inside Sundials solver
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

    ! Placeholder: if we decide to use splitting, we need to pass all the previous values of the state variables
    scalarCanairTempTrial     = realMissing
    scalarCanopyTempTrial     = scalarCanopyTempPrev
    scalarCanopyWatTrial      = realMissing
    scalarCanopyLiqTrial      = realMissing
    scalarCanopyIceTrial      = realMissing
    mLayerTempTrial           = mLayerTempPrev
    mLayerVolFracWatTrial     = realMissing
    mLayerVolFracLiqTrial     = realMissing
    mLayerVolFracIceTrial     = realMissing
    mLayerMatricHeadTrial     = mLayerMatricHeadPrev
    mLayerMatricHeadLiqTrial  = realMissing
    scalarAquiferStorageTrial = realMissing

    ! extract states from the state vector
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
                    scalarAquiferStorageTrial,& ! intent(inout):   trial value of storage of water in the aquifer (m)
                    ! output: error control
                    err,cmessage)               ! intent(out):   error control
    if(err/=0)then; message=trim(message)//trim(cmessage); return; end if  ! (check for errors)

    ! Placeholder: if we decide to use splitting, we need to pass all the previous values of the state variables
    scalarCanairTempPrime     = realMissing
    scalarCanopyTempPrime     = realMissing
    scalarCanopyWatPrime      = realMissing
    scalarCanopyLiqPrime      = realMissing
    scalarCanopyIcePrime      = realMissing
    mLayerTempPrime           = realMissing
    mLayerVolFracWatPrime     = realMissing
    mLayerVolFracLiqPrime     = realMissing
    mLayerVolFracIcePrime     = realMissing
    mLayerMatricHeadPrime     = realMissing
    mLayerMatricHeadLiqPrime  = realMissing
    scalarAquiferStoragePrime = realMissing

    call varExtract(&
                  ! input
                  stateVecPrime,            & ! intent(in):    derivative of model state vector (mixed units)
                  diag_data,                & ! intent(in):    model diagnostic variables for a local HRU
                  prog_data,                & ! intent(in):    model prognostic variables for a local HRU
                  indx_data,                & ! intent(in):    indices defining model states and layers
                  ! output: variables for the vegetation canopy
                  scalarCanairTempPrime,    & ! intent(inout): derivative of canopy air temperature (K)
                  scalarCanopyTempPrime,    & ! intent(inout): derivative of canopy temperature (K)
                  scalarCanopyWatPrime,     & ! intent(inout): derivative of canopy total water (kg m-2)
                  scalarCanopyLiqPrime,     & ! intent(inout): derivative of canopy liquid water (kg m-2)
                  ! output: variables for the snow-soil domain
                  mLayerTempPrime,          & ! intent(inout): derivative of layer temperature (K)
                  mLayerVolFracWatPrime,    & ! intent(inout): derivative of volumetric total water content (-)
                  mLayerVolFracLiqPrime,    & ! intent(inout): derivative of volumetric liquid water content (-)
                  mLayerMatricHeadPrime,    & ! intent(inout): derivative of total water matric potential (m)
                  mLayerMatricHeadLiqPrime, & ! intent(inout): derivative of liquid water matric potential (m)
                  ! output: variables for the aquifer
                  scalarAquiferStoragePrime,& ! intent(inout): derivative of storage of water in the aquifer (m)
                  ! output: error control
                  err,cmessage)               ! intent(out):   error control
    if(err/=0)then; message=trim(message)//trim(cmessage); return; end if  ! (check for errors)

    ! update diagnostic variables and derivatives
    call updateVarsWithPrime(&
                    ! input
                    .false.,                                   & ! intent(in):    logical flag if computing for Jacobian update
                    .false.,                                   & ! intent(in):    logical flag to adjust temperature to account for the energy
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
                    scalarCanopyTempPrime,                     & ! intent(inout): trial value of time derivative canopy temperature (K)
                    scalarCanopyWatPrime,                      & ! intent(inout): trial value of time derivative canopy total water (kg m-2)
                    scalarCanopyLiqPrime,                      & ! intent(inout): trial value of time derivative canopy liquid water (kg m-2)
                    scalarCanopyIcePrime,                      & ! intent(inout): trial value of time derivative canopy ice content (kg m-2)
                    ! output: variables for the snow-soil domain
                    mLayerTempTrial,                           & ! intent(inout): trial vector of layer temperature (K)
                    mLayerVolFracWatTrial,                     & ! intent(inout): trial vector of volumetric total water content (-)
                    mLayerVolFracLiqTrial,                     & ! intent(inout): trial vector of volumetric liquid water content (-)
                    mLayerVolFracIceTrial,                     & ! intent(inout): trial vector of volumetric ice water content (-)
                    mLayerMatricHeadTrial,                     & ! intent(inout): trial vector of total water matric potential (m)
                    mLayerMatricHeadLiqTrial,                  & ! intent(inout): trial vector of liquid water matric potential (m)
                    mLayerTempPrime,                           & ! intent(inout): trial vector of time derivative layer temperature (K)
                    mLayerVolFracWatPrime,                     & ! intent(inout): trial vector of time derivative volumetric total water content (-)
                    mLayerVolFracLiqPrime,                     & ! intent(inout): trial vector of time derivative volumetric liquid water content (-)
                    mLayerVolFracIcePrime,                     & ! intent(inout): trial vector of time derivative volumetric ice water content (-)
                    mLayerMatricHeadPrime,                     & ! intent(inout): trial vector of time derivative total water matric potential (m)
                    mLayerMatricHeadLiqPrime,                  & ! intent(inout): trial vector of time derivative liquid water matric potential (m)
                    ! output: error control
                    err,cmessage)                                ! intent(out):   error control
    if(err/=0)then; message=trim(message)//trim(cmessage); return; end if  ! (check for errors)

    if(updateCp)then
       ! *** compute volumetric heat capacity C_p
      if(ixHowHeatCap == enthalpyFD)then
        ! compute H_T prime without phase change
        call t2enthalpyPrime(&
                        ! input: data structures
                        diag_data,                         & ! intent(in):   model diagnostic variables for a local HRU
                        mpar_data,                         & ! intent(in):   parameter data structure
                        indx_data,                         & ! intent(in):   model indices
                        lookup_data,                       & ! intent(in):   lookup table data structure
                        ! input: state variables for the vegetation canopy   
                        scalarCanairTempPrime,             & ! intent(in):   prime value of canopy air temperature (K)
                        scalarCanopyTempTrial,             & ! intent(in):   trial value of canopy temperature (K)
                        scalarCanopyWatTrial,              & ! intent(in):   trial value of canopy total water (kg m-2)                         
                        scalarCanopyTempPrime,             & ! intent(in):   prime value of canopy temperature (K)
                        scalarCanopyWatPrime,              & ! intent(in):   prime value of canopy total water (kg m-2)
                         ! input: variables for the snow-soil domain
                        mLayerTempTrial,                   & ! intent(in):   trial vector of layer temperature (K)
                        mLayerVolFracWatTrial,             & ! intent(in):   trial vector of volumetric total water content (-)
                        mLayerMatricHeadTrial,             & ! intent(in):   trial vector of total water matric potential (m)
                        mLayerTempPrime,                   & ! intent(in):   prime vector of layer temperature (K)
                        mLayerVolFracWatPrime,             & ! intent(in):   prime vector of volumetric total water content (-)
                        mLayerMatricHeadPrime,             & ! intent(in):   prime vector of total water matric potential (m)  
                        ! input: pre-computed derivatives
                        dVolTot_dPsi0,                     & ! intent(in):   derivative in total water content w.r.t. total water matric potential (m-1)
                        d2VolTot_dPsi02,                   & ! intent(in):   second derivative in total water content w.r.t. total water matric potential (m-2)
                        ! output: enthalpy prime and derivatives
                        scalarCanairEnthalpyPrime,         & ! intent(out):  prime enthalpy of the canopy air space (J m-3)
                        scalarCanopyEnthalpyPrime,         & ! intent(out):  prime enthalpy of the vegetation canopy (J m-3)
                        mLayerEnthalpyPrime,               & ! intent(out):  prime enthalpy of each snow+soil layer (J m-3)
                        dCanEnthalpyPrime_dTk,             & ! intent(out):  derivatives in prime canopy enthalpy w.r.t. temperature
                        dCanEnthalpyPrime_dWat,            & ! intent(out):  derivatives in prime canopy enthalpy w.r.t. water state
                        dEnthalpyPrime_dTk,                & ! intent(out):  derivatives in prime layer enthalpy w.r.t. temperature
                        dEnthalpyPrime_dWat,               & ! intent(out):  derivatives in prime layer enthalpy w.r.t. water state
                        dCanEnthalpyPrime_dTkPrime,        & ! intent(out):  derivatives in prime canopy enthalpy w.r.t. prime temperature
                        dCanEnthalpyPrime_dWatPrime,       & ! intent(out):  derivatives in prime canopy enthalpy w.r.t. prime water state
                        dEnthalpyPrime_dTkPrime,           & ! intent(out):  derivatives in prime layer enthalpy w.r.t. prime temperature
                        dEnthalpyPrime_dWatPrime,          & ! intent(out):  derivatives in prime layer enthalpy w.r.t. prime water state
                        ! output: error control
                        err,cmessage)                        ! intent(out):  error control
        if(err/=0)then; message=trim(message)//trim(cmessage); return; endif

        ! *** compute volumetric heat capacity C_p = dH_T/dT
        call computHeatCapWithPrime(&
                            ! input: control variables
                            nLayers,                    & ! intent(in):    number of layers (soil+snow)
                            computeVegFlux,             & ! intent(in):    flag to denote if computing the vegetation flux
                            canopyDepth,                & ! intent(in):    canopy depth (m)
                            ! input output data structures
                            mpar_data,                  & ! intent(in):    model parameters
                            indx_data,                  & ! intent(in):    model layer indices
                            diag_data,                  & ! intent(inout): model diagnostic variables for a local HRU
                            ! input: state variables
                            scalarCanopyIceTrial,       & ! intent(in):    trial value for mass of ice on the vegetation canopy (kg m-2)
                            scalarCanopyLiqTrial,       & ! intent(in):    trial value for the liquid water on the vegetation canopy (kg m-2)
                            scalarCanopyTempTrial,      & ! intent(in):    trial value of canopy temperature (K)
                            scalarCanopyTempPrime,      & ! intent(in):    prime value of canopy temperature (K)
                            scalarCanopyEnthalpyPrime,  & ! intent(in):    prime enthalpy of the vegetation canopy (J m-3)
                            mLayerVolFracIceTrial,      & ! intent(in):    volumetric fraction of ice at the start of the sub-step (-)
                            mLayerVolFracLiqTrial,      & ! intent(in):    volumetric fraction of liquid water at the start of the sub-step (-)
                            mLayerTempTrial,            & ! intent(in):    trial temperature
                            mLayerTempPrime,            & ! intent(in):    prime temperature
                            mLayerEnthalpyPrime,        & ! intent(in):    prime enthalpy for snow and soil
                            mLayerMatricHeadTrial,      & ! intent(in):    trial total water matric potential (m)
                            ! input: pre-computed derivatives
                            dTheta_dTkCanopy,           & ! intent(in):    derivative in canopy volumetric liquid water content w.r.t. temperature (K-1)
                            scalarFracLiqVeg,           & ! intent(in):    fraction of canopy liquid water (-)
                            mLayerdTheta_dTk,           & ! intent(in):    derivative of volumetric liquid water content w.r.t. temperature (K-1)
                            mLayerFracLiqSnow,          & ! intent(in):    fraction of liquid water (-)
                            dVolTot_dPsi0,              & ! intent(in):    derivative in total water content w.r.t. total water matric potential (m-1)
                            dCanEnthalpyPrime_dTk,      & ! intent(in):    derivatives in prime canopy enthalpy w.r.t. temperature
                            dCanEnthalpyPrime_dWat,     & ! intent(in):    derivatives in prime canopy enthalpy w.r.t. water state
                            dEnthalpyPrime_dTk,         & ! intent(in):    derivatives in prime layer enthalpy w.r.t. temperature
                            dEnthalpyPrime_dWat,        & ! intent(in):    derivatives in prime layer enthalpy w.r.t. water state
                            dCanEnthalpyPrime_dTkPrime, & ! intent(in):    derivatives in prime canopy enthalpy w.r.t. prime temperature
                            dCanEnthalpyPrime_dWatPrime,& ! intent(in):    derivatives in prime canopy enthalpy w.r.t. prime water state
                            dEnthalpyPrime_dTkPrime,    & ! intent(in):    derivatives in prime layer enthalpy w.r.t. prime temperature
                            dEnthalpyPrime_dWatPrime,   & ! intent(in):    derivatives in prime layer enthalpy w.r.t. prime water state       
                            ! output
                            heatCapVegTrial,            & ! intent(out):   volumetric heat capacity of vegetation canopy
                            mLayerHeatCapTrial,         & ! intent(out):   heat capacity for snow and soil
                            dVolHtCapBulk_dPsi0,        & ! intent(out):   derivative in bulk heat capacity w.r.t. matric potential
                            dVolHtCapBulk_dTheta,       & ! intent(out):   derivative in bulk heat capacity w.r.t. volumetric water content
                            dVolHtCapBulk_dCanWat,      & ! intent(out):   derivative in bulk heat capacity w.r.t. volumetric water content
                            dVolHtCapBulk_dTk,          & ! intent(out):   derivative in bulk heat capacity w.r.t. temperature
                            dVolHtCapBulk_dTkCanopy,    & ! intent(out):   derivative in bulk heat capacity w.r.t. temperature                  
                            dVolHtCapBulk_dPsi0Prime,   & ! intent(out):   derivative in bulk heat capacity w.r.t. prime matric potential
                            dVolHtCapBulk_dThetaPrime,  & ! intent(out):   derivative in bulk heat capacity w.r.t. prime volumetric water content
                            dVolHtCapBulk_dCanWatPrime, & ! intent(out):   derivative in bulk heat capacity w.r.t. prime volumetric water content
                            dVolHtCapBulk_dTkPrime,     & ! intent(out):   derivative in bulk heat capacity w.r.t. prime temperature
                            dVolHtCapBulk_dTkCanPrime,  & ! intent(out):   derivative in bulk heat capacity w.r.t. prime temperature       
                            ! output: error control
                            err,cmessage)                    ! intent(out): error control
        if(err/=0)then; message=trim(message)//trim(cmessage); return; endif
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
                          err,cmessage)                  ! intent(out):  error control
        if(err/=0)then; message=trim(message)//trim(cmessage); return; endif
        ! set derivatives to 0 for closed form
        dVolHtCapBulk_dPsi0Prime   = 0._rkind
        dVolHtCapBulk_dThetaPrime  = 0._rkind
        dVolHtCapBulk_dCanWatPrime = 0._rkind
        dVolHtCapBulk_dTkPrime     = 0._rkind
        dVolHtCapBulk_dTkCanPrime  = 0._rkind
      endif !(choice of how compute heat capacity)

      ! compute multiplier of state vector
      call computStatMult(&
                    ! input
                    heatCapVegTrial,         & ! intent(in):  volumetric heat capacity of vegetation canopy
                    mLayerHeatCapTrial,      & ! intent(in):  volumetric heat capacity of soil and snow
                    diag_data,               & ! intent(in):  model diagnostic variables for a local HRU
                    indx_data,               & ! intent(in):  indices defining model states and layers
                    ! output
                    sMul,                    & ! intent(out): multiplier for state vector (used in the residual calculations)
                    err,cmessage)              ! intent(out): error control
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
                          ! output: error control
                          err,cmessage)                   ! intent(out): error control
      if(err/=0)then; err=55; message=trim(message)//trim(cmessage); return; end if
    else
      ! set heat capacity derivatives to 0 for constant through step
      dVolHtCapBulk_dPsi0        = 0._rkind
      dVolHtCapBulk_dTheta       = 0._rkind
      dVolHtCapBulk_dCanWat      = 0._rkind
      dVolHtCapBulk_dTk          = 0._rkind
      dVolHtCapBulk_dTkCanopy    = 0._rkind
      dVolHtCapBulk_dPsi0Prime   = 0._rkind
      dVolHtCapBulk_dThetaPrime  = 0._rkind
      dVolHtCapBulk_dCanWatPrime = 0._rkind
      dVolHtCapBulk_dTkPrime     = 0._rkind
      dVolHtCapBulk_dTkCanPrime  = 0._rkind
      dThermalC_dWatAbove        = 0._rkind
      dThermalC_dWatBelow        = 0._rkind
      dThermalC_dTempAbove       = 0._rkind
      dThermalC_dTempBelow       = 0._rkind  
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

    firstSplitOper = .false. ! after call computeFlux once in dt, no longer firstSplitOper

    ! compute soil compressibility (-) and its derivative w.r.t. matric head (m)
    ! NOTE: we already extracted trial matrix head and volumetric liquid water as part of the flux calculations
    call soilCmpresPrime(&
                    ! input:
                    ixRichards,                             & ! intent(in):    choice of option for Richards' equation
                    ixBeg,ixEnd,                            & ! intent(in):    start and end indices defining desired layers
                    mLayerMatricHeadPrime(1:nSoil),         & ! intent(in):    matric head at the start of the time step (m)
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
                      dt1,                       & ! intent(in):  length of the residual time step (seconds)
                      nSnow,                     & ! intent(in):  number of snow layers
                      nSoil,                     & ! intent(in):  number of soil layers
                      nLayers,                   & ! intent(in):  total number of layers
                      ! input: flux vectors
                      sMul,                      & ! intent(in):  state vector multiplier (used in the residual calculations)
                      fluxVec,                   & ! intent(in):  flux vector
                      ! input: state variables (already disaggregated into scalars and vectors)
                      scalarCanopyTempTrial,     & ! intent(in):  trial value for the temperature of the vegetation canopy (K)
                      mLayerTempTrial,           & ! intent(in):  trial value for the temperature of each snow and soil layer (K)
                      scalarCanairTempPrime,     & ! intent(in):  Prime value for the temperature of the canopy air space (K s-1)
                      scalarCanopyTempPrime,     & ! intent(in):  Prime value for the temperature of the vegetation canopy (K s-1)
                      scalarCanopyWatPrime,      & ! intent(in):  Prime value for the water on the vegetation canopy (kg m-2 s-1)
                      mLayerTempPrime,           & ! intent(in):  Prime value for the temperature of each snow and soil layer (K s-1)
                      scalarAquiferStoragePrime, & ! intent(in):  Prime value of storage of water in the aquifer (m s-1)
                      ! input: diagnostic variables defining the liquid water and ice content (function of state variables)
                      scalarCanopyIcePrime,      & ! intent(in):  Prime value for the ice on the vegetation canopy (kg m-2 s-1)
                      scalarCanopyLiqPrime,      & ! intent(in):  Prime value for the liq on the vegetation canopy (kg m-2 s-1)
                      mLayerVolFracIcePrime,     & ! intent(in):  Prime value for the volumetric ice in each snow and soil layer (s-1)
                      mLayerVolFracWatPrime,     & ! intent(in):  Prime value for the volumetric water in each snow and soil layer (s-1)
                      mLayerVolFracLiqPrime,     & ! intent(in):  Prime value for the volumetric liq in each snow and soil layer (s-1)
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

    else !currently not using residuals outside Sundials!
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
  use fsundials_nvector_mod
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
                ! input: previous values of variables needed in data window outside of internal IDA
                eqns_data%scalarCanopyTempPrev,          & ! intent(in):    previous value for temperature of the vegetation canopy (K)
                eqns_data%mLayerTempPrev,                & ! intent(in):    previous vector of layer temperature (K)
                eqns_data%mLayerMatricHeadPrev,          & ! intent(in):    previous value for total water matric potential (m)
                ! output: new values of variables needed in data window outside of internal IDA
                eqns_data%scalarCanopyTempTrial,         & ! intent(out):   trial value for temperature of the vegetation canopy (K)
                eqns_data%scalarCanopyWatTrial,          & ! intent(out):   trial value for mass of total water on the vegetation canopy (kg m-2)
                eqns_data%mLayerTempTrial,               & ! intent(out):   trial vector of layer temperature (K)
                eqns_data%mLayerMatricHeadTrial,         & ! intent(out):   trial value for total water matric potential (m)
                eqns_data%mLayerVolFracWatTrial,         & ! intent(out):   trial vector of volumetric total water content (-)
                ! output: new prime values of variables needed in data window outside of internal IDA
                eqns_data%scalarCanairTempPrime,         & ! intent(out):   prime value for temperature of the canopy air space (K s-1)
                eqns_data%scalarCanopyTempPrime,         & ! intent(out):   prime value for temperature of the vegetation canopy (K s-1)
                eqns_data%scalarCanopyWatPrime,          & ! intent(out):   prime value for total water content of the vegetation canopy (kg m-2 s-1)
                eqns_data%scalarCanopyIcePrime,          & ! intent(out):   prime value for mass of ice on the vegetation canopy (kg m-2 s-1)
                eqns_data%mLayerTempPrime,               & ! intent(out):   prime value for temperature of each snow and soil layer (K s-1)
                eqns_data%mLayerMatricHeadPrime,         & ! intent(out):   prime value for matric head of each snow and soil layer (m s-1)
                eqns_data%mLayerMatricHeadLiqPrime,      & ! intent(out):   prime value for liquid water matric potential (m s-1)
                eqns_data%mLayerVolFracWatPrime,         & ! intent(out):   prime value for volumetric total water content of each snow and soil layer (s-1)
                eqns_data%mLayerVolFracIcePrime,         & ! intent(out):   prime value for volumetric fraction of ice (s-1)
                ! output: enthalpy prime values    
                eqns_data%scalarCanairEnthalpyPrime,     & ! intent(out):   prime value for enthalpy of the canopy air space (J m-3 s-1)
                eqns_data%scalarCanopyEnthalpyPrime,     & ! intent(out):   prime value of enthalpy of the vegetation canopy (J m-3 s-1)
                eqns_data%mLayerEnthalpyPrime,           & ! intent(out):   prime vector of enthalpy for snow+soil layers (J m-3 s-1)
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

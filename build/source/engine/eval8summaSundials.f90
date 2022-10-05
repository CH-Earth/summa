
module eval8summaSundials_module

! data types
USE nrtype

! access missing values
USE globalData,only:integerMissing  ! missing integer
USE globalData,only:realMissing     ! missing double precision number
USE globalData,only:quadMissing     ! missing quadruple precision number

! access the global print flag
USE globalData,only:globalPrintFlag

! define access to state variables to print
USE globalData,only: iJac1          ! first layer of the Jacobian to print
USE globalData,only: iJac2          ! last layer of the Jacobian to print

! domain types
USE globalData,only:iname_veg       ! named variables for vegetation
USE globalData,only:iname_snow      ! named variables for snow
USE globalData,only:iname_soil      ! named variables for soil

! named variables to describe the state variable type
USE globalData,only:iname_nrgCanair ! named variable defining the energy of the canopy air space
USE globalData,only:iname_nrgCanopy ! named variable defining the energy of the vegetation canopy
USE globalData,only:iname_watCanopy ! named variable defining the mass of water on the vegetation canopy
USE globalData,only:iname_nrgLayer  ! named variable defining the energy state variable for snow+soil layers
USE globalData,only:iname_watLayer  ! named variable defining the total water state variable for snow+soil layers
USE globalData,only:iname_liqLayer  ! named variable defining the liquid  water state variable for snow+soil layers
USE globalData,only:iname_matLayer  ! named variable defining the matric head state variable for soil layers
USE globalData,only:iname_lmpLayer  ! named variable defining the liquid matric potential state variable for soil layers
USE globalData,only:model_decisions        ! model decision structure

! constants
USE multiconst,only:&
                    Tfreeze,      & ! temperature at freezing              (K)
                    LH_fus,       & ! latent heat of fusion                (J kg-1)
                    LH_vap,       & ! latent heat of vaporization          (J kg-1)
                    LH_sub,       & ! latent heat of sublimation           (J kg-1)
                    Cp_air,       & ! specific heat of air                 (J kg-1 K-1)
                    iden_air,     & ! intrinsic density of air             (kg m-3)
                    iden_ice,     & ! intrinsic density of ice             (kg m-3)
                    iden_water      ! intrinsic density of liquid water    (kg m-3)

! provide access to the derived types to define the data structures
USE data_types,only:&
                    var_i,        & ! data vector (i4b)
                    var_d,        & ! data vector (rkind)
                    var_ilength,  & ! data vector with variable length dimension (i4b)
                    var_dlength,  & ! data vector with variable length dimension (rkind)
                    zlookup,      &
                    model_options   ! defines the model decisions

! indices that define elements of the data structures
USE var_lookup,only:iLookDECISIONS               ! named variables for elements of the decision structure
USE var_lookup,only:iLookPARAM                   ! named variables for structure elements
USE var_lookup,only:iLookPROG                    ! named variables for structure elements
USE var_lookup,only:iLookINDEX                   ! named variables for structure elements
USE var_lookup,only:iLookDIAG                    ! named variables for structure elements
USE var_lookup,only:iLookFLUX                    ! named variables for structure elements
USE var_lookup,only:iLookDERIV                   ! named variables for structure elements

! look-up values for the choice of groundwater representation (local-column, or single-basin)
USE mDecisions_module,only:  &
 localColumn,                & ! separate groundwater representation in each local soil column
 singleBasin,                & ! single groundwater store over the entire basin
 enthalpyFD                    ! heat capacity using enthalpy

! look-up values for the choice of groundwater parameterization
USE mDecisions_module,only:  &
 qbaseTopmodel,              & ! TOPMODEL-ish baseflow parameterization
 bigBucket,                  & ! a big bucket (lumped aquifer model)
 noExplicit                    ! no explicit groundwater parameterization

! look-up values for the form of Richards' equation
USE mDecisions_module,only:  &
 moisture,                   & ! moisture-based form of Richards' equation
 mixdform                      ! mixed form of Richards' equation

implicit none
private
public::eval8summaSundials
public::eval8summa4IDA


contains

! **********************************************************************************************************
! public subroutine eval8summaSundials: compute the residual vector
! **********************************************************************************************************
subroutine eval8summaSundials(&
                      ! input: model control
                      dt_cur,                  & ! intent(in):    current stepsize
                      dt,                      & ! intent(in):    entire time step
                      nSnow,                   & ! intent(in):    number of snow layers
                      nSoil,                   & ! intent(in):    number of soil layers
                      nLayers,                 & ! intent(in):    total number of layers
                      nState,                  & ! intent(in):    total number of state variables
                      insideIDA,               & ! intent(in):    flag to indicate if we are inside Sundials solver
                      firstSubStep,            & ! intent(in):    flag to indicate if we are processing the first sub-step
                      firstFluxCall,           & ! intent(inout): flag to indicate if we are processing the first flux call
                      firstSplitOper,          & ! intent(inout): flag to indicate if we are processing the first flux call in a splitting operation
                      computeVegFlux,          & ! intent(in):    flag to indicate if we need to compute fluxes over vegetation
                      scalarSolution,          & ! intent(in):    flag to indicate the scalar solution
                      ! input: state vectors
                      stateVec,                & ! intent(in):    model state vector
                      stateVecPrime,           & ! intent(in):    derivative of model state vector
                      sMul,                    & ! intent(inout): state vector multiplier (used in the residual calculations)
                      ! input: data structures
                      model_decisions,         & ! intent(in):    model decisions
                      lookup_data,             & ! intent(in):    lookup data
                      type_data,               & ! intent(in):    type of vegetation and soil
                      attr_data,               & ! intent(in):    spatial attributes
                      mpar_data,               & ! intent(in):    model parameters
                      forc_data,               & ! intent(in):    model forcing data
                      bvar_data,               & ! intent(in):    average model variables for the entire basin
                      prog_data,               & ! intent(in):    model prognostic variables for a local HRU
                      ! input-output: data structures
                      indx_data,               & ! intent(inout): index data
                      diag_data,               & ! intent(inout)  model diagnostic variables for a local HRU
                      flux_data,               & ! intent(inout): model fluxes for a local HRU
                      deriv_data,              & ! intent(inout): derivatives in model fluxes w.r.t. relevant state variables
                      ! input-output: here we need to pass some extra variables that do not get updated in in the Sundials loops
                      scalarCanopyTempTrial,   & ! intent(out):   trial value of canopy temperature (K)
                      scalarCanopyTempPrev,    & ! intent(in):    value of canopy temperature (K)
                      scalarCanopyIceTrial,    & ! intent(out):   trial value for mass of ice on the vegetation canopy (kg m-2)
                      scalarCanopyIcePrev,     & ! intent(in):    value for mass of ice on the vegetation canopy (kg m-2)
                      scalarCanopyLiqTrial,    & ! intent(out):   trial value of canopy liquid water (kg m-2)
                      scalarCanopyLiqPrev,     & ! intent(in):    value of canopy liquid water (kg m-2)
                      scalarCanopyEnthalpyTrial,& ! intent(out):  trial value for enthalpy of the vegetation canopy (J m-3)
                      scalarCanopyEnthalpyPrev,& ! intent(in):    value for enthalpy of the vegetation canopy (J m-3)
                      mLayerTempTrial,         & ! intent(out):   trial vector of layer temperature (K)
                      mLayerTempPrev,          & ! intent(in):    vector of layer temperature (K)
                      mLayerMatricHeadLiqTrial,& ! intent(out):   trial value for liquid water matric potential (m)
                      mLayerMatricHeadTrial,   & ! intent(out):   trial value for total water matric potential (m)
                      mLayerMatricHeadPrev,    & ! intent(in):    value for total water matric potential (m)
                      mLayerVolFracWatTrial,   & ! intent(out):   trial vector of volumetric total water content (-)
                      mLayerVolFracWatPrev,    & ! intent(in):    vector of volumetric total water content (-)
                      mLayerVolFracIceTrial,   & ! intent(out):   trial vector of volumetric ice water content (-)
                      mLayerVolFracIcePrev,    & ! intent(in):    vector of volumetric ice water content (-)
                      mLayerVolFracLiqTrial,   & ! intent(out):   trial vector of volumetric liquid water content (-)
                      mLayerVolFracLiqPrev,    & ! intent(in):    vector of volumetric liquid water content (-)
                      scalarAquiferStorageTrial,& ! intent(out):  trial value of storage of water in the aquifer (m)
                      scalarAquiferStoragePrev,& ! intent(in):    value of storage of water in the aquifer (m)
                      mLayerEnthalpyPrev,      & ! intent(in):    vector of enthalpy for snow+soil layers (J m-3)
                      mLayerEnthalpyTrial,     & ! intent(out):   trial vector of enthalpy for snow+soil layers (J m-3)
                      ! input-output: baseflow
                      ixSaturation,            & ! intent(inout): index of the lowest saturated layer
                      dBaseflow_dMatric,       & ! intent(out):   derivative in baseflow w.r.t. matric head (s-1)
                      ! output: flux and residual vectors
                      feasible,                & ! intent(out):   flag to denote the feasibility of the solution
                      fluxVec,                 & ! intent(out):   flux vector
                      resSink,                 & ! intent(out):   additional (sink) terms on the RHS of the state equation
                      resVec,                  & ! intent(out):   residual vector
                      err,message)               ! intent(out):   error control
  ! --------------------------------------------------------------------------------------------------------------------------------
  ! provide access to subroutines
  USE getVectorz_module, only:varExtract                  ! extract variables from the state vector
  USE updateVarsSundials_module, only:updateVarsSundials  ! update variables
  USE t2enthalpy_module, only:t2enthalpy_T                ! compute enthalpy
  USE computFlux_module, only:soilCmpresSundials          ! compute soil compression
  USE computFlux_module, only:computFlux                  ! compute fluxes given a state vector
  USE computHeatCap_module,only:computHeatCap             ! compute heat capacity
  USE computHeatCap_module,only:computHeatCapAnalytic     ! compute heat capacity
  USE computHeatCap_module,only:computCm
  USE computHeatCap_module, only:computStatMult
  USE computResidSundials_module,only:computResidSundials ! compute residuals given a state vector
  USE computThermConduct_module,only:computThermConduct
  USE computEnthalpy_module,only:computEnthalpy
  USE computEnthalpy_module,only:computEnthalpyPrime
  implicit none
  ! --------------------------------------------------------------------------------------------------------------------------------
  ! --------------------------------------------------------------------------------------------------------------------------------
  ! input: model control
  real(rkind),intent(in)          :: dt_cur                 ! current stepsize
  real(rkind),intent(in)          :: dt                     ! entire time step
  integer(i4b),intent(in)         :: nSnow                  ! number of snow layers
  integer(i4b),intent(in)         :: nSoil                  ! number of soil layers
  integer(i4b),intent(in)         :: nLayers                ! total number of layers
  integer,intent(in)              :: nState                 ! total number of state variables
  logical(lgt),intent(in)         :: insideIDA              ! flag to indicate if we are inside Sundials solver
  logical(lgt),intent(in)         :: firstSubStep           ! flag to indicate if we are processing the first sub-step
  logical(lgt),intent(inout)      :: firstFluxCall
  logical(lgt),intent(inout)      :: firstSplitOper         ! flag to indicate if we are processing the first flux call in a splitting operation
  logical(lgt),intent(in)         :: computeVegFlux         ! flag to indicate if computing fluxes over vegetation
  logical(lgt),intent(in)         :: scalarSolution         ! flag to denote if implementing the scalar solution
  ! input: state vectors
  real(rkind),intent(in)          :: stateVec(:)            ! model state vector
  real(rkind),intent(in)          :: stateVecPrime(:)       ! model state vector
  real(qp),intent(inout)          :: sMul(:)   ! NOTE: qp   ! state vector multiplier (used in the residual calculations)
  ! input: data structures
  type(model_options),intent(in)  :: model_decisions(:)     ! model decisions
  type(zLookup),intent(in)        :: lookup_data            ! lookup tables
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
  ! input-output: here we need to pass some extra variables that do not get updated in in the Sundials loops
  real(rkind),intent(out)         :: scalarCanopyTempTrial     ! trial value for temperature of the vegetation canopy (K)
  real(rkind),intent(in)          :: scalarCanopyTempPrev      ! previous value for temperature of the vegetation canopy (K)
  real(rkind),intent(out)         :: scalarCanopyIceTrial      ! trial value for mass of ice on the vegetation canopy (kg m-2)
  real(rkind),intent(in)          :: scalarCanopyIcePrev       ! previous value for mass of ice on the vegetation canopy (kg m-2)
  real(rkind),intent(out)         :: scalarCanopyLiqTrial      ! trial value for mass of ice on the vegetation canopy (kg m-2)
  real(rkind),intent(in)          :: scalarCanopyLiqPrev       ! previous value for mass of ice on the vegetation canopy (kg m-2)
  real(rkind),intent(out)         :: scalarCanopyEnthalpyTrial ! trial value for enthalpy of the vegetation canopy (J m-3)
  real(rkind),intent(in)          :: scalarCanopyEnthalpyPrev  ! previous value of enthalpy of the vegetation canopy (J m-3)
  real(rkind),intent(out)         :: mLayerTempTrial(:)        ! trial vector of layer temperature (K)
  real(rkind),intent(in)          :: mLayerTempPrev(:)
  real(rkind),intent(out)         :: mLayerMatricHeadLiqTrial(:)  ! trial value for liquid water matric potential (m)
  real(rkind),intent(out)         :: mLayerMatricHeadTrial(:)  ! trial value for total water matric potential (m)
  real(rkind),intent(in)          :: mLayerMatricHeadPrev(:) ! value for total water matric potential (m)
  real(rkind),intent(out)         :: mLayerVolFracWatTrial(:)  ! trial vector of volumetric total water content (-)
  real(rkind),intent(in)          :: mLayerVolFracWatPrev(:) ! vector of volumetric total water content (-)
  real(rkind),intent(out)         :: mLayerVolFracIceTrial(:)  ! trial vector of volumetric ice water content (-)
  real(rkind),intent(in)          :: mLayerVolFracIcePrev(:) ! vector of volumetric ice water content (-)
  real(rkind),intent(out)         :: mLayerVolFracLiqTrial(:)  ! trial vector of volumetric liquid water content (-)
  real(rkind),intent(in)          :: mLayerVolFracLiqPrev(:) ! vector of volumetric liquid water content (-)
  real(rkind),intent(out)         :: scalarAquiferStorageTrial ! trial value of storage of water in the aquifer (m)
  real(rkind),intent(in)          :: scalarAquiferStoragePrev  ! value of storage of water in the aquifer (m)
  real(rkind),intent(in)          :: mLayerEnthalpyPrev(:)    ! vector of enthalpy for snow+soil layers (J m-3)
  real(rkind),intent(out)         :: mLayerEnthalpyTrial(:)   ! trial vector of enthalpy for snow+soil layers (J m-3)
  ! input-output: baseflow
  integer(i4b),intent(inout)      :: ixSaturation           ! index of the lowest saturated layer
  real(rkind),intent(out)         :: dBaseflow_dMatric(:,:) ! derivative in baseflow w.r.t. matric head (s-1)
  ! output: flux and residual vectors
  logical(lgt),intent(out)        :: feasible               ! flag to denote the feasibility of the solution
  real(rkind),intent(out)         :: fluxVec(:)             ! flux vector
  real(rkind),intent(out)         :: resSink(:)             ! sink terms on the RHS of the flux equation
  real(qp),intent(out)            :: resVec(:) ! NOTE: qp   ! residual vector
  ! output: error control
  integer(i4b),intent(out)        :: err                    ! error code
  character(*),intent(out)        :: message                ! error message
  ! --------------------------------------------------------------------------------------------------------------------------------
  ! local variables
  ! --------------------------------------------------------------------------------------------------------------------------------
  real(rkind)                       :: dt1                       ! residual step size
  ! state variables
  real(rkind)                       :: scalarCanairTempTrial     ! trial value for temperature of the canopy air space (K)
  real(rkind)                       :: scalarCanopyWatTrial      ! trial value for liquid water storage in the canopy (kg m-2)
  ! derivative of state variables
  real(rkind)                        :: scalarCanairTempPrime     ! derivative value for temperature of the canopy air space (K)
  real(rkind)                        :: scalarCanopyTempPrime     ! derivative value for temperature of the vegetation canopy (K)
  real(rkind)                        :: scalarCanopyWatPrime      ! derivative value for liquid water storage in the canopy (kg m-2)
  real(rkind),dimension(nLayers)     :: mLayerTempPrime           ! derivative value for temperature of layers in the snow and soil domains (K)
  real(rkind),dimension(nLayers)     :: mLayerVolFracWatPrime     ! derivative value for volumetric fraction of total water (-)
  real(rkind),dimension(nSoil)       :: mLayerMatricHeadPrime     ! derivative value for total water matric potential (m)
  real(rkind),dimension(nSoil)       :: mLayerMatricHeadLiqPrime  ! derivative value for liquid water matric potential (m)
  real(rkind)                        :: scalarAquiferStoragePrime ! derivative value of storage of water in the aquifer (m)
  ! derivative of diagnostic variables
  real(rkind)                        :: scalarCanopyLiqPrime      ! derivative value for mass of liquid water on the vegetation canopy (kg m-2)
  real(rkind)                        :: scalarCanopyIcePrime      ! derivative value for mass of ice on the vegetation canopy (kg m-2)
  real(rkind),dimension(nLayers)     :: mLayerVolFracLiqPrime     ! derivative value for volumetric fraction of liquid water (-)
  real(rkind),dimension(nLayers)     :: mLayerVolFracIcePrime     ! derivative value for volumetric fraction of ice (-)
  ! enthalpy
  real(rkind)                       :: scalarCanairEnthalpy      ! enthalpy of the canopy air space (J m-3)
  real(rkind),dimension(nLayers)    :: mLayerEnthalpyPrime       ! enthalpy of each snow+soil layer (J m-3)
  ! other local variables
  integer(i4b)                      :: iLayer                    ! index of model layer in the snow+soil domain
  integer(i4b)                      :: jState(1)                 ! index of model state for the scalar solution within the soil domain
  integer(i4b)                      :: ixBeg,ixEnd               ! index of indices for the soil compression routine
  integer(i4b),parameter            :: ixVegVolume=1             ! index of the desired vegetation control volumne (currently only one veg layer)
  real(rkind)                       :: xMin,xMax                 ! minimum and maximum values for water content
  real(rkind),parameter             :: canopyTempMax=500._rkind  ! expected maximum value for the canopy temperature (K)
  character(LEN=256)                :: cmessage                  ! error message of downwind routine
  real(rkind)                       :: scalarCanopyCmTrial       ! trial value of Cm for the canopy
  real(rkind),dimension(nLayers)    :: mLayerCmTrial             ! trial vector of Cm for snow+soil
  logical(lgt),parameter            :: updateCp=.true.           ! flag to indicate if we update Cp at each step
  logical(lgt),parameter            :: needCm=.false.            ! flag to indicate if the energy equation contains Cm = dH_T/dTheta_m

  ! --------------------------------------------------------------------------------------------------------------------------------
  ! association to variables in the data structures
  ! --------------------------------------------------------------------------------------------------------------------------------
  associate(&
    ! model decisions
    ixRichards              => model_decisions(iLookDECISIONS%f_Richards)%iDecision   ,&  ! intent(in):  [i4b]   index of the form of Richards' equation
    ! snow parameters
    snowfrz_scale           => mpar_data%var(iLookPARAM%snowfrz_scale)%dat(1)         ,&  ! intent(in):  [dp]    scaling parameter for the snow freezing curve (K-1)
    ! soil parameters
    theta_sat               => mpar_data%var(iLookPARAM%theta_sat)%dat                ,&  ! intent(in):  [dp(:)] soil porosity (-)
    specificStorage         => mpar_data%var(iLookPARAM%specificStorage)%dat(1)       ,&  ! intent(in):  [dp]    specific storage coefficient (m-1)
    theta_res               => mpar_data%var(iLookPARAM%theta_res)%dat                ,&  ! intent(in):  [dp(:)] residual volumetric water content (-)
    ! canopy and layer depth
    canopyDepth             => diag_data%var(iLookDIAG%scalarCanopyDepth)%dat(1)      ,&  ! intent(in):  [dp   ] canopy depth (m)
    mLayerDepth             => prog_data%var(iLookPROG%mLayerDepth)%dat               ,&  ! intent(in):  [dp(:)] depth of each layer in the snow-soil sub-domain (m)
    ! model state variables
    scalarSfcMeltPond       => prog_data%var(iLookPROG%scalarSfcMeltPond)%dat(1)      ,&  ! intent(in):  [dp]    ponded water caused by melt of the "snow without a layer" (kg m-2)
    mLayerVolFracIce        => prog_data%var(iLookPROG%mLayerVolFracIce)%dat          ,&  ! intent(in):  [dp(:)] volumetric fraction of ice (-)
    ! soil compression
    scalarSoilCompress      => diag_data%var(iLookDIAG%scalarSoilCompress)%dat(1)     ,&  ! intent(in): [dp]    total change in storage associated with compression of the soil matrix (kg m-2)
    mLayerCompress          => diag_data%var(iLookDIAG%mLayerCompress)%dat            ,&  ! intent(in): [dp(:)] change in storage associated with compression of the soil matrix (-)
    ! derivatives
    dVolTot_dPsi0           => deriv_data%var(iLookDERIV%dVolTot_dPsi0)%dat           ,&  ! intent(in): [dp(:)] derivative in total water content w.r.t. total water matric potential
    dCompress_dPsi          => deriv_data%var(iLookDERIV%dCompress_dPsi)%dat          ,&  ! intent(in): [dp(:)] derivative in compressibility w.r.t. matric head (m-1)
    ! mapping
    ixMapFull2Subset        => indx_data%var(iLookINDEX%ixMapFull2Subset)%dat         ,&  ! intent(in): [i4b(:)] mapping of full state vector to the state subset
    ixControlVolume         => indx_data%var(iLookINDEX%ixControlVolume)%dat          ,&  ! intent(in): [i4b(:)] index of control volume for different domains (veg, snow, soil)
    ! indices
    ixCasNrg                => indx_data%var(iLookINDEX%ixCasNrg)%dat(1)              ,&  ! intent(in): [i4b]    index of canopy air space energy state variable (nrg)
    ixVegNrg                => indx_data%var(iLookINDEX%ixVegNrg)%dat(1)              ,&  ! intent(in): [i4b]    index of canopy energy state variable (nrg)
    ixVegHyd                => indx_data%var(iLookINDEX%ixVegHyd)%dat(1)              ,&  ! intent(in): [i4b]    index of canopy hydrology state variable (mass)
    ixSnowOnlyNrg           => indx_data%var(iLookINDEX%ixSnowOnlyNrg)%dat            ,&  ! intent(in): [i4b(:)] indices for energy states in the snow subdomain
    ixSnowSoilHyd           => indx_data%var(iLookINDEX%ixSnowSoilHyd)%dat            ,&  ! intent(in): [i4b(:)] indices for hydrology states in the snow+soil subdomain
    ixStateType             => indx_data%var(iLookINDEX%ixStateType)%dat              ,&  ! intent(in): [i4b(:)] indices defining the type of the state (iname_nrgLayer...)
    ixHydCanopy             => indx_data%var(iLookINDEX%ixHydCanopy)%dat              ,&  ! intent(in): [i4b(:)] index of the hydrology states in the canopy domain
    ixHydType               => indx_data%var(iLookINDEX%ixHydType)%dat                ,&  ! intent(in): [i4b(:)] index of the type of hydrology states in snow+soil domain
    layerType               => indx_data%var(iLookINDEX%layerType)%dat                ,&  ! intent(in): [i4b(:)] layer type (iname_soil or iname_snow)
    heatCapVegTrial         =>  diag_data%var(iLookDIAG%scalarBulkVolHeatCapVeg)%dat(1),& ! intent(out): volumetric heat capacity of vegetation canopy
    mLayerHeatCapTrial      =>  diag_data%var(iLookDIAG%mLayerVolHtCapBulk)%dat        &  ! intent(out): heat capacity for snow and soil
    ) ! association to variables in the data structures
    ! --------------------------------------------------------------------------------------------------------------------------------
    ! initialize error control
    err=0; message="eval8summaSundials/"
    feasible=.true.

    ! check the feasibility of the solution only if not inside Sundials solver
    if (.not.insideIDA) then
      ! check that the canopy air space temperature is reasonable
      if(ixCasNrg/=integerMissing)then
        if(stateVec(ixCasNrg) > canopyTempMax) feasible=.false.
        if(stateVec(ixCasNrg) > canopyTempMax) message=trim(message)//'canopy air space temp high,'
        if(.not.feasible) write(*,'(a,1x,L1,1x,10(f20.10,1x))') 'feasible, max, stateVec( ixCasNrg )', feasible, canopyTempMax, stateVec(ixCasNrg)
      endif

      ! check that the canopy temperature is reasonable
      if(ixVegNrg/=integerMissing)then
        if(stateVec(ixVegNrg) > canopyTempMax) feasible=.false.
        if(stateVec(ixVegNrg) > canopyTempMax) message=trim(message)//'canopy temp high,'
        if(.not.feasible) write(*,'(a,1x,L1,1x,10(f20.10,1x))') 'feasible, max, stateVec( ixVegNrg )', feasible, canopyTempMax, stateVec(ixVegNrg)
      endif

      ! check canopy liquid water is not negative
      if(ixVegHyd/=integerMissing)then
        if(stateVec(ixVegHyd) < 0._rkind) feasible=.false.
        if(stateVec(ixVegHyd) < 0._rkind) message=trim(message)//'canopy water negative,'
        if(.not.feasible) write(*,'(a,1x,L1,1x,10(f20.10,1x))') 'feasible, min, stateVec( ixVegHyd )', feasible, 0._rkind, stateVec(ixVegHyd)
      end if

      ! check snow temperature is below freezing
      if(count(ixSnowOnlyNrg/=integerMissing)>0)then
        if(any(stateVec( pack(ixSnowOnlyNrg,ixSnowOnlyNrg/=integerMissing) ) > Tfreeze)) feasible=.false.
        if(any(stateVec( pack(ixSnowOnlyNrg,ixSnowOnlyNrg/=integerMissing) ) > Tfreeze)) message=trim(message)//'snow temp above freezing,'
        do iLayer=1,nSnow
          if(.not.feasible) write(*,'(a,1x,i4,1x,L1,1x,10(f20.10,1x))') 'iLayer, feasible, max, stateVec( ixSnowOnlyNrg(iLayer) )', iLayer, feasible, Tfreeze, stateVec( ixSnowOnlyNrg(iLayer) )
        enddo
      endif

      ! loop through non-missing hydrology state variables in the snow+soil domain
      do concurrent (iLayer=1:nLayers,ixSnowSoilHyd(iLayer)/=integerMissing)

        ! check the minimum and maximum water constraints
        if(ixHydType(iLayer)==iname_watLayer .or. ixHydType(iLayer)==iname_liqLayer)then

          ! --> minimum
          if (layerType(iLayer) == iname_soil) then
            xMin = theta_res(iLayer-nSnow)
          else
            xMin = 0._rkind
          endif

          ! --> maximum
          select case( layerType(iLayer) )
            case(iname_snow); xMax = merge(iden_ice,  1._rkind - mLayerVolFracIce(iLayer), ixHydType(iLayer)==iname_watLayer)
            case(iname_soil); xMax = merge(theta_sat(iLayer-nSnow), theta_sat(iLayer-nSnow) - mLayerVolFracIce(iLayer), ixHydType(iLayer)==iname_watLayer)
          end select

          ! --> check
          if(stateVec( ixSnowSoilHyd(iLayer) ) < xMin .or. stateVec( ixSnowSoilHyd(iLayer) ) > xMax) feasible=.false.
          if(stateVec( ixSnowSoilHyd(iLayer) ) < xMin .or. stateVec( ixSnowSoilHyd(iLayer) ) > xMax)  message=trim(message)//'layer water outside bounds,'
          if(.not.feasible) write(*,'(a,1x,i4,1x,L1,1x,10(f20.10,1x))') 'iLayer, feasible, stateVec( ixSnowSoilHyd(iLayer) ), xMin, xMax = ', iLayer, feasible, stateVec( ixSnowSoilHyd(iLayer) ), xMin, xMax

        endif  ! if water states

      end do  ! loop through non-missing hydrology state variables in the snow+soil domain

      ! early return for non-feasible solutions
      if(.not.feasible)then
        fluxVec(:) = realMissing
        resVec(:)  = quadMissing
        message=trim(message)//'non-feasible'
        err=20; return
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
    scalarCanopyTempTrial     = scalarCanopyTempPrev
    scalarCanopyLiqTrial      = scalarCanopyLiqPrev
    scalarCanopyIceTrial      = scalarCanopyIcePrev
    mLayerTempTrial           = mLayerTempPrev
    mLayerVolFracWatTrial     = mLayerVolFracWatPrev
    mLayerVolFracLiqTrial     = mLayerVolFracLiqPrev
    mLayerVolFracIceTrial     = mLayerVolFracIcePrev
    mLayerMatricHeadTrial     = mLayerMatricHeadPrev
    scalarAquiferStorageTrial = scalarAquiferStoragePrev

    ! extract states from the state vector
    call varExtract(&
                    ! input
                    stateVec,                 & ! intent(in):    model state vector (mixed units)
                    diag_data,                & ! intent(in):    model diagnostic variables for a local HRU
                    prog_data,                & ! intent(in):    model prognostic variables for a local HRU
                    indx_data,                & ! intent(in):    indices defining model states and layers
                    ! output: variables for the vegetation canopy
                    scalarCanairTempTrial,    & ! intent(inout):   trial value of canopy air temperature (K)
                    scalarCanopyTempTrial,    & ! intent(inout):   trial value of canopy temperature (K)
                    scalarCanopyWatTrial,     & ! intent(inout):   trial value of canopy total water (kg m-2)
                    scalarCanopyLiqTrial,     & ! intent(inout):   trial value of canopy liquid water (kg m-2)
                    ! output: variables for the snow-soil domain
                    mLayerTempTrial,          & ! intent(inout):   trial vector of layer temperature (K)
                    mLayerVolFracWatTrial,    & ! intent(inout):   trial vector of volumetric total water content (-)
                    mLayerVolFracLiqTrial,    & ! intent(inout):   trial vector of volumetric liquid water content (-)
                    mLayerMatricHeadTrial,    & ! intent(inout):   trial vector of total water matric potential (m)
                    mLayerMatricHeadLiqTrial, & ! intent(inout):   trial vector of liquid water matric potential (m)
                    ! output: variables for the aquifer
                    scalarAquiferStorageTrial,& ! intent(inout):   trial value of storage of water in the aquifer (m)
                    ! output: error control
                    err,cmessage)               ! intent(out):   error control
    if(err/=0)then; message=trim(message)//trim(cmessage); return; end if  ! (check for errors)

    ! initialize to state variable from the last update
    ! should all be set to previous values if splits, but for now operator splitting is not hooked up
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
                  scalarCanairTempPrime,    & ! intent(inout):   derivative of canopy air temperature (K)
                  scalarCanopyTempPrime,    & ! intent(inout):   derivative of canopy temperature (K)
                  scalarCanopyWatPrime,     & ! intent(niout):   derivative of canopy total water (kg m-2)
                  scalarCanopyLiqPrime,     & ! intent(inout):   derivative of canopy liquid water (kg m-2)
                  ! output: variables for the snow-soil domain
                  mLayerTempPrime,          & ! intent(inout):   derivative of layer temperature (K)
                  mLayerVolFracWatPrime,    & ! intent(inout):   derivative of volumetric total water content (-)
                  mLayerVolFracLiqPrime,    & ! intent(inout):   derivative of volumetric liquid water content (-)
                  mLayerMatricHeadPrime,    & ! intent(inout):   derivative of total water matric potential (m)
                  mLayerMatricHeadLiqPrime, & ! intent(inout):   derivative of liquid water matric potential (m)
                  ! output: variables for the aquifer
                  scalarAquiferStoragePrime,& ! intent(inout):   derivative of storage of water in the aquifer (m)
                  ! output: error control
                  err,cmessage)               ! intent(out):   error control
    if(err/=0)then; message=trim(message)//trim(cmessage); return; end if  ! (check for errors)

    call updateVarsSundials(&
                    ! input
                    dt_cur,                                    &
                    .false.,                                   & ! intent(in):    logical flag if computing Jacobian for Sundials solver
                    .false.,                                   & ! intent(in):    logical flag to adjust temperature to account for the energy
                    mpar_data,                                 & ! intent(in):    model parameters for a local HRU
                    indx_data,                                 & ! intent(in):    indices defining model states and layers
                    prog_data,                                 & ! intent(in):    model prognostic variables for a local HRU
                    mLayerVolFracWatPrev,                      & ! intent(in):    previous vector of total water matric potential (m)
                    mLayerMatricHeadPrev,                      & ! intent(in):    previous vector of volumetric total water content (-)
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
                    mLayerTempPrime,                           & !
                    mLayerVolFracWatPrime,                     & ! intent(inout): trial vector of time derivative volumetric total water content (-)
                    mLayerVolFracLiqPrime,                     & ! intent(inout): trial vector of time derivative volumetric liquid water content (-)
                    mLayerVolFracIcePrime,                     & ! intent(inout): trial vector of time derivative volumetric ice water content (-)
                    mLayerMatricHeadPrime,                     & ! intent(inout): trial vector of time derivative total water matric potential (m)
                    mLayerMatricHeadLiqPrime,                  & ! intent(inout): trial vector of time derivative liquid water matric potential (m)
                    ! output: error control
                    err,cmessage)                                ! intent(out):   error control
    if(err/=0)then; message=trim(message)//trim(cmessage); return; end if  ! (check for errors)

    ! print the water content
    if(globalPrintFlag)then
      if(iJac1<nSnow) write(*,'(a,10(f16.10,1x))') 'mLayerVolFracWatTrial = ', mLayerVolFracWatTrial(iJac1:min(iJac2,nSnow))
      if(iJac1<nSnow) write(*,'(a,10(f16.10,1x))') 'mLayerVolFracLiqTrial = ', mLayerVolFracLiqTrial(iJac1:min(iJac2,nSnow))
      if(iJac1<nSnow) write(*,'(a,10(f16.10,1x))') 'mLayerVolFracIceTrial = ', mLayerVolFracIceTrial(iJac1:min(iJac2,nSnow))
    endif

    if(updateCp)then
       ! *** compute volumetric heat capacity C_p
      if(model_decisions(iLookDECISIONS%howHeatCap)%iDecision == enthalpyFD)then
        ! compute H_T
        call t2enthalpy_T(&
                        ! input: data structures
                        diag_data,                   & ! intent(in):  model diagnostic variables for a local HRU
                        mpar_data,                   & ! intent(in):  parameter data structure
                        indx_data,                   & ! intent(in):  model indices
                        lookup_data,                 & ! intent(in):  lookup table data structure
                        ! input: state variables for the vegetation canopy
                        scalarCanairTempTrial,       & ! intent(in):  trial value of canopy air temperature (K)
                        scalarCanopyTempTrial,       & ! intent(in):  trial value of canopy temperature (K)
                        scalarCanopyWatTrial,        & ! intent(in):  trial value of canopy total water (kg m-2)
                        scalarCanopyIceTrial,        & ! intent(in):  trial value for canopy ice content (kg m-2)
                        ! input: variables for the snow-soil domain
                        mLayerTempTrial,             & ! intent(in):  trial vector of layer temperature (K)
                        mLayerVolFracWatTrial,       & ! intent(in):  trial vector of volumetric total water content (-)
                        mLayerMatricHeadTrial,       & ! intent(in):  trial vector of total water matric potential (m)
                        mLayerVolFracIceTrial,       & ! intent(in):  trial vector of volumetric fraction of ice (-)
                        ! output: enthalpy
                        scalarCanairEnthalpy,        & ! intent(out):  enthalpy of the canopy air space (J m-3)
                        scalarCanopyEnthalpyTrial,   & ! intent(out):  enthalpy of the vegetation canopy (J m-3)
                        mLayerEnthalpyTrial,         & ! intent(out):  enthalpy of each snow+soil layer (J m-3)
                        ! output: error control
                        err,cmessage)                  ! intent(out): error control
        if(err/=0)then; message=trim(message)//trim(cmessage); return; endif

        ! *** compute volumetric heat capacity C_p = dH_T/dT
        call computHeatCap(&
                            ! input: control variables
                            nLayers,                   & ! intent(in): number of layers (soil+snow)
                            computeVegFlux,            & ! intent(in): flag to denote if computing the vegetation flux
                            canopyDepth,               & ! intent(in): canopy depth (m)
                            ! input data structures
                            mpar_data,                 & ! intent(in): model parameters
                            indx_data,                 & ! intent(in): model layer indices
                            diag_data,                 & ! intent(in): model diagnostic variables for a local HRU
                            ! input: state variables
                            scalarCanopyIceTrial,      & ! intent(in): trial value for canopy ice content (kg m-2)
                            scalarCanopyLiqTrial,      & ! intent(in): trial value for the liquid water on the vegetation canopy (kg m-2)
                            scalarCanopyTempTrial,     & ! intent(in): trial value of canopy temperature (K)
                            scalarCanopyTempPrev,      & ! intent(in): previous value of canopy temperature (K)
                            scalarCanopyEnthalpyTrial, & ! intent(in): trial enthalpy of the vegetation canopy (J m-3)
                            scalarCanopyEnthalpyPrev,  & ! intent(in): previous enthalpy of the vegetation canopy (J m-3)
                            mLayerVolFracIceTrial,     & ! intent(in): volumetric fraction of ice at the start of the sub-step (-)
                            mLayerVolFracLiqTrial,     & ! intent(in): volumetric fraction of liquid water at the start of the sub-step (-)
                            mLayerTempTrial,           & ! intent(in): trial temperature
                            mLayerTempPrev,            & ! intent(in): previous temperature
                            mLayerEnthalpyTrial,       & ! intent(in): trial enthalpy for snow and soil
                            mLayerEnthalpyPrev,        & ! intent(in): previous enthalpy for snow and soil
                            ! output
                            heatCapVegTrial,           & ! intent(out): volumetric heat capacity of vegetation canopy
                            mLayerHeatCapTrial,        & ! intent(out): heat capacity for snow and soil
                            ! output: error control
                            err,cmessage)                    ! intent(out): error control
        if(err/=0)then; message=trim(message)//trim(cmessage); return; endif
        ! to conserve energy compute finite difference approximation of (theta_ice)'
        if(dt_cur > 1e-14_rkind) then
          scalarCanopyIcePrime = ( scalarCanopyIceTrial - scalarCanopyIcePrev ) / dt_cur
          do concurrent (iLayer=1:nLayers)
                mLayerVolFracIcePrime(iLayer) = ( mLayerVolFracIceTrial(iLayer) - mLayerVolFracIcePrev(iLayer) ) / dt_cur
          end do
        endif ! if dt_cur is not too samll
      else ! if using closed formula of heat capacity
        call computHeatCapAnalytic(&
                          ! input: control variables
                          computeVegFlux,              & ! intent(in):   flag to denote if computing the vegetation flux
                          canopyDepth,                 & ! intent(in):   canopy depth (m)
                          ! input: state variables
                          scalarCanopyIceTrial,        & ! intent(in)
                          scalarCanopyLiqTrial,        & ! intent(in)
                          mLayerVolFracIceTrial,       & ! intent(in):   volumetric fraction of ice at the start of the sub-step (-)
                          mLayerVolFracLiqTrial,       & ! intent(in):   fraction of liquid water at the start of the sub-step (-)
                          ! input data structures
                          mpar_data,                   & ! intent(in):   model parameters
                          indx_data,                   & ! intent(in):   model layer indices
                          ! output
                          heatCapVegTrial,             & ! intent(out):  volumetric heat capacity of vegetation canopy
                          mLayerHeatCapTrial,          & ! intent(out):  volumetric heat capacity of soil and snow
                          ! output: error control
                          err,cmessage)                     ! intent(out):  error control
      endif

      ! compute multiplier of state vector
      call computStatMult(&
                    ! input
                    heatCapVegTrial,                  & ! intent(in):    volumetric heat capacity of vegetation canopy
                    mLayerHeatCapTrial,               & ! intent(in):    volumetric heat capacity of soil and snow
                    diag_data,                        & ! intent(in):    model diagnostic variables for a local HRU
                    indx_data,                        & ! intent(in):    indices defining model states and layers
                    ! output
                    sMul,                             & ! intent(out):   multiplier for state vector (used in the residual calculations)
                    err,cmessage)                       ! intent(out):   error control
      if(err/=0)then; message=trim(message)//trim(cmessage); return; endif  ! (check for errors)

      ! update thermal conductivity
      call computThermConduct(&
                          ! input: control variables
                          computeVegFlux,               & ! intent(in): flag to denote if computing the vegetation flux
                          canopyDepth,                  & ! intent(in): canopy depth (m)
                          ! input: state variables
                          scalarCanopyIceTrial,         & ! intent(in)
                          scalarCanopyLiqTrial,         & ! intent(in)
                          mLayerVolFracIceTrial,        & ! intent(in): volumetric fraction of ice at the start of the sub-step (-)
                          mLayerVolFracLiqTrial,        & ! intent(in): volumetric fraction of liquid water at the start of the sub-step (-)
                            ! input/output: data structures
                          mpar_data,                    & ! intent(in):    model parameters
                          indx_data,                    & ! intent(in):    model layer indices
                          prog_data,                    & ! intent(in):    model prognostic variables for a local HRU
                          diag_data,                    & ! intent(inout): model diagnostic variables for a local HRU
                          err,cmessage)               ! intent(out): error control
      if(err/=0)then; err=55; message=trim(message)//trim(cmessage); return; end if

    endif ! updateCp

    if(needCm)then
      ! compute C_m
      call computCm(&
                      ! input: control variables
                      computeVegFlux,           & ! intent(in): flag to denote if computing the vegetation flux
                      ! input: state variables
                      scalarCanopyTempTrial,    & ! intent(in)
                      mLayerTempTrial,          & ! intent(in): volumetric fraction of liquid water at the start of the sub-step (-)
                      mLayerMatricHeadTrial,    & ! intent(in)
                      ! input data structures
                      mpar_data,                & ! intent(in):    model parameters
                      indx_data,                & ! intent(in):    model layer indices
                      ! output
                      scalarCanopyCmTrial,      & ! intent(out):   Cm for vegetation
                      mLayerCmTrial,            & ! intent(out):   Cm for soil and snow
                      err,cmessage)                ! intent(out): error control
    else
      scalarCanopyCmTrial = 0._qp
      mLayerCmTrial = 0._qp
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
                    insideIDA,                 & ! intent(in): logical flag if inside Sundials solver
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
                    dBaseflow_dMatric,         & ! intent(out):   derivative in baseflow w.r.t. matric head (s-1), we will use it later in computJacobSundials
                    fluxVec,                   & ! intent(out):   flux vector (mixed units)
                    ! output: error control
                    err,cmessage)                ! intent(out):   error code and error message
    if(err/=0)then; message=trim(message)//trim(cmessage); return; end if  ! (check for errors)

    firstSplitOper = .true.

    ! compute soil compressibility (-) and its derivative w.r.t. matric head (m)
    ! NOTE: we already extracted trial matrix head and volumetric liquid water as part of the flux calculations
    call soilCmpresSundials(&
                    ! input:
                    ixRichards,                             & ! intent(in): choice of option for Richards' equation
                    ixBeg,ixEnd,                            & ! intent(in): start and end indices defining desired layers
                    mLayerMatricHeadPrime(1:nSoil),         & ! intent(in): matric head at the start of the time step (m)
                    mLayerVolFracLiqTrial(nSnow+1:nLayers), & ! intent(in): trial value for the volumetric liquid water content in each soil layer (-)
                    mLayerVolFracIceTrial(nSnow+1:nLayers), & ! intent(in): trial value for the volumetric ice content in each soil layer (-)
                    specificStorage,                        & ! intent(in): specific storage coefficient (m-1)
                    theta_sat,                              & ! intent(in): soil porosity (-)
                    ! output:
                    mLayerCompress,                         & ! intent(inout): compressibility of the soil matrix (-)
                    dCompress_dPsi,                         & ! intent(inout): derivative in compressibility w.r.t. matric head (m-1)
                    err,cmessage)                             ! intent(out): error code and error message
    if(err/=0)then; message=trim(message)//trim(cmessage); return; end if  ! (check for errors)

    ! compute the residual vector
    if (insideIDA)then
      dt1 = 1._qp ! always 1 for sundials since using Prime derivatives

      call computResidSundials(&
                      ! input: model control
                      dt1,                       & ! intent(in):    length of the residual time step (seconds)
                      nSnow,                     & ! intent(in):    number of snow layers
                      nSoil,                     & ! intent(in):    number of soil layers
                      nLayers,                   & ! intent(in):    total number of layers
                      ! input: flux vectors
                      sMul,                      & ! intent(in):    state vector multiplier (used in the residual calculations)
                      fluxVec,                   & ! intent(in):    flux vector
                      ! input: state variables (already disaggregated into scalars and vectors)
                      scalarCanopyTempTrial,     & ! intent(in):
                      mLayerTempTrial,           & ! intent(in)
                      scalarCanairTempPrime,     & ! intent(in):    Prime value for the temperature of the canopy air space (K)
                      scalarCanopyTempPrime,     & ! intent(in):    Prime value for the temperature of the vegetation canopy (K)
                      scalarCanopyWatPrime,      &
                      mLayerTempPrime,           & ! intent(in):    Prime value for the temperature of each snow and soil layer (K)
                      scalarAquiferStoragePrime, & ! intent(in):    Prime value of storage of water in the aquifer (m)
                      ! input: diagnostic variables defining the liquid water and ice content (function of state variables)
                      scalarCanopyIcePrime,      & ! intent(in):    Prime value for the ice on the vegetation canopy (kg m-2)
                      scalarCanopyLiqPrime,      & ! intent(in):
                      mLayerVolFracIcePrime,     & ! intent(in):    Prime value for the volumetric ice in each snow and soil layer (-)
                      mLayerVolFracWatPrime,     &
                      mLayerVolFracLiqPrime,     &
                      scalarCanopyCmTrial,       & ! intent(in) Cm of vegetation canopy
                      mLayerCmTrial,             & ! intent(in) Cm of soil and snow
                      ! input: data structures
                      prog_data,                 & ! intent(in):    model prognostic variables for a local HRU
                      diag_data,                 & ! intent(in):    model diagnostic variables for a local HRU
                      flux_data,                 & ! intent(in):    model fluxes for a local HRU
                      indx_data,                 & ! intent(in):    index data
                      ! output
                      resSink,                   & ! intent(out):   additional (sink) terms on the RHS of the state equation
                      resVec,                    & ! intent(out):   residual vector
                      err,cmessage)                ! intent(out):   error control
      if(err/=0)then; message=trim(message)//trim(cmessage); return; end if  ! (check for errors)

    else !currently not using residuals outside Sundials!
      dt1 = 1._qp
    endif

  ! end association with the information in the data structures
  end associate

end subroutine eval8summaSundials


! **********************************************************************************************************
! public function eval8summa4IDA: compute the residual vector F(t,y,y') required for IDA solver
! **********************************************************************************************************
! Return values:
!    0 = success,
!    1 = recoverable error,
!   -1 = non-recoverable error
! ----------------------------------------------------------------
integer(c_int) function eval8summa4IDA(tres, sunvec_y, sunvec_yp, sunvec_r, user_data) &
      result(ierr) bind(C,name='eval8summa4IDA')

  !======= Inclusions ===========
  use, intrinsic :: iso_c_binding
  use fida_mod
  use fsundials_nvector_mod
  use fnvector_serial_mod
  use type4IDA

  !======= Declarations =========
  implicit none

  ! calling variables
  real(rkind), value          :: tres      ! current time                 t
  type(N_Vector)              :: sunvec_y  ! solution N_Vector    y
  type(N_Vector)              :: sunvec_yp ! derivative N_Vector  y'
  type(N_Vector)              :: sunvec_r  ! residual N_Vector    F(t,y,y')
  type(c_ptr), value          :: user_data ! user-defined data

  ! pointers to data in SUNDIALS vectors
  type(eqnsData), pointer     :: eqns_data ! equations data
  real(rkind), pointer        :: stateVec(:)
  real(rkind), pointer        :: stateVecPrime(:)
  real(rkind), pointer        :: rVec(:)
  logical(lgt)                :: feasible
  integer(i4b)                :: retval
  real(c_double)              :: stepsize_next(1)
  !======= Internals ============

  ! get equations data from user-defined data
  call c_f_pointer(user_data, eqns_data)

  ! get data arrays from SUNDIALS vectors
  stateVec(1:eqns_data%nState)  => FN_VGetArrayPointer(sunvec_y)
  stateVecPrime(1:eqns_data%nState) => FN_VGetArrayPointer(sunvec_yp)
  rVec(1:eqns_data%nState)  => FN_VGetArrayPointer(sunvec_r)

  retval = FIDAGetCurrentStep(eqns_data%ida_mem, stepsize_next)
  if (retval /= 0) then
    print *, 'Error in FIDAGetCurrentStep, retval = ', retval, '; halting'
    stop 1
  end if

  ! compute the flux and the residual vector for a given state vector
  call eval8summaSundials(&
                ! input: model control
                stepsize_next(1),                  & ! intent(in):    current stepsize
                eqns_data%dt,                      & ! intent(in):    data step
                eqns_data%nSnow,                   & ! intent(in):    number of snow layers
                eqns_data%nSoil,                   & ! intent(in):    number of soil layers
                eqns_data%nLayers,                 & ! intent(in):    number of layers
                eqns_data%nState,                  & ! intent(in):    number of state variables in the current subset
                .true.,                            & ! intent(in):    inside Sundials solver
                eqns_data%firstSubStep,            & ! intent(in): flag to indicate if we are processing the first sub-step
                eqns_data%firstFluxCall,           & ! intent(inout): flag to indicate if we are processing the first flux call
                eqns_data%firstSplitOper,          & ! intent(inout):    flag to indicate if we are processing the first flux call in a splitting operation
                eqns_data%computeVegFlux,          & ! intent(in):    flag to indicate if we need to compute fluxes over vegetation
                eqns_data%scalarSolution,          & ! intent(in):    flag to indicate the scalar solution
                ! input: state vectors
                stateVec,                          & ! intent(in):    model state vector
                stateVecPrime,                     & ! intent(in):    model state vector
                eqns_data%sMul,                    & ! intent(inout): state vector multiplier (used in the residual calculations)
                ! input: data structures
                model_decisions,                   & ! intent(in):    model decisions
                eqns_data%lookup_data,             & ! intent(in):    lookup data
                eqns_data%type_data,               & ! intent(in):    type of vegetation and soil
                eqns_data%attr_data,               & ! intent(in):    spatial attributes
                eqns_data%mpar_data,               & ! intent(in):    model parameters
                eqns_data%forc_data,               & ! intent(in):    model forcing data
                eqns_data%bvar_data,               & ! intent(in):    average model variables for the entire basin
                eqns_data%prog_data,               & ! intent(in):    model prognostic variables for a local HRU
                ! input-output: data structures
                eqns_data%indx_data,               & ! intent(inou):  index data
                eqns_data%diag_data,               & ! intent(inout): model diagnostic variables for a local HRU
                eqns_data%flux_data,               & ! intent(inout): model fluxes for a local HRU (initial flux structure)
                eqns_data%deriv_data,              & ! intent(inout): derivatives in model fluxes w.r.t. relevant state variables
                 ! input-output: here we need to pass some extra variables that do not get updated in in the Sundials loops
                eqns_data%scalarCanopyTempTrial,   & ! intent(in):    trial value of canopy temperature (K)
                eqns_data%scalarCanopyTempPrev,    & ! intent(in):    previous value of canopy temperature (K)
                eqns_data%scalarCanopyIceTrial,    & ! intent(out):   trial value for mass of ice on the vegetation canopy (kg m-2)
                eqns_data%scalarCanopyIcePrev,     & ! intent(in):    value for mass of ice on the vegetation canopy (kg m-2)
                eqns_data%scalarCanopyLiqTrial,    & ! intent(out):   trial value of canopy liquid water (kg m-2)
                eqns_data%scalarCanopyLiqPrev,     & ! intent(in):    value of canopy liquid water (kg m-2)
                eqns_data%scalarCanopyEnthalpyTrial,& ! intent(out):  trial value for enthalpy of the vegetation canopy (J m-3)
                eqns_data%scalarCanopyEnthalpyPrev, & ! intent(in):   value for enthalpy of the vegetation canopy (J m-3)
                eqns_data%mLayerTempTrial,         & ! intent(out):   trial vector of layer temperature (K)
                eqns_data%mLayerTempPrev,          & ! intent(in):    vector of layer temperature (K)
                eqns_data%mLayerMatricHeadLiqTrial,& ! intent(out):   trial value for liquid water matric potential (m)
                eqns_data%mLayerMatricHeadTrial,   & ! intent(out):   trial value for total water matric potential (m)
                eqns_data%mLayerMatricHeadPrev,    & ! intent(in):    value for total water matric potential (m)
                eqns_data%mLayerVolFracWatTrial,   & ! intent(out):   trial vector of volumetric total water content (-)
                eqns_data%mLayerVolFracWatPrev,    & ! intent(in):    vector of volumetric total water content (-)
                eqns_data%mLayerVolFracIceTrial,   & ! intent(out):   trial vector of volumetric ice water content (-)
                eqns_data%mLayerVolFracIcePrev,    & ! intent(in):    vector of volumetric ice water content (-)
                eqns_data%mLayerVolFracLiqTrial,   & ! intent(out):   trial vector of volumetric liquid water content (-)
                eqns_data%mLayerVolFracLiqPrev,    & ! intent(in):    vector of volumetric liquid water content (-)
                eqns_data%scalarAquiferStorageTrial,& ! intent(out):  trial value of storage of water in the aquifer (m)
                eqns_data%scalarAquiferStoragePrev, & ! intent(in):   value of storage of water in the aquifer (m)
                eqns_data%mLayerEnthalpyPrev,      & ! intent(in):    vector of enthalpy for snow+soil layers (J m-3)
                eqns_data%mLayerEnthalpyTrial,     & ! intent(out):   trial vector of enthalpy for snow+soil layers (J m-3)
                ! input-output: baseflow
                eqns_data%ixSaturation,            & ! intent(inout): index of the lowest saturated layer
                eqns_data%dBaseflow_dMatric,       & ! intent(out):   derivative in baseflow w.r.t. matric head (s-1)
                 ! output: flux and residual vectors
                feasible,                          & ! intent(out):   flag to denote the feasibility of the solution
                eqns_data%fluxVec,                 & ! intent(out):   flux vector
                eqns_data%resSink,                 & ! intent(out):   additional (sink) terms on the RHS of the state equation
                rVec,                              & ! intent(out):   residual vector
                eqns_data%err,eqns_data%message)     ! intent(out):   error control

  if(eqns_data%err > 0)then; eqns_data%message=trim(eqns_data%message); ierr=-1; return; endif
  if(eqns_data%err < 0)then; eqns_data%message=trim(eqns_data%message); ierr=1; return; endif
  if(.not.feasible)then; eqns_data%message=trim(eqns_data%message)//'state vector not feasible'; ierr = 1; return; endif

  ! return success
  ierr = 0
  return

end function eval8summa4IDA


end module eval8summaSundials_module

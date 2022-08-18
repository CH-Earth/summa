
module eval8JacDAE_module

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
USE globalData,only: ixFullMatrix   ! named variable for the full Jacobian matrix
USE globalData,only: ixBandMatrix   ! named variable for the band diagonal matrix
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
 singleBasin                   ! single groundwater store over the entire basin

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
public::eval8JacDAE

contains

 ! **********************************************************************************************************
 ! public subroutine eval8JacDAE: compute the Jacobian matrix
 ! **********************************************************************************************************
 subroutine eval8JacDAE(&
                       ! input: model control
                       cj,                      & ! intent(in):    this scalar changes whenever the step size or method order changes
                       dt,                      & ! intent(in):    time step
                       nSnow,                   & ! intent(in):    number of snow layers
                       nSoil,                   & ! intent(in):    number of soil layers
                       nLayers,                 & ! intent(in):    total number of layers
                       ixMatrix,                & ! intent(in):    form of the Jacobian matrix
                       computeVegFlux,          & ! intent(in):    flag to indicate if we need to compute fluxes over vegetation
                       scalarSolution,          & ! intent(in):    flag to indicate the scalar solution
                       ! input: state vectors
                       stateVec,                & ! intent(in):    model state vector
                       stateVecPrime,           & ! intent(in):    derivative of model state vector
                       sMul,                    & ! intent(in):    state vector multiplier (used in the residual calculations)
                       ! input: data structures
                       model_decisions,         & ! intent(in):    model decisions
                       mpar_data,               & ! intent(in):    model parameters
                       prog_data,               & ! intent(in):    model prognostic variables for a local HRU
                       ! input-output: data structures
                       indx_data,               & ! intent(inout): index data
                       diag_data,               & ! intent(inout): model diagnostic variables for a local HRU
                       deriv_data,              & ! intent(inout): derivatives in model fluxes w.r.t. relevant state variables
                       ! input: baseflow
                       dBaseflow_dMatric,       & ! intent(in):    derivative in baseflow w.r.t. matric head (s-1)
                       ! output: flux and residual vectors
                       dMat,                    & ! intent(inout): diagonal of Jacobian Matrix
                       Jac,                     & ! intent(out):   jacobian matrix
                       err,message)               ! intent(out):   error control
 ! --------------------------------------------------------------------------------------------------------------------------------
 ! provide access to subroutines
 USE varExtrSundials_module, only:varExtract2                    ! extract variables from the state vector
 USE varExtrSundials_module, only:varExtractSundials
 USE updateVars4JacDAE_module, only:updateVars4JacDAE           ! update prognostic variables
 USE computJacDAE_module,only:computJacDAE
 implicit none
 ! --------------------------------------------------------------------------------------------------------------------------------
 ! --------------------------------------------------------------------------------------------------------------------------------
 ! input: model control
 real(rkind),intent(in)             :: cj                     ! this scalar changes whenever the step size or method order changes
 real(rkind),intent(in)             :: dt                     ! time step
 integer(i4b),intent(in)            :: nSnow                  ! number of snow layers
 integer(i4b),intent(in)            :: nSoil                  ! number of soil layers
 integer(i4b),intent(in)            :: nLayers                ! total number of layers
 integer(i4b)                       :: ixMatrix               ! form of matrix (band diagonal or full matrix)
 logical(lgt),intent(in)            :: computeVegFlux         ! flag to indicate if computing fluxes over vegetation
 logical(lgt),intent(in)            :: scalarSolution         ! flag to denote if implementing the scalar solution
 ! input: state vectors
 real(rkind),intent(in)             :: stateVec(:)            ! model state vector
 real(rkind),intent(in)             :: stateVecPrime(:)       ! model state vector
 real(qp),intent(in)                :: sMul(:)   ! NOTE: qp   ! state vector multiplier (used in the residual calculations)
 ! input: data structures
 type(model_options),intent(in)     :: model_decisions(:)     ! model decisions
 type(var_dlength),  intent(in)     :: mpar_data              ! model parameters
 type(var_dlength),  intent(in)     :: prog_data              ! prognostic variables for a local HRU
 ! output: data structures
 type(var_ilength),intent(inout)    :: indx_data              ! indices defining model states and layers
 type(var_dlength),intent(inout)    :: diag_data              ! diagnostic variables for a local HRU
 type(var_dlength),intent(inout)    :: deriv_data             ! derivatives in model fluxes w.r.t. relevant state variables
 ! input-output: baseflow
 real(rkind),intent(in)             :: dBaseflow_dMatric(:,:) ! derivative in baseflow w.r.t. matric head (s-1)
 ! output: Jacobian
 real(rkind), intent(inout)         :: dMat(:)
 real(rkind), intent(out)           :: Jac(:,:)               ! jacobian matrix
 ! output: error control
 integer(i4b),intent(out)           :: err                    ! error code
 character(*),intent(out)           :: message                ! error message
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
 ! other local variables
 character(LEN=256)              :: cmessage                     ! error message of downwind routine
 real(rkind)                        :: dt1

 ! --------------------------------------------------------------------------------------------------------------------------------
 ! association to variables in the data structures
 ! --------------------------------------------------------------------------------------------------------------------------------
 associate(&
 ! model decisions
 ixRichards              => model_decisions(iLookDECISIONS%f_Richards)%iDecision   ,&  ! intent(in):  [i4b]   index of the form of Richards' equation
 ! soil parameters
 theta_sat               => mpar_data%var(iLookPARAM%theta_sat)%dat                ,&  ! intent(in):  [dp(:)] soil porosity (-)
 specificStorage         => mpar_data%var(iLookPARAM%specificStorage)%dat(1)       ,&  ! intent(in):  [dp]    specific storage coefficient (m-1)
 ! model state variables
 mLayerVolFracLiq        => prog_data%var(iLookPROG%mLayerVolFracLiq)%dat          ,&  ! intent(in):  [dp(:)] volumetric fraction of liquid water (-)
 mLayerVolFracIce        => prog_data%var(iLookPROG%mLayerVolFracIce)%dat          ,&  ! intent(in):  [dp(:)] volumetric fraction of ice (-)
 mLayerMatricHeadLiq     => diag_data%var(iLookDIAG%mLayerMatricHeadLiq)%dat       ,&  ! intent(in):  [dp(:)] liquid water matric potential (m)
 ixGroundwater           => model_decisions(iLookDECISIONS%groundwatr)%iDecision    &
 ) ! association to variables in the data structures
 ! --------------------------------------------------------------------------------------------------------------------------------
 ! initialize error control
 err=0; message="eval8JacDAE/"


 ! extract variables from the model state vector
 call varExtract2(&
                 ! input
                 stateVec,            & ! intent(in):    model state vector (mixed units)
                 diag_data,                & ! intent(in):    model diagnostic variables for a local HRU
                 prog_data,                & ! intent(in):    model prognostic variables for a local HRU
                 indx_data,                & ! intent(in):    indices defining model states and layers
                 ! output: variables for the vegetation canopy
                 scalarCanairTempTrial,    & ! intent(out):   trial value of canopy air temperature (K)
                 scalarCanopyTempTrial,    & ! intent(out):   trial value of canopy temperature (K)
                 scalarCanopyWatTrial,     & ! intent(out):   trial value of canopy total water (kg m-2)
                 scalarCanopyLiqTrial,     & ! intent(out):   trial value of canopy liquid water (kg m-2)
                 ! output: variables for the snow-soil domain
                 mLayerTempTrial,          & ! intent(out):   trial vector of layer temperature (K)
                 mLayerVolFracWatTrial,    & ! intent(out):   trial vector of volumetric total water content (-)
                 mLayerVolFracLiqTrial,    & ! intent(out):   trial vector of volumetric liquid water content (-)
                 mLayerMatricHeadTrial,    & ! intent(out):   trial vector of total water matric potential (m)
                 mLayerMatricHeadLiqTrial, & ! intent(out):   trial vector of liquid water matric potential (m)
                 ! output: variables for the aquifer
                 scalarAquiferStorageTrial,& ! intent(out):   trial value of storage of water in the aquifer (m)
                 ! output: error control
                 err,cmessage)               ! intent(out):   error control
 if(err/=0)then; message=trim(message)//trim(cmessage); return; end if  ! (check for errors)



 ! extract derivative of variables from derivative of the model state vector
 call varExtractSundials(&
                 ! input
                 stateVecPrime,            & ! intent(in):    derivative of model state vector (mixed units)
                 diag_data,                & ! intent(in):    model diagnostic variables for a local HRU
                 prog_data,                & ! intent(in):    model prognostic variables for a local HRU
                 indx_data,                & ! intent(in):    indices defining model states and layers
                 ! output: variables for the vegetation canopy
                 scalarCanairTempPrime,    & ! intent(out):   derivative of canopy air temperature (K)
                 scalarCanopyTempPrime,    & ! intent(out):   derivative of canopy temperature (K)
                 scalarCanopyWatPrime,     & ! intent(out):   derivative of canopy total water (kg m-2)
                 scalarCanopyLiqPrime,     & ! intent(out):   derivative of canopy liquid water (kg m-2)
                 ! output: variables for the snow-soil domain
                 mLayerTempPrime,          & ! intent(out):   derivative of layer temperature (K)
                 mLayerVolFracWatPrime,    & ! intent(out):   derivative of volumetric total water content (-)
                 mLayerVolFracLiqPrime,    & ! intent(out):   derivative of volumetric liquid water content (-)
                 mLayerMatricHeadPrime,    & ! intent(out):   derivative of total water matric potential (m)
                 mLayerMatricHeadLiqPrime, & ! intent(out):   derivative of liquid water matric potential (m)
                 ! output: variables for the aquifer
                 scalarAquiferStoragePrime,& ! intent(out):   derivative of storage of water in the aquifer (m)
                 ! output: error control
                 err,cmessage)               ! intent(out):   error control
 if(err/=0)then; message=trim(message)//trim(cmessage); return; end if  ! (check for errors)




 call updateVars4JacDAE(&
                 ! input
                 dt,                                        & ! intent(in):    time step
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
                 scalarCanopyTempPrime,                     & ! intent(inout): trial value of canopy temperature (K)
                 scalarCanopyWatPrime,                      & ! intent(inout): trial value of canopy total water (kg m-2)
                 scalarCanopyLiqPrime,                      & ! intent(inout): trial value of canopy liquid water (kg m-2)
                 scalarCanopyIcePrime,                      & ! intent(inout): trial value of canopy ice content (kg m-2
                 ! output: variables for the snow-soil domain
                 mLayerTempTrial,                           & ! intent(inout): trial vector of layer temperature (K)
                 mLayerVolFracWatTrial,                     & ! intent(inout): trial vector of volumetric total water content (-)
                 mLayerVolFracLiqTrial,                     & ! intent(inout): trial vector of volumetric liquid water content (-)
                 mLayerVolFracIceTrial,                     & ! intent(inout): trial vector of volumetric ice water content (-)
                 mLayerMatricHeadTrial,                     & ! intent(inout): trial vector of total water matric potential (m)
                 mLayerMatricHeadLiqTrial,                  & ! intent(inout): trial vector of liquid water matric potential (m)
                 mLayerTempPrime,                           & !
                 mLayerVolFracWatPrime,                     & ! intent(inout): Prime vector of volumetric total water content (-)
                 mLayerVolFracLiqPrime,                     & ! intent(inout): Prime vector of volumetric liquid water content (-)
                 mLayerVolFracIcePrime,                     & !
                 mLayerMatricHeadPrime,                     & ! intent(inout): Prime vector of total water matric potential (m)
                 mLayerMatricHeadLiqPrime,                  & ! intent(inout): Prime vector of liquid water matric potential (m)
                 ! output: error control
                 err,cmessage)                                ! intent(out):   error control
 if(err/=0)then; message=trim(message)//trim(cmessage); return; end if  ! (check for errors)



 ! -----
 ! * compute the Jacobian matrix...
 ! --------------------------------

 ! compute the analytical Jacobian matrix
 ! NOTE: The derivatives were computed in the previous call to computFlux
 !       This occurred either at the call to eval8DAE at the start of sysSolveSundials
 !        or in the call to eval8DAE in the previous iteration
 dt1 = 1._qp
 call computJacDAE(&
                  ! input: model control
                  cj,                             & ! intent(in):    this scalar changes whenever the step size or method order changes
                  dt1,                            & ! intent(in):    length of the time step (seconds)
                  nSnow,                          & ! intent(in):    number of snow layers
                  nSoil,                          & ! intent(in):    number of soil layers
                  nLayers,                        & ! intent(in):    total number of layers
                  computeVegFlux,                 & ! intent(in):    flag to indicate if we need to compute fluxes over vegetation
                  (ixGroundwater==qbaseTopmodel), & ! intent(in):    flag to indicate if we need to compute baseflow
                  ixMatrix,                       & ! intent(in):    form of the Jacobian matrix
                  specificStorage,                & ! intent(in):    specific storage coefficient (m-1)
                  theta_sat,                      & ! intent(in):    soil porosity (-)
                  ixRichards,                     & ! intent(in):    choice of option for Richards' equation
                  ! input: data structures
                  indx_data,                      & ! intent(in):    index data
                  prog_data,                      & ! intent(in):    model prognostic variables for a local HRU
                  diag_data,                      & ! intent(in):    model diagnostic variables for a local HRU
                  deriv_data,                     & ! intent(in):    derivatives in model fluxes w.r.t. relevant state variables
                  dBaseflow_dMatric,              & ! intent(in):    derivative in baseflow w.r.t. matric head (s-1)
                  ! input: state variables
                  mLayerTempTrial,                & ! intent(in):    trial value for the temperature of each snow and soil layer (K)
                  mLayerTempPrime,                & ! intent(in)
                  mLayerMatricHeadPrime,          & ! intent(in)
                  mLayerMatricHeadLiqPrime,       & ! intent(in)
                  mLayerVolFracWatPrime,          & ! intent(in)
                  scalarCanopyTempTrial,          & ! intent(in)
                  scalarCanopyTempPrime,          & ! intent(in) derivative value for temperature of the vegetation canopy (K)
                  scalarCanopyWatPrime,           & ! intetn(in)
                  ! input-output: Jacobian and its diagonal
                  dMat,                           & ! intent(inout): diagonal of the Jacobian matrix
                  Jac,                            & ! intent(out):   Jacobian matrix
                  ! output: error control
                  err,cmessage)                     ! intent(out):   error code and error message
   if(err/=0)then; message=trim(message)//trim(cmessage); return; end if  ! (check for errors)



 ! end association with the information in the data structures
 end associate


 end subroutine eval8JacDAE
end module eval8JacDAE_module

module type4ida

! data types
USE nrtype
USE, intrinsic :: iso_c_binding

! provide access to the derived types to define the data structures
USE data_types,only:&
                    var_i,        & ! data vector (i4b)
                    var_d,        & ! data vector (rkind)
                    var_ilength,  & ! data vector with variable length dimension (i4b)
                    var_dlength,  & ! data vector with variable length dimension (rkind)
                    zLookup,      & ! data vector with variable length dimension (rkind)
                    model_options   ! defines the model decisions

implicit none

type data4ida
  type(c_ptr)                     :: ida_mem                         ! IDA memory
  real(rkind)                     :: dt                              ! data step
  integer(i4b)                    :: nSnow                           ! number of snow layers
  integer(i4b)                    :: nSoil                           ! number of soil layers
  integer(i4b)                    :: nLayers                         ! total number of layers
  integer(i4b)                    :: nState                          ! total number of state variables
  integer(i4b)                    :: ixMatrix                        ! form of matrix (dense or banded)
  logical(lgt)                    :: firstSubStep                    ! flag to indicate if we are processing the first sub-step
  logical(lgt)                    :: firstFluxCall                   ! flag to indicate if we are processing the first flux call
  logical(lgt)                    :: firstSplitOper                  ! flag to indicate if we are processing the first flux call in a splitting operation
  logical(lgt)                    :: computeVegFlux                  ! flag to indicate if computing fluxes over vegetation
  logical(lgt)                    :: scalarSolution                  ! flag to denote if implementing the scalar solution
  type(model_options),allocatable :: model_decisions(:)              ! model decisions
  type(zLookup)                   :: lookup_data                  ! lookup tables
  type(var_i)                     :: type_data                       ! type of vegetation and soil
  type(var_d)                     :: attr_data                       ! spatial attributes
  type(var_dlength)               :: mpar_data                       ! model parameters
  type(var_d)                     :: forc_data                       ! model forcing data
  type(var_dlength)               :: bvar_data                       ! model variables for the local basin
  type(var_dlength)               :: prog_data                       ! prognostic variables for a local HRU
  type(var_ilength)               :: indx_data                       ! indices defining model states and layers
  type(var_dlength)               :: diag_data                       ! diagnostic variables for a local HRU
  type(var_dlength)               :: flux_data                       ! model fluxes for a local HRU
  type(var_dlength)               :: deriv_data                      ! derivatives in model fluxes w.r.t. relevant state variables
  real(qp), allocatable           :: sMul(:)                         ! state vector multiplier (used in the residual calculations)
  real(rkind), allocatable        :: dMat(:)                         ! diagonal of the Jacobian matrix
  real(rkind), allocatable        :: fluxVec(:)                      ! flux vector
  real(qp), allocatable           :: resVec(:)                       ! residual vector
  real(qp), allocatable           :: resSink(:)                      ! additional (sink) terms on the RHS of the state equation
  real(rkind), allocatable        :: atol(:)                         ! vector of absolute tolerances
  real(rkind), allocatable        :: rtol(:)                         ! vector of relative tolerances
  integer(i4b)                    :: ixSaturation                    ! index of the lowest saturated layer
  real(rkind), allocatable        :: dBaseflow_dMatric(:,:)          ! derivative in baseflow w.r.t. matric head (s-1)
  integer(i4b)                    :: err                             ! error code
  character(len=50)               :: message                         ! error message
  real(rkind)                     :: scalarCanopyTempPrev            ! previous value for temperature of the vegetation canopy (K)
  real(rkind), allocatable        :: mLayerTempPrev(:)               ! previous vector of layer temperature (K)
  real(rkind), allocatable        :: mLayerMatricHeadPrev(:)         ! previous value for total water matric potential (m)
  real(rkind)                     :: scalarCanopyEnthalpyTrial       ! trial value for enthalpy of the vegetation canopy (J m-2)
  real(rkind)                     :: scalarCanopyTempTrial           ! trial value for temperature of the vegetation canopy (K)
  real(rkind)                     :: scalarCanopyWatTrial            ! trial value for mass of total water on the vegetation canopy (kg m-2)
  real(rkind), allocatable        :: mLayerTempTrial(:)              ! trial vector of layer temperature (K)
  real(rkind), allocatable        :: mLayerMatricHeadTrial(:)        ! trial value for total water matric potential (m)
  real(rkind)                     :: scalarCanopyTempPrime           ! prime value for temperature of the vegetation canopy (K s-1)
  real(rkind)                     :: scalarCanopyWatPrime            ! prime value for mass of total water on the vegetation canopy (kg m-2 s-1)
  real(rkind), allocatable        :: mLayerTempPrime(:)              ! prime vector of temperature of each snow and soil layer (K s-1)
  real(rkind), allocatable        :: mLayerMatricHeadPrime(:)        ! prime vector of matric head of each snow and soil layer (m s-1)
  real(rkind), allocatable        :: mLayerVolFracWatPrime(:)        ! prime vector of volumetric total water content of each snow and soil layer (s-1)
 end type data4ida


end module type4ida






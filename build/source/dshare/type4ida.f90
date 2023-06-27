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
  type(c_ptr)                     :: ida_mem                      ! IDA memory
  real(rkind)                     :: dt                           ! data step
  integer(i4b)                    :: nSnow                        ! number of snow layers
  integer(i4b)                    :: nSoil                        ! number of soil layers
  integer(i4b)                    :: nLayers                      ! total number of layers
  integer(i4b)                    :: nState                       ! total number of state variables
  integer(i4b)                    :: ixMatrix                     ! form of matrix (dense or banded)
  logical(lgt)                    :: firstSubStep                 ! flag to indicate if we are processing the first sub-step
  logical(lgt)                    :: firstFluxCall                ! flag to indicate if we are processing the first flux call
  logical(lgt)                    :: firstSplitOper               ! flag to indicate if we are processing the first flux call in a splitting operation
  logical(lgt)                    :: computeVegFlux               ! flag to indicate if computing fluxes over vegetation
  logical(lgt)                    :: scalarSolution               ! flag to denote if implementing the scalar solution
  type(model_options),allocatable :: model_decisions(:)           ! model decisions
  type(zLookup)                   :: lookup_data                  ! lookup tables
  type(var_i)                     :: type_data                    ! type of vegetation and soil
  type(var_d)                     :: attr_data                    ! spatial attributes
  type(var_dlength)               :: mpar_data                    ! model parameters
  type(var_d)                     :: forc_data                    ! model forcing data
  type(var_dlength)               :: bvar_data                    ! model variables for the local basin
  type(var_dlength)               :: prog_data                    ! prognostic variables for a local HRU
  type(var_ilength)               :: indx_data                    ! indices defining model states and layers
  type(var_dlength)               :: diag_data                    ! diagnostic variables for a local HRU
  type(var_dlength)               :: flux_data                    ! model fluxes for a local HRU
  type(var_dlength)               :: deriv_data                   ! derivatives in model fluxes w.r.t. relevant state variables
   real(rkind)                    :: scalarCanopyTempTrial        ! trial value of canopy temperature (K)
  real(rkind)                     :: scalarCanopyTempPrev         ! previous value of canopy temperature (K)
  real(rkind)                     :: scalarCanopyIceTrial         ! trial value of canopy ice content (kg m-2)
  real(rkind)                     :: scalarCanopyIcePrev          ! value of canopy ice content (kg m-2) at previous step
  real(rkind)                     :: scalarCanopyLiqTrial         ! trial value of canopy ice content (kg m-2)
  real(rkind)                     :: scalarCanopyLiqPrev          ! value of canopy ice content (kg m-2) at previous step
  real(rkind)                     :: scalarCanopyEnthalpyTrial    ! trial enthalpy of the vegetation canopy (J m-3)
  real(rkind)                     :: scalarCanopyEnthalpyPrev     ! previous enthalpy of the vegetation canopy (J m-3)
  real(qp), allocatable           :: sMul(:)                      ! state vector multiplier (used in the residual calculations)
  real(rkind), allocatable        :: dMat(:)                      ! diagonal of the Jacobian matrix
  real(rkind), allocatable        :: fluxVec(:)                   ! flux vector
  real(qp), allocatable           :: resSink(:)                   ! additional (sink) terms on the RHS of the state equation
  real(rkind), allocatable        :: atol(:)                      ! vector of absolute tolerances
  real(rkind), allocatable        :: rtol(:)                      ! vector of relative tolerances
  real(rkind), allocatable        :: mLayerTempTrial(:)           ! trial vector of layer temperature (K)
  real(rkind), allocatable        :: mLayerTempPrev(:)            ! vector of layer temperature (K) at previous step
  real(rkind), allocatable        :: mLayerMatricHeadLiqTrial(:)  ! trial value for liquid water matric potential (m)
  real(rkind), allocatable        :: mLayerMatricHeadTrial(:)     ! trial vector of total water matric potential (m)
  real(rkind), allocatable        :: mLayerMatricHeadPrev(:)      ! vector of total water matric potential (m) at previous step
  real(rkind), allocatable        :: mLayerVolFracWatTrial(:)     ! trial value for volumetric fraction of total water (-)
  real(rkind), allocatable        :: mLayerVolFracWatPrev(:)      ! value for volumetric fraction of total water (-) at previous step
  real(rkind), allocatable        :: mLayerVolFracIceTrial(:)     ! trial value for volumetric fraction of ice (-)
  real(rkind), allocatable        :: mLayerVolFracIcePrev(:)      ! value for volumetric fraction of ice (-) at previous step
  real(rkind), allocatable        :: mLayerVolFracLiqTrial(:)     ! trial value for volumetric fraction of liquid water (-)
  real(rkind), allocatable        :: mLayerVolFracLiqPrev(:)      ! value for volumetric fraction of liquid water (-) at previous step
  real(rkind)                     :: scalarAquiferStorageTrial    ! trial value of storage of water in the aquifer (m)
  real(rkind)                     :: scalarAquiferStoragePrev     ! value of storage of water in the aquifer (m) at previous step
  real(rkind), allocatable        :: mLayerEnthalpyTrial(:)       ! trial enthalpy of snow and soil (J m-3)
  real(rkind), allocatable        :: mLayerEnthalpyPrev(:)        ! enthalpy of snow and soil (J m-3) at previous step
  real(rkind), allocatable        :: mLayerTempPrime(:)           ! derivative value for temperature of each snow and soil layer (K)
  real(rkind), allocatable        :: mLayerMatricHeadPrime(:)     ! derivative value for matric head of each snow and soil layer (m)
  real(rkind), allocatable        :: mLayerMatricHeadLiqPrime(:)  ! derivative value for liquid water matric head of each snow and soil layer (m)
  real(rkind), allocatable        :: mLayerVolFracWatPrime(:)     ! derivative value for volumetric total water content of each snow and soil layer (-)
  real(rkind)                     :: scalarCanopyTempPrime        ! derivative value for temperature of the vegetation canopy (K)
  real(rkind)                     :: scalarCanopyWatPrime         ! derivative value for total water content of the vegetation canopy (kg m-2)
  real(rkind), allocatable        :: dBaseflow_dMatric(:,:)       ! derivative in baseflow w.r.t. matric head (s-1)
  integer(i4b)                    :: ixSaturation                 ! index of the lowest saturated layer
  integer(i4b)                    :: err                          ! error code
  character(len = 50)             :: message                      ! error message
end type data4ida


end module type4ida






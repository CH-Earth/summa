module type4kinsol
  ! data types
  USE nrtype
  USE, intrinsic :: iso_c_binding      

  USE data_types,only:&
                      var_i,        & ! data vector (i4b)
                      var_d,        & ! data vector (rkind)
                      var_ilength,  & ! data vector with variable length dimension (i4b)
                      var_dlength,  & ! data vector with variable length dimension (dp)
                      model_options   ! defines the model decisions
  implicit none

type kinsol_data
  real(rkind)                          :: dt
  integer(i4b)                         :: nSnow
  integer(i4b)                         :: nSoil
  integer(i4b)                         :: nLayers
  integer(i4b)                         :: nState
  logical(lgt)                         :: computeVegFlux        ! flag to indicate if computing fluxes over vegetation
  logical(lgt)                         :: computeBaseflow
  integer(i4b)                         :: ixMatrix
  logical(lgt)                         :: firstSubStep           ! flag to indicate if we are processing the first sub-step
  logical(lgt)                         :: firstFluxCall          ! flag to indicate if we are processing the first flux call
  logical(lgt)                         :: firstSplitOper         ! flag to indicate if we are processing the first flux call in a splitting operation
  logical(lgt)                         :: scalarSolution         ! flag to denote if implementing the scalar solution
  type(var_i)                          :: type_data              ! type of vegetation and soil
  type(var_d)                          :: attr_data              ! spatial attributes
  type(var_dlength)                    :: mpar_data              ! model parameters
  type(var_d)                          :: forc_data              ! model forcing data
  type(var_dlength)                    :: bvar_data              ! model variables for the local basin
  type(var_dlength)                    :: flux_data              ! model fluxes for a local HRU
  integer(i4b)                         :: ixSaturation           ! index of the lowest saturated layer (NOTE: only computed on the first iteration)
  type(var_ilength)                    :: indx_data 
  type(var_dlength)                    :: prog_data
  type(var_dlength)                    :: diag_data
  type(var_dlength)                    :: deriv_data
  logical(lgt)                         :: feasible               ! flag to denote the feasibility of the solution
  real(rkind)                          :: fEval                  ! function evaluation
  real(rkind),allocatable              :: sMul(:)   ! NOTE: qp   ! state vector multiplier (used in the residual calculations)
  real(rkind),allocatable              :: fScale(:)              ! function scaling vector
  real(rkind),allocatable              :: fluxVec(:)             ! flux vector
  real(rkind),allocatable              :: resSink(:)             ! sink terms on the RHS of the flux equation
  type(model_options),allocatable      :: model_decisions(:)     ! model decisions
  real(rkind), allocatable             :: dBaseflow_dMatric(:,:)
  real(rkind), allocatable             :: dMat(:)
  real(rkind), allocatable             :: stateVecPrev(:)        ! state vector from the previous iteration to help with infeasibility
  logical(lgt)                         :: firstStateiteration   ! flag to denote if we computed an iteration so we know to save the state
end type kinsol_data


end module type4kinsol


module evalJac4IDA_module


  !======= Inclusions ===========
  use, intrinsic :: iso_c_binding
  use nrtype
  use type4IDA
  USE globalData,only:model_decisions        ! model decision structure
  ! provide access to the derived types to define the data structures
  USE data_types,only:&
                    var_i,        & ! data vector (i4b)
                    var_d,        & ! data vector (rkind)
                    var_ilength,  & ! data vector with variable length dimension (i4b)
                    var_dlength,  & ! data vector with variable length dimension (rkind)
                    model_options   ! defines the model decisions



  ! privacy
  implicit none
  private
  public::evalJac4IDA


contains

  ! **********************************************************************************************************
  ! public function evalJac4IDA: the interface to compute the Jacobian matrix dF/dy + c dF/dy' for IDA solver
  ! **********************************************************************************************************
  ! Return values:
  !    0 = success,
  !    1 = recoverable error,
  !   -1 = non-recoverable error
  ! ----------------------------------------------------------------
  integer(c_int) function evalJac4IDA(t, cj, sunvec_y, sunvec_yp, sunvec_r, &
                      sunmat_J, user_data, sunvec_temp1, sunvec_temp2, sunvec_temp3) &
                      result(ierr) bind(C,name='evalJac4IDA')

    !======= Inclusions ===========
    use, intrinsic :: iso_c_binding
    use fsundials_nvector_mod
    use fsundials_matrix_mod
    use fnvector_serial_mod
    use fsunmatrix_dense_mod
    use nrtype
    use type4IDA
    use eval8JacDAE_module,only:eval8JacDAE    ! compute Jacobian matrix
    !======= Declarations =========
    implicit none

    ! calling variables
    real(rkind), value            :: t              ! current time
    real(rkind), value            :: cj             ! step size scaling factor
    type(N_Vector)                :: sunvec_y       ! solution N_Vector
    type(N_Vector)                :: sunvec_yp      ! derivative N_Vector
    type(N_Vector)                :: sunvec_r       ! residual N_Vector
    type(SUNMatrix)               :: sunmat_J       ! Jacobian SUNMatrix
    type(c_ptr), value            :: user_data      ! user-defined data
    type(N_Vector)                :: sunvec_temp1   ! temporary N_Vector
    type(N_Vector)                :: sunvec_temp2   ! temporary N_Vector
    type(N_Vector)                :: sunvec_temp3   ! temporary N_Vector

    ! pointers to data in SUNDIALS vectors
    real(rkind), pointer          :: stateVec(:)    ! state vector
    real(rkind), pointer          :: stateVecPrime(:)! derivative of the state vector
    real(rkind), pointer          :: rVec(:)        ! residual vector
    real(rkind), pointer          :: Jac(:,:)       ! Jacobian matrix
    type(eqnsData), pointer       :: eqns_data      ! equations data



    !======= Internals ============

    ! get equations data from user-defined data
    call c_f_pointer(user_data, eqns_data)


    ! get data arrays from SUNDIALS vectors
    stateVec(1:eqns_data%nState)  => FN_VGetArrayPointer(sunvec_y)
    stateVecPrime(1:eqns_data%nState) => FN_VGetArrayPointer(sunvec_yp)
    rVec(1:eqns_data%nState)  => FN_VGetArrayPointer(sunvec_r)
    Jac(1:eqns_data%nState, 1:eqns_data%nState) => FSUNDenseMatrix_Data(sunmat_J)

    ! compute Jacobian matrix
    call eval8JacDAE(&
                 ! input: model control
                 cj,                                & ! intent(in):    this scalar changes whenever the step size or method order changes
                 eqns_data%dt,                      & ! intent(in):    data step
                 eqns_data%nSnow,                   & ! intent(in):    number of snow layers
                 eqns_data%nSoil,                   & ! intent(in):    number of soil layers
                 eqns_data%nLayers,                 & ! intent(in):    number of layers
                 eqns_data%ixMatrix,                & ! intent(in):    type of matrix (dense or banded)
                 eqns_data%computeVegFlux,          & ! intent(in):    flag to indicate if we need to compute fluxes over vegetation
                 eqns_data%scalarSolution,          & ! intent(in):    flag to indicate the scalar solution
                 ! input: state vectors
                 stateVec,                          & ! intent(in):    model state vector
                 stateVecPrime,                     & ! intent(in):    model state vector
                 eqns_data%sMul,                    & ! intent(in):    state vector multiplier (used in the residual calculations)
                 ! input: data structures
                 model_decisions,                   & ! intent(in):    model decisions
                 eqns_data%mpar_data,               & ! intent(in):    model parameters
                 eqns_data%prog_data,               & ! intent(in):    model prognostic variables for a local HRU
                 ! input-output: data structures
                 eqns_data%indx_data,               & ! intent(inou):  index data
                 eqns_data%diag_data,               & ! intent(inout): model diagnostic variables for a local HRU
                 eqns_data%deriv_data,              & ! intent(inout): derivatives in model fluxes w.r.t. relevant state variables
                 ! input: baseflow
                 eqns_data%dBaseflow_dMatric,       & ! intent(in):    derivative in baseflow w.r.t. matric head (s-1)
                 ! output
                 eqns_data%dMat,                    & ! intetn(inout): diagonal of the Jacobian matrix
                 Jac,                               & ! intent(out):   Jacobian matrix
                 eqns_data%err,eqns_data%message)     ! intent(out):   error control

   if(eqns_data%err > 0)then; eqns_data%message=trim(eqns_data%message); ierr=-1; return; endif
   if(eqns_data%err < 0)then; eqns_data%message=trim(eqns_data%message); ierr=1; return; endif

   ! return success
   ierr = 0
   return



 end function evalJac4IDA


end module evalJac4IDA_module

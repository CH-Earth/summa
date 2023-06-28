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

module summaSolve4kinsol_module

    !======= Inclusions ===========
USE, intrinsic :: iso_c_binding
USE nrtype
USE type4kinsol

! access the global print flag
USE globalData,only:globalPrintFlag

! access missing values
USE globalData,only:integerMissing  ! missing integer
USE globalData,only:realMissing     ! missing double precision number

! access matrix information
USE globalData,only: ixFullMatrix   ! named variable for the full Jacobian matrix
USE globalData,only: ixBandMatrix   ! named variable for the band diagonal matrix
USE globalData,only: ku             ! number of super-diagonal bands
USE globalData,only: kl             ! number of sub-diagonal bands

! global metadata
USE globalData,only:flux_meta       ! metadata on the model fluxes

! constants
USE multiconst,only: Tfreeze        ! temperature at freezing              (K)

! provide access to indices that define elements of the data structures
USE var_lookup,only:iLookPROG       ! named variables for structure elements
USE var_lookup,only:iLookDIAG       ! named variables for structure elements
USE var_lookup,only:iLookDECISIONS  ! named variables for elements of the decision structure
USE var_lookup,only:iLookDERIV     ! named variables for structure elements
USE var_lookup,only:iLookFLUX       ! named variables for structure elements
USE var_lookup,only:iLookPARAM      ! named variables for structure elements
USE var_lookup,only:iLookINDEX      ! named variables for structure elements

! provide access to the derived types to define the data structures
USE data_types,only:&
                    var_i,        & ! data vector (i4b)
                    var_d,        & ! data vector (rkind)
                    var_ilength,  & ! data vector with variable length dimension (i4b)
                    var_dlength,  & ! data vector with variable length dimension (rkind)
                    zLookup,      & ! lookup tables
                    model_options   ! defines the model decisions

! look-up values for the choice of groundwater parameterization
USE mDecisions_module,only:       &
  qbaseTopmodel,                  & ! TOPMODEL-ish baseflow parameterization
  bigBucket,                      & ! a big bucket (lumped aquifer model)
  noExplicit                         ! no explicit groundwater parameterization
                  
! look-up values for method used to compute derivative
USE mDecisions_module,only:       &
  numerical,                      & ! numerical solution
  analytical                        ! analytical solution

! privacy
 implicit none
 private::setInitialCondition
 private::setSolverParams
 private::getErrMessage
 public::summaSolve4kinsol

contains


!-------------------
! * public subroutine summaSolve4kinsol: solve F(y) = 0 by KINSOL (y is the state vector)
! ------------------
subroutine summaSolve4kinsol(&
                      dt_cur,                  & ! intent(in):    current stepsize
                      dt,                      & ! intent(in):    data time step
                      fScale,                  & ! intent(in):    characteristic scale of the function evaluations (mixed units)
                      xScale,                  & ! intent(in):    characteristic scale of the state vector (mixed units)
                      nSnow,                   & ! intent(in):    number of snow layers
                      nSoil,                   & ! intent(in):    number of soil layers
                      nLayers,                 & ! intent(in):    total number of layers
                      nStat,                   & ! intent(in):    total number of state variables
                      ixMatrix,                & ! intent(in):    type of matrix (dense or banded)
                      firstSubStep,            & ! intent(in):    flag to indicate if we are processing the first sub-step
                      computeVegFlux,          & ! intent(in):    flag to indicate if we need to compute fluxes over vegetation
                      scalarSolution,          & ! intent(in):    flag to indicate the scalar solution
                      ! input: state vectors
                      stateVecInit,            & ! intent(in):    initial state vector
                      sMul,                    & ! intent(inout): state vector multiplier (USEd in the residual calculations)
                      dMat,                    & ! intent(inout): diagonal of the Jacobian matrix (excludes fluxes)
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
                      indx_data,               & ! intent(inout): index data                      diag_data,               & ! intent(inout): model diagnostic variables for a local HRU
                      diag_data,               & ! intent(inout): model diagnostic variables for a local HRU
                      flux_data,               & ! intent(inout): model fluxes for a local HRU
                      deriv_data,              & ! intent(inout): derivatives in model fluxes w.r.t. relevant state variables
                      ! output
                      ixSaturation,            & ! intent(inout) index of the lowest saturated layer (NOTE: only computed on the first iteration)
                      kinsolSucceeds,          & ! intent(out):   flag to indicate if KINSOL successfully solved the problem in current data step
                      stateVec,                & ! intent(out):   model state vector
                      err,message)               ! intent(out):   error control

  !======= Inclusions ===========

  USE fkinsol_mod                                 ! Fortran interface to KINSOL
  USE fsundials_context_mod                       ! Fortran interface to SUNContext
  USE fnvector_serial_mod                         ! Fortran interface to serial N_Vector
  USE fsundials_nvector_mod                       ! Fortran interface to generic N_Vector
  USE fsunmatrix_dense_mod                        ! Fortran interface to dense SUNMatrix
  USE fsunmatrix_band_mod                         ! Fortran interface to banded SUNMatrix
  USE fsundials_matrix_mod                        ! Fortran interface to generic SUNMatrix
  USE fsunlinsol_dense_mod                        ! Fortran interface to dense SUNLinearSolver
  USE fsunlinsol_band_mod                         ! Fortran interface to banded SUNLinearSolver
  USE fsundials_linearsolver_mod                  ! Fortran interface to generic SUNLinearSolver
  USE allocspace_module,only:allocLocal           ! allocate local data structures
  USE getVectorz_module, only:checkFeas           ! check feasibility of state vector
  USE eval8summa_module,only:eval8summa4kinsol    ! DAE/ODE functions
  USE eval8summa_module,only:eval8summa           ! residual of DAE
  USE computJacob_module,only:computJacob4kinsol  ! system Jacobian
   
  !======= Declarations =========
  implicit none

  ! --------------------------------------------------------------------------------------------------------------------------------
  ! calling variables
  ! --------------------------------------------------------------------------------------------------------------------------------
  ! input: model control
  real(rkind),intent(in)          :: dt_cur                 ! current stepsize
  real(rkind),intent(in)          :: dt                     ! data time step
  real(rkind),intent(inout)       :: fScale(:)              ! characteristic scale of the function evaluations (mixed units)
  real(rkind),intent(inout)       :: xScale(:)              ! characteristic scale of the state vector (mixed units)
  integer(i4b),intent(in)         :: nSnow                  ! number of snow layers
  integer(i4b),intent(in)         :: nSoil                  ! number of soil layers
  integer(i4b),intent(in)         :: nLayers                ! total number of layers
  integer(i4b),intent(in)         :: nStat                  ! total number of state variables
  integer(i4b),intent(in)         :: ixMatrix               ! form of matrix (dense or banded)
  logical(lgt),intent(in)         :: firstSubStep           ! flag to indicate if we are processing the first sub-step
  logical(lgt),intent(in)         :: computeVegFlux         ! flag to indicate if computing fluxes over vegetation
  logical(lgt),intent(in)         :: scalarSolution         ! flag to denote if implementing the scalar solution
  ! input: state vectors
  real(rkind),intent(in)          :: stateVecInit(:)        ! model state vector
  real(qp),intent(in)             :: sMul(:)                ! state vector multiplier (used in the residual calculations)
  real(rkind), intent(inout)      :: dMat(:)                ! diagonal of the Jacobian matrix (excludes fluxes)
  ! input: data structures
  type(model_options),intent(in)  :: model_decisions(:)       ! model decisions
  type(zLookup),      intent(in)  :: lookup_data            ! lookup tables
  type(var_i),        intent(in)  :: type_data              ! type of vegetation and soil
  type(var_d),        intent(in)  :: attr_data              ! spatial attributes
  type(var_dlength),  intent(in)  :: mpar_data              ! model parameters
  type(var_d),        intent(in)  :: forc_data              ! model forcing data
  type(var_dlength),  intent(in)  :: bvar_data              ! model variables for the local basin
  type(var_dlength),  intent(in)  :: prog_data              ! prognostic variables for a local HRU
  ! input-output: data structures
  type(var_ilength),intent(inout) :: indx_data              ! indices defining model states and layers
  type(var_dlength),intent(inout) :: diag_data              ! diagnostic variables for a local HRU
  type(var_dlength),intent(inout) :: flux_data              ! model fluxes for a local HRU
  type(var_dlength),intent(inout) :: deriv_data             ! derivatives in model fluxes w.r.t. relevant state variables
  ! output: state vectors
  integer(i4b),intent(inout)      :: ixSaturation           ! index of the lowest saturated layer
  real(rkind),intent(inout)       :: stateVec(:)            ! model state vector (y)
  logical(lgt),intent(out)        :: kinsolSucceeds         ! flag to indicate if KINSOL is successful
  ! output: error control
  integer(i4b),intent(out)        :: err                    ! error code
  character(*),intent(out)        :: message                ! error message
 
 ! --------------------------------------------------------------------------------------------------------------------------------
  ! local variables
  ! --------------------------------------------------------------------------------------------------------------------------------
  type(N_Vector),           pointer :: sunvec_y             ! sundials solution vector
  type(SUNMatrix),          pointer :: sunmat_A             ! sundials matrix
  type(N_Vector),           pointer :: sunvec_fscale        ! vector containing diagonal elements of function scaling matrix
  type(N_Vector),           pointer :: sunvec_xscale        ! vector containing diagonal elements of state scaling matrix
  type(SUNLinearSolver),    pointer :: sunlinsol_LS         ! sundials linear solver
  type(c_ptr)                       :: kinsol_mem           ! KINSOL memory
  type(c_ptr)                       :: sunctx               ! SUNDIALS simulation context
  type(data4kinsol),        target  :: eqns_data            ! KINSOL type
  integer(i4b)                      :: retval, retvalr      ! return value
  logical(lgt)                      :: feasible             ! feasibility flag
  integer(c_long)                   :: mu, lu               ! in banded matrix mode in SUNDIALS type
  integer(c_long)                   :: nState               ! total number of state variables in SUNDIALS type
  real(rkind)                       :: rVec(nStat)          ! residual vector
  integer(i4b)                      :: iVar, i              ! indices
  character(LEN=256)                :: cmessage             ! error message of downwind routine
  logical(lgt)                      :: use_fdJac                    ! flag to use finite difference Jacobian, default false
  logical(lgt),parameter            :: offErrWarnMessage = .true.   ! flag to turn IDA warnings off, default true
 ! -----------------------------------------------------------------------------------------------------

  ! initialize error control
  err=0; message="summaSolve4kinsol/"

  ! choose Jacobian type
  select case(model_decisions(iLookDECISIONS%fDerivMeth)%iDecision) 
    case(numerical);  use_fdJac =.true.
    case(analytical); use_fdJac =.false.
    case default; err=20; message=trim(message)//'expect choice numericl or analytic to calculate derivatives for Jacobian'; return
  end select

  nState = nStat ! total number of state variables in SUNDIALS type
  kinsolSucceeds = .true.

  ! fill eqns_data which will be required later to call eval8summa
  eqns_data%dt_cur                  = dt_cur
  eqns_data%dt                      = dt
  eqns_data%nSnow                   = nSnow
  eqns_data%nSoil                   = nSoil
  eqns_data%nLayers                 = nLayers
  eqns_data%nState                  = nState
  eqns_data%ixMatrix                = ixMatrix
  eqns_data%firstFluxCall           = .false. ! already called for initial
  eqns_data%firstSplitOper          = .false. ! already called for initial and false inside solver
  eqns_data%firstSubStep            = firstSubStep
  eqns_data%computeVegFlux          = computeVegFlux
  eqns_data%scalarSolution          = scalarSolution
  eqns_data%model_decisions         = model_decisions
  eqns_data%deriv_data              = deriv_data
  eqns_data%lookup_data             = lookup_data
  eqns_data%type_data               = type_data
  eqns_data%attr_data               = attr_data
  eqns_data%mpar_data               = mpar_data
  eqns_data%forc_data               = forc_data
  eqns_data%bvar_data               = bvar_data
  eqns_data%prog_data               = prog_data
  eqns_data%indx_data               = indx_data
  eqns_data%diag_data               = diag_data
  eqns_data%flux_data               = flux_data
  eqns_data%deriv_data              = deriv_data
  eqns_data%ixSaturation            = ixSaturation

  ! allocate space and fill
  allocate( eqns_data%fScale(nState) ); eqns_data%fScale = fScale
  allocate( eqns_data%xScale(nState) ); eqns_data%xScale = xScale
  allocate( eqns_data%sMul(nState) );   eqns_data%sMul   = sMul
  allocate( eqns_data%dMat(nState) );   eqns_data%dMat   = dMat

  ! allocate space for other variables
  if(model_decisions(iLookDECISIONS%groundwatr)%iDecision==qbaseTopmodel)then
    allocate(eqns_data%dBaseflow_dMatric(nSoil,nSoil),stat=err)
  else
    allocate(eqns_data%dBaseflow_dMatric(0,0),stat=err)
  end if
  allocate( eqns_data%fluxVec(nState) )
  allocate( eqns_data%resSink(nState) )
  
  retval = FSUNContext_Create(c_null_ptr, sunctx)

  ! create serial vectors
  sunvec_y => FN_VMake_Serial(nState, stateVec, sunctx)
  if (.not. associated(sunvec_y)) then; err=20; message='summaSolve4kinsol: sunvec = NULL'; return; endif

  ! create the scaling vectors
  sunvec_fscale => FN_VMake_Serial(nState, fscale, sunctx)
  if (.not. associated(sunvec_fscale)) then; err=20; message='summaSolve4kinsol: sunvec = NULL'; return; endif
  sunvec_xscale => FN_VMake_Serial(nState, xscale, sunctx)
  if (.not. associated(sunvec_xscale)) then; err=20; message='summaSolve4kinsol: sunvec = NULL'; return; endif

  ! initialize solution vectors
  call setInitialCondition(nState, stateVecInit, sunvec_y)

  ! create memory
  kinsol_mem = FKINCreate(sunctx)
  if (.not. c_associated(kinsol_mem)) then; err=20; message='summaSolve4kinsol: kinsol_mem = NULL'; return; endif

  ! Attach user data to memory
  retval = FKINSetUserData(kinsol_mem, c_loc(eqns_data))
  if (retval /= 0) then; err=20; message='summaSolve4kinsol: error in FKINSetUserData'; return; endif

  ! Set solver parameters before calling FKINInit
   call setSolverParams(nint(mpar_data%var(iLookPARAM%maxiter)%dat(1)), kinsol_mem, retval)
   if (retval /= 0) then; err=20; message='summaSolve4kinsol: error in setSolverParams'; return; endif

  ! Set the function Kinsol will use to advance the state
  retval = FKINInit(kinsol_mem, c_funloc(eval8summa4kinsol), sunvec_y)
  if (retval /= 0) then; err=20; message='summaSolve4kinsol: error in FKINInit'; return; endif

  ! define the form of the matrix
  select case(ixMatrix)
    case(ixBandMatrix)
      mu = ku; lu = kl;
      ! Create banded SUNMatrix for use in linear solves
      sunmat_A => FSUNBandMatrix(nState, mu, lu, sunctx)
      if (.not. associated(sunmat_A)) then; err=20; message='summaSolve4kinsol: sunmat = NULL'; return; endif

      ! Create banded SUNLinearSolver object
      sunlinsol_LS => FSUNLinSol_Band(sunvec_y, sunmat_A, sunctx)
      if (.not. associated(sunlinsol_LS)) then; err=20; message='summaSolve4kinsol: sunlinsol = NULL'; return; endif

    case(ixFullMatrix)
      ! Create dense SUNMatrix for use in linear solves
      sunmat_A => FSUNDenseMatrix(nState, nState, sunctx)
      if (.not. associated(sunmat_A)) then; err=20; message='summaSolve4kinsol: sunmat = NULL'; return; endif

      ! Create dense SUNLinearSolver object
      sunlinsol_LS => FSUNLinSol_Dense(sunvec_y, sunmat_A, sunctx)
      if (.not. associated(sunlinsol_LS)) then; err=20; message='summaSolve4kinsol: sunlinsol = NULL'; return; endif

      ! check
    case default;  err=20; message='summaSolve4kinsol: error in type of matrix'; return

  end select  ! form of matrix

  ! Attach the matrix and linear solver
  retval = FKINSetLinearSolver(kinsol_mem, sunlinsol_LS, sunmat_A);
  if (retval /= 0) then; err=20; message='summaSolve4kinsol: error in FKINSetLinearSolver'; return; endif

  ! Set the user-supplied Jacobian routine
  if(.not.use_fdJac)then
    retval = FKINSetJacFn(kinsol_mem, c_funloc(computJacob4kinsol))
  if (retval /= 0) then; err=20; message='summaSolve4kinsol: error in FKINSetJacFn'; return; endif
  endif    

  ! Disable error messages and warnings
  if(offErrWarnMessage) then
    retval = FKINSetErrFile(kinsol_mem, c_null_ptr)
  endif

  !****************************** Main Solver **********************************************
  ! Call KINSol to solve problem with choice of solver, linesearch or Picard
  retval = FKINSol(kinsol_mem, sunvec_y, KIN_LINESEARCH, sunvec_xscale, sunvec_fscale)
  !retval = FKINSol(kinsol_mem, sunvec_y, KIN_PICARD, sunvec_xscale, sunvec_fscale)

  if( retvalr < 0 )then
    kinsolSucceeds = .false.
    call getErrMessage(retval,cmessage)
    message=trim(message)//trim(cmessage)
  else
    ! check the feasibility of the solution
    feasible=.true.
    call checkFeas(&
                    ! input
                    stateVec,                                  & ! intent(in):    model state vector (mixed units)
                    eqns_data%mpar_data,                       & ! intent(in):    model parameters
                    eqns_data%prog_data,                       & ! intent(in):    model prognostic variables for a local HRU
                    eqns_data%indx_data,                       & ! intent(in):    indices defining model states and layers
                    ! output: feasibility
                    feasible,                                  & ! intent(inout):   flag to denote the feasibility of the solution
                    ! output: error control
                    err,cmessage)                                 ! intent(out):   error control
  endif 
  !****************************** End of Main Solver ***************************************

  err               = eqns_data%err
  message           = eqns_data%message
  if( .not. feasible)then 
    kinsolSucceeds = .false.
    message=trim(message)//trim(cmessage)//'non-feasible'
  endif

  if(kinsolSucceeds)then
    ! copy to output data
    diag_data     = eqns_data%diag_data
    flux_data     = eqns_data%flux_data
    deriv_data    = eqns_data%deriv_data
    ixSaturation  = eqns_data%ixSaturation
    indx_data%var(iLookINDEX%numberFluxCalc)%dat(1) = eqns_data%indx_data%var(iLookINDEX%numberFluxCalc)%dat(1) !only number of Flux Calculations changes
  endif

  ! free memory
  deallocate( eqns_data%model_decisions)
  deallocate( eqns_data%fScale )
  deallocate( eqns_data%xScale )
  deallocate( eqns_data%sMul )
  deallocate( eqns_data%dMat )
  deallocate( eqns_data%dBaseflow_dMatric )
  deallocate( eqns_data%fluxVec )
  deallocate( eqns_data%resSink )

  call FKINFree(kinsol_mem)
  retval = FSUNLinSolFree(sunlinsol_LS)
  if(retval /= 0)then; err=20; message='summaSolve4kinsol: unable to free the linear solver'; return; endif
  call FSUNMatDestroy(sunmat_A)
  call FN_VDestroy(sunvec_y)
  call FN_VDestroy(sunvec_xscale)
  call FN_VDestroy(sunvec_fscale)
  retval = FSUNContext_Free(sunctx)
  if(retval /= 0)then; err=20; message='summaSolve4kinsol: unable to free the SUNDIALS context'; return; endif

end subroutine summaSolve4kinsol

! ----------------------------------------------------------------
! SetInitialCondition: routine to initialize u vector.
! ----------------------------------------------------------------
subroutine setInitialCondition(neq, y, sunvec_u)

  !======= Inclusions ===========
  USE, intrinsic :: iso_c_binding
  USE fsundials_nvector_mod
  USE fnvector_serial_mod

  !======= Declarations =========
  implicit none

  ! calling variables
  type(N_Vector)  :: sunvec_u  ! solution N_Vector
  integer(c_long) :: neq
  real(rkind)     :: y(neq)

  ! pointers to data in SUNDIALS vectors
  real(c_double), pointer :: uu(:)

  ! get data arrays from SUNDIALS vectors
  uu(1:neq) => FN_VGetArrayPointer(sunvec_u)

  uu = y

end subroutine setInitialCondition

! ----------------------------------------------------------------
! setSolverParams: private routine to set parameters in KINSOL solver
! ----------------------------------------------------------------
subroutine setSolverParams(nonlin_iter,kinsol_mem,retval)

  !======= Inclusions ===========
  USE, intrinsic :: iso_c_binding
  USE fkinsol_mod   ! Fortran interface to KINSOL

  !======= Declarations =========
  implicit none

  ! calling variables
  integer,intent(in)          :: nonlin_iter        ! maximum number of nonlinear iterations, default = 200, set in parameters
  type(c_ptr),intent(inout)   :: kinsol_mem         ! KINSOL memory
  integer(i4b),intent(out)    :: retval             ! return value

  !======= Internals ============
  integer(c_long)             :: nonlin_itr         ! maximum number of nonlinear iterations in SUNDIALS type
  integer(c_long),parameter   :: mset = 1           ! maximum number of times the solver is called without Jacobian update, pass 0 to give default of 10 times
  integer(c_long),parameter   :: msubset = 1        ! maximum number of nonlinear iterations between checks by the residual monitoring algorithm, default=5
  integer(c_long),parameter   :: maa = 0            ! maximum number of prior residuals to use acceleration, default = 0
  integer(c_long),parameter   :: beta_fail = 10     ! maximum number of beta condition failures, default = 10
  real(qp),parameter          :: fnormtol = 0.0     ! stopping tolerance on the scaled maximum norm of the system function, pass 0 to give default of unit_roundoff**(1/3)
  real(qp),parameter          :: scsteptol = 0.0    ! stopping tolerance on the minimum scaled step length, pass 0 to give default of unit_roundoff**(2/3)
      
  ! Set maximum number of times the linear solver is called without a Jacobian update
  retval = FKINSetMaxSetupCalls(kinsol_mem, mset)
  if (retval /= 0) return

  ! Every msubset iterations, test if a Jacobian evaluation is necessary
  retval = FKINSetMaxSubSetupCalls(kinsol_mem, msubset)
  if (retval /= 0) return

  ! Set maximum number of iterations   
  nonlin_itr = nonlin_iter ! maximum number of nonlinear iterations in SUNDIALS type
  retval = FKINSetNumMaxIters(kinsol_mem, nonlin_itr)
  if (retval /= 0) return

  ! Set maximum number of prior residuals to use for Anderson acceleration 
  ! ONLY in conjunction with Picard or fixed-point iteration
  retval = FKINSetMAA(kinsol_mem, maa);
  if (retval /= 0) return

  ! Set maximum number of beta condition failures in the linesearch
  retval = FKINSetMaxBetaFails(kinsol_mem, beta_fail)
  if (retval /= 0) return

  ! Set tolerances for stopping criteria: scaled maximum norm of the system function
  retval = FKINSetFuncNormTol(kinsol_mem, fnormtol)
  if (retval /= 0) return

  ! Set stopping tolerance on the scaled maximum norm of the system function
  retval = FKINSetScaledStepTol(kinsol_mem, scsteptol)
  if (retval /= 0) return

end subroutine setSolverParams

! ----------------------------------------------------------------
! getErrMessage: private routine to get error message for KINSOL solver
! ----------------------------------------------------------------
subroutine getErrMessage(retval,message)

  !======= Declarations =========
  implicit none

  ! calling variables
  integer(i4b),intent(in)    :: retval                 ! return value from KINSOL
  character(*),intent(out)   :: message                ! error message

  ! get message
  if( retval==-1 ) message = 'KIN_MEM_NULL'            ! The kin_mem argument was NULL.
  if( retval==-2 ) message = 'KIN_ILL_INPUT'           ! One of the function inputs is illegal.
  if( retval==-3 ) message = 'KIN_NO_MALLOC'           ! The KINSOL memory was not allocated by a call to KINMalloc.
  if( retval==-4 ) message = 'KIN_MEM_FAIL'            ! A memory allocation failed.
  if( retval==-5 ) message = 'KIN_LINESEARCH_NONCONV ' ! The linesearch algorithm was unable to find an iterate sufficiently distinct from the current iterate.
  if( retval==-6 ) message = 'KIN_MAXITER_REACHED'     ! The maximum number of nonlinear iterations has been reached.
  if( retval==-7 ) message = 'KIN_MXNEWT_5X_EXCEEDED'  ! Five consecutive steps have been taken that satisfy a scaled step length test.
  if( retval==-8 ) message = 'KIN_LINESEARCH_BCFAIL'   ! The linesearch algorithm was unable to satisfy the 𝛽-condition for nbcfails iterations.
  if( retval==-9 ) message = 'KIN_LINSOLV_NO_RECOVERY' ! The user-supplied routine preconditioner slve function failed recoverably, but the preconditioner is already current.
  if( retval==-10) message = 'KIN_LINIT_FAIL'          ! The linear solver’s initialization function failed.
  if( retval==-11) message = 'KIN_LSETUP_FAIL'         ! The linear solver’s setup function failed in an unrecoverable manner.
  if( retval==-12) message = 'KIN_LSOLVE_FAIL'         ! The linear solver’s solve function failed in an unrecoverable manner.
  if( retval==-13) message = 'KIN_SYSFUNC_FAIL'        ! The system function failed in an unrecoverable manner.
  if( retval==-14) message = 'KIN_FIRST_SYSFUNC_ERR'   ! The system function failed with a recoverable error at the first call.
  if( retval==-15) message = 'KIN_REPTD_SYSFUNC_ERR'   ! The system function had repeated recoverable errors.
 
end subroutine getErrMessage


end module summaSolve4kinsol_module
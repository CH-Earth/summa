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

module summaSolve4ida_module


!======= Inclusions ===========
USE, intrinsic :: iso_c_binding
USE nrtype
USE type4ida

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
 private::find_rootdir
 public::layerDisCont4ida
 private::getErrMessage
 public::summaSolve4ida

contains


!-------------------
! * public subroutine summaSolve4ida: solve F(y,y') = 0 by IDA (y is the state vector)
! ------------------
subroutine summaSolve4ida(                         &
                      dt_cur,                  & ! intent(in):    current stepsize
                      dt,                      & ! intent(in):    data time step
                      atol,                    & ! intent(in):    absolute tolerance
                      rtol,                    & ! intent(in):    relative tolerance
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
                      indx_data,               & ! intent(inout): index data
                      diag_data,               & ! intent(inout): model diagnostic variables for a local HRU
                      flux_data,               & ! intent(inout): model fluxes for a local HRU
                      flux_sum,                & ! intent(inout): sum of fluxes model fluxes for a local HRU over a dt_cur
                      deriv_data,              & ! intent(inout): derivatives in model fluxes w.r.t. relevant state variables
                      mLayerCmpress_sum,       & ! intent(inout): sum of compression of the soil matrix
                      ! output
                      ixSaturation,            & ! intent(inout)  index of the lowest saturated layer (NOTE: only computed on the first iteration)
                      idaSucceeds,             & ! intent(out):   flag to indicate if IDA successfully solved the problem in current data step
                      tooMuchMelt,             & ! intent(inout): lag to denote that there was too much melt
                      stateVec,                & ! intent(out):   model state vector
                      stateVecPrime,           & ! intent(out):   derivative of model state vector
                      err,message)               ! intent(out):   error control

  !======= Inclusions ===========
  USE fida_mod                                    ! Fortran interface to IDA
  USE fsundials_context_mod                       ! Fortran interface to SUNContext
  USE fnvector_serial_mod                         ! Fortran interface to serial N_Vector
  USE fsundials_nvector_mod                       ! Fortran interface to generic N_Vector
  USE fsunmatrix_dense_mod                        ! Fortran interface to dense SUNMatrix
  USE fsunmatrix_band_mod                         ! Fortran interface to banded SUNMatrix
  USE fsundials_matrix_mod                        ! Fortran interface to generic SUNMatrix
  USE fsunlinsol_dense_mod                        ! Fortran interface to dense SUNLinearSolver
  USE fsunlinsol_band_mod                         ! Fortran interface to banded SUNLinearSolver
  USE fsundials_linearsolver_mod                  ! Fortran interface to generic SUNLinearSolver
  USE fsunnonlinsol_newton_mod                    ! Fortran interface to Newton SUNNonlinearSolver
  USE fsundials_nonlinearsolver_mod               ! Fortran interface to generic SUNNonlinearSolver
  USE allocspace_module,only:allocLocal           ! allocate local data structures
  USE getVectorz_module, only:checkFeas           ! check feasibility of state vector
  USE eval8summaWithPrime_module,only:eval8summa4ida      ! DAE/ODE functions
  USE eval8summaWithPrime_module,only:eval8summaWithPrime ! residual of DAE
  USE computJacobWithPrime_module,only:computJacob4ida    ! system Jacobian
  USE tol4ida_module,only:computWeight4ida        ! weight required for tolerances

  !======= Declarations =========
  implicit none

  ! --------------------------------------------------------------------------------------------------------------------------------
  ! calling variables
  ! --------------------------------------------------------------------------------------------------------------------------------
  ! input: model control
  real(rkind),intent(in)          :: dt_cur                 ! current stepsize
  real(qp),intent(in)             :: dt                     ! data time step
  real(qp),intent(inout)          :: atol(:)                ! vector of absolute tolerances
  real(qp),intent(inout)          :: rtol(:)                ! vector of relative tolerances
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
  type(var_dlength),intent(inout) :: flux_sum               ! sum of fluxes model fluxes for a local HRU over a dt_cur
  type(var_dlength),intent(inout) :: deriv_data             ! derivatives in model fluxes w.r.t. relevant state variables
  real(rkind),intent(inout)       :: mLayerCmpress_sum(:)   ! sum of soil compress
  ! output: state vectors
  integer(i4b),intent(inout)      :: ixSaturation           ! index of the lowest saturated layer
  real(rkind),intent(inout)       :: stateVec(:)            ! model state vector (y)
  real(rkind),intent(inout)       :: stateVecPrime(:)       ! model state vector (y')
  logical(lgt),intent(out)        :: idaSucceeds            ! flag to indicate if IDA is successful
  logical(lgt),intent(inout)      :: tooMuchMelt            ! flag to denote that there was too much melt
  ! output: error control
  integer(i4b),intent(out)        :: err                    ! error code
  character(*),intent(out)        :: message                ! error message

  ! --------------------------------------------------------------------------------------------------------------------------------
  ! local variables
  ! --------------------------------------------------------------------------------------------------------------------------------
  type(N_Vector),           pointer :: sunvec_y             ! sundials solution vector
  type(N_Vector),           pointer :: sunvec_yp            ! sundials derivative vector
  type(N_Vector),           pointer :: sunvec_av            ! sundials tolerance vector
  type(SUNMatrix),          pointer :: sunmat_A             ! sundials matrix
  type(SUNLinearSolver),    pointer :: sunlinsol_LS         ! sundials linear solver
  type(SUNNonLinearSolver), pointer :: sunnonlin_NLS        ! sundials nonlinear solver
  type(c_ptr)                       :: ida_mem              ! IDA memory
  type(c_ptr)                       :: sunctx               ! SUNDIALS simulation context
  type(data4ida),           target  :: eqns_data            ! IDA type
  integer(i4b)                      :: retval, retvalr      ! return value
  logical(lgt)                      :: feasible             ! feasibility flag
  real(qp)                          :: t0                   ! starting time
  real(qp)                          :: dt_last(1)           ! last time step
  real(qp)                          :: dt_diff              ! difference from previous timestep
  integer(c_long)                   :: mu, lu               ! in banded matrix mode in SUNDIALS type
  integer(c_long)                   :: nState               ! total number of state variables in SUNDIALS type
  real(rkind)                       :: rVec(nStat)          ! residual vector
  integer(i4b)                      :: iVar, i              ! indices
  integer(i4b)                      :: nRoot                ! total number of roots (events) to find
  real(qp)                          :: tret(1)              ! time in data window
  real(qp)                          :: tretPrev             ! previous time in data window
  integer(i4b),allocatable          :: rootsfound(:)        ! crossing direction of discontinuities
  integer(i4b),allocatable          :: rootdir(:)           ! forced crossing direction of discontinuities
  logical(lgt)                      :: tinystep             ! if step goes below small size
  type(var_dlength)                 :: flux_prev            ! previous model fluxes for a local HRU
  character(LEN=256)                :: cmessage             ! error message of downwind routine
  real(rkind),allocatable           :: mLayerMatricHeadPrimePrev(:) ! previous derivative value for total water matric potential (m s-1)
  real(rkind),allocatable           :: dCompress_dPsiPrev(:)        ! previous derivative value soil compression
  logical(lgt)                      :: use_fdJac                    ! flag to use finite difference Jacobian, default false
  logical(lgt),parameter            :: offErrWarnMessage = .true.   ! flag to turn IDA warnings off, default true
  logical(lgt),parameter            :: detect_events = .true.       ! flag to do event detection and restarting, default true

  ! -----------------------------------------------------------------------------------------------------

  ! initialize error control
  err=0; message="summaSolve4ida/"

  ! choose Jacobian type
  select case(model_decisions(iLookDECISIONS%fDerivMeth)%iDecision) 
    case(numerical);  use_fdJac =.true.
    case(analytical); use_fdJac =.false.
    case default; err=20; message=trim(message)//'expect choice numericl or analytic to calculate derivatives for Jacobian'; return
  end select

  nState = nStat ! total number of state variables in SUNDIALS type
  idaSucceeds = .true.

  ! fill eqns_data which will be required later to call eval8summaWithPrime
  eqns_data%dt                      = dt
  eqns_data%nSnow                   = nSnow
  eqns_data%nSoil                   = nSoil
  eqns_data%nLayers                 = nLayers
  eqns_data%nState                  = nState
  eqns_data%ixMatrix                = ixMatrix
  eqns_data%firstSubStep            = firstSubStep
  eqns_data%computeVegFlux          = computeVegFlux
  eqns_data%scalarSolution          = scalarSolution
  eqns_data%model_decisions         = model_decisions
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
  allocate( eqns_data%atol(nState) ); eqns_data%atol = atol
  allocate( eqns_data%rtol(nState) ); eqns_data%rtol = rtol
  allocate( eqns_data%sMul(nState) ); eqns_data%sMul = sMul
  allocate( eqns_data%dMat(nState) ); eqns_data%dMat = dMat

  ! allocate space for the to save previous fluxes
  call allocLocal(flux_meta(:),flux_prev,nSnow,nSoil,err,message)
  if(err/=0)then; err=20; message=trim(message)//trim(message); return; endif
  flux_prev                         = eqns_data%flux_data

  ! allocate space for other variables
  if(model_decisions(iLookDECISIONS%groundwatr)%iDecision==qbaseTopmodel)then
    allocate(eqns_data%dBaseflow_dMatric(nSoil,nSoil),stat=err)
  else
    allocate(eqns_data%dBaseflow_dMatric(0,0),stat=err)
  end if
  allocate( eqns_data%mLayerMatricHeadLiqTrial(nSoil) )
  allocate( eqns_data%mLayerMatricHeadTrial(nSoil) )
  allocate( eqns_data%mLayerMatricHeadPrev(nSoil) )
  allocate( eqns_data%mLayerVolFracWatTrial(nLayers) )
  allocate( eqns_data%mLayerVolFracWatPrev(nLayers) )
  allocate( eqns_data%mLayerTempTrial(nLayers) )
  allocate( eqns_data%mLayerTempPrev(nLayers) )
  allocate( eqns_data%mLayerVolFracIceTrial(nLayers) )
  allocate( eqns_data%mLayerVolFracIcePrev(nLayers) )
  allocate( eqns_data%mLayerVolFracLiqTrial(nLayers) )
  allocate( eqns_data%mLayerVolFracLiqPrev(nLayers) )
  allocate( eqns_data%mLayerEnthalpyTrial(nLayers) )
  allocate( eqns_data%mLayerEnthalpyPrev(nLayers) )
  allocate( eqns_data%mLayerTempPrime(nLayers) )
  allocate( eqns_data%mLayerMatricHeadPrime(nSoil) )
  allocate( eqns_data%mLayerMatricHeadLiqPrime(nSoil) )
  allocate( eqns_data%mLayerVolFracWatPrime(nLayers) )
  allocate( mLayerMatricHeadPrimePrev(nSoil) )
  allocate( dCompress_dPsiPrev(nSoil) )
  allocate( eqns_data%fluxVec(nState) )
  allocate( eqns_data%resSink(nState) )

  ! need the following values for the first substep
  eqns_data%scalarCanopyTempPrev     = prog_data%var(iLookPROG%scalarCanopyTemp)%dat(1)
  eqns_data%scalarCanopyIcePrev      = prog_data%var(iLookPROG%scalarCanopyIce)%dat(1)
  eqns_data%scalarCanopyLiqPrev      = prog_data%var(iLookPROG%scalarCanopyLiq)%dat(1)
  eqns_data%scalarCanopyEnthalpyPrev = diag_data%var(iLookDIAG%scalarCanopyEnthalpy)%dat(1)
  eqns_data%mLayerTempPrev(:)        = prog_data%var(iLookPROG%mLayerTemp)%dat(:)
  eqns_data%mLayerMatricHeadPrev(:)  = prog_data%var(iLookPROG%mLayerMatricHead)%dat(:)
  eqns_data%mLayerVolFracWatPrev(:)  = prog_data%var(iLookPROG%mLayerVolFracWat)%dat(:)
  eqns_data%mLayerVolFracIcePrev(:)  = prog_data%var(iLookPROG%mLayerVolFracIce)%dat(:)
  eqns_data%mLayerVolFracLiqPrev(:)  = prog_data%var(iLookPROG%mLayerVolFracLiq)%dat(:)
  eqns_data%mLayerEnthalpyPrev(:)    = diag_data%var(iLookDIAG%mLayerEnthalpy)%dat(:)
  eqns_data%scalarAquiferStoragePrev = prog_data%var(iLookPROG%scalarAquiferStorage)%dat(1)
  mLayerMatricHeadPrimePrev(:)       = 0._rkind
  dCompress_dPsiPrev(:)              = 0._rkind

  retval = FSUNContext_Create(c_null_ptr, sunctx)

  ! create serial vectors
  sunvec_y => FN_VMake_Serial(nState, stateVec, sunctx)
  if (.not. associated(sunvec_y)) then; err=20; message='summaSolve4ida: sunvec = NULL'; return; endif
  sunvec_yp => FN_VMake_Serial(nState, stateVecPrime, sunctx)
  if (.not. associated(sunvec_yp)) then; err=20; message='summaSolve4ida: sunvec = NULL'; return; endif

  ! initialize solution vectors
  call setInitialCondition(nState, stateVecInit, sunvec_y, sunvec_yp)

  ! create memory
  ida_mem = FIDACreate(sunctx)
  if (.not. c_associated(ida_mem)) then; err=20; message='summaSolve4ida: ida_mem = NULL'; return; endif

  ! Attach user data to memory
  retval = FIDASetUserData(ida_mem, c_loc(eqns_data))
  if (retval /= 0) then; err=20; message='summaSolve4ida: error in FIDASetUserData'; return; endif

  ! Set the function IDA will use to advance the state
  t0 = 0._rkind
  retval = FIDAInit(ida_mem, c_funloc(eval8summa4ida), t0, sunvec_y, sunvec_yp)
  if (retval /= 0) then; err=20; message='summaSolve4ida: error in FIDAInit'; return; endif

  ! set tolerances
  retval = FIDAWFtolerances(ida_mem, c_funloc(computWeight4ida))
  if (retval /= 0) then; err=20; message='summaSolve4ida: error in FIDAWFtolerances'; return; endif

  ! initialize rootfinding problem and allocate space, counting roots
  if(detect_events)then
    nRoot = 0
    if(eqns_data%indx_data%var(iLookINDEX%ixVegNrg)%dat(1)/=integerMissing) nRoot = nRoot+1
    if(nSnow>0)then
      do i = 1,nSnow
        if(eqns_data%indx_data%var(iLookINDEX%ixSnowOnlyNrg)%dat(i)/=integerMissing) nRoot = nRoot+1
      enddo
    endif
    if(nSoil>0)then
      do i = 1,nSoil
        if(eqns_data%indx_data%var(iLookINDEX%ixSoilOnlyHyd)%dat(i)/=integerMissing) nRoot = nRoot+1
        if(eqns_data%indx_data%var(iLookINDEX%ixSoilOnlyNrg)%dat(i)/=integerMissing) nRoot = nRoot+1
      enddo
    endif
    allocate( rootsfound(nRoot) )
    allocate( rootdir(nRoot) )
    rootdir = 0
    retval = FIDARootInit(ida_mem, nRoot, c_funloc(layerDisCont4ida))
    if (retval /= 0) then; err=20; message='summaSolve4ida: error in FIDARootInit'; return; endif
  else ! will not use, allocate at something
    nRoot = 1
    allocate( rootsfound(nRoot) )
    allocate( rootdir(nRoot) )
  endif

  ! define the form of the matrix
  select case(ixMatrix)
    case(ixBandMatrix)
      mu = ku; lu = kl;
      ! Create banded SUNMatrix for use in linear solves
      sunmat_A => FSUNBandMatrix(nState, mu, lu, sunctx)
      if (.not. associated(sunmat_A)) then; err=20; message='summaSolve4ida: sunmat = NULL'; return; endif

      ! Create banded SUNLinearSolver object
      sunlinsol_LS => FSUNLinSol_Band(sunvec_y, sunmat_A, sunctx)
      if (.not. associated(sunlinsol_LS)) then; err=20; message='summaSolve4ida: sunlinsol = NULL'; return; endif

    case(ixFullMatrix)
      ! Create dense SUNMatrix for use in linear solves
      sunmat_A => FSUNDenseMatrix(nState, nState, sunctx)
      if (.not. associated(sunmat_A)) then; err=20; message='summaSolve4ida: sunmat = NULL'; return; endif

      ! Create dense SUNLinearSolver object
      sunlinsol_LS => FSUNLinSol_Dense(sunvec_y, sunmat_A, sunctx)
      if (.not. associated(sunlinsol_LS)) then; err=20; message='summaSolve4ida: sunlinsol = NULL'; return; endif

      ! check
    case default;  err=20; message='summaSolve4ida: error in type of matrix'; return

  end select  ! form of matrix

  ! Attach the matrix and linear solver
  ! For the nonlinear solver, IDA uses a Newton SUNNonlinearSolver-- it is not necessary to create and attach it
  retval = FIDASetLinearSolver(ida_mem, sunlinsol_LS, sunmat_A);
  if (retval /= 0) then; err=20; message='summaSolve4ida: error in FIDASetLinearSolver'; return; endif

  ! Set the user-supplied Jacobian routine
  if(.not.use_fdJac)then
    retval = FIDASetJacFn(ida_mem, c_funloc(computJacob4ida))
    if (retval /= 0) then; err=20; message='summaSolve4ida: error in FIDASetJacFn'; return; endif
  endif

  ! Enforce the solver to stop at end of the time step
  retval = FIDASetStopTime(ida_mem, dt_cur)
  if (retval /= 0) then; err=20; message='summaSolve4ida: error in FIDASetStopTime'; return; endif

  ! Set solver parameters at end of setup
  call setSolverParams(dt_cur, nint(mpar_data%var(iLookPARAM%maxiter)%dat(1)), ida_mem, retval)
  if (retval /= 0) then; err=20; message='summaSolve4ida: error in setSolverParams'; return; endif

  ! Disable error messages and warnings
  if(offErrWarnMessage) then
    retval = FIDASetErrFile(ida_mem, c_null_ptr)
    retval = FIDASetNoInactiveRootWarn(ida_mem)
  endif

  !*********************** Main Solver * loop on one_step mode *****************************
  tinystep = .false.
  tret(1) = t0           ! initial time
  tretPrev = tret(1)
  do while(tret(1) < dt_cur)

    ! call this at beginning of step to reduce root bouncing (only looking in one direction)
    if(detect_events .and. .not.tinystep)then
      call find_rootdir(eqns_data, rootdir)
      retval = FIDASetRootDirection(ida_mem, rootdir)
      if (retval /= 0) then; err=20; message='summaSolve4ida: error in FIDASetRootDirection'; return; endif
    endif

    eqns_data%firstFluxCall = .false. ! already called for initial
    eqns_data%firstSplitOper = .true. ! always true at start of dt_cur since no splitting

    ! call IDASolve, advance solver just one internal step
    retvalr = FIDASolve(ida_mem, dt_cur, tret, sunvec_y, sunvec_yp, IDA_ONE_STEP)
    if( retvalr < 0 )then
      idaSucceeds = .false.
      call getErrMessage(retvalr,cmessage)
      message=trim(message)//trim(cmessage)
      exit
    end if

    tooMuchMelt = .false.
    ! loop through non-missing energy state variables in the snow domain to see if need to merge
    do concurrent (i=1:nSnow,indx_data%var(iLookINDEX%ixSnowOnlyNrg)%dat(i)/=integerMissing)
      if (stateVec(indx_data%var(iLookINDEX%ixSnowOnlyNrg)%dat(i)) > Tfreeze) tooMuchMelt = .true. !need to merge
    enddo
    if(tooMuchMelt)exit

    ! get the last stepsize and difference from previous end time, not necessarily the same
    retval = FIDAGetLastStep(ida_mem, dt_last)
    dt_diff = tret(1) - tretPrev

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

    ! early return for non-feasible solutions, will fail in current IDA formulation
    if(.not.feasible)then
      eqns_data%fluxVec(:) = realMissing
      message=trim(message)//trim(cmessage)//'non-feasible'
      return
    end if

    ! sum of fluxes smoothed over the time step, average from instantaneous values
    do iVar=1,size(flux_meta)
      flux_sum%var(iVar)%dat(:) = flux_sum%var(iVar)%dat(:) + ( eqns_data%flux_data%var(iVar)%dat(:) &
                                                                + flux_prev%var(iVar)%dat(:) ) *dt_diff/2._rkind
    end do
    mLayerCmpress_sum(:) = mLayerCmpress_sum(:) + ( eqns_data%deriv_data%var(iLookDERIV%dCompress_dPsi)%dat(:) * eqns_data%mLayerMatricHeadPrime(:) &
                                                    + dCompress_dPsiPrev(:)  * mLayerMatricHeadPrimePrev(:) ) * dt_diff/2._rkind

    ! save required quantities for next step
    eqns_data%scalarCanopyTempPrev     = eqns_data%scalarCanopyTempTrial
    eqns_data%scalarCanopyIcePrev      = eqns_data%scalarCanopyIceTrial
    eqns_data%scalarCanopyLiqPrev      = eqns_data%scalarCanopyLiqTrial
    eqns_data%scalarCanopyEnthalpyPrev = eqns_data%scalarCanopyEnthalpyTrial
    eqns_data%mLayerTempPrev(:)        = eqns_data%mLayerTempTrial(:)
    eqns_data%mLayerMatricHeadPrev(:)  = eqns_data%mLayerMatricHeadTrial(:)
    eqns_data%mLayerVolFracWatPrev(:)  = eqns_data%mLayerVolFracWatTrial(:)
    eqns_data%mLayerVolFracIcePrev(:)  = eqns_data%mLayerVolFracIceTrial(:)
    eqns_data%mLayerVolFracLiqPrev(:)  = eqns_data%mLayerVolFracLiqTrial(:)
    eqns_data%mLayerEnthalpyPrev(:)    = eqns_data%mLayerEnthalpyTrial(:)
    eqns_data%scalarAquiferStoragePrev = eqns_data%scalarAquiferStorageTrial
    mLayerMatricHeadPrimePrev(:)       = eqns_data%mLayerMatricHeadPrime(:)
    dCompress_dPsiPrev(:)              = eqns_data%deriv_data%var(iLookDERIV%dCompress_dPsi)%dat(:)
    tretPrev                           = tret(1)
    flux_prev                          = eqns_data%flux_data

    ! Restart for where vegetation and layers cross freezing point
    if(detect_events)then
      if (retvalr .eq. IDA_ROOT_RETURN) then !IDASolve succeeded and found one or more roots at tret(1)
        ! rootsfound[i]= +1 indicates that gi is increasing, -1 g[i] decreasing, 0 no root
        !retval = FIDAGetRootInfo(ida_mem, rootsfound)
        !if (retval < 0) then; err=20; message='summaSolve4ida: error in FIDAGetRootInfo'; return; endif
        !print '(a,f15.7,2x,17(i2,2x))', "time, rootsfound[] = ", tret(1), rootsfound
        ! Reininitialize solver for running after discontinuity and restart
        retval = FIDAReInit(ida_mem, tret(1), sunvec_y, sunvec_yp)
        if (retval /= 0) then; err=20; message='summaSolve4ida: error in FIDAReInit'; return; endif
        if(dt_last(1) < 0.1_rkind)then ! don't keep calling if step is small (more accurate with this tiny but getting hung up)
          retval = FIDARootInit(ida_mem, 0, c_funloc(layerDisCont4ida))
          tinystep = .true.
        else
          retval = FIDARootInit(ida_mem, nRoot, c_funloc(layerDisCont4ida))
          tinystep = .false.
        endif
        if (retval /= 0) then; err=20; message='summaSolve4ida: error in FIDARootInit'; return; endif
      endif
    endif

  enddo ! while loop on one_step mode until time dt_cur
  !****************************** End of Main Solver ***************************************

  err               = eqns_data%err
  message           = eqns_data%message
  if( .not. feasible)then 
    idaSucceeds = .false.
    message=trim(message)//trim(cmessage)//'non-feasible'
  endif

  if(idaSucceeds)then
    ! copy to output data
    diag_data     = eqns_data%diag_data
    flux_data     = eqns_data%flux_data
    deriv_data    = eqns_data%deriv_data
    ixSaturation  = eqns_data%ixSaturation
    indx_data%var(iLookINDEX%numberFluxCalc)%dat(1) = eqns_data%indx_data%var(iLookINDEX%numberFluxCalc)%dat(1) !only number of Flux Calculations changes
  endif

  ! free memory
  deallocate( eqns_data%model_decisions)
  deallocate( eqns_data%sMul )
  deallocate( eqns_data%dMat )
  deallocate( eqns_data%dBaseflow_dMatric )
  deallocate( eqns_data%mLayerMatricHeadLiqTrial )
  deallocate( eqns_data%mLayerMatricHeadTrial )
  deallocate( eqns_data%mLayerMatricHeadPrev )
  deallocate( eqns_data%mLayerVolFracWatTrial )
  deallocate( eqns_data%mLayerVolFracWatPrev )
  deallocate( eqns_data%mLayerVolFracIceTrial )
  deallocate( eqns_data%mLayerTempPrev )
  deallocate( eqns_data%mLayerTempTrial )
  deallocate( eqns_data%mLayerVolFracIcePrev )
  deallocate( eqns_data%mLayerVolFracLiqPrev )
  deallocate( eqns_data%mLayerEnthalpyTrial )
  deallocate( eqns_data%mLayerEnthalpyPrev )
  deallocate( eqns_data%mLayerTempPrime )
  deallocate( eqns_data%mLayerMatricHeadPrime )
  deallocate( eqns_data%mLayerMatricHeadLiqPrime )
  deallocate( eqns_data%mLayerVolFracWatPrime )

  deallocate( mLayerMatricHeadPrimePrev )
  deallocate( dCompress_dPsiPrev )
  deallocate( eqns_data%fluxVec )
  deallocate( eqns_data%resSink )
  deallocate( rootsfound )
  deallocate( rootdir )

  call FIDAFree(ida_mem)
  retval = FSUNLinSolFree(sunlinsol_LS)
  if(retval /= 0)then; err=20; message='summaSolve4ida: unable to free the linear solver'; return; endif
  call FSUNMatDestroy(sunmat_A)
  call FN_VDestroy(sunvec_y)
  call FN_VDestroy(sunvec_yp)
  retval = FSUNContext_Free(sunctx)
  if(retval /= 0)then; err=20; message='summaSolve4ida: unable to free the SUNDIALS context'; return; endif

end subroutine summaSolve4ida

! ----------------------------------------------------------------
! SetInitialCondition: routine to initialize u and up vectors.
! ----------------------------------------------------------------
subroutine setInitialCondition(neq, y, sunvec_u, sunvec_up)

  !======= Inclusions ===========
  USE, intrinsic :: iso_c_binding
  USE fsundials_nvector_mod
  USE fnvector_serial_mod

  !======= Declarations =========
  implicit none

  ! calling variables
  type(N_Vector)  :: sunvec_u  ! solution N_Vector
  type(N_Vector)  :: sunvec_up ! derivative N_Vector
  integer(c_long) :: neq
  real(rkind)     :: y(neq)

  ! pointers to data in SUNDIALS vectors
  real(c_double), pointer :: uu(:)
  real(c_double), pointer :: up(:)

  ! get data arrays from SUNDIALS vectors
  uu(1:neq) => FN_VGetArrayPointer(sunvec_u)
  up(1:neq) => FN_VGetArrayPointer(sunvec_up)

  uu = y
  up = 0._rkind

end subroutine setInitialCondition

! ----------------------------------------------------------------
! setSolverParams: private routine to set parameters in IDA solver
! ----------------------------------------------------------------
subroutine setSolverParams(dt_cur,nonlin_iter,ida_mem,retval)

  !======= Inclusions ===========
  USE, intrinsic :: iso_c_binding
  USE fida_mod   ! Fortran interface to IDA

  !======= Declarations =========
  implicit none

  ! calling variables
  real(rkind),intent(in)      :: dt_cur             ! current whole time step
  integer,intent(in)          :: nonlin_iter        ! maximum number of nonlinear iterations, default = 4, set in parameters
  type(c_ptr),intent(inout)   :: ida_mem            ! IDA memory
  integer(i4b),intent(out)    :: retval             ! return value

  !======= Internals ============
  integer,parameter           :: max_order = 5      ! maximum BDF order,  default = 5
  real(qp),parameter          :: coef_nonlin = 0.33 ! Coeff. in the nonlinear convergence test, default = 0.33
  integer,parameter           :: acurtest_fail = 50 ! maximum number of error test failures, default = 10
  integer,parameter           :: convtest_fail = 50 ! maximum number of convergence test failures, default = 10
  integer(c_long),parameter   :: max_step = 999999  ! maximum number of steps,  default = 500
  real(qp)                    :: h_max              ! maximum stepsize,  default = infinity
  real(qp),parameter          :: h_init = 0         ! initial stepsize
 
  ! Set the maximum BDF order
  retval = FIDASetMaxOrd(ida_mem, max_order)
  if (retval /= 0) return

  ! Set Coeff. in the nonlinear convergence test
  retval = FIDASetNonlinConvCoef(ida_mem, coef_nonlin)
  if (retval /= 0) return

  ! Set maximun number of nonliear iterations
   retval = FIDASetMaxNonlinIters(ida_mem, nonlin_iter)
  if (retval /= 0) return

  !  Set maximum number of convergence test failures
  retval = FIDASetMaxConvFails(ida_mem, convtest_fail)
  if (retval /= 0) return

  !  Set maximum number of error test failures
  retval = FIDASetMaxErrTestFails(ida_mem, acurtest_fail)
  if (retval /= 0) return

  ! Set maximum number of steps
  retval = FIDASetMaxNumSteps(ida_mem, max_step)
  if (retval /= 0) return

  ! Set maximum stepsize
  h_max = dt_cur
  retval = FIDASetMaxStep(ida_mem, h_max)
  if (retval /= 0) return

  ! Set initial stepsize
  retval = FIDASetInitStep(ida_mem, h_init)
  if (retval /= 0) return

end subroutine setSolverParams

! ----------------------------------------------------------------------------------------
! find_rootdir: private routine to determine which direction to look for the root, by
!  determining if the variable is greater or less than the root. Need to do this to prevent
!  bouncing around solution
! ----------------------------------------------------------------------------------------
subroutine find_rootdir(eqns_data,rootdir)

  !======= Inclusions ===========
  use, intrinsic :: iso_c_binding
  use fsundials_nvector_mod
  use fnvector_serial_mod
  use soil_utils_module,only:crit_soilT  ! compute the critical temperature below which ice exists
  use globalData,only:integerMissing     ! missing integer
  use var_lookup,only:iLookINDEX         ! named variables for structure elements
  use multiconst,only:Tfreeze            ! freezing point of pure water (K)

  !======= Declarations =========
  implicit none

  ! calling variables
  type(data4ida),intent(in)  :: eqns_data  ! equations data
  integer(i4b),intent(inout) :: rootdir(:) ! root function directions to search

  ! local variables
  integer(i4b)               :: i,ind     ! indices
  integer(i4b)               :: nState    ! number of states
  integer(i4b)               :: nSnow     ! number of snow layers
  integer(i4b)               :: nSoil     ! number of soil layers
  real(rkind)                :: xPsi      ! matric head at layer (m)
  real(rkind)                :: TcSoil    ! critical point when soil begins to freeze (K)

  ! get equations data variables
  nState = eqns_data%nState
  nSnow = eqns_data%nSnow
  nSoil = eqns_data%nSoil

  ! initialize
  ind = 0

  ! identify the critical point when vegetation begins to freeze
  if(eqns_data%indx_data%var(iLookINDEX%ixVegNrg)%dat(1)/=integerMissing)then
    ind = ind+1
    rootdir(ind) = 1
    if(eqns_data%scalarCanopyTempPrev > Tfreeze) rootdir(ind) = -1
  endif

  if(nSnow>0)then
    do i = 1,nSnow
      ! identify the critical point when the snow layer begins to freeze
      if(eqns_data%indx_data%var(iLookINDEX%ixSnowOnlyNrg)%dat(i)/=integerMissing)then
        ind = ind+1
        rootdir(ind) = 1
        if(eqns_data%mLayerTempPrev(i) > Tfreeze) rootdir(ind) = -1
      endif
    end do
  endif

  if(nSoil>0)then
    do i = 1,nSoil
      xPsi = eqns_data%mLayerMatricHeadPrev(i)
      ! identify the critical point when soil matrix potential goes below 0 and Tfreeze depends only on temp
      if (eqns_data%indx_data%var(iLookINDEX%ixSoilOnlyHyd)%dat(i)/=integerMissing)then
        ind = ind+1
        rootdir(ind) = 1
        if(xPsi > 0._rkind ) rootdir(ind) = -1
      endif
      ! identify the critical point when the soil layer begins to freeze
      if(eqns_data%indx_data%var(iLookINDEX%ixSoilOnlyNrg)%dat(i)/=integerMissing)then
        ind = ind+1
        TcSoil = crit_soilT(xPsi)
        rootdir(ind) = 1
        if(eqns_data%mLayerTempPrev(i+nSnow) > TcSoil) rootdir(ind) = -1
      endif
    end do
  endif

end subroutine find_rootdir

! ----------------------------------------------------------------------------------------
! layerDisCont4ida: The root function routine to find soil matrix potential = 0,
!  soil temp = critical frozen point, and snow and veg temp = Tfreeze
! ----------------------------------------------------------------------------------------
! Return values:
!    0 = success,
!    1 = recoverable error,
!   -1 = non-recoverable error
! ----------------------------------------------------------------------------------------
integer(c_int) function layerDisCont4ida(t, sunvec_u, sunvec_up, gout, user_data) &
      result(ierr) bind(C,name='layerDisCont4ida')

  !======= Inclusions ===========
  use, intrinsic :: iso_c_binding
  use fsundials_nvector_mod
  use fnvector_serial_mod
  use soil_utils_module,only:crit_soilT  ! compute the critical temperature below which ice exists
  use globalData,only:integerMissing     ! missing integer
  use var_lookup,only:iLookINDEX         ! named variables for structure elements
  use multiconst,only:Tfreeze            ! freezing point of pure water (K)

  !======= Declarations =========
  implicit none

  ! calling variables
  real(c_double), value      :: t         ! current time
  type(N_Vector)             :: sunvec_u  ! solution N_Vector
  type(N_Vector)             :: sunvec_up ! derivative N_Vector
  real(c_double)             :: gout(999) ! root function values, if (nVeg + nSnow + 2*nSoil)>999, problem
  type(c_ptr),    value      :: user_data ! user-defined data

  ! local variables
  integer(i4b)               :: i,ind     ! indices
  integer(i4b)               :: nState    ! number of states
  integer(i4b)               :: nSnow     ! number of snow layers
  integer(i4b)               :: nSoil     ! number of soil layers
  real(rkind)                :: xPsi      ! matric head at layer (m)
  real(rkind)                :: TcSoil    ! critical point when soil begins to freeze (K)

  ! pointers to data in SUNDIALS vectors
  real(c_double), pointer :: uu(:)
  type(data4ida), pointer :: eqns_data      ! equations data

  !======= Internals ============
  ! get equations data from user-defined data
  call c_f_pointer(user_data, eqns_data)
  nState = eqns_data%nState
  nSnow = eqns_data%nSnow
  nSoil = eqns_data%nSoil

  ! get data array from SUNDIALS vector
  uu(1:nState) => FN_VGetArrayPointer(sunvec_u)

  ! initialize
  ind = 0

  ! identify the critical point when vegetation begins to freeze
  if(eqns_data%indx_data%var(iLookINDEX%ixVegNrg)%dat(1)/=integerMissing)then
    ind = ind+1
    gout(ind) = uu(eqns_data%indx_data%var(iLookINDEX%ixVegNrg)%dat(1)) - Tfreeze
  endif

  if(nSnow>0)then
    do i = 1,nSnow
      ! identify the critical point when the snow layer begins to freeze
      if(eqns_data%indx_data%var(iLookINDEX%ixSnowOnlyNrg)%dat(i)/=integerMissing)then
        ind = ind+1
        gout(ind) = uu(eqns_data%indx_data%var(iLookINDEX%ixSnowOnlyNrg)%dat(i)) - Tfreeze
      endif
    end do
  endif

  if(nSoil>0)then
    do i = 1,nSoil
      ! identify the critical point when soil matrix potential goes below 0 and Tfreeze depends only on temp
      if (eqns_data%indx_data%var(iLookINDEX%ixSoilOnlyHyd)%dat(i)/=integerMissing)then
        ind = ind+1
        xPsi = uu(eqns_data%indx_data%var(iLookINDEX%ixSoilOnlyHyd)%dat(i))
        gout(ind) = uu(eqns_data%indx_data%var(iLookINDEX%ixSoilOnlyHyd)%dat(i))
      else
        xPsi = eqns_data%prog_data%var(iLookPROG%mLayerMatricHead)%dat(i)
      endif
      ! identify the critical point when the soil layer begins to freeze
      if(eqns_data%indx_data%var(iLookINDEX%ixSoilOnlyNrg)%dat(i)/=integerMissing)then
        ind = ind+1
        TcSoil = crit_soilT(xPsi)
        gout(ind) = uu(eqns_data%indx_data%var(iLookINDEX%ixSoilOnlyNrg)%dat(i)) - TcSoil
      endif
    end do
  endif

  ! return success
  ierr = 0
  return

end function layerDisCont4ida

! ----------------------------------------------------------------
! getErrMessage: private routine to get error message for IDA solver
! ----------------------------------------------------------------
subroutine getErrMessage(retval,message)

  !======= Declarations =========
  implicit none

  ! calling variables
  integer(i4b),intent(in)    :: retval              ! return value from IDA
  character(*),intent(out)   :: message             ! error message

  ! get message
   if( retval==-1 ) message = 'IDA_TOO_MUCH_WORK'   ! The solver took mxstep internal steps but could not reach tout.
   if( retval==-2 ) message = 'IDA_TOO_MUCH_ACC'    ! The solver could not satisfy the accuracy demanded by the user for some internal step.
   if( retval==-3 ) message = 'IDA_ERR_FAIL'        ! Error test failures occurred too many times during one internal timestep or minimum step size was reached.
   if( retval==-4 ) message = 'IDA_CONV_FAIL'       ! Convergence test failures occurred too many times during one internal time step or minimum step size was reached.
   if( retval==-5 ) message = 'IDA_LINIT_FAIL'      ! The linear solver’s initialization function failed.
   if( retval==-6 ) message = 'IDA_LSETUP_FAIL'     ! The linear solver’s setup function failed in an unrecoverable manner.
   if( retval==-7 ) message = 'IDA_LSOLVE_FAIL'     ! The linear solver’s solve function failed in an unrecoverable manner.
   if( retval==-8 ) message = 'IDA_RES_FAIL'        ! The user-provided residual function failed in an unrecoverable manner.
   if( retval==-9 ) message = 'IDA_REP_RES_FAIL'    ! The user-provided residual function repeatedly returned a recoverable error flag, but the solver was unable to recover.
   if( retval==-10) message = 'IDA_RTFUNC_FAIL'     ! The rootfinding function failed in an unrecoverable manner.
   if( retval==-11) message = 'IDA_CONSTR_FAIL'     ! The inequality constraints were violated and the solver was unable to recover.
   if( retval==-12) message = 'IDA_FIRST_RES_FAIL'  ! The user-provided residual function failed recoverably on the first call.
   if( retval==-13) message = 'IDA_LINESEARCH_FAIL' ! The line search failed.
   if( retval==-14) message = 'IDA_NO_RECOVERY'     ! The residual function, linear solver setup function, or linear solver solve function had a recoverable failure, but IDACalcIC could not recover.
   if( retval==-15) message = 'IDA_NLS_INIT_FAIL'   ! The nonlinear solver’s init routine failed.
   if( retval==-16) message = 'IDA_NLS_SETUP_FAIL'  ! The nonlinear solver’s setup routine failed.
   if( retval==-20) message = 'IDA_MEM_NULL'        ! The ida_mem argument was NULL.
   if( retval==-21) message = 'IDA_MEM_FAIL'        ! A memory allocation failed.
   if( retval==-22) message = 'IDA_ILL_INPUT'       ! One of the function inputs is illegal.
   if( retval==-23) message = 'IDA_NO_MALLOC'       ! The IDA memory was not allocated by a call to IDAInit.
   if( retval==-24) message = 'IDA_BAD_EWT'         ! Zero value of some error weight component.
   if( retval==-25) message = 'IDA_BAD_K'           ! The k-th derivative is not available.
   if( retval==-26) message = 'IDA_BAD_T'           ! The time t is outside the last step taken.
   if( retval==-27) message = 'IDA_BAD_DKY'         ! The vector argument where derivative should be stored is NULL.

end subroutine getErrMessage


end module summaSolve4ida_module



module fidaSolver_module


  !======= Inclusions ===========
  use, intrinsic :: iso_c_binding
  use nrtype
  use fida_datatypes

! access the global print flag
USE globalData,only:globalPrintFlag

! access missing values
USE globalData,only:integerMissing  ! missing integer
USE globalData,only:realMissing     ! missing double precision number
USE globalData,only:quadMissing     ! missing quadruple precision number

! access matrix information
USE globalData,only: nBands         ! length of the leading dimension of the band diagonal matrix
USE globalData,only: ixFullMatrix   ! named variable for the full Jacobian matrix
USE globalData,only: ixBandMatrix   ! named variable for the band diagonal matrix
USE globalData,only: iJac1          ! first layer of the Jacobian to print
USE globalData,only: iJac2          ! last layer of the Jacobian to print
USE globalData,only: ku             ! number of super-diagonal bands
USE globalData,only: kl             ! number of sub-diagonal bands

! domain types
USE globalData,only:iname_veg       ! named variables for vegetation
USE globalData,only:iname_snow      ! named variables for snow
USE globalData,only:iname_soil      ! named variables for soil

! state variable type
USE globalData,only:iname_nrgCanair ! named variable defining the energy of the canopy air space
USE globalData,only:iname_nrgCanopy ! named variable defining the energy of the vegetation canopy
USE globalData,only:iname_watCanopy ! named variable defining the mass of total water on the vegetation canopy
USE globalData,only:iname_liqCanopy ! named variable defining the mass of liquid water on the vegetation canopy
USE globalData,only:iname_nrgLayer  ! named variable defining the energy state variable for snow+soil layers
USE globalData,only:iname_watLayer  ! named variable defining the total water state variable for snow+soil layers
USE globalData,only:iname_liqLayer  ! named variable defining the liquid  water state variable for snow+soil layers
USE globalData,only:iname_matLayer  ! named variable defining the matric head state variable for soil layers
USE globalData,only:iname_lmpLayer  ! named variable defining the liquid matric potential state variable for soil layers
USE globalData,only:model_decisions ! model decision structure

! global metadata
USE globalData,only:flux_meta                        ! metadata on the model fluxes
USE globalData,only:diag_meta                        ! metadata on the model diagnostic variables
USE globalData,only:prog_meta                        ! metadata on the model prognostic variables
USE globalData,only:deriv_meta                       ! metadata on the model derivatives

! constants
USE multiconst,only:&
                    LH_fus,       & ! latent heat of fusion                (J K-1)
                    Tfreeze,      & ! temperature at freezing              (K)
                    iden_ice,     & ! intrinsic density of ice             (kg m-3)
                    iden_water      ! intrinsic density of liquid water    (kg m-3)

! provide access to indices that define elements of the data structures
USE var_lookup,only:iLookPROG       ! named variables for structure elements
USE var_lookup,only:iLookDIAG       ! named variables for structure elements
USE var_lookup,only:iLookFLUX       ! named variables for structure elements
USE var_lookup,only:iLookFORCE      ! named variables for structure elements
USE var_lookup,only:iLookPARAM      ! named variables for structure elements
USE var_lookup,only:iLookINDEX      ! named variables for structure elements
USE var_lookup,only:iLookDECISIONS  ! named variables for elements of the decision structure
 USE var_lookup,only:iLookDERIV     ! named variables for structure elements

! provide access to the derived types to define the data structures
USE data_types,only:&
                    var_i,        & ! data vector (i4b)
                    var_d,        & ! data vector (dp)
                    var_ilength,  & ! data vector with variable length dimension (i4b)
                    var_dlength,  & ! data vector with variable length dimension (dp)
                    zLookup,      & 
                    model_options   ! defines the model decisions

! look-up values for the choice of groundwater representation (local-column, or single-basin)
USE mDecisions_module,only:       &
 localColumn,                     & ! separate groundwater representation in each local soil column
 singleBasin                        ! single groundwater store over the entire basin

! look-up values for the choice of groundwater parameterization
USE mDecisions_module,only:      &
 qbaseTopmodel,                  & ! TOPMODEL-ish baseflow parameterization
 bigBucket,                      & ! a big bucket (lumped aquifer model)
 noExplicit                        ! no explicit groundwater parameterization

  

  ! privacy
  implicit none
  private
  public::fidaSolver


contains

 !-------------------
 ! * subroutine fidaSolver: solve F(y,y') = 0 by FIDA. Here, y is the state vector
 ! ------------------

 subroutine fidaSolver(                         &  
                       dt,                      & ! end of the current time step
                       h_init,                  & ! intent(in) initial stepsize
                       atol,                    & ! absolute telerance
                       rtol,                    & ! relative tolerance                                              
                       nSnow,                   & ! intent(in):    number of snow layers
                       nSoil,                   & ! intent(in):    number of soil layers
                       nLayers,                 & ! intent(in):    total number of layers
                       nState,                  & ! intent(in):    total number of state variables
                       ixMatrix,                & ! intent(in):    type of matrix (dense or banded)
                       ixQuadrature,            & ! intent(in):    type of quadrature method for approximating average flux
                       firstSubStep,            & ! intent(in):    flag to indicate if we are processing the first sub-step
                       firstFluxCall,           & ! intent(inout): flag to indicate if we are processing the first flux call
                       computeVegFlux,          & ! intent(in):    flag to indicate if we need to compute fluxes over vegetation
                       scalarSolution,          & ! intent(in):    flag to indicate the scalar solution
                       ! input: state vectors
                       stateVecInit,            & ! intent(in):    initial state vector              
                       fScale,                  & ! intent(in):    function scaling vector
                       sMul,                    & ! intent(inout):    state vector multiplier (used in the residual calculations)
                       dMat,                    & ! intent(inout)
                       numDiscon,               & ! intent(in)
                       ! input: data structures
                       lookup_data,             &
                       type_data,               & ! intent(in):    type of vegetation and soil
                       attr_data,               & ! intent(in):    spatial attributes
                       mpar_data,               & ! intent(in):    model parameters
                       forc_data,               & ! intent(in):    model forcing data
                       bvar_data,               & ! intent(in):    average model variables for the entire basin
                       prog_data,               & ! intent(in):    model prognostic variables for a local HRU
                       ! input-output: data structures
                       indx_data,               & ! intent(in): index data
                       diag_data,               & ! intent(inout): model diagnostic variables for a local HRU
                       flux_temp,               & ! intent(inout)
                       flux_data,               & ! intent(inout): model fluxes for a local HRU
                       flux_sum,                & ! intent(inout)
                       deriv_data,              & ! intent(inout): derivatives in model fluxes w.r.t. relevant state variables
                       ! input-output: baseflow
                       ixSaturation,            & ! intent(inout): index of the lowest saturated layer (NOTE: only computed on the first iteration)
                       dBaseflow_dMatric,       & ! intent(out):   derivative in baseflow w.r.t. matric head (s-1)
                       mLayerCmpress_sum,       &
                       ! output 
                       tret,                    & ! time which the solution is returned, if successfull tret = tend
                       dt_last,                 &
                       t_last,                  &
                       stepsize_past,           &
                       stateVec,                & ! intent(out):    model state vector
                       stateVecPrime,           & ! intent(out):    derivative of model state vector   
                       fluxVec,                 & ! intent(out):   flux vector
                       resSink,                 & ! intent(out):   additional (sink) terms on the RHS of the state equation
                       rVec,                    &
                       err,message              & ! intent(out):   error control
                      )

  !======= Inclusions ===========
  use, intrinsic :: iso_c_binding
  use nrtype
  use fida_datatypes
  use evalEqnsFida_module,only:evalEqnsFida
  use allocspace_module,only:allocLocal                ! allocate local data structures


  use fida_mod                      ! Fortran interface to IDA
  use fnvector_serial_mod           ! Fortran interface to serial N_Vector
  use fsunmatrix_dense_mod          ! Fortran interface to dense SUNMatrix
  use fsunlinsol_dense_mod          ! Fortran interface to dense SUNLinearSolver
  use fsunmatrix_band_mod           ! Fortran interface to banded SUNMatrix
  use fsunlinsol_band_mod           ! Fortran interface to banded SUNLinearSolver
  use fsunnonlinsol_newton_mod      ! Fortran interface to Newton SUNNonlinearSolver
  use fsundials_matrix_mod          ! Fortran interface to generic SUNMatrix
  use fsundials_nvector_mod         ! Fortran interface to generic N_Vector
  use fsundials_linearsolver_mod    ! Fortran interface to generic SUNLinearSolver
  use fsundials_nonlinearsolver_mod ! Fortran interface to generic SUNNonlinearSolver
  use evalEqnsFida_module,only:evalEqnsFida           ! ODE functions
  use evalJacFida_module,only:evalJacFida   
  use tolFida_module,only:computWeightFida
  use convTestFida_module,only:convTestFida
  use eval8summaFida_module,only:eval8summaFida
  USE computEnthalpy_module,only:computEnthalpy
!  use findDiscontinuity_module,only:findDiscontinuity

  !======= Declarations =========
  implicit none
  
 ! --------------------------------------------------------------------------------------------------------------------------------
 ! calling variables
 ! --------------------------------------------------------------------------------------------------------------------------------
 ! input: model control
 real(qp),intent(in)             :: dt
 real(qp),intent(in)             :: h_init
 real(qp),intent(inout)          :: atol(:)
 real(qp),intent(inout)          :: rtol(:)
 integer(i4b),intent(in)         :: nSnow                  ! number of snow layers
 integer(i4b),intent(in)         :: nSoil                  ! number of soil layers
 integer(i4b),intent(in)         :: nLayers                ! total number of layers
 integer(c_long),intent(in)      :: nState                 ! total number of state variables
 integer(i4b)                    :: ixMatrix               ! form of matrix (dense or banded)
 integer(i4b)                    :: ixQuadrature           ! type of quadrature method for approximating average flux
 logical(lgt),intent(in)         :: firstSubStep           ! flag to indicate if we are processing the first sub-step
 logical(lgt),intent(inout)      :: firstFluxCall          ! flag to indicate if we are processing the first flux call
 logical(lgt),intent(in)         :: computeVegFlux         ! flag to indicate if computing fluxes over vegetation
 logical(lgt),intent(in)         :: scalarSolution         ! flag to denote if implementing the scalar solution
 ! input: state vectors
 real(dp),intent(in)             :: stateVecInit(:)        ! model state vector
 real(dp),intent(in)             :: fScale(:)              ! function scaling vector
 real(qp),intent(inout)          :: sMul(:)   ! NOTE: qp   ! state vector multiplier (used in the residual calculations)
 real(dp), intent(inout)         :: dMat(:)
 integer(i4b), intent(in)        :: numDiscon
 ! input: data structures
 type(zLookup),intent(in)        :: lookup_data            ! lookup tables
 type(var_i),        intent(in)  :: type_data              ! type of vegetation and soil
 type(var_d),        intent(in)  :: attr_data              ! spatial attributes
 type(var_dlength),  intent(in)  :: mpar_data              ! model parameters
 type(var_d),        intent(in)  :: forc_data              ! model forcing data
 type(var_dlength),  intent(in)  :: bvar_data              ! model variables for the local basin
 type(var_dlength),  intent(in)  :: prog_data              ! prognostic variables for a local HRU
 ! output: data structures
 type(var_ilength),intent(in)    :: indx_data              ! indices defining model states and layers
 type(var_dlength),intent(inout) :: diag_data              ! diagnostic variables for a local HRU
 type(var_dlength),intent(inout) :: flux_temp             ! model fluxes for a local HRU
  type(var_dlength),intent(inout) :: flux_data              ! model fluxes for a local HRU
 type(var_dlength),intent(inout) :: flux_sum
 type(var_dlength),intent(inout) :: deriv_data             ! derivatives in model fluxes w.r.t. relevant state variables
 ! input-output: baseflow
 integer(i4b),intent(inout)      :: ixSaturation           ! index of the lowest saturated layer (NOTE: only computed on the first iteration)
 real(dp),intent(out)            :: dBaseflow_dMatric(:,:) ! derivative in baseflow w.r.t. matric head (s-1)
 real(dp),intent(inout)          :: mLayerCmpress_sum(:)
 ! output: flux and residual vectors
  real(dp),intent(inout)         :: stateVec(:)       ! model state vector
 real(dp),intent(inout)          :: stateVecPrime(:)       ! model state vector
 real(dp),intent(out)            :: fluxVec(:)             ! flux vector
 real(dp),intent(out)            :: resSink(:)             ! sink terms on the RHS of the flux equation
 real(dp),intent(out)            :: rVec(:)
 ! output: error control
 integer(i4b),intent(out)        :: err                    ! error code
 character(*),intent(out)        :: message                ! error message
 real(qp),intent(out)            :: tret(1)
 real(qp),intent(out)            :: dt_last(1)
 real(qp),intent(out)            :: t_last(1)
 real(qp),intent(out)            :: stepsize_past
 logical(lgt)                    :: scalling_on
  
 ! --------------------------------------------------------------------------------------------------------------------------------
 ! local variables
 ! --------------------------------------------------------------------------------------------------------------------------------
  type(N_Vector),           pointer :: sunvec_y      ! sundials solution vector
  type(N_Vector),           pointer :: sunvec_yp     ! sundials derivative vector
  type(N_Vector),           pointer :: sunvec_av     ! sundials tolerance vector
  type(SUNMatrix),          pointer :: sunmat_A      ! sundials matrix
  type(SUNLinearSolver),    pointer :: sunlinsol_LS  ! sundials linear solver
  type(SUNNonLinearSolver), pointer :: sunnonlin_NLS ! sundials nonlinear solver
  type(c_ptr)                       :: ida_mem       ! IDA memory
  type(eqnsData),           target  :: eqns_data
  integer(i4b)                      :: retval
  logical(lgt)                      :: feasible                      ! feasibility flag
  real(qp)                          :: t0
  integer (kind = 8) 				:: maxstep
  integer(kind = 8) 				:: mu, lu
  real(qp)            				:: h_max
  real(qp)            				:: coef_nonlin
  integer(i4b) 						:: lsflag
  integer(i4b)                   	:: iLayer                        ! index of layer in the snow+soil domain
  integer(i4b)                    	:: iState                        ! index of model state
  real(dp) :: idenIW
  real(c_double)             		:: stepsize_cur(1)
  integer(i4b),parameter     		:: ixRectangular=1    ! reza: should put them in a data type later
  integer(i4b),parameter     		:: ixTrapezoidal=2
  integer(i4b)               		:: iVar  
  logical(lgt)               		:: startQuadrature
  real(dp)  				 		:: mLayerMatricHeadLiqPrev(nSoil)
 globalVars: associate(& 
 nSnowSoilNrg            => indx_data%var(iLookINDEX%nSnowSoilNrg )%dat(1)         ,& ! intent(in): 
 ixSnowSoilNrg           => indx_data%var(iLookINDEX%ixSnowSoilNrg)%dat            ,& ! intent(in):
 layerType               => indx_data%var(iLookINDEX%layerType)%dat                ,&
 nSnowSoilHyd            => indx_data%var(iLookINDEX%nSnowSoilHyd )%dat(1)         ,& ! intent(in):
 ixSnowSoilHyd           => indx_data%var(iLookINDEX%ixSnowSoilHyd)%dat             & ! intent(in):
 )
  
  
  
  
  !======= Internals ============
  
  
  
  ! fill eqns_data which will be required later to call eval8summaFida 
  eqns_data%dt                      = dt
  eqns_data%nSnow                   = nSnow       
  eqns_data%nSoil                   = nSoil
  eqns_data%nLayers                 = nLayers
  eqns_data%nState                  = nState   
  eqns_data%ixMatrix                = ixMatrix           
  eqns_data%ixQuadrature            = ixQuadrature   
  eqns_data%firstSubStep            = firstSubStep
  eqns_data%firstFluxCall           = firstFluxCall 
  eqns_data%computeVegFlux          = computeVegFlux
  eqns_data%scalarSolution          = scalarSolution

  allocate( eqns_data%atol(nState) )
  eqns_data%atol = atol
  
  allocate( eqns_data%rtol(nState) )
  eqns_data%rtol = rtol
    
  allocate( eqns_data%fScale(nState) )
  eqns_data%fScale = fScale
  
  
  allocate( eqns_data%sMul(nState) )
  eqns_data%sMul                    = sMul
  
  allocate( eqns_data%dMat(nState) )
  eqns_data%dMat                   = dMat
  
 ! allocate space for the temporary prognostic variable structure
 call allocLocal(prog_meta(:),eqns_data%prog_data,nSnow,nSoil,err,message)
 if(err/=0)then; err=20; message=trim(message)//trim(message); return; endif
 eqns_data%prog_data               = prog_data

 ! allocate space for the temporary diagnostic variable structure
 call allocLocal(diag_meta(:),eqns_data%diag_data,nSnow,nSoil,err,message)
 if(err/=0)then; err=20; message=trim(message)//trim(message); return; endif
 eqns_data%diag_data               = diag_data

 ! allocate space for the temporary flux variable structure
 call allocLocal(flux_meta(:),eqns_data%flux_temp,nSnow,nSoil,err,message)
 if(err/=0)then; err=20; message=trim(message)//trim(message); return; endif
 eqns_data%flux_temp               = flux_temp

 ! allocate space for the temporary flux variable structure
 call allocLocal(flux_meta(:),eqns_data%flux_data,nSnow,nSoil,err,message)
 if(err/=0)then; err=20; message=trim(message)//trim(message); return; endif
 eqns_data%flux_data               = flux_data
 
 ! allocate space for the temporary flux variable structure
 call allocLocal(flux_meta(:),eqns_data%flux_sum,nSnow,nSoil,err,message)
 if(err/=0)then; err=20; message=trim(message)//trim(message); return; endif
 eqns_data%flux_sum               = flux_sum

 ! allocate space for the derivative structure
 call allocLocal(deriv_meta(:),eqns_data%deriv_data,nSnow,nSoil,err,message)
 if(err/=0)then; err=20; message=trim(message)//trim(message); return; end if
 eqns_data%deriv_data              = deriv_data
  
  eqns_data%lookup_data             = lookup_data
  eqns_data%type_data               = type_data
  eqns_data%attr_data               = attr_data
  eqns_data%mpar_data               = mpar_data    
  eqns_data%forc_data               = forc_data 
  eqns_data%bvar_data               = bvar_data
  eqns_data%indx_data               = indx_data
  eqns_data%ixSaturation            = ixSaturation
  
 ! allocate space for the baseflow derivatives
 ! NOTE: needs allocation because only used when baseflow sinks are active
 if(model_decisions(iLookDECISIONS%groundwatr)%iDecision==qbaseTopmodel)then
  allocate(eqns_data%dBaseflow_dMatric(nSoil,nSoil),stat=err)  ! baseflow depends on total storage in the soil column
 else
  allocate(eqns_data%dBaseflow_dMatric(0,0),stat=err)          ! allocate zero-length dimnensions to avoid passing around an unallocated matrix
 end if
 
  eqns_data%dBaseflow_dMatric       = dBaseflow_dMatric
  
  allocate( eqns_data%mLayerCmpress_sum(nSoil) )
  eqns_data%mLayerCmpress_sum = mLayerCmpress_sum
  
  allocate( eqns_data%mLayerMatricHeadLiqTrial(nSoil) )

  allocate( eqns_data%mLayerMatricHeadTrial(nSoil) )
  allocate( eqns_data%mLayerMatricHeadPrev(nSoil) )
    
  allocate( eqns_data%mLayerVolFracWatTrial(nLayers) )
  allocate( eqns_data%mLayerVolFracWatPrev(nLayers) )
  
  allocate( eqns_data%mLayerTempTrial(nLayers) )
  allocate( eqns_data%mLayerTempPrev(nLayers) )
  
  allocate( eqns_data%mLayerVolFracIceTrial(nLayers) )
  allocate( eqns_data%mLayerVolFracIcePrev(nLayers) )
  
  allocate( eqns_data%mLayerEnthalpyTrial(nLayers) )
  allocate( eqns_data%mLayerEnthalpyPrev(nLayers) )
    
  allocate( eqns_data%fluxVec(nState) )
  eqns_data%fluxVec                 = fluxVec
  
  allocate( eqns_data%resSink(nState) )
  eqns_data%resSink                 = resSink
  
  
  allocate( eqns_data%resVec(nState) )
  
  eqns_data%err                     = err
  eqns_data%message                 = message
  eqns_data%startQuadrature         = .true.
  
  

  ! create serial vectors
  sunvec_y => FN_VMake_Serial(nState, stateVec)
  if (.not. associated(sunvec_y)) then
     print *, 'ERROR: sunvec = NULL'
     stop 1
  end if

  sunvec_yp => FN_VMake_Serial(nState, stateVecPrime)
  if (.not. associated(sunvec_yp)) then
     print *, 'ERROR: sunvec = NULL'
     stop 1
  end if
  
  ! Initialize solution vectors   reza: we should provide and pass stateVecPrimeInit later
  call setInitialCondition(nState, stateVecInit, sunvec_y, sunvec_yp)


  ! Call FIDACreate and FIDAInit to initialize IDA memory
  ida_mem = FIDACreate()
  if (.not. c_associated(ida_mem)) then
     print *, 'ERROR: ida_mem = NULL'
     stop 1
  end if
  
  eqns_data%ida_mem = ida_mem


  retval = FIDASetUserData(ida_mem, c_loc(eqns_data))
  if (retval /= 0) then
     print *, 'Error in FIDASetUserData, retval = ', retval, '; halting'
     stop 1
  end if
  
  t0 = 0._dp
  retval = FIDAInit(ida_mem, c_funloc(evalEqnsFida), t0, sunvec_y, sunvec_yp)
  if (retval /= 0) then
     print *, 'Error in FIDAInit, retval = ', retval, '; halting'
     stop 1
  end if

  ! set tolerances
  retval = FIDAWFtolerances(ida_mem, c_funloc(computWeightFida))
  if (retval /= 0) then
     print *, 'Error in FIDAWFtolerances, retval = ', retval, '; halting'
     stop 1
  end if
  
  ! specify the root function
!  retval = FIDARootInit(ida_mem, numDiscon, c_funloc(findDiscontinuity))
!  if (retval /= 0) then
!     print *, 'Error in FIDARootInit, retval = ', retval, '; halting'
!     stop 1
!  end if
  
 ! define the form of the matrix
 select case(ixMatrix)
  case(ixBandMatrix)
  mu = ku; lu = kl;
  ! Create banded SUNMatrix for use in linear solves
     sunmat_A => FSUNBandMatrix(nState, mu, lu)
     if (.not. associated(sunmat_A)) then
        print *, 'ERROR: sunmat = NULL'
        stop 1
     end if

     ! Create banded SUNLinearSolver object
     sunlinsol_LS => FSUNLinSol_Band(sunvec_y, sunmat_A)
     if (.not. associated(sunlinsol_LS)) then
        print *, 'ERROR: sunlinsol = NULL'
        stop 1
     end if
  
  case(ixFullMatrix)
    ! Create dense SUNMatrix for use in linear solves
     sunmat_A => FSUNDenseMatrix(nState, nState)
     if (.not. associated(sunmat_A)) then
        print *, 'ERROR: sunmat = NULL'
        stop 1
     end if

     ! Create dense SUNLinearSolver object
     sunlinsol_LS => FSUNDenseLinearSolver(sunvec_y, sunmat_A)
     if (.not. associated(sunlinsol_LS)) then
        print *, 'ERROR: sunlinsol = NULL'
        stop 1
     end if
   
  ! check
  case default;  print *, 'unable to identify option for the type of matrix'; stop 1
  
 end select  ! form of matrix

  ! Attach the matrix and linear solver
  retval = FIDASetLinearSolver(ida_mem, sunlinsol_LS, sunmat_A);
  if (retval /= 0) then
     print *, 'Error in FIDASetLinearSolver, retval = ', retval, '; halting'
     stop 1
  end if
  
  
  
  if(ixMatrix == ixFullMatrix)then
     ! Set the user-supplied Jacobian routine    
    retval = FIDASetJacFn(ida_mem, c_funloc(evalJacFida))
   if (retval /= 0) then
      print *, 'Error in FIDASetJacFn, retval = ', retval, '; halting'
      stop 1
   end if
  
  endif

  ! Create Newton SUNNonlinearSolver object. IDA uses a
  ! Newton SUNNonlinearSolver by default, so it is not necessary
  ! to create it and attach it. It is done here for demonstration purposes.
  sunnonlin_NLS => FSUNNonlinSol_Newton(sunvec_y)
  if (.not. associated(sunnonlin_NLS)) then
     print *, 'ERROR: sunnonlinsol = NULL'
     stop 1
  end if
  

  ! Attach the nonlinear solver
  retval = FIDASetNonlinearSolver(ida_mem, sunnonlin_NLS)
  if (retval /= 0) then
     print *, 'Error in FIDASetNonlinearSolver, retval = ', retval, '; halting'
     stop 1
  end if
  
!  Set the user-supplied nonlinear solver convergence test function 
!    retval = FSUNNonlinSolSetConvTestFn(sunnonlin_NLS, c_funloc(convTestFida), c_loc(eqns_data));
!  if (retval /= 0) then
!     print *, 'Error in FSUNNonlinSolSetConvTestFn, retval = ', retval, '; halting'
!     stop 1
!  end if
  
  
  ! Set the maximum BDF order,  default = 5
  retval = FIDASetMaxOrd(ida_mem, 5)
  if (retval /= 0) then
     print *, 'Error in FIDASetMaxOrd, retval = ', retval, '; halting'
     stop 1
  end if
  
  
  ! Set Coeff. in the nonlinear convergence test, default = 0.33
  coef_nonlin = 0.33
  retval = FIDASetNonlinConvCoef(ida_mem, coef_nonlin)
  if (retval /= 0) then
     print *, 'Error in FIDASetNonlinConvCoef, retval = ', retval, '; halting'
     stop 1
  end if
  
   
  ! Set maximun number of nonliear iterations, default = 4
  retval = FIDASetMaxNonlinIters(ida_mem, 100)
  if (retval /= 0) then
     print *, 'Error in IDASetMaxNonlinIters, retval = ', retval, '; halting'
     stop 1
  end if
  
  !  Set maximum number of convergence test failures, default = 10
  retval = FIDASetMaxConvFails(ida_mem, 50)
  if (retval /= 0) then
     print *, 'Error in FIDASetMaxConvFails, retval = ', retval, '; halting'
     stop 1
  end if
  
  !  Set maximum number of error test failures, default = 10
  retval = FIDASetMaxErrTestFails(ida_mem, 50)
  if (retval /= 0) then
     print *, 'Error in FIDASetMaxErrTestFails, retval = ', retval, '; halting'
     stop 1
  end if
  
  ! Set maximum number of steps,  dafault = 500
   maxstep = 999999
  retval = FIDASetMaxNumSteps(ida_mem, maxstep)
  if (retval /= 0) then
     print *, 'Error in FIDASetMaxNumSteps, retval = ', retval, '; halting'
     stop 1
  end if
  
  ! Set maximum stepsize,  dafault = infinity
 ! h_max = 1
 ! retval = FIDASetMaxStep(ida_mem, h_max)
 ! if (retval /= 0) then
 !    print *, 'Error in FIDASetMaxSteps, retval = ', retval, '; halting'
 !    stop 1
 ! end if
  
  ! Set initial stepsize
  retval = FIDASetInitStep(ida_mem, h_init)
  if (retval /= 0) then
     print *, 'Error in FIDASetInitStep, retval = ', retval, '; halting'
     stop 1
  end if
  
  retval = FIDASetStopTime(ida_mem, dt)
  if (retval /= 0) then
     print *, 'Error in FIDASetStopTime, retval = ', retval, '; halting'
     stop 1
  end if
  
 
 ! retval = FIDASetLinearSolutionScaling(ida_mem, 0)  ! scalling_on 1 and off 0
 
 tret(1) = t0
 ! need the following values for the first substep
 eqns_data%scalarCanopyTempPrev		= prog_data%var(iLookPROG%scalarCanopyTemp)%dat(1)
 eqns_data%scalarCanopyIcePrev      = prog_data%var(iLookPROG%scalarCanopyIce)%dat(1) 
 eqns_data%mLayerVolFracWatPrev(:) 	= prog_data%var(iLookPROG%mLayerVolFracWat)%dat(:)
 eqns_data%mLayerTempPrev(:) 		= prog_data%var(iLookPROG%mLayerTemp)%dat(:)
 eqns_data%mLayerVolFracIcePrev(:) 	= prog_data%var(iLookPROG%mLayerVolFracIce)%dat(:)   
 eqns_data%mLayerMatricHeadPrev(:) 	= prog_data%var(iLookPROG%mLayerMatricHead)%dat(:) 
 eqns_data%mLayerEnthalpyPrev(:) 	= diag_data%var(iLookDIAG%mLayerEnthalpy)%dat(:)
 eqns_data%scalarCanopyEnthalpyPrev = diag_data%var(iLookDIAG%scalarCanopyEnthalpy)%dat(1)
 mLayerMatricHeadLiqPrev(:) 		= diag_data%var(iLookDIAG%mLayerMatricHeadLiq)%dat(:)
 
 
 !**********************************************************************************
 !****************************** Main Solver ***************************************
 !************************* loop on one_step mode **********************************
 !**********************************************************************************                                 
 do while(tret(1) < dt) 
  ! call IDASolve
  retval = FIDASolve(ida_mem, dt, tret, sunvec_y, sunvec_yp, IDA_ONE_STEP) 
  
  ! get the last stepsize  
  retval = FIDAGetLastStep(ida_mem, stepsize_cur)
  if (retval /= 0) then
     print *, 'Error in FIDAGetLastStep, retval = ', retval, '; halting'
     stop 1
  end if  
    
    ! compute the flux and the residual vector for a given state vector
  call eval8summaFida(&
                 ! input: model control
                 stepsize_cur(1),                   &
                 eqns_data%dt,                      &
                 eqns_data%nSnow,                   & ! intent(in):    number of snow layers
                 eqns_data%nSoil,                   & ! intent(in):    number of soil layers
                 eqns_data%nLayers,                 & ! intent(in):    number of layers
                 eqns_data%nState,                  & ! intent(in):    number of state variables in the current subset
                 eqns_data%firstSubStep,            & ! intent(in):    flag to indicate if we are processing the first sub-step
                 eqns_data%firstFluxCall,           & ! intent(inout): flag to indicate if we are processing the first flux call
                 eqns_data%computeVegFlux,          & ! intent(in):    flag to indicate if we need to compute fluxes over vegetation
                 eqns_data%scalarSolution,          & ! intent(in):    flag to indicate the scalar solution
                 ! input: state vectors
                 stateVec,                          & ! intent(in):    model state vector
                 stateVecPrime,                     & ! intent(in):    model state vector
                 eqns_data%sMul,                    & ! intent(inout):    state vector multiplier (used in the residual calculations)
                 ! input: data structures
                 model_decisions,                   & ! intent(in):    model decisions
                 eqns_data%lookup_data,             &
                 eqns_data%type_data,               & ! intent(in):    type of vegetation and soil
                 eqns_data%attr_data,               & ! intent(in):    spatial attributes
                 eqns_data%mpar_data,               & ! intent(in):    model parameters
                 eqns_data%forc_data,               & ! intent(in):    model forcing data
                 eqns_data%bvar_data,               & ! intent(in):    average model variables for the entire basin
                 eqns_data%prog_data,               & ! intent(in):    model prognostic variables for a local HRU
                 ! input-output: data structures
                 eqns_data%indx_data,               & ! intent(inou):    index data
                 eqns_data%diag_data,               & ! intent(inout): model diagnostic variables for a local HRU
                 eqns_data%flux_data,               & ! intent(inout): model fluxes for a local HRU (initial flux structure)
                 eqns_data%deriv_data,              & ! intent(inout): derivatives in model fluxes w.r.t. relevant state variables
                 ! input-output: baseflow
                 eqns_data%ixSaturation,             & ! intent(inout): index of the lowest saturated layer (NOTE: only computed on the first iteration)
                 eqns_data%dBaseflow_dMatric,        & ! intent(out):   derivative in baseflow w.r.t. matric head (s-1), we will use it later for Jacobian
                 eqns_data%scalarCanopyTempTrial,    & ! intent(in):  trial value of canopy temperature (K)
                 eqns_data%scalarCanopyTempPrev,     & ! intent(in):  previous value of canopy temperature (K)
                 eqns_data%scalarCanopyIceTrial,	 &
                 eqns_data%scalarCanopyIcePrev,		 &
                 eqns_data%scalarCanopyEnthalpyTrial,& ! intent(in):  trial enthalpy of the vegetation canopy (J m-3)
                 eqns_data%scalarCanopyEnthalpyPrev, & ! intent(in):  previous enthalpy of the vegetation canopy (J m-3)
                 eqns_data%mLayerTempTrial,          &
                 eqns_data%mLayerTempPrev,           &
                 eqns_data%mLayerMatricHeadLiqTrial, &
                 eqns_data%mLayerMatricHeadTrial,    &
                 eqns_data%mLayerMatricHeadPrev,     &
                 eqns_data%mLayerVolFracWatTrial,    &
                 eqns_data%mLayerVolFracWatPrev,     &
                 eqns_data%mLayerVolFracIceTrial,    &
                 eqns_data%mLayerVolFracIcePrev,     &
                 eqns_data%mLayerEnthalpyPrev,       & ! intent(in)
                 eqns_data%mLayerEnthalpyTrial,      & ! intent(out)
                 ! output
                 feasible,                          & ! intent(out):   flag to denote the feasibility of the solution
                 eqns_data%fluxVec,                 & ! intent(out):   flux vector
                 eqns_data%resSink,                 & ! intent(out):   additional (sink) terms on the RHS of the state equation
                 rVec,                  			& ! intent(out):   residual vector
                 eqns_data%err,eqns_data%message)     ! intent(out):   error control 
                 
  
  select case(ixQuadrature)
       ! sum of flux
       case(ixRectangular)
            do  iVar=1,size(flux_meta) 
              flux_sum%var(iVar)%dat(:) = flux_sum%var(iVar)%dat(:) + eqns_data%flux_data%var(iVar)%dat(:) *  stepsize_cur(1) 
            end do
       case(ixTrapezoidal)
            if(startQuadrature)then
                 do  iVar=1,size(flux_meta) 
                     flux_sum%var(iVar)%dat(:) =  flux_temp%var(iVar)%dat(:) *  stepsize_cur(1) 
                 end do 
                 startQuadrature = .false.
            else
                 do  iVar=1,size(flux_meta) 
                     flux_sum%var(iVar)%dat(:) = flux_sum%var(iVar)%dat(:) + flux_temp%var(iVar)%dat(:) & 
                                                                           *  ( stepsize_past + stepsize_cur(1) )
                 end do        
           endif 
           do  iVar=1,size(flux_meta) 
              flux_temp%var(iVar)%dat(:) = flux_data%var(iVar)%dat(:) 
           end do
           stepsize_past = stepsize_cur(1)
       case default; err=20; message=trim(message)//'expect case to be ixRecangular, ixTrapezoidal'; return
  end select
  
   ! sum of mLayerCmpress
   mLayerCmpress_sum(:) = mLayerCmpress_sum(:) + eqns_data%deriv_data%var(iLookDERIV%dCompress_dPsi)%dat(:) &
                                    * ( eqns_data%mLayerMatricHeadLiqTrial(:) - mLayerMatricHeadLiqPrev(:) )
   ! save values of some quantities for next step
   eqns_data%scalarCanopyTempPrev		= eqns_data%scalarCanopyTempTrial
   eqns_data%scalarCanopyIcePrev		= eqns_data%scalarCanopyIceTrial
   eqns_data%mLayerTempPrev(:) 			= eqns_data%mLayerTempTrial(:)
   mLayerMatricHeadLiqPrev(:) 			= eqns_data%mLayerMatricHeadLiqTrial(:)
   eqns_data%mLayerMatricHeadPrev(:) 	= eqns_data%mLayerMatricHeadTrial(:)
   eqns_data%mLayerVolFracWatPrev(:) 	= eqns_data%mLayerVolFracWatTrial(:)
   eqns_data%mLayerVolFracIcePrev(:) 	= eqns_data%mLayerVolFracIceTrial(:)
   eqns_data%mLayerEnthalpyPrev(:) 		= eqns_data%mLayerEnthalpyTrial(:)
   eqns_data%scalarCanopyEnthalpyPrev 	= eqns_data%scalarCanopyEnthalpyTrial

 end do ! while loop on one_step mode
 
 !****************************** End of Main Solver ***************************************
 
  firstFluxCall 	= eqns_data%firstFluxCall
  fluxVec 			= eqns_data%fluxVec         
  diag_data 		= eqns_data%diag_data 
  flux_temp 		= eqns_data%flux_temp             
  flux_data 		= eqns_data%flux_data             
  deriv_data 		= eqns_data%deriv_data  
  dBaseflow_dMatric = eqns_data%dBaseflow_dMatric  
  ixSaturation 		= eqns_data%ixSaturation   
  resSink 			= eqns_data%resSink 
  stepsize_past 	= eqns_data%stepsize_past
  err 				= eqns_data%err
  message 			= eqns_data%message        


  retval = FIDAGetLastStep(ida_mem, dt_last)
  if (retval /= 0) then
     print *, 'Error in FIDAGetLastStep, retval = ', retval, '; halting'
     stop 1
  end if
  
  retval = FIDAGetCurrentTime(ida_mem, t_last)
  if (retval /= 0) then
     print *, 'Error in FIDAGetCurrentTime, retval = ', retval, '; halting'
     stop 1
  end if
  

  ! free memory  
  deallocate(eqns_data%fScale)
  deallocate(eqns_data%sMul)
  deallocate(eqns_data%dMat)
  deallocate(eqns_data%dBaseflow_dMatric)
  deallocate(eqns_data%mLayerCmpress_sum)
  deallocate(eqns_data%mLayerMatricHeadLiqTrial)
  deallocate(eqns_data%mLayerMatricHeadTrial)
  deallocate(eqns_data%mLayerMatricHeadPrev)
  deallocate( eqns_data%fluxVec )  
  deallocate( eqns_data%resSink )
  deallocate( eqns_data%resVec )
  deallocate( eqns_data%mLayerVolFracWatTrial )
  deallocate( eqns_data%mLayerVolFracWatPrev )
  deallocate( eqns_data%mLayerVolFracIceTrial )
  deallocate( eqns_data%mLayerTempPrev )
  deallocate( eqns_data%mLayerTempTrial )
  deallocate( eqns_data%mLayerVolFracIcePrev )
  deallocate( eqns_data%mLayerEnthalpyTrial )
  deallocate( eqns_data%mLayerEnthalpyPrev )
  
  call FIDAFree(ida_mem)
  retval = FSUNNonlinSolFree(sunnonlin_NLS)
  retval = FSUNLinSolFree(sunlinsol_LS)
  call FSUNMatDestroy(sunmat_A)
  call FN_VDestroy(sunvec_y)
  call FN_VDestroy(sunvec_yp)
  
   ! end associate statements
   end associate globalVars

  

 end subroutine fidaSolver
 

! ----------------------------------------------------------------
! SetInitialCondition: routine to initialize u and up vectors.
! ----------------------------------------------------------------
subroutine setInitialCondition(neq, y, sunvec_u, sunvec_up)

  !======= Inclusions ===========
  use, intrinsic :: iso_c_binding
  use fsundials_nvector_mod
  use fnvector_serial_mod

  !======= Declarations =========
  implicit none

  ! calling variables
  type(N_Vector)  :: sunvec_u  ! solution N_Vector
  type(N_Vector)  :: sunvec_up ! derivative N_Vector
  integer(c_long) :: neq
  real(dp)        :: y(neq)

  ! pointers to data in SUNDIALS vectors
  real(c_double), pointer :: uu(:)
  real(c_double), pointer :: up(:)


  !======= Internals ============

  ! get data arrays from SUNDIALS vectors
  uu(1:neq) => FN_VGetArrayPointer(sunvec_u)
  up(1:neq) => FN_VGetArrayPointer(sunvec_up)


  uu = y
  up = 0._dp


end subroutine setInitialCondition

end module fidaSolver_module








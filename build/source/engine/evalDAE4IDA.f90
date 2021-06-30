

module evalDAE4IDA_module


  !======= Inclusions ===========
  use, intrinsic :: iso_c_binding
  use nrtype
  use type4IDA
  USE globalData,only:model_decisions        ! model decision structure
  USE globalData,only:flux_meta                        ! metadata on the model fluxes
  ! provide access to the derived types to define the data structures
  USE data_types,only:&
                    var_i,        & ! data vector (i4b)
                    var_d,        & ! data vector (rkind)
                    var_ilength,  & ! data vector with variable length dimension (i4b)
                    var_dlength,  & ! data vector with variable length dimension (rkind)
                    model_options   ! defines the model decisions
 USE multiconst,only:iden_water      ! intrinsic density of liquid water    (kg m-3)
 USE var_lookup,only:iLookDIAG 
 USE var_lookup,only:iLookPROG 
 USE var_lookup,only:iLookINDEX      ! named variables for structure elements
 USE var_lookup,only:iLookDERIV                   ! named variables for structure elements
  

  ! privacy
  implicit none
  private
  public::evalDAE4IDA


contains

  ! **********************************************************************************************************
  ! public function evalDAE4IDA: compute the residual vector F(t,y,y') required for FIDA solver
  ! **********************************************************************************************************
  ! Return values:
  !    0 = success,
  !    1 = recoverable error,
  !   -1 = non-recoverable error
  ! ----------------------------------------------------------------
  integer(c_int) function evalDAE4IDA(tres, sunvec_y, sunvec_yp, sunvec_r, user_data) &
       result(ierr) bind(C,name='evalDAE4IDA')

    !======= Inclusions ===========
    use, intrinsic :: iso_c_binding
    use fida_mod 
    use fsundials_nvector_mod
    use fnvector_serial_mod
    use nrtype
    use type4IDA
    use eval8DAE_module,only:eval8DAE

    !======= Declarations =========
    implicit none

    ! calling variables
    real(rkind), value         :: tres      ! current time                 t
    type(N_Vector)          :: sunvec_y  ! solution N_Vector    y
    type(N_Vector)          :: sunvec_yp ! derivative N_Vector  y'
    type(N_Vector)          :: sunvec_r  ! residual N_Vector    F(t,y,y')
    type(c_ptr), value      :: user_data ! user-defined data  


    ! pointers to data in SUNDIALS vectors
    type(eqnsData), pointer    :: eqns_data ! equations data
    real(rkind), pointer          :: stateVec(:)
    real(rkind), pointer          :: stateVecPrime(:)
    real(rkind), pointer          :: rVec(:)
    logical(lgt)               :: feasible   
    integer(i4b)               :: retval  
    real(c_double)             :: stepsize_next(1)



    
    

    !======= Internals ============
    
    ! get equations data from user-defined data
    call c_f_pointer(user_data, eqns_data)

     
 
    ! get data arrays from SUNDIALS vectors
    stateVec  => FN_VGetArrayPointer(sunvec_y)
    stateVecPrime => FN_VGetArrayPointer(sunvec_yp)
    rVec  => FN_VGetArrayPointer(sunvec_r) 
    
   retval = FIDAGetCurrentStep(eqns_data%ida_mem, stepsize_next)
   if (retval /= 0) then
     print *, 'Error in FIDAGetCurrentStep, retval = ', retval, '; halting'
     stop 1
  end if 
 

    ! compute the flux and the residual vector for a given state vector
    call eval8DAE(&
                 ! input: model control
                 stepsize_next(1),                  &
                 eqns_data%dt,                      &
                 eqns_data%nSnow,                   & ! intent(in):    number of snow layers
                 eqns_data%nSoil,                   & ! intent(in):    number of soil layers
                 eqns_data%nLayers,                 & ! intent(in):    number of layers
                 eqns_data%nState,                  & ! intent(in):    number of state variables in the current subset
                 eqns_data%firstSubStep,            & ! intent(in):    flag to indicate if we are processing the first sub-step
                 eqns_data%firstFluxCall,           & ! intent(inout)
                 eqns_data%firstSplitOper,			& ! intent(in)
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
                 eqns_data%dBaseflow_dMatric,       & ! intent(out):   derivative in baseflow w.r.t. matric head (s-1), we will use it later for Jacobian
                 eqns_data%scalarCanopyTempTrial,    & ! intent(in):  trial value of canopy temperature (K)
                 eqns_data%scalarCanopyTempPrev,     & ! intent(in):  previous value of canopy temperature (K)
                 eqns_data%scalarCanopyIceTrial,	&
                 eqns_data%scalarCanopyIcePrev,		&
                 eqns_data%scalarCanopyLiqTrial,	 &
                 eqns_data%scalarCanopyLiqPrev,		 &
                 eqns_data%scalarCanopyEnthalpyTrial,& ! intent(in):  trial enthalpy of the vegetation canopy (J m-3)
                 eqns_data%scalarCanopyEnthalpyPrev, & ! intent(in):  previous enthalpy of the vegetation canopy (J m-3)
                 eqns_data%mLayerTempTrial,         &
                 eqns_data%mLayerTempPrev,          &
                 eqns_data%mLayerMatricHeadLiqTrial,&
                 eqns_data%mLayerMatricHeadTrial, 	&
                 eqns_data%mLayerMatricHeadPrev, 	&
                 eqns_data%mLayerVolFracWatTrial,   &
                 eqns_data%mLayerVolFracWatPrev,    &
                 eqns_data%mLayerVolFracIceTrial,   &
                 eqns_data%mLayerVolFracIcePrev,   	&
                 eqns_data%mLayerVolFracLiqTrial,   &
                 eqns_data%mLayerVolFracLiqPrev,   	&
                 eqns_data%scalarAquiferStorageTrial, &
                 eqns_data%scalarAquiferStoragePrev, &   
                 eqns_data%mLayerEnthalpyPrev,      & ! intent(in)
                 eqns_data%mLayerEnthalpyTrial,     & ! intent(out)   
                 eqns_data%ixSaturation,			& ! intent(inout)             
                 ! output
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

 end function evalDAE4IDA


end module evalDAE4IDA_module

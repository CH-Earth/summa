

module findDiscontinuity_module


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
  

  ! privacy
  implicit none
  private
  public::findDiscontinuity


contains

  ! **********************************************************************************************************
  ! public function findDiscontinuity: the root function 
  ! **********************************************************************************************************
  ! Return values:
  !    0 = success,
  !    1 = recoverable error,
  !   -1 = non-recoverable error
  ! ----------------------------------------------------------------
  integer(c_int) function findDiscontinuity(tres, sunvec_y, sunvec_yp, fval, user_data) &
       result(ierr) bind(C,name='findDiscontinuity')

    !======= Inclusions ===========
    use, intrinsic :: iso_c_binding
    use fida_mod 
    use fsundials_nvector_mod
    use fnvector_serial_mod
    use nrtype
    use type4IDA
    use varExtrFida_module, only:residDiscontinuity

    !======= Declarations =========
    implicit none

    ! calling variables
    real(rkind), value         :: tres      ! current time                 
    type(N_Vector)          :: sunvec_y  ! solution N_Vector    y
    type(N_Vector)          :: sunvec_yp ! derivative N_Vector  y'
    real(c_double)          :: fval(1000)   ! root function values         
    type(c_ptr), value      :: user_data ! user-defined data  


    ! pointers to data in SUNDIALS vectors
    type(eqnsData), pointer    :: eqns_data ! equations data
    real(rkind), pointer          :: stateVec(:)
    real(rkind), pointer          :: stateVecPrime(:) 
    integer(i4b)               :: retval  

    
    

    !======= Internals ============
    
    ! get equations data from user-defined data
    call c_f_pointer(user_data, eqns_data)
 
    ! get data arrays from SUNDIALS vectors
    stateVec  => FN_VGetArrayPointer(sunvec_y)
  !  stateVecPrime => FN_VGetArrayPointer(sunvec_yp) 
    
   call residDiscontinuity(&
                       ! input
                       stateVec,                                  & ! intent(in):    model state vector (mixed units)
                       eqns_data%diag_data,                       & ! intent(in):    model diagnostic variables for a local HRU
                       eqns_data%prog_data,                       & ! intent(in):    model prognostic variables for a local HRU
                       eqns_data%indx_data,                       & ! intent(in):    indices defining model states and layers
                       ! output
                       fval,                                      & ! intent(out) 
                       eqns_data%err,eqns_data%message)                                 ! intent(out):   error control

  
   ! return success
   ierr = 0
   return

 end function findDiscontinuity


end module findDiscontinuity_module

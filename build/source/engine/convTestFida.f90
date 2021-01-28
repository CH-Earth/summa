

module convTestFida_module


  !======= Inclusions ===========
  use, intrinsic :: iso_c_binding
  use nrtype
  use fida_datatypes
  USE globalData,only:model_decisions        ! model decision structure
  USE globalData,only:flux_meta                        ! metadata on the model fluxes
  ! access missing values
  USE globalData,only:integerMissing  ! missing integer
  USE globalData,only:realMissing     ! missing double precision number
  
   ! indices that define elements of the data structures
   USE var_lookup,only:iLookPARAM                   ! named variables for structure elements
   USE var_lookup,only:iLookPROG                    ! named variables for structure elements
   USE var_lookup,only:iLookINDEX                   ! named variables for structure elements

   USE multiconst,only:&
                    Cp_air,       & ! specific heat of air                 (J kg-1 K-1)
                    iden_air,     & ! intrinsic density of air             (kg m-3)
                    iden_ice,     & ! intrinsic density of ice             (kg m-3)
                    iden_water      ! intrinsic density of liquid water    (kg m-3)

  ! provide access to the derived types to define the data structures
  USE data_types,only:&
                    var_i,        & ! data vector (i4b)
                    var_d,        & ! data vector (dp)
                    var_ilength,  & ! data vector with variable length dimension (i4b)
                    var_dlength,  & ! data vector with variable length dimension (dp)
                    model_options   ! defines the model decisions

  

  ! privacy
  implicit none
  private::checkConvFida
  public::convTestFida


contains

  ! **********************************************************************************************************
  ! public function convTestFida: 
  ! **********************************************************************************************************
  ! Return values:
  !    0 = success,
  !   -1 = fail
  ! ----------------------------------------------------------------
  integer(c_int) function convTestFida(NLS, sunvec_ycor, sunvec_del, tol, sunvec_ewt, user_data) &
       result(ierr) bind(C,name='convTestFida')

    !======= Inclusions ===========
    use, intrinsic :: iso_c_binding
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
    use nrtype
    use fida_datatypes


    !======= Declarations =========
    implicit none

    ! calling variables
    type(SUNNonLinearSolver)            ::  NLS
    type(N_Vector)                      :: sunvec_ycor  
    type(N_Vector)                      :: sunvec_del 
    real(dp), value                     :: tol
    type(N_Vector)                      :: sunvec_ewt  
    type(c_ptr), value                  :: user_data   
    


    ! pointers to data in SUNDIALS vectors
    type(eqnsData), pointer   :: eqns_data ! equations data
    real(dp), pointer         :: xVec(:)
    real(dp), pointer         :: xInc(:)
    logical(lgt)              :: isConv
    integer(i4b)              :: mSoil                    ! number of soil layers in solution vector
    
    ! get equations data from user-defined data
    call c_f_pointer(user_data, eqns_data)
    xVec  => FN_VGetArrayPointer(sunvec_ycor)
    xInc  => FN_VGetArrayPointer(sunvec_del)
    
    ! get the number of soil layers in the solution vector
   mSoil = size(eqns_data%indx_data%var(iLookINDEX%ixMatOnly)%dat)
    
   isConv = .true.
  
   
   call checkConvFida(xVec,xInc,eqns_data%nSnow,eqns_data%resVec,eqns_data%mpar_data,eqns_data%prog_data,eqns_data%indx_data,eqns_data%scalarSolution,msoil,isConv)

  ! convergence check
  if(isConv)then
   ierr = 0
  else 
   ierr = SUN_NLS_CONTINUE
  endif

  return

 end function convTestFida
 
  ! *********************************************************************************************************
  ! private subroutine checkConvFida: check convergence based on the residual vector
  ! *********************************************************************************************************
 
 subroutine checkConvFida(xVec,xInc,nSnow,resVec,mpar_data,prog_data,indx_data,scalarSolution,msoil,isConv)
  
    implicit none
    
    real(dp),intent(in)             :: xVec(:)
    real(dp),intent(in)             :: xInc(:)
    integer(i4b),intent(in)         :: nSnow
    real(qp),intent(in)             :: resVec(:)              ! residual vector
    type(var_dlength),intent(in)    :: mpar_data              ! model parameters
    type(var_dlength),intent(in)    :: prog_data              ! prognostic variables for a local HRU
    type(var_ilength),intent(in)    :: indx_data              ! indices defining model states and layers  
    integer(i4b),intent(in)         :: mSoil                    ! number of soil layers in solution vector
    logical(lgt),intent(in)         :: scalarSolution         ! flag to denote if implementing the scalar solution
    logical(lgt),intent(out)        :: isConv
    ! locals
    real(dp),dimension(mSoil) :: psiScale               ! scaling factor for matric head
    real(dp),parameter        :: xSmall=1.e-0_dp        ! a small offset
    real(dp),parameter        :: scalarTighten=0.1_dp   ! scaling factor for the scalar solution
    real(dp)                  :: soilWatbalErr          ! error in the soil water balance
    real(dp)                  :: canopy_max             ! absolute value of the residual in canopy water (kg m-2)
    real(dp),dimension(1)     :: energy_max             ! maximum absolute value of the energy residual (J m-3)
    real(dp),dimension(1)     :: liquid_max             ! maximum absolute value of the volumetric liquid water content residual (-)
    real(dp),dimension(1)     :: matric_max             ! maximum absolute value of the matric head iteration increment (m)
    real(dp)                  :: aquifer_max            ! absolute value of the residual in aquifer water (m)
    logical(lgt)              :: canopyConv             ! flag for canopy water balance convergence
    logical(lgt)              :: watbalConv             ! flag for soil water balance convergence
    logical(lgt)              :: liquidConv             ! flag for residual convergence
    logical(lgt)              :: matricConv             ! flag for matric head convergence
    logical(lgt)              :: energyConv             ! flag for energy convergence
    logical(lgt)              :: aquiferConv            ! flag for aquifer water balance convergence

    real(dp)                        :: absTolWatVeg = 1e-6
    real(dp)                        :: absTolWatSnowSoil = 1e-6
    real(dp)                        :: relTolMatric = 1e-7
    real(dp)                        :: absTolAquifr = 1e-6
    real(dp)                        :: absTolWatBal = 1e-6
    real(dp)                        :: absTolEerngy = 1e-2
  ! -------------------------------------------------------------------------------------------------------------------------------------------------
  ! association to variables in the data structures

  ! -------------------------------------------------------------------------------------------------------------------------------------------------

  ! check convergence based on the canopy water balance
  if(indx_data%var(iLookINDEX%ixVegHyd)%dat(1)/=integerMissing)then
   canopy_max = real(abs(resVec(indx_data%var(iLookINDEX%ixVegHyd)%dat(1))), dp)*iden_water
   canopyConv = (canopy_max    < absTolWatVeg)  ! absolute error in canopy water balance (mm)
  else
   canopy_max = realMissing
   canopyConv = .true.
  endif

  ! check convergence based on the residuals for energy (J m-3)
  if(size(indx_data%var(iLookINDEX%ixNrgOnly)%dat)>0)then
   energy_max = real(maxval(abs( resVec(indx_data%var(iLookINDEX%ixNrgOnly)%dat) )), dp)
   energyConv = (energy_max(1) < absTolEerngy)  ! (based on the residual)
  else
   energy_max = realMissing
   energyConv = .true.
  endif

  ! check convergence based on the residuals for volumetric liquid water content (-)
  if(size(indx_data%var(iLookINDEX%ixHydOnly)%dat)>0)then
   liquid_max = real(maxval(abs( resVec(indx_data%var(iLookINDEX%ixHydOnly)%dat) ) ), dp)
   ! (tighter convergence for the scalar solution)
   if(scalarSolution)then
    liquidConv = (liquid_max(1) < absTolWatSnowSoil * scalarTighten)   ! (based on the residual)
   else
    liquidConv = (liquid_max(1) < absTolWatSnowSoil)                 ! (based on the residual)
   endif
  else
   liquid_max = realMissing
   liquidConv = .true.
  endif

  ! check convergence based on the iteration increment for matric head
  ! NOTE: scale by matric head to avoid unnecessairly tight convergence when there is no water
  if(size(indx_data%var(iLookINDEX%ixMatOnly)%dat)>0)then
   psiScale   = abs( xVec(indx_data%var(iLookINDEX%ixMatOnly)%dat) ) + xSmall ! avoid divide by zero
   matric_max = maxval(abs( xInc(indx_data%var(iLookINDEX%ixMatOnly)%dat)/psiScale ) )
   matricConv = (matric_max(1) < relTolMatric)  ! NOTE: based on iteration increment
  else
   matric_max = realMissing
   matricConv = .true.
  endif
  
  ! check convergence based on the soil water balance error (m)
  if(size(indx_data%var(iLookINDEX%ixMatOnly)%dat)>0)then
   soilWatBalErr = sum( real(resVec(indx_data%var(iLookINDEX%ixMatOnly)%dat), dp)*prog_data%var(iLookPROG%mLayerDepth)%dat(nSnow + indx_data%var(iLookINDEX%ixMatricHead)%dat ) )
   watbalConv    = (abs(soilWatbalErr) < absTolWatBal)  ! absolute error in total soil water balance (m)
  else
   soilWatbalErr = realMissing
   watbalConv    = .true.
  endif

  ! check convergence based on the aquifer storage
  if(indx_data%var(iLookINDEX%ixAqWat)%dat(1)/=integerMissing)then
   aquifer_max = real(abs(resVec(indx_data%var(iLookINDEX%ixAqWat)%dat(1))), dp)*iden_water
   aquiferConv = (aquifer_max    < absTolAquifr)  ! absolute error in aquifer water balance (mm)
  else
   aquifer_max = realMissing
   aquiferConv = .true.
  endif
  
!  print *, 'canopyConv = ', canopyConv
!  print *, 'matricConv = ', matricConv
!  print *, 'liquidConv = ', liquidConv
!  print *, 'energyConv = ', energyConv
!  print *, 'aquiferConv = ', aquiferConv
!  print *, '-----------------------------------'

  isConv = (watbalConv .and. canopyConv .and. matricConv .and. liquidConv .and. energyConv .and. aquiferConv)

  
  end subroutine checkConvFida


end module convTestFida_module

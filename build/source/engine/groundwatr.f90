module groundwatr_module
USE nrtype
! named variables for model decisions
USE mDecisions_module,only:  &                      ! look-up values for method used to compute derivative
 numerical,   & ! numerical solution
 analytical     ! analytical solution
implicit none
! constant parameters
real(dp),parameter     :: valueMissing=-9999._dp    ! missing value parameter
real(dp),parameter     :: verySmall=epsilon(1.0_dp) ! a very small number (used to avoid divide by zero)
real(dp),parameter     :: dx=1.e-8_dp               ! finite difference increment
private
public::q_baseflow
public::groundwatr
contains

 ! ************************************************************************************************
 ! new subroutine: compute water balance for the aquifer
 ! ************************************************************************************************
 subroutine groundwatr(&
                       ! input: model control
                       dt,                   & ! intent(in): time step (s)
                       ixDerivMethod,        & ! intent(in): method used to calculate derivatives
                       ! input: effective parameters
                       aquiferHydCond,       & ! intent(in): effective hydraulic conductivity (m s-1)
                       aquiferScaleFactor,   & ! intent(in): scaling factor for aquifer storage (m)
                       aquiferBaseflowExp,   & ! intent(in): exponent in bucket baseflow parameterization (-)
                       ! input: aquifer fluxes
                       aquiferRecharge,      & ! intent(in): aquifer recharge (m s-1)
                       aquiferTranspire,     & ! intent(in): aquifer transpiration (m s-1)
                       ! input-output
                       aquiferStorage,       & ! intent(inout): aquifer storage (m)
                       ! output
                       aquiferBaseflow,      & ! intent(out): aquifer baseflow (m s-1)
                       err,message)            ! intent(out): error control
 implicit none
 ! -----------------------------------------------------------------------------------------------------------------------------------------
 ! input: model control
 real(dp),intent(in)              :: dt                      ! time step (s)
 integer(i4b),intent(in)          :: ixDerivMethod           ! method used to calculate derivatives
 ! input: effective parameters
 real(dp),intent(in)              :: aquiferHydCond          ! effective hydraulic conductivity (m s-1)
 real(dp),intent(in)              :: aquiferScaleFactor      ! scaling factor for aquifer storage (m)
 real(dp),intent(in)              :: aquiferBaseflowExp      ! exponent in bucket baseflow parameterization (-)
 ! input: aquifer fluxes
 real(dp),intent(in)              :: aquiferRecharge         ! aquifer recharge (m s-1)
 real(dp),intent(in)              :: aquiferTranspire        ! aquifer transpiration (m s-1)
 ! input-output
 real(dp),intent(inout)           :: aquiferStorage          ! aquifer storage (m)
 ! output
 real(dp),intent(out)             :: aquiferBaseflow         ! aquifer baseflow (m s-1)
 integer(i4b),intent(out)         :: err                     ! error code
 character(*),intent(out)         :: message                 ! error message
 ! -----------------------------------------------------------------------------------------------------------------------------------------
 ! local variables
 character(LEN=256)               :: cmessage                ! error message of downwind routine
 integer(i4b)                     :: iter                    ! iteration index
 integer(i4b),parameter           :: maxiter=20              ! maximum number of iterations
 real(dp)                         :: aquiferStorageTrial     ! trial value of aquifer storage
 real(dp)                         :: scalarBaseflow          ! baseflow (m s-1)
 real(dp)                         :: scalarBaseflowDeriv     ! derivative in baseflow w.r.t. aquifer storage (s-1)
 real(dp)                         :: res                     ! residual in water balance (m)
 real(dp)                         :: aquiferIncr             ! iteration increment (m)
 real(dp),parameter               :: tolRes=1.e-8_dp         ! convergence tolerance for the residual
 real(dp),parameter               :: tolInc=1.e-10_dp        ! convergence tolerance for the iteration increment
 ! -----------------------------------------------------------------------------------------------------------------------------------------
 ! initialize error control
 err=0; message='groundwatr/'

 ! initialize state variables
 aquiferStorageTrial   = aquiferStorage ! initialize as start-of-step value

 ! iterate
 do iter=1,maxiter

  ! compute baseflow
  call q_baseflow(&
                  ! input: model decisions
                  .true.,               & ! intent(in): flag indicating if derivatives are desired
                  ixDerivMethod,        & ! intent(in): method used to calculate derivatives
                  ! input: effective parameters
                  aquiferHydCond,       & ! intent(in): effective hydraulic conductivity (m s-1)
                  aquiferScaleFactor,   & ! intent(in): scaling factor for aquifer storage (m)
                  aquiferBaseflowExp,   & ! intent(in): exponent in bucket baseflow parameterization (-)
                  ! input: state variables
                  aquiferStorageTrial,  & ! intent(in): choice of groundwater parameterization
                  ! output
                  scalarBaseflow,       & ! intent(out): total baseflow (m s-1)
                  scalarBaseflowDeriv,  & ! intent(out): derivative in baseflow flux w.r.t. water table depth (m s-1)
                  err,message)            ! intent(out): error control
  if(err/=0)then; message=trim(message)//trim(cmessage); return; endif

  ! compute the residual
  res = (aquiferRecharge + aquiferTranspire - scalarBaseflow)*dt - (aquiferStorageTrial - aquiferStorage)

  ! compute the iteration increment
  aquiferIncr = res/(1._dp + scalarBaseflowDeriv*dt)

  ! print progress
  !print*, 'scalarBaseflowDeriv = ', scalarBaseflowDeriv
  !write(*,'(a,i4,1x,2(f20.10,1x),2(e20.10,1x))') 'iter, aquiferStorageTrial, scalarBaseflow, res, aquiferIncr = ', &
  !                                                iter, aquiferStorageTrial, scalarBaseflow, res, aquiferIncr

  ! update the aquifer
  aquiferStorageTrial = aquiferStorageTrial + aquiferIncr
  !if(iter > 10) pause

  ! check convergence
  if(res < tolRes .or. aquiferIncr < tolInc) exit

  ! check that we converged
  if(iter==maxiter)then; err=20; message=trim(message)//'failed to converge'; return; endif

 end do ! iterating

 ! return aquifer storaage and baseflow
 aquiferStorage  = aquiferStorageTrial
 aquiferBaseflow = scalarBaseflow

 end subroutine groundwatr

 ! ************************************************************************************************
 ! new subroutine: compute baseflow
 ! ************************************************************************************************
 subroutine q_baseflow(&
                       ! input: model decisions
                       deriv_desired,        & ! intent(in): flag indicating if derivatives are desired
                       ixDerivMethod,        & ! intent(in): method used to calculate derivatives
                       ! input: effective parameters
                       aquiferHydCond,       & ! intent(in): effective hydraulic conductivity (m s-1)
                       aquiferScaleFactor,   & ! intent(in): scaling factor for aquifer storage (m)
                       aquiferBaseflowExp,   & ! intent(in): exponent in bucket baseflow parameterization (-)
                       ! input: state variables
                       scalarAquiferStorage, & ! intent(in): aquifer storage (m)
                       ! output
                       scalarBaseflow,       & ! intent(out): total baseflow (m s-1)
                       scalarBaseflowDeriv,  & ! intent(out): derivative in baseflow flux w.r.t. water table depth (m s-1)
                       err,message)            ! intent(out): error control
 implicit none
 ! -----------------------------------------------------------------------------------------------------------------------------------------
 ! input: model decisions
 logical(lgt),intent(in)          :: deriv_desired           ! flag indicating if derivatives are desired
 integer(i4b),intent(in)          :: ixDerivMethod           ! method used to calculate derivatives
 ! input: effective parameters
 real(dp),intent(in)              :: aquiferHydCond          ! effective hydraulic conductivity (m s-1)
 real(dp),intent(in)              :: aquiferScaleFactor      ! scaling factor for aquifer storage (m)
 real(dp),intent(in)              :: aquiferBaseflowExp      ! exponent in bucket baseflow parameterization (-)
 ! input: state variables
 real(dp),intent(in)              :: scalarAquiferStorage    ! trial value of aquifer strorage (m)
 ! output
 real(dp),intent(out)             :: scalarBaseflow          ! baseflow (m s-1)
 real(dp),intent(out)             :: scalarBaseflowDeriv     ! derivative in baseflow flux w.r.t. water table depth (m s-1)
 integer(i4b),intent(out)         :: err                     ! error code
 character(*),intent(out)         :: message                 ! error message
 ! -----------------------------------------------------------------------------------------------------------------------------------------
 ! local variables
 real(dp)                         :: aquiferStorageTrial     ! trial value of aquifer strorage (m)
 real(dp)                         :: scaledStorage           ! scaled storage (-)
 integer(i4b)                     :: itry                    ! index of different flux calculations
 integer(i4b)                     :: nFlux                   ! number of flux calculations required (>1 = numerical derivatives)
 integer(i4b),parameter           :: unperturbed=0           ! named variable to identify the case of unperturbed state variables
 integer(i4b),parameter           :: perturbState=1          ! named variable to identify the case where we perturb the state in the current layer
 integer(i4b),parameter           :: perturbStateAbove=2     ! named variable to identify the case where we perturb the state layer above
 integer(i4b),parameter           :: perturbStateBelow=3     ! named variable to identify the case where we perturb the state layer below
 real(dp)                         :: scalarFlux              ! baseflow flux (m s-1)
 real(dp)                         :: scalarFlux_dStateAbove  ! baseflow flux with perturbation to the state above (m s-1)
 ! -----------------------------------------------------------------------------------------------------------------------------------------
 ! initialize error control
 err=0; message='q_baseflow/'
 
 ! identify number of ADDITIONAL flux evaluations (used when computing numerical derivatives)
 if(ixDerivMethod==numerical .and. deriv_desired)then
  nFlux=3   ! NOTE: we cycle through undesired perturbations, so only actually do one additional flux eval
 else
  nFlux=0
 endif

 ! *** loop to compute numerical derivatives
 ! either one or multiple flux calls, depending on if using analytical or numerical derivatives
 do itry=nFlux,0,-1  ! (work backwards to ensure all computed fluxes come from the un-perturbed case)

  ! skip undesired perturbations
  if(itry==perturbState .or. itry==perturbStateBelow) cycle

  ! define the trial value for aquifer storage
  select case(itry)
   case(unperturbed);       aquiferStorageTrial  = scalarAquiferStorage
   case(perturbStateAbove); aquiferStorageTrial  = scalarAquiferStorage + dx
    case(perturbStateBelow,perturbState); err=10; message=trim(message)//'only perturb aquifer storage when computing baseflow flux -- should not get here'; return
   case default; err=10; message=trim(message)//"unknown perturbation"; return
  end select ! (type of perturbation)

  ! compute scaled storage (-)
  scaledStorage  = aquiferStorageTrial/aquiferScaleFactor

  ! compute baseflow
  scalarBaseflow = aquiferHydCond*(scaledStorage**aquiferBaseflowExp)

  ! compute derivative in baseflow
  if(ixDerivMethod==analytical .and. deriv_desired)then
   scalarBaseflowDeriv = (aquiferHydCond/aquiferScaleFactor)*aquiferBaseflowExp*scaledStorage**(aquiferBaseflowExp - 1._dp)
  else
   scalarBaseflowDeriv = valueMissing
  endif

  ! get copies of baseflow flux to compute derivatives
  if(deriv_desired .and. ixDerivMethod==numerical)then
   select case(itry)
    case(unperturbed);       scalarFlux             = scalarBaseflow
    case(perturbStateAbove); scalarFlux_dStateAbove = scalarBaseflow
    case(perturbStateBelow,perturbState); err=10; message=trim(message)//'only perturb aquifer storage when computing baseflow flux -- should not get here'; return
    case default; err=10; message=trim(message)//'unknown perturbation'; return
   end select
  endif

 end do  ! (multiple flux calls for computing  numerical derivatives

 ! * compute derivatives
 ! NOTE: baseflow derivatives w.r.t. state below are *actually* w.r.t. water table depth, so need to be corrected for aquifer storage
 if(deriv_desired)then
  if(ixDerivMethod==numerical) scalarBaseflowDeriv = (scalarFlux_dStateAbove - scalarFlux)/dx
 else
  scalarBaseflowDeriv = valueMissing
 endif

 end subroutine q_baseflow

end module groundwatr_module

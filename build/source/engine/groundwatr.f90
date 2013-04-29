module groundwatr_module
USE nrtype
! provide access to look-up values for model decisions
USE mDecisions_module,only:  &
 ! look-up values for method used to compute derivative
 numerical,                  & ! numerical solution
 analytical,                 & ! analytical solution
 ! look-up values for the type of hydraulic conductivity profile
 constant,                   & ! constant hydraulic conductivity with depth
 exp_profile,                & ! exponential profile
 powerLaw_profile,           & ! power-law profile
 linear_profile,             & ! linear profile
 ! look-up values for the choice of groundwater parameterization
 equilWaterTable,            & ! equilibrium water table
 pseudoWaterTable,           & ! pseudo water table
 bigBucket,                  & ! a big bucket (lumped aquifer model)
 noExplicit                    ! no explicit groundwater parameterization
implicit none
! constant parameters
real(dp),parameter     :: valueMissing=-9999._dp    ! missing value parameter
real(dp),parameter     :: verySmall=epsilon(1.0_dp) ! a very small number (used to avoid divide by zero)
real(dp),parameter     :: dx=1.e-8_dp               ! finite difference increment
private
public::q_baseflow
public::basinAqfr
contains

 ! ************************************************************************************************
 ! new subroutine: compute water balance for the aquifer
 ! ************************************************************************************************
 subroutine basinAqfr(&
                      ! input
                      dt,              &  ! time step (s)
                      aquiferRecharge, &  ! aquifer recharge (m s-1)
                      aquiferTranspire,&  ! aquifer transpiration (m s-1)
                      ! input-output
                      aquiferStorage,  &  ! aquifer storage (m)
                      ! output
                      aquiferBaseflow, &  ! aquifer baseflow (m s-1)
                      err,message)        ! error control
 ! model decision structures
 USE data_struc,only:model_decisions      ! model decision structure
 USE var_lookup,only:iLookDECISIONS       ! named variables for elements of the decision structure
 ! model variables, parameters, forcing data, etc.
 USE data_struc,only:bpar_data            ! basin pameters
 USE var_lookup,only:iLookBPAR            ! named variables denoting position of parameter in data structure
 implicit none
 ! -----------------------------------------------------------------------------------------------------------------------------------------
 ! input
 real(dp),intent(in)              :: dt                      ! time step (s)
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
 real(dp)                         :: scalarWaterTableDepth   ! water table depth (not used)
 real(dp)                         :: scalarBaseflow          ! baseflow (m s-1)
 real(dp)                         :: scalarBaseflowDeriv     ! derivative in baseflow w.r.t. aquifer storage (s-1)
 real(dp)                         :: res                     ! residual in water balance (m)
 real(dp)                         :: aquiferIncr             ! iteration increment (m)
 real(dp),parameter               :: tolRes=1.e-8_dp         ! convergence tolerance for the residual
 real(dp),parameter               :: tolInc=1.e-10_dp        ! convergence tolerance for the iteration increment
 ! -----------------------------------------------------------------------------------------------------------------------------------------
 ! initialize error control
 err=0; message='basinAqfr/'

 ! basic checks -- only implement basin aquifer for the big bucket
 if(model_decisions(iLookDECISIONS%groundwatr)%iDecision /= bigBucket)then
  message=trim(message)//'the only supported option for basin aquifer is the big bucket'
  err=20; return
 endif

 ! initialize state variables
 aquiferStorageTrial   = aquiferStorage ! initialize as start-of-step value
 scalarWaterTableDepth = valueMissing   ! water table depth should not be used here, so set it to missing in order to cause problems

 ! iterate
 do iter=1,maxiter

  ! compute baseflow
  call q_baseflow(&
                  ! input: model decisions
                  .true.,                                               & ! intent(in): flag indicating if derivatives are desired
                  model_decisions(iLookDECISIONS%fDerivMeth)%iDecision, & ! intent(in): method used to calculate flux derivatives 
                  model_decisions(iLookDECISIONS%groundwatr)%iDecision, & ! intent(in): groundwater parameterization
                  model_decisions(iLookDECISIONS%hc_Profile)%iDecision, & ! intent(in): option for the hydraulic conductivity profile
                  ! input: state and diagnostic variables
                  scalarWaterTableDepth,                                & ! intent(in): water table depth (m)
                  aquiferStorageTrial,                                  & ! intent(in): choice of groundwater parameterization
                  bpar_data%var(iLookBPAR%basin__hydCond),              & ! intent(in): saturated hydraulic conductivity at the surface (m s-1)
                  ! input: parameters
                  bpar_data%var(iLookBPAR%basin__kAnisotropic),         & ! intent(in): anisotropy factor for lateral hydraulic conductivity (-) 
                  bpar_data%var(iLookBPAR%basin__zScale_TOPMODEL),      & ! intent(in): scale factor for TOPMODEL-ish baseflow parameterization (m)
                  bpar_data%var(iLookBPAR%basin__aquiferScaleFactor),   & ! intent(in): scaling factor for aquifer storage in the big bucket (m)
                  bpar_data%var(iLookBPAR%basin__bucketBaseflowExp),    & ! intent(in): exponent in bucket baseflow parameterization (-)
                  ! output
                  scalarBaseflow,                                       & ! intent(out): total baseflow (m s-1)
                  scalarBaseflowDeriv,                                  & ! intent(out): derivative in baseflow flux w.r.t. water table depth (m s-1)
                  err,message)                                            ! intent(out): error control
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

 end subroutine basinAqfr

 ! ************************************************************************************************
 ! new subroutine: compute baseflow
 ! ************************************************************************************************
 subroutine q_baseflow(&
                       ! input: model decisions
                       deriv_desired,               & ! intent(in): flag indicating if derivatives are desired
                       ixDerivMethod,               & ! intent(in): choice of method used to compute derivative
                       ixGroundwater,               & ! intent(in): choice of groundwater parameterization
                       ixHydcondProfile,            & ! intent(in): choice of hydraulic conductivity profile
                       ! input: state and diagnostic variables
                       scalarWaterTableDepth,       & ! intent(in): water table depth (m)
                       scalarAquiferStorage,        & ! intent(in): aquifer storage (m)
                       surfaceSatHydCond,           & ! intent(in): saturated hydraulic conductivity at the surface (m s-1)
                       ! input: parameters
                       kAnisotropic,                & ! intent(in): anisotropy factor for lateral hydraulic conductivity (-)
                       zScale_TOPMODEL,             & ! intent(in): scale factor for TOPMODEL-ish baseflow parameterization (m)
                       aquiferScaleFactor,          & ! intent(in): scaling factor for aquifer storage in the big bucket (m)
                       bucketBaseflowExp,           & ! intent(in): exponent in bucket baseflow parameterization
                       ! output
                       scalarBaseflow,              & ! intent(out): total baseflow (m s-1)
                       scalarBaseflowDeriv,         & ! intent(out): derivative in baseflow flux w.r.t. water table depth (m s-1)
                       err,message)                   ! intent(out): error control
 implicit none
 ! -----------------------------------------------------------------------------------------------------------------------------------------
 ! input: model decisions
 logical(lgt),intent(in)          :: deriv_desired                 ! flag indicating if derivatives are desired
 integer(i4b),intent(in)          :: ixDerivMethod                 ! index defining option for computing derivatives (analytical or numerical)
 integer(i4b),intent(in)          :: ixGroundwater                 ! index defining choice of groundwater parameterization
 integer(i4b),intent(in)          :: ixHydcondProfile              ! index defining hydraulic conductivity profile
 ! input: state and diagnostic variables
 real(dp),intent(in)              :: scalarWaterTableDepth         ! trial value of water table depth (m)
 real(dp),intent(in)              :: scalarAquiferStorage          ! trial value of aquifer strorage (m)
 real(dp),intent(in)              :: surfaceSatHydCond             ! saturated hydraulic conductivity at the surface (m s-1) 
 ! input: parameters
 real(dp),intent(in)              :: kAnisotropic                  ! anisotropy factor for lateral hydraulic conductivity (-)
 real(dp),intent(in)              :: zScale_TOPMODEL               ! scale factor for TOPMODEL-ish baseflow parameterization (m)
 real(dp),intent(in)              :: aquiferScaleFactor            ! scaling factor for aquifer storage in the big bucket (m)
 real(dp),intent(in)              :: bucketBaseflowExp             ! exponent in bucket baseflow parameterization
 ! output
 real(dp),intent(out)             :: scalarBaseflow                ! baseflow (m s-1)
 real(dp),intent(out)             :: scalarBaseflowDeriv           ! derivative in baseflow flux w.r.t. water table depth (m s-1)
 integer(i4b),intent(out)         :: err                           ! error code
 character(*),intent(out)         :: message                       ! error message
 ! -----------------------------------------------------------------------------------------------------------------------------------------
 ! local variables
 character(LEN=256)               :: cmessage                      ! error message of downwind routine
 real(dp)                         :: scalarWaterTableDepthTrial    ! trial value of water table depth (m)
 real(dp)                         :: scalarAquiferStorageTrial     ! trial value of aquifer strorage (m)
 integer(i4b)                     :: itry                          ! index of different flux calculations
 integer(i4b)                     :: nFlux                         ! number of flux calculations required (>1 = numerical derivatives)
 integer(i4b),parameter           :: unperturbed=0                 ! named variable to identify the case of unperturbed state variables
 integer(i4b),parameter           :: perturbState=1                ! named variable to identify the case where we perturb the state in the current layer
 integer(i4b),parameter           :: perturbStateAbove=2           ! named variable to identify the case where we perturb the state layer above
 integer(i4b),parameter           :: perturbStateBelow=3           ! named variable to identify the case where we perturb the state layer below
 real(dp)                         :: scalarFlux                    ! baseflow flux (m s-1)
 real(dp)                         :: scalarFlux_dStateAbove        ! baseflow flux with perturbation to the state above (m s-1)
 ! -----------------------------------------------------------------------------------------------------------------------------------------
 ! initialize error control
 err=0; message='q_baseflow/'
 
 ! identify number of ADDITIONAL flux evaluations
 if(ixDerivMethod==numerical .and. deriv_desired)then
  nFlux=3   ! NOTE: cycle through undesired perturbations, so only actually do one additional flux eval
 else
  nFlux=0
 endif

 ! ------------------------------------------------------------------------------------------------------------------------------------------------
 ! *****
 ! case of no explicit deep groundwater
 if(ixGroundwater == noExplicit)then
  scalarBaseflow      = 0._dp  ! baseflow from the aquifer (m s-1)
  scalarBaseflowDeriv = 0._dp  ! derivative in baseflow w.r.t. aquifer storage (s-1)
  return
 endif

 ! *** loop to compute numerical derivatives
 ! either one or multiple flux calls, depending on if using analytical or numerical derivatives
 do itry=nFlux,0,-1  ! (work backwards to ensure all computed fluxes come from the un-perturbed case)

  ! only perturb states for the pseudoWaterTable
  if(ixGroundwater == equilWaterTable)then
   if(itry /= unperturbed) cycle
  endif

  ! * identify the type of perturbation
  select case(itry)
   ! (skip undersired perturbations)
   case(perturbStateBelow); cycle   ! assume baseflow at "bottom" of aquifer, so perturb w.r.t. state above
   case(perturbState); cycle        ! assume baseflow at "bottom" of aquifer, so perturb w.r.t. state above
   ! (unperturbed)
   case(unperturbed)
    select case(ixGroundwater)
     case(pseudoWaterTable,equilWaterTable); scalarWaterTableDepthTrial = scalarWaterTableDepth
     case(bigBucket);                        scalarAquiferStorageTrial  = scalarAquiferStorage
     case default; err=20; message=trim(message)//'unable to identify choice of groundwater representation'; return
    end select
   ! (perturbed)
   case(perturbStateAbove)
    select case(ixGroundwater)
     case(pseudoWaterTable,equilWaterTable); scalarWaterTableDepthTrial = scalarWaterTableDepth + dx
     case(bigBucket);                        scalarAquiferStorageTrial  = scalarAquiferStorage + dx
     case default; err=20; message=trim(message)//'unable to identify choice of groundwater representation'; return
    end select
   case default; err=10; message=trim(message)//"unknown perturbation"; return
  end select ! (type of perturbation)

  ! select groundwater option
  select case(ixGroundwater)

   ! * compute baseflow using topmodel methods
   case(pseudoWaterTable,equilWaterTable)
    call QBtopmodel(&
                    (ixDerivMethod==analytical .and. deriv_desired), & ! input: flag indicating if derivatives are desired
                    ixHydcondProfile,            & ! input: index defining hydraulic conductivity profile
                    scalarWaterTableDepthTrial,  & ! input: water table depth (m)
                    surfaceSatHydCond,           & ! input: saturated hydraulic conductivity at the surface (m s-1)
                    kAnisotropic,                & ! input: anisotropy factor for lateral hydraulic conductivity (-)
                    zScale_TOPMODEL,             & ! input: scale factor for TOPMODEL-ish baseflow parameterization (m)
                    scalarBaseflow,              & ! output: total baseflow (m s-1)
                    scalarBaseflowDeriv,         & ! output: derivative in baseflow flux w.r.t. water table depth (m s-1)
                    err,cmessage)                  ! output: error control
    if(err/=0)then; message=trim(message)//trim(cmessage); return; endif

   ! * compute baseflow using the conceptual big-bucket
   case(bigBucket)
    call QBbigbuckt(&
                    (ixDerivMethod==analytical .and. deriv_desired), & ! input: flag indicating if derivatives are desired
                    scalarAquiferStorageTrial,   & ! input: trial value of aquifer storage (m)
                    aquiferScaleFactor,          & ! input: scaling factor for aquifer storage in the big bucket (m)
                    surfaceSatHydCond,           & ! input: saturated hydraulic conductivity at the surface (m s-1)
                    kAnisotropic,                & ! input: anisotropy factor for lateral hydraulic conductivity (-)
                    bucketBaseflowExp,           & ! input: exponent in bucket baseflow parameterization
                    scalarBaseflow,              & ! output: total baseflow (m s-1)
                    scalarBaseflowDeriv,         & ! output: derivative in baseflow flux w.r.t. aquifer storage (m s-1)
                    err,cmessage)                  ! output: error control
    if(err/=0)then; message=trim(message)//trim(cmessage); return; endif
   
   ! * check that we identified the groundwater option
   case default
    message=trim(message)//'unable to identify choice of groundwater representation'
    err=20; return

  end select  ! choice of groundwater option

  ! * get copies of baseflow flux to compute derivatives
  if(deriv_desired .and. ixDerivMethod==numerical)then
   select case(itry)
    case(unperturbed);       scalarFlux             = scalarBaseflow
    case(perturbStateAbove); scalarFlux_dStateAbove = scalarBaseflow
    case(perturbStateBelow,perturbState); err=10; message=trim(message)//'only perturb water table depth when computing baseflow flux -- should not get here'; return
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





 ! ************************************************************************************************
 ! ************************************************************************************************
 ! ************************************************************************************************
 ! ************************************************************************************************
 ! ***** PRIVATE SUBROUTINES **********************************************************************
 ! ************************************************************************************************
 ! ************************************************************************************************
 ! ************************************************************************************************

 ! ************************************************************************************************
 ! private subroutine: compute baseflow flux using a topmodel-type approach
 ! ************************************************************************************************
 subroutine QBtopmodel(&
                       deriv_desired,               & ! input: flag indicating if derivatives are desired
                       ixHydcondProfile,            & ! input: index defining hydraulic conductivity profile
                       scalarWaterTableDepth,       & ! input: water table depth (m)
                       k_surf,                      & ! input: saturated hydraulic conductivity at the surface (m s-1)
                       kAnisotropic,                & ! input: anisotropy factor for lateral hydraulic conductivity (-)
                       zScale_TOPMODEL,             & ! input: scale factor for TOPMODEL-ish baseflow parameterization (m)
                       scalarBaseflow,              & ! output: baseflow from each soil layer (m s-1)
                       scalarBaseflowDeriv,         & ! output: derivative in baseflow flux w.r.t. water table depth (m s-1)
                       err,message)                   ! output: error control
 implicit none
 ! -----------------------------------------------------------------------------------------------------------------------------------------
 ! input
 logical(lgt),intent(in)   :: deriv_desired           ! flag to indicate if derivatives are desired
 integer(i4b),intent(in)   :: ixHydcondProfile        ! index defining hydraulic conductivity profile
 real(dp),intent(in)       :: scalarWaterTableDepth   ! trial value of water table depth (m)
 real(dp),intent(in)       :: k_surf                  ! saturated hydraulic conductivity at the surface (m s-1) 
 real(dp),intent(in)       :: kAnisotropic            ! anisotropy factor for lateral hydraulic conductivity (-)
 real(dp),intent(in)       :: zScale_TOPMODEL         ! scale factor for TOPMODEL-ish baseflow parameterization (m)
 ! output
 real(dp),intent(out)      :: scalarBaseflow          ! baseflow (m s-1)
 real(dp),intent(out)      :: scalarBaseflowDeriv     ! derivative in baseflow flux w.r.t. water table depth (m s-1)
 integer(i4b),intent(out)  :: err                     ! error code
 character(*),intent(out)  :: message                 ! error message
 ! -----------------------------------------------------------------------------------------------------------------------------------------
 ! initialize error control
 err=0; message='QBtopmodel/'
 ! compute the baseflow (m s-1)
 select case(ixHydcondProfile)
  ! (exponential transmissivity profile)
  case(exp_profile)
   scalarBaseflow = k_surf*kAnisotropic*exp(-scalarWaterTableDepth/zScale_TOPMODEL)
   if(deriv_desired)then
    scalarBaseflowDeriv = -scalarBaseflow/zScale_TOPMODEL
   else
    scalarBaseflowDeriv = valueMissing
   endif
  ! (linear transmissivity profile)
  case(linear_profile)
   scalarBaseflow = k_surf*kAnisotropic*(1._dp - scalarWaterTableDepth/zScale_TOPMODEL)
   if(deriv_desired)then
    scalarBaseflowDeriv = -k_surf*kAnisotropic/zScale_TOPMODEL
   else
    scalarBaseflowDeriv = valueMissing
   endif
  ! (constant transmissivity profile)
  case(constant)
   scalarBaseflow = k_surf*kAnisotropic*exp(-scalarWaterTableDepth/zScale_TOPMODEL)
   if(deriv_desired)then
    scalarBaseflowDeriv = -scalarBaseflow/zScale_TOPMODEL
   else
    scalarBaseflowDeriv = valueMissing
   endif
  ! (power-law transmissivity profile)
  case(powerLaw_profile)
   message=trim(message)//"power-law hydraulic conductivity profile not implemented yet"
   err=20; return
  ! (unknown transmissivity profile)
  case default
   message=trim(message)//"unknown hydraulic conductivity profile"
   err=20; return
 end select
 end subroutine QBtopmodel


 ! ************************************************************************************************
 ! private subroutine: compute baseflow flux using a conceptual bog bucket
 ! ************************************************************************************************
 subroutine QBbigbuckt(&
                       deriv_desired,               &    ! input: flag indicating if derivatives are desired
                       scalarAquiferStorageTrial,   &    ! input: trial value of aquifer storage (m)
                       aquiferScaleFactor,          &    ! input: scaling factor for aquifer storage in the big bucket (m)
                       k_surf,                      &    ! input: saturated hydraulic conductivity at the surface (m s-1)
                       kAnisotropic,                &    ! input: anisotropy factor for lateral hydraulic conductivity (-)
                       bucketBaseflowExp,           &    ! input: exponent in bucket baseflow parameterization
                       scalarAquiferBaseflow,       &    ! output: total baseflow (m s-1)
                       scalarAquiferBaseflowDeriv,  &    ! output: derivative in baseflow flux w.r.t. aquifer storage (s-1)
                       err,message)                      ! output: error control
 implicit none
 ! -----------------------------------------------------------------------------------------------------------------------------------------
 ! input
 logical(lgt),intent(in)   :: deriv_desired              ! flag to indicate if derivatives are desired
 real(dp),intent(in)       :: scalarAquiferStorageTrial  ! trial value of water table depth (m)
 real(dp),intent(in)       :: aquiferScaleFactor         ! scaling factor for aquifer storage in the big bucket (m)
 real(dp),intent(in)       :: k_surf                     ! saturated hydraulic conductivity at the surface (m s-1) 
 real(dp),intent(in)       :: kAnisotropic               ! anisotropy factor for lateral hydraulic conductivity (-)
 real(dp),intent(in)       :: bucketBaseflowExp          ! exponent in bucket baseflow parameterization
 ! output
 real(dp),intent(out)      :: scalarAquiferBaseflow      ! baseflow flux (m s-1)
 real(dp),intent(out)      :: scalarAquiferBaseflowDeriv ! derivative in baseflow flux w.r.t. water table depth (s-1)
 integer(i4b),intent(out)  :: err                        ! error code
 character(*),intent(out)  :: message                    ! error message
 ! local
 real(dp)                  :: scaledStorage              ! scaled storage (-)
 real(dp)                  :: maxBaseflowRate            ! maximum baseflow rate (m s-1)
 ! -----------------------------------------------------------------------------------------------------------------------------------------
 ! initialize error control
 err=0; message='QBbigbuckt/'
 ! get temporary variables
 scaledStorage    = scalarAquiferStorageTrial/aquiferScaleFactor
 maxBaseflowRate  = kAnisotropic*k_surf
 ! compute baseflow flux (m s-1)
 scalarAquiferBaseflow = maxBaseflowRate*scaledStorage**bucketBaseflowExp
 ! compute derivative in baseflow flux w.r.t. aquifer storage (s-1)
 if(deriv_desired)then
  scalarAquiferBaseflowDeriv = (maxBaseflowRate/aquiferScaleFactor)*bucketBaseflowExp*scaledStorage**(bucketBaseflowExp - 1._dp)
 else
  scalarAquiferBaseflowDeriv = valueMissing
 endif
 end subroutine QBbigbuckt







end module groundwatr_module

module canopySnow_module
USE nrtype
implicit none
private
public::canopySnow

contains

 ! ************************************************************************************************
 ! new subroutine: compute change in snow stored on the vegetation canopy
 ! ************************************************************************************************
 subroutine canopySnow(&
                       dt,                          & ! intent(in): time step (seconds)
                       computeVegFlux,              & ! intent(in): flag to denote if computing energy flux over vegetation
                       err,message)                   ! intent(out): error control
 ! ------------------------------------------------------------------------------------------------
 ! model variables, parameters, forcing data, etc.
 USE data_struc,only:attr_data,type_data,mpar_data,forc_data,mvar_data,indx_data    ! data structures
 USE var_lookup,only:iLookATTR,iLookTYPE,iLookPARAM,iLookFORCE,iLookMVAR,iLookINDEX ! named variables for structure elements
 ! model decision structures
 USE data_struc,only:model_decisions    ! model decision structure
 USE var_lookup,only:iLookDECISIONS     ! named variables for elements of the decision structure
 implicit none
 ! dummy variables
 real(dp),intent(in)           :: dt                  ! time step (seconds)
 logical(lgt),intent(in)       :: computeVegFlux      ! flag to indicate if we are computing fluxes over vegetation (.false. means veg is buried with snow)
 integer(i4b),intent(out)      :: err                 ! error code
 character(*),intent(out)      :: message             ! error message
 ! local variables
 character(LEN=256)            :: cmessage            ! error message of downwind routine
 real(dp)                      :: exposedVAI          ! exposed vegetation area index (-)
 ! initialize error control
 err=0; message='canopySnow/'

 ! compute exposed vegetation area index
 exposedVAI = mvar_data%var(iLookMVAR%scalarExposedLAI)%dat(1) + mvar_data%var(iLookMVAR%scalarExposedSAI)%dat(1)

 ! compute change in snow stored on the vegetation canopy
 call canopySnow_muster(&
                        ! input: model control
                        dt,                                                          & ! intent(in): time step (seconds)
                        computeVegFlux,                                              & ! intent(in): flag to indicate if we are computing fluxes over vegetation (.false. means veg is buried with snow)
                        model_decisions(iLookDECISIONS%snowIncept)%iDecision,        & ! intent(in): choice of option to determine maximum snow interception capacity
                        ! input-output: state variables
                        mvar_data%var(iLookMVAR%scalarCanopyIce)%dat(1),             & ! intent(inout): storage of ice on the vegetation canopy (kg m-2)
                        ! input: diagnostic variables
                        exposedVAI,                                                  & ! intent(in): exposed vegetation area index (m2 m-2)
                        forc_data%var(iLookFORCE%airtemp),                           & ! intent(in): air temperature (K)
                        mvar_data%var(iLookMVAR%scalarSnowfall)%dat(1),              & ! intent(in): computed snowfall rate (kg m-2 s-1)
                        mvar_data%var(iLookMVAR%scalarNewSnowDensity)%dat(1),        & ! intent(in): density of new snow (kg m-3)
                        mvar_data%var(iLookMVAR%scalarCanopyLiqDrainage)%dat(1),     & ! intent(in): liquid drainage from the vegetation canopy (kg m-2 s-1)
                        ! input: parameters
                        mpar_data%var(iLookPARAM%refInterceptCapSnow),               & ! intent(in): reference canopy interception capacity for snow per unit leaf area (kg m-2)
                        mpar_data%var(iLookPARAM%ratioDrip2Unloading),               & ! intent(in): ratio of canopy drip to snow unloading (-)
                        mpar_data%var(iLookPARAM%snowUnloadingCoeff),                & ! intent(in): time constant for unloading of snow from the forest canopy (s-1)
                        ! output: diagnostic variables
                        mvar_data%var(iLookMVAR%scalarThroughfallSnow)%dat(1),       & ! intent(out): snow that reaches the ground without ever touching the canopy (kg m-2 s-1)
                        mvar_data%var(iLookMVAR%scalarCanopySnowUnloading)%dat(1),   & ! intent(out): unloading of snow from the vegetion canopy (kg m-2 s-1)
                        ! output: error control
                        err,cmessage                                                 ) ! intent(out): error control
 if(err/=0)then; err=20; message=trim(message)//trim(cmessage); return; endif

 end subroutine canopySnow



 ! *****************************************************************************************************************************************************************
 ! *****************************************************************************************************************************************************************
 ! *****************************************************************************************************************************************************************
 ! ***** PRIVATE SUBROUTINES ***************************************************************************************************************************************
 ! *****************************************************************************************************************************************************************
 ! *****************************************************************************************************************************************************************
 ! *****************************************************************************************************************************************************************

 ! ************************************************************************************************
 ! new subroutine: compute change in snow stored on the vegetation canopy
 ! ************************************************************************************************
 subroutine canopySnow_muster(&
                              ! input: model control
                              dt,                          & ! intent(in): time step (seconds)
                              computeVegFlux,              & ! intent(in): flag to denote if computing energy flux over vegetation
                              ixSnowInterception,          & ! intent(in): choice of option to determine maximum snow interception capacity
                              ! input-output: state variables
                              scalarCanopyIce,             & ! intent(inout): storage of ice on the vegetation canopy (kg m-2)
                              ! input: forcing and diagnostic variables
                              exposedVAI,                  & ! intent(in): exposed vegetation area index (m2 m-2)
                              scalarAirtemp,               & ! intent(in): air temperature (K)
                              scalarSnowfall,              & ! intent(in): computed snowfall rate (kg m-2 s-1)
                              scalarNewSnowDensity,        & ! intent(in): density of new snow (kg m-3)
                              scalarCanopyLiqDrainage,     & ! intent(in): liquid drainage from the vegetation canopy (kg m-2 s-1)
                              ! input: parameters
                              refInterceptCapSnow,         & ! intent(in): reference canopy interception capacity for snow per unit leaf area (kg m-2)
                              ratioDrip2Unloading,         & ! intent(in): ratio of canopy drip to snow unloading (-)
                              snowUnloadingCoeff,          & ! intent(in): time constant for unloading of snow from the forest canopy (s-1)
                              ! output: diagnostic variables
                              scalarThroughfallSnow,       & ! intent(out): snow that reaches the ground without ever touching the canopy (kg m-2 s-1)
                              scalarCanopySnowUnloading,   & ! intent(out): unloading of snow from the vegetion canopy (kg m-2 s-1)
                              ! output: error control
                              err,message)                   ! intent(out): error control
 ! -----------------------------------------------------------------------------------------------------------------------------------------------------
 ! physical constants
 USE multiconst,only:Tfreeze         ! freezing point of pure water (K)
 ! model decisions
 USE mDecisions_module,only:       &
                       stickySnow, & ! maximum interception capacity an increasing function of temerature
                       lightSnow     ! maximum interception capacity an inverse function of new snow densit
 implicit none
 ! input: control
 real(dp),intent(in)           :: dt                         ! time step (seconds)
 logical(lgt),intent(in)       :: computeVegFlux             ! flag to indicate if we are computing fluxes over vegetation (.false. means veg is buried with snow)
 integer(i4b),intent(in)       :: ixSnowInterception         ! choice of option to determine maximum snow interception capacity
 ! input-output: state variables
 real(dp),intent(inout)        :: scalarCanopyIce            ! mass of ice on the vegetation canopy (kg m-2)
 ! input: diagnostic variables
 real(dp),intent(in)           :: exposedVAI                 ! exposed vegetation area index -- leaf + stem -- after burial by snow (m2 m-2)
 real(dp),intent(in)           :: scalarAirtemp              ! air temperature (K)
 real(dp),intent(in)           :: scalarSnowfall             ! computed snowfall rate (kg m-2 s-1)
 real(dp),intent(in)           :: scalarNewSnowDensity       ! density of new snow (kg m-3)
 real(dp),intent(in)           :: scalarCanopyLiqDrainage    ! liquid drainage from the vegetation canopy (kg m-2 s-1)
 ! input: parameters
 real(dp),intent(in)           :: refInterceptCapSnow        ! reference canopy interception capacity for snow per unit leaf area (kg m-2)
 real(dp),intent(in)           :: ratioDrip2Unloading        ! ratio of canopy drip to snow unloading (-)
 real(dp),intent(in)           :: snowUnloadingCoeff         ! time constant for unloading of snow from the forest canopy (s-1)
 ! output: diagnostic variables
 real(dp),intent(out)          :: scalarThroughfallSnow      ! snow that reaches the ground without ever touching the canopy (kg m-2 s-1)
 real(dp),intent(out)          :: scalarCanopySnowUnloading  ! unloading of snow from the vegetion canopy (kg m-2 s-1)
 ! output: error control
 integer(i4b),intent(out)      :: err                        ! error code
 character(*),intent(out)      :: message                    ! error message
 ! -------------------------------------------------------------------------------------------------------------------------------
 ! local variables
 real(dp),parameter            :: valueMissing=-9999._dp     ! missing value
 integer(i4b)                  :: iter                       ! iteration index
 integer(i4b),parameter        :: maxiter=50                 ! maximum number of iterations
 integer(i4b)                  :: itry                       ! index of loop used for testing
 real(dp)                      :: unloading_melt             ! unloading associated with canopy drip (kg m-2 s-1)
 real(dp)                      :: airtemp_degC               ! value of air temperature in degrees Celcius
 real(dp)                      :: leafScaleFactor            ! scaling factor for interception based on temperature (-)
 real(dp)                      :: leafInterceptCapSnow       ! storage capacity for snow per unit leaf area (kg m-2)
 real(dp)                      :: canopyIceScaleFactor       ! capacity scaling factor for throughfall (kg m-2)
 real(dp)                      :: throughfallDeriv           ! derivative in throughfall flux w.r.t. canopy storage (s-1)
 real(dp)                      :: unloadingDeriv             ! derivative in unloading flux w.r.t. canopy storage (s-1)
 real(dp)                      :: scalarCanopyIceIter        ! trial value for mass of ice on the vegetation canopy (kg m-2) (kg m-2)
 real(dp)                      :: flux                       ! net flux (kg m-2 s-1)
 real(dp)                      :: delS                       ! change in storage (kg m-2)
 real(dp)                      :: resMass                    ! residual in mass equation (kg m-2) 
 real(dp),parameter            :: convTolerMass=0.0001_dp    ! convergence tolerance for mass (kg m-2)
 ! -------------------------------------------------------------------------------------------------------------------------------
 ! initialize error control
 err=0; message='canopySnow_muster/'

 ! *****
 ! compute unloading due to melt drip...
 ! *************************************

 if(computeVegFlux)then
  unloading_melt = min(ratioDrip2Unloading*scalarCanopyLiqDrainage, scalarCanopyIce/dt)  ! kg m-2 s-1
 else
  unloading_melt = 0._dp
 endif
 scalarCanopyIce = scalarCanopyIce - unloading_melt

 ! *****
 ! compute the ice balance due to snowfall and unloading...
 ! ********************************************************

 ! check for early returns
 if(.not.computeVegFlux .or. (scalarSnowfall<tiny(dt) .and. scalarCanopyIce<tiny(dt)))then
  scalarThroughfallSnow     = scalarSnowfall    ! throughfall of snow through the canopy (kg m-2 s-1)
  scalarCanopySnowUnloading = unloading_melt    ! unloading of snow from the canopy (kg m-2 s-1)
  return
 endif

 ! get a trial value for canopy storage
 scalarCanopyIceIter = scalarCanopyIce

 ! iterate
 do iter=1,maxiter

  ! ** compute throughfall

  ! no snowfall
  if(scalarSnowfall<tiny(dt))then ! no snow
   ! compute throughfall -- note this is effectively zero (no snow case)
   scalarThroughfallSnow = scalarSnowfall  ! throughfall (kg m-2 s-1)
   throughfallDeriv      = 0._dp

  ! snowfall: compute interception
  else
 
   ! ** process different options for maximum branch snow interception
   select case(ixSnowInterception)

    ! * option 1: maximum interception capacity an inverse function of new snow density (e.g., Mahat and Tarboton, HydProc 2013)
    case(lightSnow)  
     ! (check new snow density is valid)
     if(scalarNewSnowDensity < 0._dp)then; err=20; message=trim(message)//'invalid new snow density'; return; endif
     ! (compute storage capacity of new snow)
     leafInterceptCapSnow  = refInterceptCapSnow*(0.27_dp + 46._dp/scalarNewSnowDensity)  ! per unit leaf area (kg m-2)

    ! * option 2: maximum interception capacity an increasing function of air temerature
    case(stickySnow)
     airtemp_degC = scalarAirtemp - Tfreeze
     if    (airtemp_degC > -1._dp)then; leafScaleFactor = 4.0_dp
     elseif(airtemp_degC > -3._dp)then; leafScaleFactor = 1.5_dp*airtemp_degC + 5.5_dp
                                  else; leafScaleFactor = 1.0_dp
     endif
     leafInterceptCapSnow = refInterceptCapSnow*leafScaleFactor
     !write(*,'(a,1x,2(f20.10,1x))') 'airtemp_degC, leafInterceptCapSnow = ', airtemp_degC, leafInterceptCapSnow
     !pause 'in stickysnow'
 
    ! check we found the case
    case default
     message=trim(message)//'unable to identify option for maximum branch interception capacity'
     err=20; return

   end select ! identifying option for maximum branch interception capacity

   ! compute maximum interception capacity for the canopy
   canopyIceScaleFactor = leafInterceptCapSnow*exposedVAI

   ! (compute throughfall)
   scalarThroughfallSnow = scalarSnowfall*(scalarCanopyIceIter/canopyIceScaleFactor)
   throughfallDeriv      = scalarSnowfall/canopyIceScaleFactor

  endif  ! (if snow is falling)

  ! ** compute unloading
  scalarCanopySnowUnloading = snowUnloadingCoeff*scalarCanopyIceIter
  unloadingDeriv            = snowUnloadingCoeff

  ! ** compute iteration increment  
  flux = scalarSnowfall - scalarThroughfallSnow - scalarCanopySnowUnloading  ! net flux (kg m-2 s-1)
  delS = (flux*dt - (scalarCanopyIceIter - scalarCanopyIce))/(1._dp + (throughfallDeriv + unloadingDeriv)*dt)
  !print*, 'scalarCanopyIceIter, flux, delS = ', scalarCanopyIceIter, flux, delS  

  ! ** check for convergence
  resMass = scalarCanopyIceIter - (scalarCanopyIce + flux*dt)
  if(abs(resMass) < convTolerMass)exit

  ! ** check for non-convengence
  if(iter==maxiter)then; err=20; message=trim(message)//'failed to converge [mass]'; return; endif

  ! ** update value  
  scalarCanopyIceIter = scalarCanopyIceIter + delS

 end do  ! iterating

 ! add the unloading associated with melt drip
 scalarCanopySnowUnloading = scalarCanopySnowUnloading + unloading_melt

 ! update mass of ice on the canopy (kg m-2)
 scalarCanopyIce = scalarCanopyIceIter

 end subroutine canopySnow_muster


end module canopySnow_module

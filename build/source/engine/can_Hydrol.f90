module can_Hydrol_module
USE nrtype
implicit none
private
public::can_Hydrol
contains

 ! ************************************************************************************************
 ! new subroutine: compute water balance for the vegetation canopy
 ! ************************************************************************************************
 subroutine can_Hydrol(&
                       ! input: control
                       dt,                       & ! intent(in): time step (seconds)
                       iter,                     & ! intent(in): iteration index
                       computeVegFlux,           & ! intent(in): flag to denote if computing energy flux over vegetation
                       ! input: state variables at the current iteration
                       scalarCanopyIceIter,      & ! intent(in): trial mass of ice on the vegetation canopy at the current iteration (kg m-2)
                       scalarCanopyLiqIter,      & ! intent(in): trial mass of liquid water on the vegetation canopy at the current iteration (kg m-2)
                       ! input: canopy evaporation (from energy routine)
                       scalarCanopyEvaporation,  & ! intent(in): canopy evaporation/condensation (kg m-2 s-1)
                       ! output: updated state variables
                       scalarCanopyLiqNew,       & ! intent(out): updated mass of liquid water on the vegetation canopy at the current iteration (kg m-2)
                       ! output: error control
                       err,message)                ! intent(out): error control
 ! model variables, parameters, forcing data, etc.
 USE data_struc,only:mpar_data,mvar_data           ! data structures
 USE var_lookup,only:iLookPARAM,iLookMVAR          ! named variables for structure elements
 implicit none
 ! input: control
 real(dp),intent(in)           :: dt                         ! time step (seconds)
 integer(i4b),intent(in)       :: iter                       ! current iteration count
 logical(lgt),intent(in)       :: computeVegFlux             ! flag to indicate if we are computing fluxes over vegetation (.false. means veg is buried with snow)
 ! input: state variables at the current iteration
 real(dp),intent(in)           :: scalarCanopyIceIter        ! trial mass of ice on the vegetation canopy at the current iteration (kg m-2)
 real(dp),intent(in)           :: scalarCanopyLiqIter        ! trial mass of liquid water on the vegetation canopy at the current iteration (kg m-2)
 ! input: canopy evaporation (from energy routine)
 real(dp),intent(in)           :: scalarCanopyEvaporation    ! canopy evaporation (kg m-2 s-1)
 ! output: updated state variables
 real(dp),intent(out)          :: scalarCanopyLiqNew         ! updated mass of liquid water on the vegetation canopy at the current iteration (kg m-2)
 ! output: error control
 integer(i4b),intent(out)      :: err                        ! error code
 character(*),intent(out)      :: message                    ! error message
 ! ----------------------------------------------------------------------------------------------------
 ! define local variables
 character(LEN=256)            :: cmessage                   ! error message of downwind routine

 ! initialize error control
 err=0; message="can_Hydrol/"

 ! *****
 ! wrapper for the canopy hydrology sub-routine...
 ! ***********************************************
 ! used as a trick to both
 !  1) save typing; and
 !  2) "protect" variables through the use of the intent attribute
 call can_Hydrol_muster(&
                        ! input: control
                        dt,                                                        & ! intent(in): time step (seconds)
                        iter,                                                      & ! intent(in): iteration index
                        computeVegFlux,                                            & ! intent(in): flag to denote if computing energy flux over vegetation
                        ! input: state variables at the current iteration
                        scalarCanopyIceIter,                                       & ! intent(in): trial mass of ice on the vegetation canopy at the current iteration (kg m-2)
                        scalarCanopyLiqIter,                                       & ! intent(in): trial mass of liquid water on the vegetation canopy at the current iteration (kg m-2)
                        ! input: canopy evaporation (from energy routine)
                        scalarCanopyEvaporation,                                   & ! intent(in): canopy evaporation/condensation (kg m-2 s-1)
                        ! input: state variables at the start of the sub-step
                        mvar_data%var(iLookMVAR%scalarCanopyIce)%dat(1),           & ! intent(in): mass of ice on the vegetation canopy at the start of the sub-step (kg m-2)
                        mvar_data%var(iLookMVAR%scalarCanopyLiq)%dat(1),           & ! intent(in): mass of liquid water on the vegetation canopy at the start of the sub-step (kg m-2)
                        ! input: forcing and diagnostic variables
                        mvar_data%var(iLookMVAR%scalarRainfall)%dat(1),            & ! intent(in): computed rainfall rate (kg m-2 s-1)
                        mvar_data%var(iLookMVAR%scalarExposedLAI)%dat(1),          & ! intent(in): exposed leaf area index after burial by snow (m2 m-2)
                        mvar_data%var(iLookMVAR%scalarExposedSAI)%dat(1),          & ! intent(in): exposed stem area index after burial by snow (m2 m-2)                       
                        ! input: parameters
                        mpar_data%var(iLookPARAM%refInterceptCapRain),             & ! intent(in): reference canopy interception capacity for rain per unit leaf area (kg m-2)
                        mpar_data%var(iLookPARAM%throughfallScaleRain),            & ! intent(in): scaling factor for throughfall (rain) (-)
                        mpar_data%var(iLookPARAM%canopyDrainageCoeff),             & ! intent(in): time constant for drainage of liquid water from the forest canopy (s-1)
                        ! output: diagnostic variables
                        mvar_data%var(iLookMVAR%scalarCanopyLiqMax)%dat(1),        & ! intent(out): maximum interception storage capacity for liquid water (kg m-2)
                        mvar_data%var(iLookMVAR%scalarThroughfallRain)%dat(1),     & ! intent(out): rain that reaches the ground without ever touching the canopy (kg m-2 s-1)
                        mvar_data%var(iLookMVAR%scalarCanopyLiqDrainage)%dat(1),   & ! intent(out): drainage of liquid water from the vegetation canopy (kg m-2 s-1)
                        ! output: updated state variables
                        scalarCanopyLiqNew,                                        & ! intent(out): updated mass of liquid water on the vegetation canopy at the current iteration (kg m-2)
                        ! output: error control
                        err,cmessage)                                                ! intent(out): error control

 ! check for errors
 if(err/=0)then; message=trim(message)//trim(cmessage); return; endif

 end subroutine can_Hydrol





 ! ************************************************************************************************
 ! private subroutine: compute water balance for the vegetation canopy
 ! ************************************************************************************************
 subroutine can_Hydrol_muster(&
                              ! input: control
                              dt,                          & ! intent(in): time step (seconds)
                              iter,                        & ! intent(in): iteration index
                              computeVegFlux,              & ! intent(in): flag to denote if computing energy flux over vegetation
                              ! input: state variables at the current iteration
                              scalarCanopyIceIter,         & ! intent(in): trial mass of ice on the vegetation canopy at the current iteration (kg m-2)
                              scalarCanopyLiqIter,         & ! intent(in): trial mass of liquid water on the vegetation canopy at the current iteration (kg m-2)
                              ! input: canopy evaporation (from energy routine)
                              scalarCanopyEvaporation,     & ! intent(in): canopy evaporation (kg m-2 s-1)
                              ! input: state variables at the start of the sub-step
                              scalarCanopyIce,             & ! intent(in): mass of ice on the vegetation canopy at the start of the sub-step (kg m-2)
                              scalarCanopyLiq,             & ! intent(in): mass of liquid water on the vegetation canopy at the start of the sub-step (kg m-2)
                              ! input: forcing and diagnostic variables
                              scalarRainfall,              & ! intent(in): computed rainfall rate (kg m-2 s-1)
                              scalarExposedLAI,            & ! intent(in): exposed leaf area index after burial by snow (m2 m-2)
                              scalarExposedSAI,            & ! intent(in): exposed stem area index after burial by snow (m2 m-2)
                              ! input: parameters
                              refInterceptCapRain,         & ! intent(in): reference canopy interception capacity for rain per unit leaf area (kg m-2)
                              throughfallScaleRain,        & ! intent(in): scaling factor for throughfall (rain) (-)
                              canopyDrainageCoeff,         & ! intent(in): time constant for drainage of liquid water from the forest canopy (s-1)
                              ! output: diagnostic variables
                              scalarCanopyLiqMax,          & ! intent(out): maximum interception storage capacity for liquid water (kg m-2)
                              scalarThroughfallRain,       & ! intent(out): rain that reaches the ground without ever touching the canopy (kg m-2 s-1)
                              scalarCanopyLiqDrainage,     & ! intent(out): drainage of liquid water from the vegetation canopy (kg m-2 s-1)
                              ! output: updated state variables
                              scalarCanopyLiqNew,          & ! intent(out): updated mass of liquid water on the vegetation canopy at the current iteration (kg m-2)
                              ! output: error control
                              err,message)                   ! intent(out): error control
 implicit none
 ! input: control
 real(dp),intent(in)           :: dt                         ! time step (seconds)
 integer(i4b),intent(in)       :: iter                       ! current iteration count
 logical(lgt),intent(in)       :: computeVegFlux             ! flag to indicate if we are computing fluxes over vegetation (.false. means veg is buried with snow)
 ! input: state variables at the current iteration
 real(dp),intent(in)           :: scalarCanopyIceIter        ! trial mass of ice on the vegetation canopy at the current iteration (kg m-2)
 real(dp),intent(in)           :: scalarCanopyLiqIter        ! trial mass of liquid water on the vegetation canopy at the current iteration (kg m-2)
 ! input: canopy evaporation (from energy routine)
 real(dp),intent(in)           :: scalarCanopyEvaporation    ! canopy evaporation (kg m-2 s-1)
 ! input: state variables at the start of the sub-step
 real(dp),intent(in)           :: scalarCanopyIce            ! mass of ice on the vegetation canopy at the start of the sub-step (kg m-2)
 real(dp),intent(in)           :: scalarCanopyLiq            ! mass of liquid water on the vegetation canopy at the start of the sub-step (kg m-2)
 ! input: forcing and diagnostic variables
 real(dp),intent(in)           :: scalarRainfall             ! computed rainfall rate (kg m-2 s-1)
 real(dp),intent(in)           :: scalarExposedLAI           ! exposed leaf area index after burial by snow (m2 m-2)
 real(dp),intent(in)           :: scalarExposedSAI           ! exposed stem area index after burial by snow (m2 m-2)
 ! input: parameters
 real(dp),intent(in)           :: refInterceptCapRain        ! reference canopy interception capacity for rain per unit leaf area (kg m-2)
 real(dp),intent(in)           :: throughfallScaleRain       ! scaling factor for throughfall (rain) (-)
 real(dp),intent(in)           :: canopyDrainageCoeff        ! time constant for drainage of liquid water from the forest canopy (s-1)
 ! output: diagnostic variables
 real(dp),intent(out)          :: scalarCanopyLiqMax         ! maximum interception storage capacity for liquid water (kg m-2)
 real(dp),intent(out)          :: scalarThroughfallRain      ! rain that reaches the ground without ever touching the canopy (kg m-2 s-1)
 real(dp),intent(out)          :: scalarCanopyLiqDrainage    ! drainage of liquid water from the vegetation canopy (kg m-2 s-1)
 ! output: updated state variables
 real(dp),intent(out)          :: scalarCanopyLiqNew         ! updated mass of liquid water on the vegetation canopy at the current iteration (kg m-2)
 ! output: error control
 integer(i4b),intent(out)      :: err                        ! error code
 character(*),intent(out)      :: message                    ! error message
 ! ----------------------------------------------------------------------------------------------------
 ! define local variables
 integer(i4b)                  :: jiter                      ! local iteration
 integer(i4b),parameter        :: maxiter=10                 ! maximum number of iterations
 real(dp)                      :: exposedVAI                 ! exposed vegetation area index (LAI + SAI)
 real(dp)                      :: scalarCanopyLiqTrial       ! trial value of canopy liquid water (kg m-2)
 real(dp)                      :: canopyLiqDrainageDeriv     ! derivative in the drainage function for canopy liquid water (s-1)
 real(dp)                      :: flux                       ! net flux (kg m-2 s-1)
 real(dp)                      :: delS                       ! change in storage (kg m-2)
 real(dp)                      :: res                        ! residual (kg m-2)
 real(dp),parameter            :: convToler=0.0001_dp        ! convergence tolerance (kg m-2)
 ! ----------------------------------------------------------------------------------------------------
 ! initialize error control
 err=0; message="can_Hydrol_muster/"

 ! set throughfall to inputs if vegetation is completely buried with snow
 if(.not.computeVegFlux)then
  scalarThroughfallRain     = scalarRainfall
  scalarCanopyLiqDrainage   = 0._dp
  return
 endif

 ! initialize canopy liquid water
 scalarCanopyLiqTrial = scalarCanopyLiqIter

 ! compute the maximum storage of liquid water (kg m-2)
 scalarCanopyLiqMax = refInterceptCapRain*(scalarExposedSAI + scalarExposedLAI)
 !print*, 'refInterceptCapRain, (scalarExposedSAI + scalarExposedLAI), scalarCanopyLiqMax = ', &
 !         refInterceptCapRain, (scalarExposedSAI + scalarExposedLAI), scalarCanopyLiqMax

 ! set throughfall to zero (kg m-2)
 scalarThroughfallRain = 0._dp

 ! test
 !print*, 'in canHydrol: scalarCanopyEvaporation = ', scalarCanopyEvaporation

 ! ***** estimate the updated value for liquid water
 do jiter=1,maxiter

  ! compute drainage of liquid water from the canopy (kg m-2 s-1)
  if(scalarCanopyLiqTrial > scalarCanopyLiqMax)then
   scalarCanopyLiqDrainage = canopyDrainageCoeff*(scalarCanopyLiqTrial - scalarCanopyLiqMax)
   canopyLiqDrainageDeriv  = canopyDrainageCoeff
  else
   scalarCanopyLiqDrainage = 0._dp
   canopyLiqDrainageDeriv  = 0._dp
  endif

  ! ** compute iteration increment  
  flux = scalarRainfall - scalarThroughfallRain + scalarCanopyEvaporation - scalarCanopyLiqDrainage  ! net flux (kg m-2 s-1)
  delS = (flux*dt - (scalarCanopyLiqTrial - scalarCanopyLiq) - (scalarCanopyIceIter - scalarCanopyIce) ) / (1._dp + canopyLiqDrainageDeriv*dt)
  !print*, 'scalarRainfall, scalarThroughfallRain, scalarCanopyEvaporation, scalarCanopyLiqDrainage = ', &
  !         scalarRainfall, scalarThroughfallRain, scalarCanopyEvaporation, scalarCanopyLiqDrainage
  !print*, 'jiter, scalarCanopyLiqTrial, scalarCanopyLiq, scalarCanopyIceIter, scalarCanopyIce = ', &
  !         jiter, scalarCanopyLiqTrial, scalarCanopyLiq, scalarCanopyIceIter, scalarCanopyIce
  !print*, '-(scalarCanopyLiqTrial - scalarCanopyLiq) = ', -(scalarCanopyLiqTrial - scalarCanopyLiq)
  !print*, '(scalarCanopyIceIter - scalarCanopyIce) = ', (scalarCanopyIceIter - scalarCanopyIce)
  !print*, '-(scalarCanopyLiqTrial - scalarCanopyLiq) - (scalarCanopyIceIter - scalarCanopyIce) = ', &
  !         -(scalarCanopyLiqTrial - scalarCanopyLiq) - (scalarCanopyIceIter - scalarCanopyIce)

  ! ** check for convergence
  res = scalarCanopyLiqTrial - (scalarCanopyLiq + flux*dt) + (scalarCanopyIceIter - scalarCanopyIce)
  !print*, 'res, delS = ', res, delS
  if(abs(res) < convToler)exit

  ! ** check for non-convengence
  if(jiter==maxiter)then; err=20; message=trim(message)//'failed to converge'; return; endif

  ! ** update value  
  scalarCanopyLiqTrial = scalarCanopyLiqTrial + delS

 end do  ! iterating

 !write(*,'(a,1x,10(e20.10,1x))') 'scalarCanopyLiqTrial, scalarCanopyLiqDrainage = ', scalarCanopyLiqTrial, scalarCanopyLiqDrainage

 ! update state variables
 scalarCanopyLiqNew = scalarCanopyLiqTrial

 end subroutine can_Hydrol_muster

end module can_Hydrol_module

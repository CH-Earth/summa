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
                       ! input: state variables
                       scalarCanopyIceIter,      & ! intent(in): trial mass of ice on the vegetation canopy at the current iteration (kg m-2)
                       scalarCanopyLiqIter,      & ! intent(in): trial mass of liquid water on the vegetation canopy at the current iteration (kg m-2)
                       ! output: updated state variables
                       scalarCanopyIceNew,       & ! intent(out): updated mass of ice on the vegetation canopy at the current iteration (kg m-2)
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
 ! input: state variables
 real(dp),intent(in)           :: scalarCanopyLiqIter        ! trial mass of liquid water on the vegetation canopy at the current iteration (kg m-2)
 real(dp),intent(in)           :: scalarCanopyIceIter        ! trial mass of ice on the vegetation canopy at the current iteration (kg m-2)
 ! output: updated state variables
 real(dp),intent(out)          :: scalarCanopyLiqNew         ! updated mass of liquid water on the vegetation canopy at the current iteration (kg m-2)
 real(dp),intent(out)          :: scalarCanopyIceNew         ! updated mass of ice on the vegetation canopy at the current iteration (kg m-2)
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
                        ! input: state variables
                        scalarCanopyIceIter,                                       & ! intent(in): trial mass of ice on the vegetation canopy at the current iteration (kg m-2)
                        scalarCanopyLiqIter,                                       & ! intent(in): trial mass of liquid water on the vegetation canopy at the current iteration (kg m-2)
                        ! input: forcing variables
                        mvar_data%var(iLookMVAR%scalarSnowfall)%dat(1),            & ! intent(in): computed snowfall rate (kg m-2 s-1)
                        mvar_data%var(iLookMVAR%scalarRainfall)%dat(1),            & ! intent(in): computed rainfall rate (kg m-2 s-1)
                        ! input: diagnostic variables
                        mvar_data%var(iLookMVAR%scalarExposedLAI)%dat(1),          & ! intent(in): exposed leaf area index after burial by snow (m2 m-2)
                        mvar_data%var(iLookMVAR%scalarExposedSAI)%dat(1),          & ! intent(in): exposed stem area index after burial by snow (m2 m-2)                       
                        mvar_data%var(iLookMVAR%scalarCanopyIceMax)%dat(1),        & ! intent(in): maximum interception storage capacity for ice (kg m-2)
                        mvar_data%var(iLookMVAR%scalarCanopyLiqMax)%dat(1),        & ! intent(in): maximum interception storage capacity for liquid water (kg m-2)
                        ! input: parameters
                        mpar_data%var(iLookPARAM%throughfallScaleSnow),            & ! intent(in): scaling factor for throughfall (snow) (-)
                        mpar_data%var(iLookPARAM%throughfallScaleRain),            & ! intent(in): scaling factor for throughfall (rain) (-)
                        mpar_data%var(iLookPARAM%snowUnloadingCoeff),              & ! intent(in): time constant for unloading of snow from the forest canopy (s-1)
                        mpar_data%var(iLookPARAM%canopyDrainageCoeff),             & ! intent(in): time constant for drainage of liquid water from the forest canopy (s-1)
                        ! output: diagnostic variables
                        mvar_data%var(iLookMVAR%scalarThroughfallSnow)%dat(1),     & ! intent(out): snow that reaches the ground without ever touching the canopy (kg m-2 s-1)
                        mvar_data%var(iLookMVAR%scalarThroughfallRain)%dat(1),     & ! intent(out): rain that reaches the ground without ever touching the canopy (kg m-2 s-1)
                        mvar_data%var(iLookMVAR%scalarCanopySnowUnloading)%dat(1), & ! intent(out): unloading of snow from the vegetion canopy (kg m-2 s-1)
                        mvar_data%var(iLookMVAR%scalarCanopyLiqDrainage)%dat(1),   & ! intent(out): drainage of liquid water from the vegetation canopy (kg m-2 s-1)
                        ! output: updated state variables
                        scalarCanopyIceNew,                                        & ! intent(out): updated mass of ice on the vegetation canopy at the current iteration (kg m-2)
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
                              ! input: state variables
                              scalarCanopyIceIter,         & ! intent(in): trial mass of ice on the vegetation canopy at the current iteration (kg m-2)
                              scalarCanopyLiqIter,         & ! intent(in): trial mass of liquid water on the vegetation canopy at the current iteration (kg m-2)
                              ! input: forcing variables
                              scalarSnowfall,              & ! intent(in): computed snowfall rate (kg m-2 s-1)
                              scalarRainfall,              & ! intent(in): computed rainfall rate (kg m-2 s-1)
                              ! input: diagnostic variables
                              scalarExposedLAI,            & ! intent(in): exposed leaf area index after burial by snow (m2 m-2)
                              scalarExposedSAI,            & ! intent(in): exposed stem area index after burial by snow (m2 m-2)
                              scalarCanopyIceMax,          & ! intent(in): maximum interception storage capacity for ice (kg m-2)
                              scalarCanopyLiqMax,          & ! intent(in): maximum interception storage capacity for liquid water (kg m-2)
                              ! input: parameters
                              throughfallScaleSnow,        & ! intent(in): scaling factor for throughfall (snow) (-)
                              throughfallScaleRain,        & ! intent(in): scaling factor for throughfall (rain) (-)
                              snowUnloadingCoeff,          & ! intent(in): time constant for unloading of snow from the forest canopy (s-1)
                              canopyDrainageCoeff,         & ! intent(in): time constant for drainage of liquid water from the forest canopy (s-1)
                              ! output: diagnostic variables
                              scalarThroughfallSnow,       & ! intent(out): snow that reaches the ground without ever touching the canopy (kg m-2 s-1)
                              scalarThroughfallRain,       & ! intent(out): rain that reaches the ground without ever touching the canopy (kg m-2 s-1)
                              scalarCanopySnowUnloading,   & ! intent(out): unloading of snow from the vegetion canopy (kg m-2 s-1)
                              scalarCanopyLiqDrainage,     & ! intent(out): drainage of liquid water from the vegetation canopy (kg m-2 s-1)
                              ! output: updated state variables
                              scalarCanopyIceNew,          & ! intent(out): updated mass of ice on the vegetation canopy at the current iteration (kg m-2)
                              scalarCanopyLiqNew,          & ! intent(out): updated mass of liquid water on the vegetation canopy at the current iteration (kg m-2)
                              ! output: error control
                              err,message)                   ! intent(out): error control
 implicit none
 ! input: control
 real(dp),intent(in)           :: dt                         ! time step (seconds)
 integer(i4b),intent(in)       :: iter                       ! current iteration count
 logical(lgt),intent(in)       :: computeVegFlux             ! flag to indicate if we are computing fluxes over vegetation (.false. means veg is buried with snow)
 ! input: state variables
 real(dp),intent(in)           :: scalarCanopyIceIter        ! trial mass of ice on the vegetation canopy at the current iteration (kg m-2)
 real(dp),intent(in)           :: scalarCanopyLiqIter        ! trial mass of liquid water on the vegetation canopy at the current iteration (kg m-2)
 ! input: forcing variables
 real(dp),intent(in)           :: scalarSnowfall             ! computed snowfall rate (kg m-2 s-1)
 real(dp),intent(in)           :: scalarRainfall             ! computed rainfall rate (kg m-2 s-1)
 ! input: diagnostic variables
 real(dp),intent(in)           :: scalarExposedLAI           ! exposed leaf area index after burial by snow (m2 m-2)
 real(dp),intent(in)           :: scalarExposedSAI           ! exposed stem area index after burial by snow (m2 m-2)
 real(dp),intent(in)           :: scalarCanopyIceMax         ! maximum interception storage capacity for ice (kg m-2)
 real(dp),intent(in)           :: scalarCanopyLiqMax         ! maximum interception storage capacity for liquid water (kg m-2)
 ! input: parameters
 real(dp),intent(in)           :: throughfallScaleSnow       ! scaling factor for throughfall (snow) (-)
 real(dp),intent(in)           :: throughfallScaleRain       ! scaling factor for throughfall (rain) (-)
 real(dp),intent(in)           :: snowUnloadingCoeff         ! time constant for unloading of snow from the forest canopy (s-1)
 real(dp),intent(in)           :: canopyDrainageCoeff        ! time constant for drainage of liquid water from the forest canopy (s-1)
 ! output: diagnostic variables
 real(dp),intent(out)          :: scalarThroughfallSnow      ! snow that reaches the ground without ever touching the canopy (kg m-2 s-1)
 real(dp),intent(out)          :: scalarThroughfallRain      ! rain that reaches the ground without ever touching the canopy (kg m-2 s-1)
 real(dp),intent(out)          :: scalarCanopySnowUnloading  ! unloading of snow from the vegetion canopy (kg m-2 s-1)
 real(dp),intent(out)          :: scalarCanopyLiqDrainage    ! drainage of liquid water from the vegetation canopy (kg m-2 s-1)
 ! output: updated state variables
 real(dp),intent(out)          :: scalarCanopyIceNew         ! updated mass of ice on the vegetation canopy at the current iteration (kg m-2)
 real(dp),intent(out)          :: scalarCanopyLiqNew         ! updated mass of liquid water on the vegetation canopy at the current iteration (kg m-2)
 ! output: error control
 integer(i4b),intent(out)      :: err                        ! error code
 character(*),intent(out)      :: message                    ! error message
 ! ----------------------------------------------------------------------------------------------------
 ! define local variables
 integer(i4b)                  :: jiter                      ! local iteration
 integer(i4b),parameter        :: maxiter=1                  ! maximum number of iterations
 real(dp)                      :: exposedVAI                 ! exposed vegetation area index (LAI + SAI)
 real(dp)                      :: canopyLiqDrainageDeriv     ! derivative in the drainage function for canopy liquid water (s-1)

 ! initialize error control
 err=0; message="can_Hydrol_muster/"

 ! set throughfall to inputs if vegetation is completely buried with snow
 if(.not.computeVegFlux)then
  scalarThroughfallSnow     = scalarSnowfall
  scalarThroughfallRain     = scalarRainfall
  scalarCanopySnowUnloading = 0._dp
  scalarCanopyLiqDrainage   = 0._dp
  return
 endif

 print*, 'scalarCanopyLiqIter, scalarCanopyLiqMax = ', scalarCanopyLiqIter, scalarCanopyLiqMax

 ! compute the exposed vegetation area index (m2 m-2)
 exposedVAI = scalarExposedSAI + scalarExposedLAI

 ! compute throughfall (kg m-2 s-1)
 scalarThroughfallSnow = scalarSnowfall*exp(-throughfallScaleSnow*exposedVAI)
 scalarThroughfallRain = scalarRainfall*exp(-throughfallScaleRain*exposedVAI)

 ! compute unloading of snow from the canopy (kg m-2 s-1)
 scalarCanopySnowUnloading = snowUnloadingCoeff*scalarCanopyIceIter

 ! ***** estimate the updated value for liquid water
 do jiter=1,maxiter

  ! compute drainage of liquid water from the canopy (kg m-2 s-1)
  if(scalarCanopyLiqIter > scalarCanopyLiqMax)then
   scalarCanopyLiqDrainage = canopyDrainageCoeff*(scalarCanopyLiqIter - scalarCanopyLiqMax)
   canopyLiqDrainageDeriv  = canopyDrainageCoeff
  else
   scalarCanopyLiqDrainage = 0._dp
   canopyLiqDrainageDeriv  = 0._dp
  endif

  print*, 'scalarCanopyLiqDrainage = ', scalarCanopyLiqDrainage  

 end do  ! iterations for canopy liquid water 





 ! check parameters
 print*, 'throughfallScaleSnow  = ', throughfallScaleSnow     ! scaling factor for throughfall (snow) (-)
 print*, 'throughfallScaleRain  = ', throughfallScaleRain     ! scaling factor for throughfall (rain) (-)
 print*, 'snowUnloadingCoeff    = ', snowUnloadingCoeff       ! time constant for unloading of snow from the forest canopy (s-1)
 print*, 'canopyDrainageCoeff   = ', canopyDrainageCoeff      ! time constant for drainage of liquid water from the forest canopy (s-1)
 pause 'in can_Hydrol'


 ! update state variables
 scalarCanopyIceNew = scalarCanopyIceIter
 scalarCanopyLiqNew = scalarCanopyLiqIter

 end subroutine can_Hydrol_muster

end module can_Hydrol_module

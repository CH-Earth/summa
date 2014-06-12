module vegLiqFlux_module
USE nrtype
implicit none
private
public::vegLiqFlux
contains

 ! ************************************************************************************************
 ! new subroutine: compute water balance for the vegetation canopy
 ! ************************************************************************************************
 subroutine vegLiqFlux(&
                       ! input
                       computeVegFlux,               & ! intent(in): flag to denote if computing energy flux over vegetation
                       scalarCanopyLiqTrial,         & ! intent(in): trial mass of liquid water on the vegetation canopy at the current iteration (kg m-2)
                       scalarRainfall,               & ! intent(in): rainfall (kg m-2 s-1)
                       scalarCanopyLiqMax,           & ! intent(in): maximum storage before canopy drainage begins (kg m-2 s-1)
                       scalarCanopyDrainageCoeff,    & ! intent(in): canopy drainage coefficient (s-1)
                       ! output
                       scalarThroughfallRain,        & ! intent(out): rain that reaches the ground without ever touching the canopy (kg m-2 s-1)
                       scalarCanopyLiqDrainage,      & ! intent(out): drainage of liquid water from the vegetation canopy (kg m-2 s-1)
                       scalarCanopyLiqDrainageDeriv, & ! intent(out): derivative in canopy drainage w.r.t. canopy liquid water (s-1)
                       err,message)                    ! intent(out): error control
 implicit none
 ! input
 logical(lgt),intent(in)       :: computeVegFlux               ! flag to indicate if we are computing fluxes over vegetation (.false. means veg is buried with snow)
 real(dp),intent(in)           :: scalarCanopyLiqTrial         ! trial mass of liquid water on the vegetation canopy at the current iteration (kg m-2)
 real(dp),intent(in)           :: scalarRainfall               ! rainfall (kg m-2 s-1)
 real(dp),intent(in)           :: scalarCanopyLiqMax           ! maximum storage before canopy drainage begins (kg m-2 s-1)
 real(dp),intent(in)           :: scalarCanopyDrainageCoeff    ! canopy drainage coefficient (s-1)
 ! output
 real(dp),intent(out)          :: scalarThroughfallRain        ! rain that reaches the ground without ever touching the canopy (kg m-2 s-1)
 real(dp),intent(out)          :: scalarCanopyLiqDrainage      ! drainage of liquid water from the vegetation canopy (kg m-2 s-1)
 real(dp),intent(out)          :: scalarCanopyLiqDrainageDeriv ! derivative in canopy drainage w.r.t. canopy liquid water (s-1)
 integer(i4b),intent(out)      :: err                          ! error code
 character(*),intent(out)      :: message                      ! error message
 ! ----------------------------------------------------------------------------------------------------
 ! initialize error control
 err=0; message="vegLiqFlux/"

 ! set throughfall to inputs if vegetation is completely buried with snow
 if(.not.computeVegFlux)then
  scalarThroughfallRain        = scalarRainfall
  scalarCanopyLiqDrainage      = 0._dp
  scalarCanopyLiqDrainageDeriv = 0._dp
  return
 endif

 ! compute canopy drainage
 if(scalarCanopyLiqTrial > scalarCanopyLiqMax)then
  scalarCanopyLiqDrainage       = scalarCanopyDrainageCoeff*(scalarCanopyLiqTrial - scalarCanopyLiqMax)
  scalarCanopyLiqDrainageDeriv  = scalarCanopyDrainageCoeff
 else
  scalarCanopyLiqDrainage       = 0._dp
  scalarCanopyLiqDrainageDeriv  = 0._dp
 endif

 end subroutine vegLiqFlux

end module vegLiqFlux_module

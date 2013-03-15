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
                       ! input
                       dt,                       & ! intent(in): time step (seconds)
                       iter,                     & ! intent(in): iteration index
                       scalarCanopyIceIter,      & ! intent(in): trial mass of ice on the vegetation canopy at the current iteration (kg m-2)
                       scalarCanopyLiqIter,      & ! intent(in): trial mass of liquid water on the vegetation canopy at the current iteration (kg m-2)
                       ! output
                       scalarCanopyIceNew,       & ! intent(out): updated mass of ice on the vegetation canopy at the current iteration (kg m-2)
                       scalarCanopyLiqNew,       & ! intent(out): updated mass of liquid water on the vegetation canopy at the current iteration (kg m-2)
                       err,message)                ! intent(out): error control
 implicit none
 ! input variables
 real(dp),intent(in)           :: dt                         ! time step (seconds)
 integer(i4b),intent(in)       :: iter                       ! current iteration count
 real(dp),intent(in)           :: scalarCanopyIceIter        ! trial mass of ice on the vegetation canopy at the current iteration (kg m-2)
 real(dp),intent(in)           :: scalarCanopyLiqIter        ! trial mass of liquid water on the vegetation canopy at the current iteration (kg m-2)
 ! output variables
 real(dp),intent(out)          :: scalarCanopyIceNew         ! updated mass of ice on the vegetation canopy at the current iteration (kg m-2)
 real(dp),intent(out)          :: scalarCanopyLiqNew         ! updated mass of liquid water on the vegetation canopy at the current iteration (kg m-2)
 integer(i4b),intent(out)      :: err                        ! error code
 character(*),intent(out)      :: message                    ! error message
 ! ----------------------------------------------------------------------------------------------------
 ! define local variables
 real(dp),parameter            :: verySmall=epsilon(dt)      ! a very small number

 ! initialize error control
 err=0; message="can_Hydrol/"

 ! placeholder for actual routine
 if(abs(scalarCanopyIceIter)>verySmall .or. abs(scalarCanopyLiqIter)>verySmall)then
  message=trim(message)//'canopy water balance not implemented yet -- canopy liquid and ice content must be zero'
  err=20; return
 endif

 ! update state variables
 scalarCanopyIceNew = scalarCanopyIceIter
 scalarCanopyLiqNew = scalarCanopyIceIter

 end subroutine can_Hydrol

end module can_Hydrol_module

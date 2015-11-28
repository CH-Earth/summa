! SUMMA - Structure for Unifying Multiple Modeling Alternatives
! Copyright (C) 2014-2015 NCAR/RAL
!
! This file is part of SUMMA
!
! For more information see: http://www.ral.ucar.edu/projects/summa
!
! This program is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! This program is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License
! along with this program.  If not, see <http://www.gnu.org/licenses/>.

module vegLiqFlux_module
USE nrtype
! look-up values for the choice of canopy shortwave radiation method
USE mDecisions_module,only:         &
                      unDefined,    & ! original model (no flexibility in canopy interception): 100% of rainfall is intercepted by the vegetation canopy
                      sparseCanopy, & ! fraction of rainfall that never hits the canopy (throughfall); drainage above threshold
                      storageFunc     ! throughfall a function of canopy storage; 100% throughfall when canopy is at capacity
implicit none
private
public::vegLiqFlux
contains


 ! ************************************************************************************************
 ! public subroutine vegLiqFlux: compute water balance for the vegetation canopy
 ! ************************************************************************************************
 subroutine vegLiqFlux(&
                       ! input
                       computeVegFlux,               & ! intent(in): flag to denote if computing energy flux over vegetation
                       scalarCanopyLiqTrial,         & ! intent(in): trial mass of liquid water on the vegetation canopy at the current iteration (kg m-2)
                       scalarRainfall,               & ! intent(in): rainfall rate (kg m-2 s-1)
                       ! output
                       scalarThroughfallRain,        & ! intent(out): rain that reaches the ground without ever touching the canopy (kg m-2 s-1)
                       scalarCanopyLiqDrainage,      & ! intent(out): drainage of liquid water from the vegetation canopy (kg m-2 s-1)
                       scalarThroughfallRainDeriv,   & ! intent(out): derivative in throughfall w.r.t. canopy liquid water (s-1)
                       scalarCanopyLiqDrainageDeriv, & ! intent(out): derivative in canopy drainage w.r.t. canopy liquid water (s-1)
                       err,message)                    ! intent(out): error control
 ! model decisions
 USE data_struc,only:model_decisions                              ! model decision structure
 USE var_lookup,only:iLookDECISIONS                               ! named variables for elements of the decision structure
 ! model variables, parameters, forcing data, etc.
 USE data_struc,only:mpar_data,mvar_data    ! data structures
 USE var_lookup,only:iLookATTR,iLookTYPE,iLookPARAM,iLookFORCE,iLookMVAR,iLookINDEX ! named variables for structure elements
 implicit none
 ! input
 logical(lgt),intent(in)       :: computeVegFlux               ! flag to indicate if we are computing fluxes over vegetation (.false. means veg is buried with snow)
 real(dp),intent(in)           :: scalarCanopyLiqTrial         ! trial mass of liquid water on the vegetation canopy at the current iteration (kg m-2)
 real(dp),intent(in)           :: scalarRainfall               ! rainfall (kg m-2 s-1)
 ! output
 real(dp),intent(out)          :: scalarThroughfallRain        ! rain that reaches the ground without ever touching the canopy (kg m-2 s-1)
 real(dp),intent(out)          :: scalarCanopyLiqDrainage      ! drainage of liquid water from the vegetation canopy (kg m-2 s-1)
 real(dp),intent(out)          :: scalarThroughfallRainDeriv   ! derivative in throughfall w.r.t. canopy liquid water (s-1)
 real(dp),intent(out)          :: scalarCanopyLiqDrainageDeriv ! derivative in canopy drainage w.r.t. canopy liquid water (s-1)
 integer(i4b),intent(out)      :: err                          ! error code
 character(*),intent(out)      :: message                      ! error message
 ! ------------------------------------------------------------------------------------------------------------------------------------------------------
 ! local variables
 character(LEN=256)            :: cmessage                     ! error message of downwind routine
 ! ----------------------------------------------------------------------------------------------------
 ! initialize error control
 err=0; message="vegLiqFlux/"

 ! wrapper routine (makes use of data structures and protects variables with the intent attribute)
 call vegLiqFlux_muster(&
                        ! input
                        computeVegFlux,                                       & ! intent(in): flag to denote if computing energy flux over vegetation
                        model_decisions(iLookDECISIONS%cIntercept)%iDecision, & ! intent(in): index defining choice of parameterization for canopy interception
                        scalarCanopyLiqTrial,                                 & ! intent(in): trial mass of liquid water on the vegetation canopy at the current iteration (kg m-2)
                        scalarRainfall,                                       & ! intent(in): rainfall rate (kg m-2 s-1)
                        ! input: forcing and parameters from data structures
                        mvar_data%var(iLookMVAR%scalarCanopyLiqMax)%dat(1),   & ! intent(in): maximum storage before canopy drainage begins (kg m-2 s-1)
                        mpar_data%var(iLookPARAM%throughfallScaleRain),       & ! intent(in): fraction of rain that hits the ground without touching the canopy (-)
                        mpar_data%var(iLookPARAM%canopyDrainageCoeff),        & ! intent(in): canopy drainage coefficient (s-1)
                        ! output
                        scalarThroughfallRain,                                & ! intent(out): rain that falls through the canopy (kg m-2 s-1)
                        scalarCanopyLiqDrainage,                              & ! intent(out): drainage of liquid water from the vegetation canopy (kg m-2 s-1)
                        scalarThroughfallRainDeriv,                           & ! intent(out): derivative in throughfall w.r.t. canopy liquid water (s-1)
                        scalarCanopyLiqDrainageDeriv,                         & ! intent(out): derivative in canopy drainage w.r.t. canopy liquid water (s-1)
                        err,cmessage)                                           ! intent(out): error control
 if(err/=0)then; message=trim(message)//trim(cmessage); return; endif


 end subroutine vegLiqFlux


 ! ************************************************************************************************
 ! private subroutine vegLiqFlux_muster: compute water balance for the vegetation canopy
 ! ************************************************************************************************
 subroutine vegLiqFlux_muster(&
                              ! input
                              computeVegFlux,               & ! intent(in): flag to denote if computing energy flux over vegetation
                              ixCanopyInterception,         & ! intent(in): index defining choice of parameterization for canopy interception
                              scalarCanopyLiqTrial,         & ! intent(in): trial mass of liquid water on the vegetation canopy at the current iteration (kg m-2)
                              scalarRainfall,               & ! intent(in): rainfall (kg m-2 s-1)
                              scalarCanopyLiqMax,           & ! intent(in): maximum storage before canopy drainage begins (kg m-2 s-1)
                              scalarThroughfallScaleRain,   & ! intent(in): fraction of rain that hits the ground without touching the canopy (-)
                              scalarCanopyDrainageCoeff,    & ! intent(in): canopy drainage coefficient (s-1)
                              ! output
                              scalarThroughfallRain,        & ! intent(out): rain that reaches the ground without ever touching the canopy (kg m-2 s-1)
                              scalarCanopyLiqDrainage,      & ! intent(out): drainage of liquid water from the vegetation canopy (kg m-2 s-1)
                              scalarThroughfallRainDeriv,   & ! intent(out): derivative in throughfall w.r.t. canopy liquid water (s-1)
                              scalarCanopyLiqDrainageDeriv, & ! intent(out): derivative in canopy drainage w.r.t. canopy liquid water (s-1)
                              err,message)                    ! intent(out): error control
 implicit none
 ! input
 logical(lgt),intent(in)       :: computeVegFlux               ! flag to indicate if we are computing fluxes over vegetation (.false. means veg is buried with snow)
 integer(i4b),intent(in)       :: ixCanopyInterception         ! index defining choice of parameterization for canopy interception
 real(dp),intent(in)           :: scalarCanopyLiqTrial         ! trial mass of liquid water on the vegetation canopy at the current iteration (kg m-2)
 real(dp),intent(in)           :: scalarRainfall               ! rainfall (kg m-2 s-1)
 real(dp),intent(in)           :: scalarCanopyLiqMax           ! maximum storage before canopy drainage begins (kg m-2 s-1)
 real(dp),intent(in)           :: scalarThroughfallScaleRain   ! fraction of rain that hits the ground without touching the canopy (-)
 real(dp),intent(in)           :: scalarCanopyDrainageCoeff    ! canopy drainage coefficient (s-1)
 ! output
 real(dp),intent(out)          :: scalarThroughfallRain        ! rain that reaches the ground without ever touching the canopy (kg m-2 s-1)
 real(dp),intent(out)          :: scalarCanopyLiqDrainage      ! drainage of liquid water from the vegetation canopy (kg m-2 s-1)
 real(dp),intent(out)          :: scalarThroughfallRainDeriv   ! derivative in throughfall w.r.t. canopy liquid water (s-1)
 real(dp),intent(out)          :: scalarCanopyLiqDrainageDeriv ! derivative in canopy drainage w.r.t. canopy liquid water (s-1)
 integer(i4b),intent(out)      :: err                          ! error code
 character(*),intent(out)      :: message                      ! error message
 ! ----------------------------------------------------------------------------------------------------
 ! initialize error control
 err=0; message="vegLiqFlux_muster/"

 ! set throughfall to inputs if vegetation is completely buried with snow
 if(.not.computeVegFlux)then
  scalarThroughfallRain        = scalarRainfall
  scalarCanopyLiqDrainage      = 0._dp
  scalarThroughfallRainDeriv   = 0._dp
  scalarCanopyLiqDrainageDeriv = 0._dp
  return
 endif

 ! compute throughfall
 select case(ixCanopyInterception)

  ! original model (no flexibility in canopy interception): 100% of rainfall is intercepted by the vegetation canopy
  ! NOTE: this could be done with scalarThroughfallScaleRain=0, though requires setting scalarThroughfallScaleRain in all test cases
  case(unDefined)
   scalarThroughfallRain      = 0._dp
   scalarThroughfallRainDeriv = 0._dp

  ! fraction of rainfall hits the ground without ever touching the canopy
  case(sparseCanopy)
   scalarThroughfallRain      = scalarThroughfallScaleRain*scalarRainfall
   scalarThroughfallRainDeriv = 0._dp

  ! throughfall a function of canopy storage
  case(storageFunc)

   ! throughfall during wetting-up phase
   if(scalarCanopyLiqTrial < scalarCanopyLiqMax)then
    scalarThroughfallRain      = scalarRainfall*(scalarCanopyLiqTrial/scalarCanopyLiqMax)
    scalarThroughfallRainDeriv = scalarRainfall/scalarCanopyLiqMax

   ! all rain falls through the canopy when the canopy is at capacity
   else
    scalarThroughfallRain      = scalarRainfall
    scalarThroughfallRainDeriv = 0._dp
   endif

  case default; err=20; message=trim(message)//'unable to identify option for canopy interception'; return

 end select ! (option for canopy interception)

 ! compute canopy drainage
 if(scalarCanopyLiqTrial > scalarCanopyLiqMax)then
  scalarCanopyLiqDrainage       = scalarCanopyDrainageCoeff*(scalarCanopyLiqTrial - scalarCanopyLiqMax)
  scalarCanopyLiqDrainageDeriv  = scalarCanopyDrainageCoeff
 else
  scalarCanopyLiqDrainage       = 0._dp
  scalarCanopyLiqDrainageDeriv  = 0._dp
 endif

 end subroutine vegLiqFlux_muster


end module vegLiqFlux_module

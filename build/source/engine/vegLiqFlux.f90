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
                       scalarCanopyLiqDrainageDeriv, & ! intent(out): derivative in canopy drainage w.r.t. canopy liquid water (s-1)
                       err,message)                    ! intent(out): error control
 ! model variables, parameters, forcing data, etc.
 USE data_struc,only:attr_data,type_data,mpar_data,forc_data,mvar_data,indx_data    ! data structures
 USE var_lookup,only:iLookATTR,iLookTYPE,iLookPARAM,iLookFORCE,iLookMVAR,iLookINDEX ! named variables for structure elements
 implicit none
 ! input
 logical(lgt),intent(in)       :: computeVegFlux               ! flag to indicate if we are computing fluxes over vegetation (.false. means veg is buried with snow)
 real(dp),intent(in)           :: scalarCanopyLiqTrial         ! trial mass of liquid water on the vegetation canopy at the current iteration (kg m-2)
 real(dp),intent(in)           :: scalarRainfall               ! rainfall (kg m-2 s-1)
 ! output
 real(dp),intent(out)          :: scalarThroughfallRain        ! rain that reaches the ground without ever touching the canopy (kg m-2 s-1)
 real(dp),intent(out)          :: scalarCanopyLiqDrainage      ! drainage of liquid water from the vegetation canopy (kg m-2 s-1)
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
                        computeVegFlux,                                     & ! intent(in): flag to denote if computing energy flux over vegetation
                        scalarCanopyLiqTrial,                               & ! intent(in): trial mass of liquid water on the vegetation canopy at the current iteration (kg m-2)
                        scalarRainfall,                                     & ! intent(in): rainfall rate (kg m-2 s-1)
                        ! input: forcing and parameters from data structures
                        mvar_data%var(iLookMVAR%scalarCanopyLiqMax)%dat(1), & ! intent(in): maximum storage before canopy drainage begins (kg m-2 s-1)
                        mpar_data%var(iLookPARAM%canopyDrainageCoeff),      & ! intent(in): canopy drainage coefficient (s-1)
                        ! output
                        scalarThroughfallRain,                              & ! intent(out): rain that reaches the ground without ever touching the canopy (kg m-2 s-1)
                        scalarCanopyLiqDrainage,                            & ! intent(out): drainage of liquid water from the vegetation canopy (kg m-2 s-1)
                        scalarCanopyLiqDrainageDeriv,                       & ! intent(out): derivative in canopy drainage w.r.t. canopy liquid water (s-1)
                        err,cmessage)                                         ! intent(out): error control
 if(err/=0)then; message=trim(message)//trim(cmessage); return; endif


 end subroutine vegLiqFlux


 ! ************************************************************************************************
 ! private subroutine vegLiqFlux_muster: compute water balance for the vegetation canopy
 ! ************************************************************************************************
 subroutine vegLiqFlux_muster(&
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
 err=0; message="vegLiqFlux_muster/"

 ! set throughfall to inputs if vegetation is completely buried with snow
 if(.not.computeVegFlux)then
  scalarThroughfallRain        = scalarRainfall
  scalarCanopyLiqDrainage      = 0._dp
  scalarCanopyLiqDrainageDeriv = 0._dp
  return
 endif

 ! set throughfall to zero (throughfall only used where there is no canopy)
 scalarThroughfallRain = 0._dp

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

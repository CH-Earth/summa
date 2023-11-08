! SUMMA - Structure for Unifying Multiple Modeling Alternatives
! Copyright (C) 2014-2020 NCAR/RAL; University of Saskatchewan; University of Washington
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

! data types
USE nrtype

! data types
USE data_types,only:var_d                ! x%var(:)       (rkind)
USE data_types,only:var_dlength          ! x%var(:)%dat   (rkind)
USE data_types,only:in_type_vegLiqFlux   ! class type for intent(in) arguments
USE data_types,only:out_type_vegLiqFlux  ! class type for intent(out) arguments

! named variables
USE var_lookup,only:iLookPARAM,iLookDIAG ! named variables for structure elements

! model decisions
USE globalData,only:model_decisions      ! model decision structure
USE var_lookup,only:iLookDECISIONS       ! named variables for elements of the decision structure

! decisions on canopy interception parameterization
USE mDecisions_module,only:         &
                      unDefined,    &    ! original model (no flexibility in canopy interception): 100% of rainfall is intercepted by the vegetation canopy
                      sparseCanopy, &    ! fraction of rainfall that never hits the canopy (throughfall); drainage above threshold
                      storageFunc        ! throughfall a function of canopy storage; 100% throughfall when canopy is at capacity

! privacy
implicit none
private
public :: vegLiqFlux
contains
! ************************************************************************************************
! public subroutine vegLiqFlux: compute water balance for the vegetation canopy
! ************************************************************************************************
subroutine vegLiqFlux(&
                    ! input:
                    in_vegLiqFlux,                & ! intent(in): model control, trial value, and rainfall rate
                    ! input-output: data structures
                    mpar_data,                    & ! intent(in): model parameters
                    diag_data,                    & ! intent(in): local HRU model diagnostic variables
                    ! output
                    out_vegLiqFlux)                 ! intent(out): output rates, derivatives, and error control
  implicit none
  ! input
  type(in_type_vegLiqFlux)        :: in_vegLiqFlux                  ! model control, trial value, and rainfall rate
  ! input-output: data structures
  type(var_dlength),intent(in)    :: mpar_data                      ! model parameters
  type(var_dlength),intent(inout) :: diag_data                      ! model diagnostic variables for the local basin
  ! output
  type(out_type_vegLiqFlux)       :: out_vegLiqFlux                 ! output rates, derivatives, and error control
  ! ------------------------------------------------------------------------------------------------------------------------------------------------------
  ! make association of local variables with information in the data structures
  associate(&
    computeVegFlux       => in_vegLiqFlux % computeVegFlux,       & ! intent(in): flag to indicate if we are computing fluxes over vegetation (.false. means veg is buried with snow)
    scalarCanopyLiqTrial => in_vegLiqFlux % scalarCanopyLiqTrial, & ! intent(in): trial mass of liquid water on the vegetation canopy at the current iteration (kg m-2)
    scalarRainfall       => in_vegLiqFlux % scalarRainfall,       & ! intent(in): rainfall (kg m-2 s-1)
    ixCanopyInterception       => model_decisions(iLookDECISIONS%cIntercept)%iDecision, & ! intent(in): index defining choice of parameterization for canopy interception
    scalarCanopyLiqMax         => diag_data%var(iLookDIAG%scalarCanopyLiqMax)%dat(1),   & ! intent(in): maximum storage before canopy drainage begins (kg m-2 s-1)
    scalarThroughfallScaleRain => mpar_data%var(iLookPARAM%throughfallScaleRain)%dat(1),& ! intent(in): fraction of rain that hits the ground without touching the canopy (-)
    scalarCanopyDrainageCoeff  => mpar_data%var(iLookPARAM%canopyDrainageCoeff)%dat(1),  & ! intent(in): canopy drainage coefficient (s-1)
    scalarThroughfallRain        => out_vegLiqFlux % scalarThroughfallRain,        & ! intent(out): rain that reaches the ground without ever touching the canopy (kg m-2 s-1)
    scalarCanopyLiqDrainage      => out_vegLiqFlux % scalarCanopyLiqDrainage,      & ! intent(out): drainage of liquid water from the vegetation canopy (kg m-2 s-1)
    scalarThroughfallRainDeriv   => out_vegLiqFlux % scalarThroughfallRainDeriv,   & ! intent(out): derivative in throughfall w.r.t. canopy liquid water (s-1)
    scalarCanopyLiqDrainageDeriv => out_vegLiqFlux % scalarCanopyLiqDrainageDeriv, & ! intent(out): derivative in canopy drainage w.r.t. canopy liquid water (s-1)
    err                          => out_vegLiqFlux % err,                          & ! intent(out): error code
    message                      => out_vegLiqFlux % cmessage                      & ! intent(out): error message
    )
    ! ------------------------------------------------------------------------------------------------------------------------------------------------------
    ! initialize error control
    err=0; message="vegLiqFlux/"

    ! set throughfall to inputs if vegetation is completely buried with snow
    if (.not.computeVegFlux) then
      scalarThroughfallRain        = scalarRainfall
      scalarCanopyLiqDrainage      = 0._rkind
      scalarThroughfallRainDeriv   = 0._rkind
      scalarCanopyLiqDrainageDeriv = 0._rkind
      return
    end if

    ! compute throughfall
    select case(ixCanopyInterception)
      ! original model (no flexibility in canopy interception): 100% of rainfall is intercepted by the vegetation canopy
      ! NOTE: this could be done with scalarThroughfallScaleRain=0, though requires setting scalarThroughfallScaleRain in all test cases
      case(unDefined)
        scalarThroughfallRain      = 0._rkind
        scalarThroughfallRainDeriv = 0._rkind
      ! fraction of rainfall hits the ground without ever touching the canopy
      case(sparseCanopy)
        scalarThroughfallRain      = scalarThroughfallScaleRain*scalarRainfall
        scalarThroughfallRainDeriv = 0._rkind
      ! throughfall a function of canopy storage
      case(storageFunc)
        ! throughfall during wetting-up phase
        if(scalarCanopyLiqTrial < scalarCanopyLiqMax)then
          scalarThroughfallRain      = scalarRainfall*(scalarCanopyLiqTrial/scalarCanopyLiqMax)
          scalarThroughfallRainDeriv = scalarRainfall/scalarCanopyLiqMax
        ! all rain falls through the canopy when the canopy is at capacity
        else
          scalarThroughfallRain      = scalarRainfall
          scalarThroughfallRainDeriv = 0._rkind
        end if
      case default; err=20; message=trim(message)//'unable to identify option for canopy interception'; return
    end select ! (option for canopy interception)

    ! compute canopy drainage
    if(scalarCanopyLiqTrial > scalarCanopyLiqMax)then
      scalarCanopyLiqDrainage       = scalarCanopyDrainageCoeff*(scalarCanopyLiqTrial - scalarCanopyLiqMax)
      scalarCanopyLiqDrainageDeriv  = scalarCanopyDrainageCoeff
    else
      scalarCanopyLiqDrainage       = 0._rkind
      scalarCanopyLiqDrainageDeriv  = 0._rkind
    end if

  end associate ! end association of local variables with information in the data structures

end subroutine vegLiqFlux

end module vegLiqFlux_module

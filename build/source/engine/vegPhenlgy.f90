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

module vegPhenlgy_module

! data types
USE nr_type

! global variables
USE globalData,only:&
                    realMissing,        & ! missing value for real numbers
                    urbanVegCategory,   & ! vegetation category for urban areas
                    minExpLogHgtFac       ! factor for minimum height of transition from the exponential to the logarithmic wind profile


! provide access to the derived types to define the data structures
USE data_types,only:&
                    var_i,            & ! data vector (i4b)
                    var_d,            & ! data vector (rkind)
                    var_dlength,      & ! data vector with variable length dimension (rkind)
                    model_options       ! defines the model decisions

! named variables defining elements in the data structures
USE var_lookup,only:iLookTYPE,iLookATTR,iLookPARAM,iLookDIAG,iLookPROG  ! named variables for structure elements
USE var_lookup,only:iLookDECISIONS                                      ! named variables for elements of the decision structure

! look-up values for the boundary conditions
USE mDecisions_module,only:      &
 prescribedHead,                 &      ! prescribed head (volumetric liquid water content for mixed form of Richards' eqn)
 prescribedTemp,                 &      ! prescribed temperature
 zeroFlux                               ! zero flux

! look-up values for the choice of canopy shortwave radiation method
USE mDecisions_module,only:      &
 noah_mp,                        &      ! full Noah-MP implementation (including albedo)
 CLM_2stream,                    &      ! CLM 2-stream model (see CLM documentation)
 UEB_2stream,                    &      ! UEB 2-stream model (Mahat and Tarboton, WRR 2011)
 NL_scatter,                     &      ! Simplified method Nijssen and Lettenmaier (JGR 1999)
 BeersLaw                               ! Beer's Law (as implemented in VIC)

! privacy
implicit none
private
public::vegPhenlgy
contains


 ! ************************************************************************************************
 ! public subroutine vegPhenlgy: compute vegetation phenology
 ! ************************************************************************************************
 subroutine vegPhenlgy(&
                       ! model control
                       nSnow,                       & ! intent(in):    number of snow layers
                       model_decisions,             & ! intent(in):    model decisions
                       fracJulDay,                  & ! intent(in):    fractional julian days since the start of year
                       yearLength,                  & ! intent(in):    number of days in the current year
                       noVeg,                       & ! intent(in):    flag to indicate if there is no vegetation (lake or glacier)
                       ! input/output: data structures
                       type_data,                   & ! intent(in):    type of vegetation and soil
                       attr_data,                   & ! intent(in):    spatial attributes
                       mpar_data,                   & ! intent(in):    model parameters
                       prog_data,                   & ! intent(inout): prognostic variables for a local HRU
                       diag_data,                   & ! intent(inout): diagnostic variables for a local HRU
                       ! output
                       computeVegFlux,              & ! intent(out): flag to indicate if we are computing fluxes over vegetation (.false. means veg is buried with snow)
                       canopyDepth,                 & ! intent(out): canopy depth (m)
                       exposedVAI,                  & ! intent(out): exposed vegetation area index (LAI + SAI)
                       err,message)                   ! intent(out): error control

 ! -------------------------------------------------------------------------------------------------
 ! modules
 USE NOAHMP_ROUTINES,only:phenology         ! determine vegetation phenology
 implicit none
 ! -------------------------------------------------------------------------------------------------
 ! input/output
 integer(i4b),intent(in)         :: nSnow               ! number of snow layers
 type(model_options),intent(in)  :: model_decisions(:)  ! model decisions
 real(rkind),intent(in)          :: fracJulDay          ! fractional julian days since the start of year
 integer(i4b),intent(in)         :: yearLength          ! number of days in the current year
 logical(lgt),intent(in)         :: noVeg               ! flag to indicate if there is no vegetation (lake or glacier)
 type(var_i),intent(in)          :: type_data           ! type of vegetation and soil
 type(var_d),intent(in)          :: attr_data           ! spatial attributes
 type(var_dlength),intent(in)    :: mpar_data           ! model parameters
 type(var_dlength),intent(inout) :: prog_data           ! prognostic variables for a local HRU
 type(var_dlength),intent(inout) :: diag_data           ! diagnostic variables for a local HRU
 ! output
 logical(lgt),intent(out)        :: computeVegFlux      ! flag to indicate if we are computing fluxes over vegetation (.false. means veg is buried with snow)
 real(rkind),intent(out)         :: canopyDepth         ! canopy depth (m)
 real(rkind),intent(out)         :: exposedVAI          ! exposed vegetation area index (LAI + SAI)
 integer(i4b),intent(out)        :: err                 ! error code
 character(*),intent(out)        :: message             ! error message
 ! -------------------------------------------------------------------------------------------------
 ! local
 real(rkind)                     :: z0Ground                   ! roughness length of the ground (ground below the canopy or non-vegetated surface) (m)
 real(rkind)                     :: notUsed_heightCanopyTop    ! height of the top of the canopy layer (m)
 real(rkind)                     :: heightAboveSnow            ! height top of canopy is above the snow surface (m)
 real(rkind)                     :: minExpLogHgt               ! minimum height above ground for logarithmic wind profile (m)

 ! initialize error control
 err=0; message="vegPhenlgy/"
 ! ----------------------------------------------------------------------------------------------------------------------------------
 ! associate variables in the data structure
 associate(&
 ! input: model decisions
 ix_bcUpprTdyn                   => model_decisions(iLookDECISIONS%bcUpprTdyn)%iDecision,      & ! intent(in): [i4b] choice of upper boundary condition for thermodynamics
 ix_bcUpprSoiH                   => model_decisions(iLookDECISIONS%bcUpprSoiH)%iDecision,      & ! intent(in): [i4b] index of method used for the upper boundary condition for soil hydrology
 ! local attributes
 vegTypeIndex                    => type_data%var(iLookTYPE%vegTypeIndex),                     & ! intent(in): [i4b] vegetation type index
 latitude                        => attr_data%var(iLookATTR%latitude),                         & ! intent(in): [dp] latitude
 ! model state variables
 scalarSnowDepth                 => prog_data%var(iLookPROG%scalarSnowDepth)%dat(1),           & ! intent(in):    [dp] snow depth on the ground surface (m)
 scalarCanopyTemp                => prog_data%var(iLookPROG%scalarCanopyTemp)%dat(1),          & ! intent(in):    [dp] temperature of the vegetation canopy at the start of the sub-step (K)
 ! diagnostic variables and parameters (input)
 z0Snow                          => mpar_data%var(iLookPARAM%z0Snow)%dat(1),                   & ! intent(in): [dp] roughness length of snow (m)
 z0Soil                          => mpar_data%var(iLookPARAM%z0Soil)%dat(1),                   & ! intent(in): [dp] roughness length of soil (m)
 heightCanopyTop                 => mpar_data%var(iLookPARAM%heightCanopyTop)%dat(1),          & ! intent(in): [dp] height of the top of the canopy layer (m)
 heightCanopyBottom              => mpar_data%var(iLookPARAM%heightCanopyBottom)%dat(1),       & ! intent(in): [dp] height of the bottom of the canopy layer (m)
 ! diagnostic variables and parameters (input/output)
 scalarLAI                       => diag_data%var(iLookDIAG%scalarLAI)%dat(1),                 & ! intent(inout): [dp] one-sided leaf area index (m2 m-2)
 scalarSAI                       => diag_data%var(iLookDIAG%scalarSAI)%dat(1),                 & ! intent(inout): [dp] one-sided stem area index (m2 m-2)
 ! diagnostic variables and parameters (output)
 scalarExposedLAI                => diag_data%var(iLookDIAG%scalarExposedLAI)%dat(1),          & ! intent(out): [dp] exposed leaf area index after burial by snow (m2 m-2)
 scalarExposedSAI                => diag_data%var(iLookDIAG%scalarExposedSAI)%dat(1),          & ! intent(out): [dp] exposed stem area index after burial by snow (m2 m-2)
 scalarGrowingSeasonIndex        => diag_data%var(iLookDIAG%scalarGrowingSeasonIndex)%dat(1),  & ! intent(out): [dp] growing season index (0=off, 1=on)
 scalarGroundSnowFraction        => diag_data%var(iLookDIAG%scalarGroundSnowFraction)%dat(1)   & ! intent(out): [dp] fraction of ground covered with snow (-)

 ) ! associate variables in data structure
 ! ----------------------------------------------------------------------------------------------------------------------------------
  if (nSnow>0) then ! case when there is snow on the ground (EXCLUDE "snow without a layer" -- in this case, evaporate from the soil)
    scalarGroundSnowFraction  = 1._rkind
  else ! case when the ground is less than a layer of snow (e.g., bare soil or snow without a layer)
    scalarGroundSnowFraction  = 0._rkind
  end if  ! (there is snow enough for a layer on the ground)
  
 ! check if we are on non-upland domain (noVeg)
 if(noVeg)then

  ! we are on non-upland domain (noVeg), no vegetation: do not compute fluxes over vegetation
   computeVegFlux           = .false. 

   ! set vegetation phenology variables to zero (no vegetation)
   scalarLAI                = 0._rkind    ! one-sided leaf area index (m2 m-2)
   scalarSAI                = 0._rkind    ! one-sided stem area index (m2 m-2)
   scalarExposedLAI         = 0._rkind    ! exposed leaf area index after burial by snow (m2 m-2)
   scalarExposedSAI         = 0._rkind    ! exposed stem area index after burial by snow (m2 m-2)
   scalarGrowingSeasonIndex = 0._rkind    ! growing season index (0=off, 1=on)
   exposedVAI               = 0._rkind    ! exposed vegetation area index (m2 m-2)
   canopyDepth              = 0._rkind    ! canopy depth (m)
   heightAboveSnow          = 0._rkind    ! height top of canopy is above the snow surface (m)

 ! check if we have isolated the snow-soil domain (used in test cases)
 elseif(ix_bcUpprTdyn == prescribedTemp .or. ix_bcUpprTdyn == zeroFlux .or. ix_bcUpprSoiH == prescribedHead) then

   ! isolated snow-soil domain: do not compute fluxes over vegetation
   computeVegFlux = .false.

   ! set vegetation phenology variables to missing
   scalarLAI                = realMissing    ! one-sided leaf area index (m2 m-2)
   scalarSAI                = realMissing    ! one-sided stem area index (m2 m-2)
   scalarExposedLAI         = realMissing    ! exposed leaf area index after burial by snow (m2 m-2)
   scalarExposedSAI         = realMissing    ! exposed stem area index after burial by snow (m2 m-2)
   scalarGrowingSeasonIndex = realMissing    ! growing season index (0=off, 1=on)
   exposedVAI               = realMissing    ! exposed vegetation area index (m2 m-2)
   canopyDepth              = realMissing    ! canopy depth (m)
   heightAboveSnow          = realMissing    ! height top of canopy is above the snow surface (m)

 ! determine vegetation phenology
 ! NOTE: recomputing phenology every sub-step accounts for changes in exposed vegetation associated with changes in snow depth
 else
   call phenology(&
                 ! input
                 vegTypeIndex,                & ! intent(in): vegetation type index
                 urbanVegCategory,            & ! intent(in): vegetation category for urban areas
                 scalarSnowDepth,             & ! intent(in): snow depth on the ground surface (m)
                 scalarCanopyTemp,            & ! intent(in): temperature of the vegetation canopy at the start of the sub-step (K)
                 latitude,                    & ! intent(in): latitude
                 yearLength,                  & ! intent(in): number of days in the current year
                 fracJulDay,                  & ! intent(in): fractional julian days since the start of year
                 scalarLAI,                   & ! intent(inout): one-sided leaf area index (m2 m-2)
                 scalarSAI,                   & ! intent(inout): one-sided stem area index (m2 m-2)
                 ! output
                 notUsed_heightCanopyTop,     & ! intent(out): height of the top of the canopy layer (m)
                 scalarExposedLAI,            & ! intent(out): exposed leaf area index after burial by snow (m2 m-2)
                 scalarExposedSAI,            & ! intent(out): exposed stem area index after burial by snow (m2 m-2)
                 scalarGrowingSeasonIndex     ) ! intent(out): growing season index (0=off, 1=on)

  ! determine additional phenological variables
  exposedVAI      = scalarExposedLAI + scalarExposedSAI   ! exposed vegetation area index (m2 m-2)
  canopyDepth     = heightCanopyTop - heightCanopyBottom  ! canopy depth (m)
  heightAboveSnow = heightCanopyTop - scalarSnowDepth     ! height top of canopy is above the snow surface (m)

  ! compute the roughness length of the ground (ground below the canopy or non-vegetated surface)
  z0Ground = z0Soil*(1._rkind - scalarGroundSnowFraction) + z0Snow*scalarGroundSnowFraction     ! roughness length (m)

  ! determine if need to include vegetation in the energy flux routines
  minExpLogHgt = minExpLogHgtFac*sqrt(heightCanopyTop) ! minimum height above ground for logarithmic wind profile (m)
  computeVegFlux = (exposedVAI > 0.05_rkind .and. heightAboveSnow > z0Ground + minExpLogHgt) ! check for complete burial of vegetatio

 end if  ! (check if the snow-soil column is isolated)

 ! end association to variables in the data structure
 end associate

 end subroutine vegPhenlgy


end module vegPhenlgy_module

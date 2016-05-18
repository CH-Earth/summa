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

module indexState_module
! data types
USE nrtype
! named variables that define the layer type
USE globalData,only:ix_soil         ! soil
USE globalData,only:ix_snow         ! snow
! named variables to describe the state variable type
USE globalData,only:iname_nrgCanair ! named variable defining the energy of the canopy air space
USE globalData,only:iname_nrgCanopy ! named variable defining the energy of the vegetation canopy
USE globalData,only:iname_watCanopy ! named variable defining the mass of water on the vegetation canopy
USE globalData,only:iname_nrgLayer  ! named variable defining the energy state variable for snow+soil layers
USE globalData,only:iname_watLayer  ! named variable defining the total water state variable for snow+soil layers
USE globalData,only:iname_matLayer  ! named variable defining the matric head state variable for soil layers
! provide access to the derived types to define the data structures
USE data_types,only:var_ilength     ! data vector with variable length dimension (i4b)
! provide access to the metadata
USE globalData,only:indx_meta       ! metadata for the variables in the index structure
! provide access to indices that define elements of the data structures
USE var_lookup,only:iLookINDEX      ! named variables for structure elements
implicit none
private
public::indexState
! control parameters
integer(i4b),parameter :: missingInteger=-9999
contains


 ! **********************************************************************************************************
 ! public subroutine indexState: define list of indices for each state variable 
 ! **********************************************************************************************************
 subroutine indexstate(computeVegFlux,          & ! intent(in):    flag to denote if computing the vegetation flux
                       indx_data,               & ! intent(inout): indices defining model states and layers
                       err,message)               ! intent(out):   error control
 ! --------------------------------------------------------------------------------------------------------------------------------
 ! provide access to the numerical recipes utility modules
 USE nr_utility_module,only:arth                           ! creates a sequence of numbers (start, incr, n)
 ! --------------------------------------------------------------------------------------------------------------------------------
 ! input
 logical(lgt),intent(in)         :: computeVegFlux         ! flag to denote if computing the vegetation flux
 type(var_ilength),intent(inout) :: indx_data              ! indices defining model states and layers
 ! output: error control
 integer(i4b),intent(out)        :: err                    ! error code
 character(*),intent(out)        :: message                ! error message
 ! --------------------------------------------------------------------------------------------------------------------------------
 ! local variables
 integer(i4b),parameter          :: nVarSnowSoil=2         ! number of state variables in the snow and soil domain (energy and total water/matric head)
 integer(i4b)                    :: iVar                   ! index of variable within the index data structure
 integer(i4b)                    :: nVec                   ! number of elements in the index vector
 ! --------------------------------------------------------------------------------------------------------------------------------
 ! make association with variables in the data structures
 associate(&
 ! number of state variables of different type
 nVegNrg       => indx_data%var(iLookINDEX%nVegNrg)%dat(1)   , & ! number of energy state variables for vegetation
 nVegMass      => indx_data%var(iLookINDEX%nVegMass)%dat(1)  , & ! number of hydrology states for vegetation (mass of water)
 nVegState     => indx_data%var(iLookINDEX%nVegState)%dat(1) , & ! number of vegetation state variables
 nNrgState     => indx_data%var(iLookINDEX%nNrgState)%dat(1) , & ! number of energy state variables
 nWatState     => indx_data%var(iLookINDEX%nWatState)%dat(1) , & ! number of "total water" states (vol. total water content)
 nMatState     => indx_data%var(iLookINDEX%nMatState)%dat(1) , & ! number of matric head state variables
 nMassState    => indx_data%var(iLookINDEX%nMassState)%dat(1), & ! number of hydrology state variables (mass of water)
 nState        => indx_data%var(iLookINDEX%nState)%dat(1)    , & ! total number of model state variables
 ! number of model layers, and layer indices
 nSnow         => indx_data%var(iLookINDEX%nSnow)%dat(1)     , & ! number of snow layers
 nSoil         => indx_data%var(iLookINDEX%nSoil)%dat(1)     , & ! number of soil layers
 nLayers       => indx_data%var(iLookINDEX%nLayers)%dat(1)   , & ! total number of layers
 layerType     => indx_data%var(iLookINDEX%layerType)%dat    , & ! index defining type of layer (soil or snow)
 ! indices of model state variables
 ixCasNrg      => indx_data%var(iLookINDEX%ixCasNrg)%dat(1)  , & ! index of canopy air space energy state variable
 ixVegNrg      => indx_data%var(iLookINDEX%ixVegNrg)%dat(1)  , & ! index of canopy energy state variable
 ixVegWat      => indx_data%var(iLookINDEX%ixVegWat)%dat(1)  , & ! index of canopy hydrology state variable (mass)
 ixTopNrg      => indx_data%var(iLookINDEX%ixTopNrg)%dat(1)  , & ! index of upper-most energy state in the snow-soil subdomain
 ixTopWat      => indx_data%var(iLookINDEX%ixTopWat)%dat(1)  , & ! index of upper-most total water state in the snow-soil subdomain
 ixTopMat      => indx_data%var(iLookINDEX%ixTopMat)%dat(1)  , & ! index of upper-most matric head state in the soil subdomain
 ixSnowSoilNrg => indx_data%var(iLookINDEX%ixSnowSoilNrg)%dat, & ! indices for energy states in the snow-soil subdomain
 ixSnowSoilWat => indx_data%var(iLookINDEX%ixSnowSoilWat)%dat, & ! indices for total water states in the snow-soil subdomain
 ixSnowOnlyNrg => indx_data%var(iLookINDEX%ixSnowOnlyNrg)%dat, & ! indices for energy states in the snow subdomain
 ixSnowOnlyWat => indx_data%var(iLookINDEX%ixSnowOnlyWat)%dat, & ! indices for total water states in the snow subdomain
 ixSoilOnlyNrg => indx_data%var(iLookINDEX%ixSoilOnlyNrg)%dat, & ! indices for energy states in the soil subdomain
 ixSoilOnlyHyd => indx_data%var(iLookINDEX%ixSoilOnlyHyd)%dat, & ! indices for hydrology states in the soil subdomain
 ! type of model state variables
 ixSoilState   => indx_data%var(iLookINDEX%ixSoilState)%dat  , & ! list of indices for all soil layers
 ixLayerState  => indx_data%var(iLookINDEX%ixLayerState)%dat   & ! list of indices for all model layers
 ) ! association to variables in the data structures
 ! --------------------------------------------------------------------------------------------------------------------------------
 ! initialize error control
 err=0; message='indexState/'

 ! -----
 ! * define the number of state variables...
 ! -----------------------------------------

 ! define the number of vegetation state variables (defines position of snow-soil states in the state vector)
 if(computeVegFlux)then
  nVegNrg   = 2
  nVegMass  = 1
  nVegState = nVegNrg + nVegMass
 else
  nVegNrg   = 0
  nVegMass  = 0
  nVegState = 0
 endif

 ! define the number state variables of different type
 nNrgState  = nVegNrg + nLayers  ! number of energy state variables
 nWatState  = nSnow              ! number of "total water" state variables -- will be modified later if using primary variable switching 
 nMatState  = nSoil              ! number of matric head state variables -- will be modified later if using primary variable switching
 nMassState = nVegMass           ! number of mass state variables -- currently restricted to canopy water

 ! define the number of model state variables
 nState = nVegState + nLayers*nVarSnowSoil   ! *nVarSnowSoil (both energy and total water)

 ! -----
 ! * define the indices of state variables...
 ! ------------------------------------------

 if(computeVegFlux)then
  ixCasNrg  = 1
  ixVegNrg  = 2
  ixVegWat  = 3
 else
  ixCasNrg  = missingInteger
  ixVegNrg  = missingInteger
  ixVegWat  = missingInteger
 endif

 ! define the index of the top layer
 ixTopNrg = nVegState + 1                       ! energy
 ixTopWat = nVegState + 2                       ! total water (only snow)
 ixTopMat = nVegState + nSnow*nVarSnowSoil + 2  ! matric head (only soil)

 ! define the indices within the snow+soil domain
 ixSnowSoilNrg = arth(ixTopNrg,nVarSnowSoil,nLayers)  ! energy
 ixSnowSoilWat = arth(ixTopWat,nVarSnowSoil,nLayers)  ! total water

 ! define indices just for the soil domains
 ixSoilOnlyNrg = arth(ixTopNrg + nSnow*nVarSnowSoil,nVarSnowSoil,nSoil)    ! energy
 ixSoilOnlyHyd = arth(ixTopMat,nVarSnowSoil,nSoil)                         ! soil hydrology states (matric head or total water)

 ! define indices just for the snow domain
 if(nSnow>0)then  ! (total water in snow only defined if snow layers exist)
  ixSnowOnlyNrg = arth(ixTopNrg,nVarSnowSoil,nSnow)    ! energy
  ixSnowOnlyWat = arth(ixTopWat,nVarSnowSoil,nSnow)    ! total water
 endif

 ! -----
 ! * re-allocate index vectors (if needed)...
 ! ------------------------------------------

 ! loop through index variables
 do iVar=1,size(indx_data%var)

  ! get the length of the desired variable
  select case(iVar)
   case(iLookINDEX%ixStateType); nVec=nState
   case(iLookINDEX%ixAllState);  nVec=nState
   case(iLookINDEX%ixNrgOnly);   nVec=nNrgState
   case(iLookINDEX%ixWatOnly);   nVec=nWatState
   case(iLookINDEX%ixMatOnly);   nVec=nMatState
   case(iLookINDEX%ixMassOnly);  nVec=nMassState
   case default; cycle ! only need to process the above variables
  end select  ! iVar

  ! check if we need to reallocate space
  if(size(indx_data%var(iVar)%dat)==nVec) cycle

  ! deallocate space
  deallocate(indx_data%var(iVar)%dat,stat=err)
  if(err/=0)then
   message=trim(message)//'unable to deallocate space for variable '//trim(indx_meta(ivar)%varname)
   err=20; return
  endif

  ! allocate space
  allocate(indx_data%var(iVar)%dat(nVec),stat=err)
  if(err/=0)then
   message=trim(message)//'unable to allocate space for variable '//trim(indx_meta(ivar)%varname)
   err=20; return
  endif

 end do  ! looping through variables

 ! -----
 ! * define the type of model states...
 ! ------------------------------------

 ! make an association to the ALLOCATABLE variables in the data structures
 ! NOTE: we need to do this here since the size may have changed above
 associate(&
 ixStateType => indx_data%var(iLookINDEX%ixStateType)%dat , & ! indices defining the type of the state (iname_nrgLayer...)
 ixAllState  => indx_data%var(iLookINDEX%ixAllState)%dat  , & ! list of indices for all model state variables
 ixNrgOnly   => indx_data%var(iLookINDEX%ixNrgOnly)%dat   , & ! list of indices for all energy states
 ixWatOnly   => indx_data%var(iLookINDEX%ixWatOnly)%dat   , & ! list of indices for all "total water" states
 ixMatOnly   => indx_data%var(iLookINDEX%ixMatOnly)%dat   , & ! list of indices for matric head state variables
 ixMassOnly  => indx_data%var(iLookINDEX%ixMassOnly)%dat    & ! list of indices for hydrology states (mass of water)
 )  ! making an association to variables in the data structures

 ! define the state type for the vegetation canopy
 if(computeVegFlux)then
  ixStateType(ixCasNrg) = iname_nrgCanair
  ixStateType(ixVegNrg) = iname_nrgCanopy
  ixStateType(ixVegWat) = iname_watCanopy
 endif

 ! define the state type for the snow-soil domain
 ixStateType(ixSnowSoilNrg) = iname_nrgLayer
 ixStateType(ixSoilOnlyHyd) = iname_matLayer ! refine later to be either iname_watLayer or iname_matLayer
 if(nSnow>0) ixStateType(ixSnowOnlyWat) = iname_watLayer

 ! define indices for state variables
 ixAllState   = arth(1,1,nState)
 ixSoilState  = arth(1,1,nSoil)
 ixLayerState = arth(1,1,nLayers)

 ! define vector of indices for specific states
                  ixNrgOnly  = pack(ixAllState, ixStateType==iname_nrgCanair .or. ixStateType==iname_nrgCanopy .or. ixStateType==iname_nrgLayer)
 if(nWatState>0)  ixWatOnly  = pack(ixAllState, ixStateType==iname_watLayer)
 if(nMatState>0)  ixMatOnly  = pack(ixAllState, ixStateType==iname_matLayer)
 if(nMassState>0) ixMassOnly = pack(ixAllState, ixStateType==iname_watCanopy)

 ! end association to the ALLOCATABLE variables in the data structures
 end associate 

 ! --------------------------------------------------------------------------------------------------------------------------------
 ! --------------------------------------------------------------------------------------------------------------------------------

 end associate  ! end association to variables in the data structures
 end subroutine indexState


end module indexState_module

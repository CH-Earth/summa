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
! missing data
USE globalData,only:integerMissing  ! missing integer
! domain types
USE globalData,only:iname_cas       ! canopy air space
USE globalData,only:iname_veg       ! vegetation
USE globalData,only:iname_snow      ! snow
USE globalData,only:iname_soil      ! soil
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
public::resizeIndx
! control parameters
integer(i4b),parameter :: missingInteger=-9999
contains


 ! **********************************************************************************************************
 ! public subroutine indexState: define list of indices for each state variable 
 ! **********************************************************************************************************
 subroutine indexState(computeVegFlux,          & ! intent(in):    flag to denote if computing the vegetation flux
                       nSnow,nSoil,nLayers,     & ! intent(in):    number of snow and soil layers, and total number of layers
                       indx_data,               & ! intent(inout): indices defining model states and layers
                       err,message)               ! intent(out):   error control
 ! --------------------------------------------------------------------------------------------------------------------------------
 ! provide access to the numerical recipes utility modules
 USE nr_utility_module,only:arth                           ! creates a sequence of numbers (start, incr, n)
 ! --------------------------------------------------------------------------------------------------------------------------------
 ! input
 logical(lgt),intent(in)         :: computeVegFlux         ! flag to denote if computing the vegetation flux
 integer(i4b),intent(in)         :: nSnow,nSoil,nLayers    ! number of snow and soil layers, and total number of layers
 type(var_ilength),intent(inout) :: indx_data              ! indices defining model states and layers
 ! output: error control
 integer(i4b),intent(out)        :: err                    ! error code
 character(*),intent(out)        :: message                ! error message
 ! --------------------------------------------------------------------------------------------------------------------------------
 ! general local variables
 character(len=256)              :: cmessage               ! message of downwind routine
 integer(i4b),parameter          :: nVarSnowSoil=2         ! number of state variables in the snow and soil domain (energy and total water/matric head)
 ! indices of model state variables
 integer(i4b)                    :: ixCasNrg               ! index of canopy air space energy state variable
 integer(i4b)                    :: ixVegNrg               ! index of canopy energy state variable
 integer(i4b)                    :: ixVegWat               ! index of canopy hydrology state variable (mass)
 integer(i4b)                    :: ixTopNrg               ! index of upper-most energy state in the snow-soil subdomain
 integer(i4b)                    :: ixTopWat               ! index of upper-most total water state in the snow-soil subdomain
 integer(i4b)                    :: ixTopMat               ! index of upper-most matric head state in the soil subdomain
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
 ! vectors of indices for specfic state types within specific sub-domains IN THE FULL STATE VECTOR
 ixNrgLayer    => indx_data%var(iLookINDEX%ixNrgLayer)%dat   , & ! indices IN THE FULL VECTOR for energy states in the snow+soil domain
 ixHydLayer    => indx_data%var(iLookINDEX%ixHydLayer)%dat   , & ! indices IN THE FULL VECTOR for hydrology states in the snow+soil domain
 ! indices for model state variables
 ixSoilState   => indx_data%var(iLookINDEX%ixSoilState)%dat  , & ! list of indices for all soil layers
 ixLayerState  => indx_data%var(iLookINDEX%ixLayerState)%dat , & ! list of indices for all model layers
 ! type of hydrology state variables in the snow+soil domains
 ixHydType     => indx_data%var(iLookINDEX%ixHydType)%dat      & ! index of the type of hydrology states in snow+soil domain
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
 end if

 ! define the number state variables of different type
 nNrgState  = nVegNrg + nLayers  ! number of energy state variables
 nWatState  = nSnow              ! number of "total water" state variables -- will be modified later if using primary variable switching 
 nMatState  = nSoil              ! number of matric head state variables -- will be modified later if using primary variable switching
 nMassState = nVegMass           ! number of mass state variables -- currently restricted to canopy water

 ! define the number of model state variables
 nState = nVegState + nLayers*nVarSnowSoil   ! *nVarSnowSoil (both energy and total water)

 ! -----
 ! * define the indices of state variables WITHIN THE FULL STATE VECTOR...
 ! -----------------------------------------------------------------------

 ! NOTE: Local variables are used here, since the actual indices (in the data structures) depend on the selection of state variables

 ! define indices in the vegetation domain
 if(computeVegFlux)then
  ixCasNrg  = 1
  ixVegNrg  = 2
  ixVegWat  = 3
 else
  ixCasNrg  = missingInteger
  ixVegNrg  = missingInteger
  ixVegWat  = missingInteger
 end if

 ! define the index of the top layer
 ixTopNrg = nVegState + 1                       ! energy
 ixTopWat = nVegState + 2                       ! total water (only snow)
 ixTopMat = nVegState + nSnow*nVarSnowSoil + 2  ! matric head (only soil)

 ! define the indices within the snow+soil domain
 ixNrgLayer = arth(ixTopNrg,nVarSnowSoil,nLayers)  ! energy
 ixHydLayer = arth(ixTopWat,nVarSnowSoil,nLayers)  ! total water

 ! re-allocate index vectors (if needed)...
 call resizeIndx( (/iLookINDEX%ixMapFull2Subset, iLookINDEX%ixControlVolume, iLookINDEX%ixDomainType, iLookINDEX%ixStateType, iLookINDEX%ixAllState/), & ! desired variables
                  indx_data,  & ! data structure
                  nState,     & ! vector length
                  err,cmessage) ! error control
 if(err/=0)then; message=trim(message)//trim(cmessage); return; endif

 ! -----
 ! * define the type of model states...
 ! ------------------------------------

 ! make an association to the ALLOCATABLE variables in the data structures
 ! NOTE: we need to do this here since the size may have changed above
 associate(&
 ixControlVolume => indx_data%var(iLookINDEX%ixControlVolume)%dat , & ! index of control volume for different domains (veg, snow, soil)
 ixDomainType    => indx_data%var(iLookINDEX%ixDomainType)%dat    , & ! indices defining the type of the domain (iname_veg, iname_snow, iname_soil)
 ixStateType     => indx_data%var(iLookINDEX%ixStateType)%dat     , & ! indices defining the type of the state (iname_nrgLayer...)
 ixAllState      => indx_data%var(iLookINDEX%ixAllState)%dat        & ! list of indices for all model state variables
 )  ! making an association to variables in the data structures

 ! define indices for state variables
 ixAllState   = arth(1,1,nState)
 ixSoilState  = arth(1,1,nSoil)
 ixLayerState = arth(1,1,nLayers)

 ! define the state type for the vegetation canopy
 if(computeVegFlux)then
  ixStateType(ixCasNrg) = iname_nrgCanair
  ixStateType(ixVegNrg) = iname_nrgCanopy
  ixStateType(ixVegWat) = iname_watCanopy
 endif

 ! define the state type for the snow+soil domain (energy)
 ixStateType(ixNrgLayer) = iname_nrgLayer

 ! define the state type for the snow+soil domain (hydrology)
 if(nSnow>0) ixStateType( ixHydLayer(      1:nSnow)   ) = iname_watLayer
             ixStateType( ixHydLayer(nSnow+1:nLayers) ) = iname_matLayer ! refine later to be either iname_watLayer or iname_matLayer

 ! define the domain type for vegetation
 if(computeVegFlux)then
  ixDomainType(ixCasNrg) = iname_cas
  ixDomainType(ixVegNrg) = iname_veg
  ixDomainType(ixVegWat) = iname_veg
 endif

 ! define the domain type for snow
 if(nSnow>0)then
  ixDomainType( ixNrgLayer(1:nSnow) ) = iname_snow
  ixDomainType( ixHydLayer(1:nSnow) ) = iname_snow
 endif

 ! define the domain type for soil
 ixDomainType( ixNrgLayer(nSnow+1:nLayers) ) = iname_soil
 ixDomainType( ixHydLayer(nSnow+1:nLayers) ) = iname_soil

 ! define the index of each control volume in the vegetation domains
 if(computeVegFlux)then
  ixControlVolume(ixCasNrg) = 1
  ixControlVolume(ixVegNrg) = 1
  ixControlVolume(ixVegWat) = 1
 endif

 ! define the index of the each control volume in the snow domain
 if(nSnow>0)then
  ixControlVolume( ixNrgLayer(1:nSnow) ) = ixLayerState(1:nSnow)
  ixControlVolume( ixHydLayer(1:nSnow) ) = ixLayerState(1:nSnow)
 endif

 ! define the index of the each control volume in the soil domain
 ixControlVolume( ixNrgLayer(nSnow+1:nLayers) ) = ixSoilState(1:nSoil)
 ixControlVolume( ixHydLayer(nSnow+1:nLayers) ) = ixSoilState(1:nSoil)

 ! define the type of variable in the snow+soil domain
 ixHydType(1:nLayers) = ixStateType( ixHydLayer(1:nLayers) )

 ! end association to the ALLOCATABLE variables in the data structures
 end associate 

 ! --------------------------------------------------------------------------------------------------------------------------------
 ! --------------------------------------------------------------------------------------------------------------------------------

 end associate  ! end association to variables in the data structures
 end subroutine indexState

 ! **********************************************************************************************************
 ! public subroutine resizeIndx: re-size specific index vectors 
 ! **********************************************************************************************************
 subroutine resizeIndx(ixDesire,indx_data,nVec,err,message)
 ! input
 integer(i4b)     ,intent(in)    :: ixDesire(:)            ! variables needing to be re-sized
 type(var_ilength),intent(inout) :: indx_data              ! indices defining model states and layers
 integer(i4b)     ,intent(in)    :: nVec                   ! desired vector length 
 ! output
 integer(i4b)     ,intent(out)   :: err                    ! error code
 character(*)     ,intent(out)   :: message                ! error message 
 ! local variables
 integer(i4b)                    :: jVar,iVar              ! vatiable index
 ! initialize error control
 err=0; message='resizeIndx/'

 ! loop through variables
 do jVar=1,size(ixDesire)

  ! define index in index array
  iVar = ixDesire(jVar)

  ! check iVar is within range
  if(iVar<1 .or. iVar>size(indx_data%var))then
   message=trim(message)//'desired variable is out of range'
   err=20; return
  endif

  ! check if we need to reallocate space
  if(size(indx_data%var(iVar)%dat) == nVec) cycle

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

  ! set to missing
  indx_data%var(iVar)%dat = integerMissing

 end do  ! looping through variables

 end subroutine resizeIndx

end module indexState_module

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

module getVectorz_module

! data types
USE nrtype

! missing values
USE globalData,only:integerMissing  ! missing integer
USE globalData,only:realMissing     ! missing real number

! layer types
USE globalData,only:ix_soil,ix_snow ! named variables for snow and soil

! named variables to describe the state variable type
USE globalData,only:iname_nrgCanair ! named variable defining the energy of the canopy air space
USE globalData,only:iname_nrgCanopy ! named variable defining the energy of the vegetation canopy
USE globalData,only:iname_watCanopy ! named variable defining the mass of water on the vegetation canopy
USE globalData,only:iname_nrgLayer  ! named variable defining the energy state variable for snow+soil layers
USE globalData,only:iname_watLayer  ! named variable defining the total water state variable for snow+soil layers
USE globalData,only:iname_matLayer  ! named variable defining the matric head state variable for soil layers

! constants
USE multiconst,only:&
                    Cp_air,       & ! specific heat of air                 (J kg-1 K-1)
                    iden_air,     & ! intrinsic density of air             (kg m-3)
                    iden_ice,     & ! intrinsic density of ice             (kg m-3)
                    iden_water      ! intrinsic density of liquid water    (kg m-3)

! provide access to the derived types to define the data structures
USE data_types,only:&
                    var_i,        & ! data vector (i4b)
                    var_d,        & ! data vector (dp)
                    var_ilength,  & ! data vector with variable length dimension (i4b)
                    var_dlength     ! data vector with variable length dimension (dp)

! provide access to indices that define elements of the data structures
USE var_lookup,only:iLookPROG       ! named variables for structure elements
USE var_lookup,only:iLookDIAG       ! named variables for structure elements
USE var_lookup,only:iLookPARAM      ! named variables for structure elements
USE var_lookup,only:iLookINDEX      ! named variables for structure elements

implicit none
private
public::popStateVec
public::varExtract

! common variables
real(dp),parameter :: valueMissing=-9999._dp ! missing value

contains


 ! **********************************************************************************************************
 ! public subroutine popStateVec: populate model state vectors 
 ! **********************************************************************************************************
 subroutine popStateVec(&
                        ! input: model control
                        computeVegFlux,          & ! intent(in):    flag to denote if computing energy flux over vegetation
                        stateType,               & ! intent(in):    type of desired model state variables
                        nState,                  & ! intent(in):    number of desired state variables
                        ! input-output: data structures
                        prog_data,               & ! intent(in):    model prognostic variables for a local HRU
                        diag_data,               & ! intent(in):    model diagnostic variables for a local HRU
                        indx_data,               & ! intent(inout): indices defining model states and layers
                        ! output
                        stateVec,                & ! intent(out):   model state vector
                        fScale,                  & ! intent(out):   function scaling vector (mixed units)
                        xScale,                  & ! intent(out):   variable scaling vector (mixed units)
                        sMul,                    & ! intent(out):   multiplier for state vector (used in the residual calculations)
                        dMat,                    & ! intent(out):   diagonal of the Jacobian matrix (excludes fluxes) 
                        err,message)               ! intent(out):   error control
 ! --------------------------------------------------------------------------------------------------------------------------------
 ! --------------------------------------------------------------------------------------------------------------------------------
 ! input: model control
 logical(lgt),intent(in)         :: computeVegFlux         ! flag to indicate if we are computing fluxes over vegetation (.false. means veg is buried with snow)
 integer(i4b),intent(in)         :: stateType(:)           ! type of desired model state variables
 integer(i4b),intent(in)         :: nState                 ! number of desired state variables
 ! input-output: data structures
 type(var_dlength),intent(in)    :: prog_data              ! prognostic variables for a local HRU
 type(var_dlength),intent(in)    :: diag_data              ! diagnostic variables for a local HRU
 type(var_ilength),intent(inout) :: indx_data              ! indices defining model states and layers
 ! output: state vectors
 real(dp),intent(out)            :: stateVec(:)            ! model state vector (mixed units)
 real(dp),intent(out)            :: fScale(:)              ! function scaling vector (mixed units)
 real(dp),intent(out)            :: xScale(:)              ! variable scaling vector (mixed units)
 real(qp),intent(out)            :: sMul(:)    ! NOTE: qp  ! multiplier for state vector (used in the residual calculations)
 real(dp),intent(out)            :: dMat(:)                ! diagonal of the Jacobian matrix (excludes fluxes)
 ! output: error control
 integer(i4b),intent(out)        :: err                    ! error code
 character(*),intent(out)        :: message                ! error message
 ! --------------------------------------------------------------------------------------------------------------------------------
 ! local variables
 ! --------------------------------------------------------------------------------------------------------------------------------
 ! scaling parameters
 real(dp),parameter              :: fScaleLiq=0.01_dp      ! func eval: characteristic scale for volumetric liquid water content (-)
 real(dp),parameter              :: fScaleMat=10._dp       ! func eval: characteristic scale for matric head (m)
 real(dp),parameter              :: fScaleNrg=1000000._dp  ! func eval: characteristic scale for energy (J m-3)
 real(dp),parameter              :: xScaleLiq=0.1_dp       ! state var: characteristic scale for volumetric liquid water content (-)
 real(dp),parameter              :: xScaleMat=10._dp       ! state var: characteristic scale for matric head (m)
 real(dp),parameter              :: xScaleTemp=1._dp       ! state var: characteristic scale for temperature (K)
 ! temporary vectors
 integer(i4b),dimension(nState)  :: ixSequence             ! sequential index in model state vector
 ! --------------------------------------------------------------------------------------------------------------------------------
 ! --------------------------------------------------------------------------------------------------------------------------------
 ! make association with variables in the data structures
 associate(&
 ! model states for the vegetation canopy
 scalarCanairTemp  => prog_data%var(iLookPROG%scalarCanairTemp)%dat(1)       ,& ! intent(in): [dp] temperature of the canopy air space (K)
 scalarCanopyTemp  => prog_data%var(iLookPROG%scalarCanopyTemp)%dat(1)       ,& ! intent(in): [dp] temperature of the vegetation canopy (K)
 scalarCanopyWat   => prog_data%var(iLookPROG%scalarCanopyWat)%dat(1)        ,& ! intent(in): [dp] mass of total water on the vegetation canopy (kg m-2)
 ! model state variable vectors for the snow-soil layers
 mLayerTemp        => prog_data%var(iLookPROG%mLayerTemp)%dat                ,& ! intent(in): [dp(:)] temperature of each snow/soil layer (K)
 mLayerVolFracWat  => prog_data%var(iLookPROG%mLayerVolFracWat)%dat          ,& ! intent(in): [dp(:)] volumetric fraction of total water (-)
 mLayerMatricHead  => prog_data%var(iLookPROG%mLayerMatricHead)%dat          ,& ! intent(in): [dp(:)] matric head (m)
 ! model diagnostic variables
 canopyDepth       => diag_data%var(iLookDIAG%scalarCanopyDepth)%dat(1)      ,&  ! intent(in): [dp   ] canopy depth (m)
 volHeatCapVeg     => diag_data%var(iLookDIAG%scalarBulkVolHeatCapVeg)%dat(1),& ! intent(in): [dp   ] bulk volumetric heat capacity of vegetation (J m-3 K-1)
 mLayerVolHeatCap  => diag_data%var(iLookDIAG%mLayerVolHtCapBulk)%dat        ,& ! intent(in): [dp(:)] bulk volumetric heat capacity in each snow and soil layer (J m-3 K-1)
 ! indices defining model states and layers
 ixCasNrg          => indx_data%var(iLookINDEX%ixCasNrg)%dat                 ,& ! intent(in): [i4b] index of canopy air space energy state variable
 ixVegNrg          => indx_data%var(iLookINDEX%ixVegNrg)%dat                 ,& ! intent(in): [i4b] index of canopy energy state variable
 ixVegWat          => indx_data%var(iLookINDEX%ixVegWat)%dat                 ,& ! intent(in): [i4b] index of canopy hydrology state variable (mass)
 ixTopNrg          => indx_data%var(iLookINDEX%ixTopNrg)%dat                 ,& ! intent(in): [i4b] index of upper-most energy state in the snow-soil subdomain
 ixTopWat          => indx_data%var(iLookINDEX%ixTopWat)%dat                 ,& ! intent(in): [i4b] index of upper-most total water state in the snow-soil subdomain
 ixTopMat          => indx_data%var(iLookINDEX%ixTopMat)%dat                 ,& ! intent(in): [i4b] index of upper-most matric head state in the soil subdomain



 ixSnowSoilNrg     => indx_data%var(iLookINDEX%ixSnowSoilNrg)%dat            ,& ! intent(in): [i4b(:)] indices for energy states in the snow-soil subdomain
 ixSnowOnlyWat     => indx_data%var(iLookINDEX%ixSnowOnlyWat)%dat            ,& ! intent(in): [i4b(:)] indices for total water states in the snow subdomain
 ixSnowSoilWat     => indx_data%var(iLookINDEX%ixSnowSoilWat)%dat            ,& ! intent(in): [i4b(:)] indices for total water states in the snow-soil subdomain
 ixSoilOnlyHyd     => indx_data%var(iLookINDEX%ixSoilOnlyHyd)%dat            ,& ! intent(in): [i4b(:)] indices for hydrology states in the soil subdomain
 ixNrgOnly         => indx_data%var(iLookINDEX%ixNrgOnly)%dat                ,& ! intent(in): [i4b(:)] list of indices for all energy states
 ixWatOnly         => indx_data%var(iLookINDEX%ixWatOnly)%dat                ,& ! intent(in): [i4b(:)] list of indices for all "total water" states
 ixMatOnly         => indx_data%var(iLookINDEX%ixMatOnly)%dat                ,& ! intent(in): [i4b(:)] list of indices for matric head state variables
 ixMassOnly        => indx_data%var(iLookINDEX%ixMassOnly)%dat               ,& ! intent(in): [i4b(:)] list of indices for hydrology states (mass of water)
 ! number of different states
 nSnow             => indx_data%var(iLookINDEX%nSnow)%dat(1)                 ,& ! intent(in): [i4b] number of snow layers
 nSoil             => indx_data%var(iLookINDEX%nSoil)%dat(1)                 ,& ! intent(in): [i4b] number of soil layers
 nLayers           => indx_data%var(iLookINDEX%nLayers)%dat(1)               ,& ! intent(in): [i4b] total number of layers
 nNrgState         => indx_data%var(iLookINDEX%nNrgState)%dat(1)             ,& ! intent(in): [i4b] number of energy state variables
 nWatState         => indx_data%var(iLookINDEX%nWatState)%dat(1)             ,& ! intent(in): [i4b] number of "total water" states (vol. total water content)
 nMatState         => indx_data%var(iLookINDEX%nMatState)%dat(1)             ,& ! intent(in): [i4b] number of matric head state variables
 nMassState        => indx_data%var(iLookINDEX%nMassState)%dat(1)             & ! intent(in): [i4b] number of hydrology state variables (mass of water)
 )  ! end association with variables in the data structures
 ! --------------------------------------------------------------------------------------------------------------------------------
 ! initialize error control
 err=0; message='popStateVec/'

 ! -----
 ! * identify position in the state vector for specific variables...
 ! -----------------------------------------------------------------

 print*, 'stateType = ', stateType 
 print*, 'merge = ', merge(1, 0, stateType==iname_watCanopy)

 ! initialize indices
 ixCasNrg = integerMissing ! energy in the canopy air space
 ixVegNrg = integerMissing ! energy of the vegetation canopy
 ixVegWat = integerMissing ! mass of water in the vegetation canopy
 ixTopNrg = integerMissing ! upper-most energy state variable in the snow+soil domain
 ixTopWat = integerMissing ! upper-most hydrology state variable in the snow+soil domain
 ixTopMat = integerMissing ! upper-most matric head state variable in the snow+soil domain

 ! get the number of state variables
 nCasNrg = count(stateType==iname_nrgCanair)
 nVegNrg = count(stateType==iname_nrgCanopy)
 nVegWat = count(stateType==iname_watCanopy)
 nLayWat = count(stateType==iname_watLayer)
 nLayMat = count(stateType==iname_matLayer)

 ! define indices for vegetation variables
 ! NOTE: workaround for (not-yet-implemented) f2008 intrinsic findloc -- merge provides a vector with 1s where mask is true and 0s otherwise
 if(nCasNrg==1) ixCasNrg = maxloc( merge(1, 0, stateType==iname_nrgCanair) )
 if(nVegNrg==1) ixVegNrg = maxloc( merge(1, 0, stateType==iname_nrgCanopy) )
 if(nVegWat==1) ixVegWat = maxloc( merge(1, 0, stateType==iname_watCanopy) )

 ! define index for the upper-most energy state variable in the snow+soil domain
 ! NOTE: workaround for (not-yet-implemented) f2008 intrinsic findloc -- merge provides a vector with 1s where mask is true and 0s otherwise
 ixTopNrg = maxloc( merge(1, 0, stateType==iname_nrgLayer) )

 ! define index for the upper-most matric head state variable
 ! NOTE: workaround for (not-yet-implemented) f2008 intrinsic findloc -- merge provides a vector with 1s where mask is true and 0s otherwise
 if(nLayMat>0) ixTopMat = maxloc( merge(1, 0, stateType==iname_matLayer) )

 ! define index for the upper-most water state variable
 ! NOTE: workaround for (not-yet-implemented) f2008 intrinsic findloc -- merge provides a vector with 1s where mask is true and 0s otherwise
 if(nLayWat>0)then
  ixTopWat = maxloc( merge(1, 0, stateType==iname_watLayer) )
 else
  ixTopWat = ixTopMat ! no water state variables, then upper-most water variable is the upper-most matric head variable
 endif

 print*, 'ixCasNrg = ', ixCasNrg
 print*, 'ixVegNrg = ', ixVegNrg
 print*, 'ixVegWat = ', ixVegWat

 ! -----
 ! * initialize state vectors...
 ! -----------------------------

 ! define index in full state vector
 ixSequence = arth(1,1,nState)

 ! build elements of the state vector for the vegetation canopy
 if(ixCasNrg(1)/=integerMissing) stateVec(ixCasNrg) = scalarCanairTemp
 if(ixVegNrg(1)/=integerMissing) stateVec(ixVegNrg) = scalarCanopyTemp
 if(ixVegWat(1)/=integerMissing) stateVec(ixVegWat) = scalarCanopyWat  ! kg m-2

 ! build the state vector for the snow and soil domain
 if(nLayNrg>0) stateVec( pack(ixSequence, stateType==iname_nrgLayer) ) = mLayerTemp(1:nLayers)
 if(nLayWat>0) stateVec( pack(ixSequence, stateType==iname_watLayer) ) = mLayerVolFracWat(1:nSnow)
 if(nLayMat>0) stateVec( pack(ixSequence, stateType==iname_matLayer) ) = mLayerMatricHead(1:nSoil)

 ! -----
 ! * define scaling vectors...
 ! ---------------------------

 ! define the function and variable scaling factors for energy
 where(stateType==iname_nrgCanair .or. stateType==iname_nrgCanopy .or. stateType==iname_nrgLayer)
  fScale = 1._dp / fScaleNrg  ! 1/(J m-3)
  xScale = 1._dp  ! K
 endwhere

 ! define the function and variable scaling factors for water on the vegetation canopy
 where(stateType==iname_watCanopy)
  fScale = 1._dp / (fScaleLiq*canopyDepth*iden_water)  ! 1/(kg m-2)
  xScale = 1._dp  ! (kg m-2)
 endwhere

 ! define the function and variable scaling factors for water in the snow+soil domain
 where(stateType==iname_watLayer)
  fScale = 1._dp / fScaleLiq  ! (-)
  xScale = 1._dp  ! (-)
 end where
 
 ! define the function and variable scaling factors for water in the snow+soil domain
 where(stateType==iname_matLayer)
  fScale = 1._dp / fScaleLiq  ! (-)
  xScale = 1._dp  ! (m)
 end where

 ! -----
 ! * define components of derivative matrices that are constant over a time step (substep)...
 ! ------------------------------------------------------------------------------------------

 ! define additional vectors used in the residual calculations
 sMul(:) = 1._dp         ! multiplier for the state vector
 dMat(:) = valueMissing  ! diagonal of the Jacobian matrix (excludes fluxes) 

 ! define the multiplier for the state vector for residual calculations (vegetation canopy)
 if(computeVegFlux)then
  sMul(ixCasNrg) = Cp_air*iden_air   ! volumetric heat capacity of air (J m-3 K-1)
  sMul(ixVegNrg) = volHeatCapVeg     ! volumetric heat capacity of the vegetation (J m-3 K-1)
  sMul(ixVegWat) = 1._dp             ! nothing else on the left hand side
 endif

 ! define the multiplier for the state vector for residual calculations (snow-soil domain)
 sMul(ixSnowSoilNrg) = mLayerVolHeatCap(1:nLayers)
 sMul(ixSnowSoilWat) = 1._dp

 ! compute terms in the Jacobian for vegetation (excluding fluxes)
 ! NOTE: this is computed outside the iteration loop because it does not depend on state variables
 ! NOTE: energy for vegetation is computed *within* the iteration loop as it includes phase change
 if(computeVegFlux)then
  dMat(ixCasNrg) = Cp_air*iden_air          ! volumetric heat capacity of air (J m-3 K-1)
  dMat(ixVegWat) = 1._dp                    ! nothing else on the left hand side
 endif

 ! compute terms in the Jacobian for the snow domain (excluding fluxes)
 ! NOTE: this is computed outside the iteration loop because it does not depend on state variables
 if(nSnow>0)&  ! (liquid water in snow only defined if snow layers exist)
 dMat(ixSnowOnlyWat) = 1._dp

 end associate  ! end association to variables in the data structures
 end subroutine popStateVec




 ! **********************************************************************************************************
 ! public subroutine varExtract: extract variables from the state vector and compute diagnostic variables
 ! **********************************************************************************************************
 subroutine varExtract(&
                       ! input
                       stateVec,                                  & ! intent(in):    model state vector (mixed units)
                       indx_data,                                 & ! intent(in):    indices defining model states and layers
                       snowfrz_scale,                             & ! intent(in):    scaling parameter for the snow freezing curve (K-1)
                       vGn_alpha,vGn_n,theta_sat,theta_res,vGn_m, & ! intent(in):    van Genutchen soil parameters
                       ! output: variables for the vegetation canopy
                       fracLiqVeg,                                & ! intent(out):   fraction of liquid water on the vegetation canopy (-)
                       scalarCanairTempTrial,                     & ! intent(out):   trial value of canopy air temperature (K)
                       scalarCanopyTempTrial,                     & ! intent(out):   trial value of canopy temperature (K)
                       scalarCanopyWatTrial,                      & ! intent(out):   trial value of canopy total water (kg m-2)
                       scalarCanopyLiqTrial,                      & ! intent(out):   trial value of canopy liquid water (kg m-2)
                       scalarCanopyIceTrial,                      & ! intent(out):   trial value of canopy ice content (kg m-2)
                       ! output: variables for the snow-soil domain
                       fracLiqSnow,                               & ! intent(out):   volumetric fraction of water in each snow layer (-)
                       mLayerTempTrial,                           & ! intent(out):   trial vector of layer temperature (K)
                       mLayerVolFracWatTrial,                     & ! intent(out):   trial vector of volumetric total water content (-) 
                       mLayerVolFracLiqTrial,                     & ! intent(out):   trial vector of volumetric liquid water content (-) 
                       mLayerVolFracIceTrial,                     & ! intent(out):   trial vector of volumetric ice water content (-) 
                       mLayerMatricHeadTrial,                     & ! intent(out):   trial vector of matric head (m)
                       ! output: error control 
                       err,message)                                 ! intent(out):   error control
 ! --------------------------------------------------------------------------------------------------------------------------------
 USE snow_utils_module,only:fracliquid                         ! compute the fraction of liquid water at a given temperature (snow)
 USE updatState_module,only:updateSnow                         ! update snow states
 USE updatState_module,only:updateSoil                         ! update soil states
 implicit none
 ! --------------------------------------------------------------------------------------------------------------------------------
 ! input
 real(dp),intent(in)             :: stateVec(:)                ! model state vector (mixed units)
 type(var_ilength),intent(in)    :: indx_data                  ! indices defining model states and layers                 
 real(dp),intent(in)             :: snowfrz_scale              ! scaling parameter for the snow freezing curve (K-1)
 real(dp),intent(in)             :: vGn_alpha,vGn_n,theta_sat,theta_res,vGn_m  ! van Genutchen soil parameters
 ! output: variables for the vegetation canopy
 real(dp),intent(out)            :: fracLiqVeg                 ! fraction of liquid water on the vegetation canopy (-)
 real(dp),intent(out)            :: scalarCanairTempTrial      ! trial value of canopy air temperature (K)
 real(dp),intent(out)            :: scalarCanopyTempTrial      ! trial value of canopy temperature (K)
 real(dp),intent(out)            :: scalarCanopyWatTrial       ! trial value of canopy total water (kg m-2)
 real(dp),intent(out)            :: scalarCanopyLiqTrial       ! trial value of canopy liquid water (kg m-2)
 real(dp),intent(out)            :: scalarCanopyIceTrial       ! trial value of canopy ice content (kg m-2)
 ! output: variables for the snow-soil domain
 real(dp),intent(out)            :: fracLiqSnow(:)             ! volumetric fraction of water in each snow layer (-)
 real(dp),intent(out)            :: mLayerTempTrial(:)         ! trial vector of layer temperature (K)
 real(dp),intent(out)            :: mLayerVolFracWatTrial(:)   ! trial vector of volumetric total water content (-)
 real(dp),intent(out)            :: mLayerVolFracLiqTrial(:)   ! trial vector of volumetric liquid water content (-)
 real(dp),intent(out)            :: mLayerVolFracIceTrial(:)   ! trial vector of volumetric ice water content (-)
 real(dp),intent(out)            :: mLayerMatricHeadTrial(:)   ! trial vector of matric head (m)
 ! output: error control 
 integer(i4b),intent(out)        :: err                        ! error code
 character(*),intent(out)        :: message                    ! error message
 ! --------------------------------------------------------------------------------------------------------------------------------
 ! local variables
 integer(i4b)                    :: iLayer                     ! index of layer in the snow and soil sub-domains
 character(len=256)              :: cMessage                   ! error message of downwind routine
 ! --------------------------------------------------------------------------------------------------------------------------------
 ! make association with variables in the data structures
 associate(&
 ! number of model layers, and layer type
 nSnow             => indx_data%var(iLookINDEX%nSnow)%dat(1)                 ,& ! intent(in): [i4b]    total number of snow layers
 nSoil             => indx_data%var(iLookINDEX%nSoil)%dat(1)                 ,& ! intent(in): [i4b]    total number of soil layers
 nLayers           => indx_data%var(iLookINDEX%nLayers)%dat(1)               ,& ! intent(in): [i4b]    total number of snow and soil layers
 nVegState         => indx_data%var(iLookINDEX%nVegState)%dat(1)             ,& ! intent(in): [i4b]    number of vegetation state variables
 layerType         => indx_data%var(iLookINDEX%layerType)%dat                ,& ! intent(in): [i4b(:)] index defining type of layer (soil or snow)
 ! indices defining model states and layers
 ixCasNrg          => indx_data%var(iLookINDEX%ixCasNrg)%dat(1)              ,& ! intent(in): [i4b] index of canopy air space energy state variable
 ixVegNrg          => indx_data%var(iLookINDEX%ixVegNrg)%dat(1)              ,& ! intent(in): [i4b] index of canopy energy state variable
 ixVegWat          => indx_data%var(iLookINDEX%ixVegWat)%dat(1)              ,& ! intent(in): [i4b] index of canopy hydrology state variable (mass)
 ixSnowSoilNrg     => indx_data%var(iLookINDEX%ixSnowSoilNrg)%dat            ,& ! intent(in): [i4b(:)] indices for energy states in the snow-soil subdomain
 ixSnowOnlyWat     => indx_data%var(iLookINDEX%ixSnowOnlyWat)%dat            ,& ! intent(in): [i4b(:)] indices for total water states in the snow subdomain
 ixSoilOnlyHyd     => indx_data%var(iLookINDEX%ixSoilOnlyHyd)%dat             & ! intent(in): [i4b(:)] indices for hydrology states in the soil subdomain
 ) ! association with variables in the data structures

 ! --------------------------------------------------------------------------------------------------------------------------------
 ! --------------------------------------------------------------------------------------------------------------------------------

 ! initialize canopy variables
 ! NOTE: Currently we expect three vegetation state variables (canopy air space temperature, canopy temperature, and intercepted water).
 !       The code here assumes that there are no other vegetation states, hence the check.
 if(nVegState==3)then
  ! (extract state variables)
  scalarCanairTempTrial = stateVec(ixCasNrg)
  scalarCanopyTempTrial = stateVec(ixVegNrg)
  scalarCanopyWatTrial  = stateVec(ixVegWat)                    ! total water storage on the canopy (kg m-2)
  ! (compute liquid water content)
  fracLiqVeg    = fracliquid(scalarCanopyTempTrial,snowfrz_scale)    ! fraction of liquid water (-)
  scalarCanopyLiqTrial = fracLiqVeg*scalarCanopyWatTrial             ! mass of liquid water on the canopy (kg m-2)
  scalarCanopyIceTrial = (1._dp - fracLiqVeg)*scalarCanopyWatTrial   ! mass of ice on the canopy (kg m-2)
 elseif(nVegState==0)then ! vegetation buried by snow
  ! (state variables)
  scalarCanairTempTrial = realMissing
  scalarCanopyTempTrial = realMissing
  scalarCanopyWatTrial  = realMissing
  ! (diagnostic variables)
  fracLiqVeg            = realMissing
  scalarCanopyLiqTrial  = realMissing
  scalarCanopyIceTrial  = realMissing
 else  ! unexpected
  message=trim(message)//'unexpected number of vegetation state variables'
  err=20; return
 endif

 ! extract state variables for layers in the snow-soil system
 mLayerTempTrial(1:nLayers)     = stateVec(ixSnowSoilNrg)
 mLayerVolFracWatTrial(1:nSnow) = stateVec(ixSnowOnlyWat)
 mLayerMatricHeadTrial(1:nSoil) = stateVec(ixSoilOnlyHyd)

 ! compute diagnostic variables in the snow and soil sub-domains
 do iLayer=1,nLayers
  select case(layerType(iLayer))

   !** snow
   case(ix_snow)
    call updateSnow(&
                    ! input
                    mLayerTempTrial(iLayer),                   & ! intent(in): layer temperature (K)
                    mLayerVolFracWatTrial(iLayer),             & ! intent(in): volumetric fraction of total water (-)
                    snowfrz_scale,                             & ! intent(in): scaling parameter for the snow freezing curve (K-1)
                    ! output
                    mLayerVolFracLiqTrial(iLayer),             & ! intent(out): volumetric fraction of liquid water (-)
                    mLayerVolFracIceTrial(iLayer),             & ! intent(out): volumetric fraction of ice (-)
                    fracLiqSnow(iLayer),                       & ! intent(out): fraction of liquid water in each snow layer (-)
                    err,cmessage)                                ! intent(out): error control
    if(err/=0)then; message=trim(message)//trim(cmessage); return; endif  ! (check for errors)

   !** soil
   case(ix_soil)
    call updateSoil(&
                    ! input
                    mLayerTempTrial(iLayer),                   & ! intent(in): layer temperature (K)
                    mLayerMatricHeadTrial(iLayer-nSnow),       & ! intent(in): matric head (m)
                    vGn_alpha,vGn_n,theta_sat,theta_res,vGn_m, & ! intent(in): van Genutchen soil parameters
                    ! output
                    mLayerVolFracWatTrial(iLayer),             & ! intent(out): volumetric fraction of total water (-)
                    mLayerVolFracLiqTrial(iLayer),             & ! intent(out): volumetric fraction of liquid water (-)
                    mLayerVolFracIceTrial(iLayer),             & ! intent(out): volumetric fraction of ice (-)
                    err,cmessage)                                ! intent(out): error control
    if(err/=0)then; message=trim(message)//trim(cmessage); return; endif  ! (check for errors)

   !** check errors
   case default; err=40; message=trim(message)//"cannot identify the layer as snow or soil"; return

  endselect  ! identifying type of layer

 end do   ! looping through layers

 ! end association to the variables in the data structures
 end associate

 end subroutine varExtract

end module getVectorz_module

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
USE multiconst,only:integerMissing  ! missing integer
USE multiconst,only:realMissing     ! missing real number

! domain types
USE globalData,only:iname_cas       ! named variables for canopy air space
USE globalData,only:iname_veg       ! named variables for vegetation canopy
USE globalData,only:iname_snow      ! named variables for snow
USE globalData,only:iname_soil      ! named variables for soil

! named variables to describe the state variable type
USE globalData,only:iname_nrgCanair ! named variable defining the energy of the canopy air space
USE globalData,only:iname_nrgCanopy ! named variable defining the energy of the vegetation canopy
USE globalData,only:iname_watCanopy ! named variable defining the mass of water on the vegetation canopy
USE globalData,only:iname_nrgLayer  ! named variable defining the energy state variable for snow+soil layers
USE globalData,only:iname_watLayer  ! named variable defining the total water state variable for snow+soil layers
USE globalData,only:iname_matLayer  ! named variable defining the matric head state variable for soil layers

! metadata for information in the data structures
USE globalData,only:indx_meta       ! metadata for the variables in the index structure

! constants
USE multiconst,only:&
                    Cp_air,       & ! specific heat of air                 (J kg-1 K-1)
                    LH_fus,       & ! latent heat of fusion                (J kg-1)
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
USE var_lookup,only:iLookDIAG       ! named variables for structure elements
USE var_lookup,only:iLookPROG       ! named variables for structure elements
USE var_lookup,only:iLookDERIV      ! named variables for structure elements
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
                        domainType,              & ! intent(in):    id of domain for desired model state variables
                        stateType,               & ! intent(in):    type of desired model state variables
                        stateIndx,               & ! intent(in):    indices for the desired state variables
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
 USE nr_utility_module,only:arth                   ! get a sequence of numbers arth(start, incr, count)
 ! --------------------------------------------------------------------------------------------------------------------------------
 ! input: model control
 logical(lgt),intent(in)         :: computeVegFlux         ! flag to indicate if we are computing fluxes over vegetation (.false. means veg is buried with snow)
 integer(i4b),intent(in)         :: domainType(:)          ! id of domain for desired model state variables
 integer(i4b),intent(in)         :: stateType(:)           ! type of desired model state variables
 integer(i4b),intent(in)         :: stateIndx(:)           ! indices for the desired model state variables
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
 ! number of state variables
 integer(i4b)                    :: nCasNrg                ! number of variables defining energy in the canopy air space (currently expect 1)
 integer(i4b)                    :: nVegNrg                ! number of variables defining energy in the vegetation canopy (currently expect 1)
 integer(i4b)                    :: nVegWat                ! number of variables defining water in the vegetation canopy (currently expect 1)
 integer(i4b)                    :: nLayNrg                ! number of variables defining energy in the snow+soil domain
 integer(i4b)                    :: nLayWat                ! number of variables defining water in the snow+soil domain
 integer(i4b)                    :: nLayMat                ! number of variables defining matric head in the snow+soil domain
 ! state subsets
 integer(i4b)                    :: iVar                   ! variable index 
 integer(i4b),dimension(nState)  :: ixSequence             ! sequential index in model state vector
 logical(lgt),dimension(nState)  :: stateMask              ! mask of state vector for specific state subsets
 character(len=256)              :: cmessage               ! error message of downwind routine
 ! --------------------------------------------------------------------------------------------------------------------------------
 ! --------------------------------------------------------------------------------------------------------------------------------
 ! make association with variables in the data structures
 fixedLength: associate(&
 ! model states for the vegetation canopy
 scalarCanairTemp  => prog_data%var(iLookPROG%scalarCanairTemp)%dat(1)       ,& ! intent(in): [dp] temperature of the canopy air space (K)
 scalarCanopyTemp  => prog_data%var(iLookPROG%scalarCanopyTemp)%dat(1)       ,& ! intent(in): [dp] temperature of the vegetation canopy (K)
 scalarCanopyWat   => prog_data%var(iLookPROG%scalarCanopyWat)%dat(1)        ,& ! intent(in): [dp] mass of total water on the vegetation canopy (kg m-2)
 ! model state variable vectors for the snow-soil layers
 mLayerTemp        => prog_data%var(iLookPROG%mLayerTemp)%dat                ,& ! intent(in): [dp(:)] temperature of each snow/soil layer (K)
 mLayerVolFracWat  => prog_data%var(iLookPROG%mLayerVolFracWat)%dat          ,& ! intent(in): [dp(:)] volumetric fraction of total water (-)
 mLayerMatricHead  => prog_data%var(iLookPROG%mLayerMatricHead)%dat          ,& ! intent(in): [dp(:)] matric head (m)
 ! model diagnostic variables
 canopyDepth       => diag_data%var(iLookDIAG%scalarCanopyDepth)%dat(1)      ,& ! intent(in): [dp   ] canopy depth (m)
 volHeatCapVeg     => diag_data%var(iLookDIAG%scalarBulkVolHeatCapVeg)%dat(1),& ! intent(in): [dp   ] bulk volumetric heat capacity of vegetation (J m-3 K-1)
 mLayerVolHeatCap  => diag_data%var(iLookDIAG%mLayerVolHtCapBulk)%dat        ,& ! intent(in): [dp(:)] bulk volumetric heat capacity in each snow and soil layer (J m-3 K-1)
 ! indices defining specific model states
 ixCasNrg          => indx_data%var(iLookINDEX%ixCasNrg)%dat                 ,& ! intent(in): [i4b] index of canopy air space energy state variable
 ixVegNrg          => indx_data%var(iLookINDEX%ixVegNrg)%dat                 ,& ! intent(in): [i4b] index of canopy energy state variable
 ixVegWat          => indx_data%var(iLookINDEX%ixVegWat)%dat                 ,& ! intent(in): [i4b] index of canopy hydrology state variable (mass)
 ixTopNrg          => indx_data%var(iLookINDEX%ixTopNrg)%dat                 ,& ! intent(in): [i4b] index of upper-most energy state in the snow-soil subdomain
 ixTopWat          => indx_data%var(iLookINDEX%ixTopWat)%dat                 ,& ! intent(in): [i4b] index of upper-most total water state in the snow-soil subdomain
 ixTopMat          => indx_data%var(iLookINDEX%ixTopMat)%dat                 ,& ! intent(in): [i4b] index of upper-most matric head state in the soil subdomain
 ! indices of the entire state vector, all model layers, and soil layers
 ixSoilState       => indx_data%var(iLookINDEX%ixSoilState)%dat              ,& ! intent(in): [i4b(:)] list of indices for all soil layers
 ixLayerState      => indx_data%var(iLookINDEX%ixLayerState)%dat             ,& ! intent(in): [i4b(:)] list of indices for all model layers
 ixMapSubset2Full  => indx_data%var(iLookINDEX%ixMapSubset2Full)%dat         ,& ! intent(in): [i4b(:)] list of indices in the full state vector that are in the state subset
 ! type of hydrology states in the snow+soil domain
 ixHydType         => indx_data%var(iLookINDEX%ixHydType)%dat                ,& ! intent(in): [i4b(:)] index of the type of hydrology states in snow+soil domain
 ! number of layers
 nSnow             => indx_data%var(iLookINDEX%nSnow)%dat(1)                 ,& ! intent(in): [i4b] number of snow layers
 nSoil             => indx_data%var(iLookINDEX%nSoil)%dat(1)                 ,& ! intent(in): [i4b] number of soil layers
 nLayers           => indx_data%var(iLookINDEX%nLayers)%dat(1)                & ! intent(in): [i4b] total number of layers
 )  ! end association with variables in the data structures
 ! --------------------------------------------------------------------------------------------------------------------------------
 ! initialize error control
 err=0; message='popStateVec/'

 ! -----
 ! * identify position in the state vector for specific variables...
 ! -----------------------------------------------------------------

 ! initialize indices
 ixCasNrg = integerMissing ! energy in the canopy air space
 ixVegNrg = integerMissing ! energy of the vegetation canopy
 ixVegWat = integerMissing ! mass of water in the vegetation canopy
 ixTopNrg = integerMissing ! upper-most energy state variable in the snow+soil domain
 ixTopWat = integerMissing ! upper-most hydrology state variable in the snow+soil domain
 ixTopMat = integerMissing ! upper-most matric head state variable in the snow+soil domain

 ! get the number of state variables in the vegetation canopy
 nCasNrg = count(stateType==iname_nrgCanair)
 nVegNrg = count(stateType==iname_nrgCanopy)
 nVegWat = count(stateType==iname_watCanopy)

 ! get the number of state variables in the snow+soil domains
 nLayNrg = count(stateType==iname_nrgLayer)
 nLayWat = count(stateType==iname_watLayer)
 nLayMat = count(stateType==iname_matLayer)

 ! define indices for vegetation variables
 ! NOTE: workaround for (not-yet-implemented) f2008 intrinsic findloc -- merge provides a vector with 1s where mask is true and 0s otherwise
 if(nCasNrg==1) ixCasNrg = maxloc( merge(1, 0, stateType==iname_nrgCanair) )
 if(nVegNrg==1) ixVegNrg = maxloc( merge(1, 0, stateType==iname_nrgCanopy) )
 if(nVegWat==1) ixVegWat = maxloc( merge(1, 0, stateType==iname_watCanopy) )

 ! define index for the upper-most energy state variable in the snow+soil domain
 ! NOTE: workaround for (not-yet-implemented) f2008 intrinsic findloc -- merge provides a vector with 1s where mask is true and 0s otherwise
 if(nLayNrg>0) ixTopNrg = maxloc( merge(1, 0, stateType==iname_nrgLayer) )

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

 ! -----
 ! * define vectors of specific variable types...
 ! ----------------------------------------------

 ! define index in full state vector
 ixSequence = arth(1,1,nState)

 ! loop through index variables
 do iVar=1,size(indx_data%var)

  ! define the mask
  select case(iVar)
   ! snow+soil domain
   case(iLookINDEX%ixSnowSoilNrg); stateMask = (stateType==iname_nrgLayer)                                 ! indices for energy states in the snow+soil subdomain
   case(iLookINDEX%ixSnowSoilWat); stateMask = (stateType==iname_watLayer)                                 ! indices for total water states in the snow+soil subdomain
   ! snow only domain
   case(iLookINDEX%ixSnowOnlyNrg); stateMask = (stateType==iname_nrgLayer .and. domainType==iname_snow)    ! indices for energy states in the snow subdomain
   case(iLookINDEX%ixSnowOnlyWat); stateMask = (stateType==iname_watLayer .and. domainType==iname_snow)    ! indices for total water states in the snow subdomain
   ! soil only domain
   case(iLookINDEX%ixSoilOnlyNrg); stateMask = (stateType==iname_nrgLayer .and. domainType==iname_soil)    ! indices for energy states in the soil subdomain
   case(iLookINDEX%ixSoilOnlyHyd); stateMask = (stateType==iname_matLayer .and. domainType==iname_soil)    ! indices for hydrology states in the soil subdomain
   ! vectors defining specific model states
   case(iLookINDEX%ixNrgOnly);     stateMask = (stateType==iname_nrgCanair .or. stateType==iname_nrgCanopy .or. stateType==iname_nrgLayer)  ! list of indices for all energy states
   case(iLookINDEX%ixWatOnly);     stateMask = (stateType==iname_watLayer )   ! list of indices for all "total water" states
   case(iLookINDEX%ixMatOnly);     stateMask = (stateType==iname_matLayer )   ! list of indices for matric head state variables
   case(iLookINDEX%ixMassOnly);    stateMask = (stateType==iname_watCanopy)   ! list of indices for hydrology states (mass of water)
   ! ** ignore all other variables
   case default; cycle ! only need to process the above variables
  end select  ! iVar

  ! get the subset of indices
  ! NOTE: indxSubset(subset, fullVector, mask), provides subset of fullVector where mask==.true.
  call indxSubset(indx_data%var(iVar)%dat,ixSequence,stateMask,err,cmessage)
  if(err/=0)then; message=trim(message)//trim(cmessage)//'[varname='//trim(indx_meta(ivar)%varname)//']'; return; endif

 end do  ! looping through variables in the data structure

 ! get the indices of the state subset in the full state vector
 ! NOTE: indxSubset(subset, fullVector, mask), provides subset of fullVector where mask==.true.
 call indxSubset(indx_data%var(iLookINDEX%ixMapSubset2Full)%dat, stateIndx, stateIndx/=integerMissing, err, cmessage)
 if(err/=0)then; message=trim(message)//trim(cmessage)//'[varname='//trim(indx_meta(iLookINDEX%ixVolFracWat)%varname)//']'; return; endif

 ! get the subset of indices in the snow+soil domain layers where the state variable is volumetric total water
 ! NOTE: indxSubset(subset, fullVector, mask), provides subset of fullVector where mask==.true.
 call indxSubset(indx_data%var(iLookINDEX%ixVolFracWat)%dat, ixLayerState, ixHydType==iname_watLayer, err, cmessage)
 if(err/=0)then; message=trim(message)//trim(cmessage)//'[varname='//trim(indx_meta(iLookINDEX%ixVolFracWat)%varname)//']'; return; endif

 ! get the subset of indices in the soil domain layers where the state variable is matric head
 ! NOTE: indxSubset(subset, fullVector, mask), provides subset of fullVector where mask==.true.
 call indxSubset(indx_data%var(iLookINDEX%ixMatricHead)%dat, ixSoilState, ixHydType(nSnow+1:nLayers)==iname_matLayer, err, cmessage)
 if(err/=0)then; message=trim(message)//trim(cmessage)//'[varname='//trim(indx_meta(iLookINDEX%ixVolFracWat)%varname)//']'; return; endif

 !  make association with indices defining variable types for specific sub-domains
 variableLength: associate(&
 ixSnowSoilNrg     => indx_data%var(iLookINDEX%ixSnowSoilNrg)%dat            ,& ! intent(in): [i4b] indices IN THE FULL VECTOR for energy states in the snow-soil subdomain
 ixSnowSoilWat     => indx_data%var(iLookINDEX%ixSnowSoilWat)%dat            ,& ! intent(in): [i4b] indices IN THE FULL VECTOR for total water states in the snow-soil subdomain
 ixSnowOnlyNrg     => indx_data%var(iLookINDEX%ixSnowOnlyNrg)%dat            ,& ! intent(in): [i4b] indices IN THE FULL VECTOR for energy states in the snow subdomain
 ixSnowOnlyWat     => indx_data%var(iLookINDEX%ixSnowOnlyWat)%dat            ,& ! intent(in): [i4b] indices IN THE FULL VECTOR for total water states in the snow subdomain
 ixSoilOnlyNrg     => indx_data%var(iLookINDEX%ixSoilOnlyNrg)%dat            ,& ! intent(in): [i4b] indices IN THE FULL VECTOR for energy states in the soil subdomain
 ixSoilOnlyHyd     => indx_data%var(iLookINDEX%ixSoilOnlyHyd)%dat            ,& ! intent(in): [i4b] indices IN THE FULL VECTOR for hydrology states in the soil subdomain
 ixVolFracWat      => indx_data%var(iLookINDEX%ixVolFracWat)%dat             ,& ! intent(in): [i4b] indices IN THE SNOW+SOIL VECTOR for hydrology states
 ixMatricHead      => indx_data%var(iLookINDEX%ixMatricHead)%dat              & ! intent(in): [i4b] indices IN THE SOIL VECTOR for hydrology states
 )

 ! -----
 ! * initialize state vectors...
 ! -----------------------------

 ! build elements of the state vector for the vegetation canopy
 if(ixCasNrg(1)/=integerMissing) stateVec(ixCasNrg) = scalarCanairTemp
 if(ixVegNrg(1)/=integerMissing) stateVec(ixVegNrg) = scalarCanopyTemp
 if(ixVegWat(1)/=integerMissing) stateVec(ixVegWat) = scalarCanopyWat  ! kg m-2

 ! build the energy state vector for the snow and soil domain
 if(nLayNrg>0) stateVec(ixSnowSoilNrg) = mLayerTemp(1:nLayers)

 ! build the hydrology state vector for the snow+soil domains
 ! NOTE: This enables primary variable switching and can use volFracWat as the state variable for soil layers
 if(nLayWat>0) stateVec(ixSnowSoilWat) = mLayerVolFracWat(ixVolFracWat)
 if(nLayMat>0) stateVec(ixSoilOnlyHyd) = mLayerMatricHead(ixMatricHead)

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
  if(ixCasNrg(1)/=integerMissing) sMul(ixCasNrg) = Cp_air*iden_air   ! volumetric heat capacity of air (J m-3 K-1)
  if(ixVegNrg(1)/=integerMissing) sMul(ixVegNrg) = volHeatCapVeg     ! volumetric heat capacity of the vegetation (J m-3 K-1)
  if(ixVegWat(1)/=integerMissing) sMul(ixVegWat) = 1._dp             ! nothing else on the left hand side
 endif

 ! define the multiplier for the state vector for residual calculations (snow-soil domain)
 if(nLayNrg>0) sMul(ixSnowSoilNrg) = mLayerVolHeatCap(1:nLayers)
 if(nLayWat>0) sMul(ixSnowSoilWat) = 1._dp 
 if(nLayMat>0) sMul(ixSoilOnlyHyd) = 1._dp 

 ! compute terms in the Jacobian for vegetation (excluding fluxes)
 ! NOTE: this is computed outside the iteration loop because it does not depend on state variables
 ! NOTE: energy for vegetation is computed *within* the iteration loop as it includes phase change
 if(computeVegFlux)then
  if(ixCasNrg(1)/=integerMissing) dMat(ixCasNrg) = Cp_air*iden_air          ! volumetric heat capacity of air (J m-3 K-1)
  if(ixVegWat(1)/=integerMissing) dMat(ixVegWat) = 1._dp                    ! nothing else on the left hand side
 endif

 ! compute terms in the Jacobian for the snow domain (excluding fluxes)
 ! NOTE: this is computed outside the iteration loop because it does not depend on state variables
 if(size(ixSnowOnlyWat)>0) dMat(ixSnowOnlyWat) = 1._dp

 ! ------------------------------------------------------------------------------------------
 ! ------------------------------------------------------------------------------------------
 
 end associate variableLength   ! end association to variables in the data structure where vector length may change
 end associate fixedLength      ! end association to variables in the data structure where vector length does not change
 end subroutine popStateVec




 ! **********************************************************************************************************
 ! public subroutine varExtract: extract variables from the state vector and compute diagnostic variables
 ! **********************************************************************************************************
 subroutine varExtract(&
                       ! input
                       do_adjustTemp,                             & ! intent(in):    logical flag to adjust temperature to account for the energy used in melt+freeze
                       fluxVec,         & ! OPTIONAL              & ! intent(in):    model flux vector (mixed units)
                       stateVec,                                  & ! intent(in):    model state vector (mixed units)
                       mpar_data,                                 & ! intent(in):    model parameters for a local HRU
                       diag_data,                                 & ! intent(in):    model diagnostic variables for a local HRU
                       prog_data,                                 & ! intent(in):    model prognostic variables for a local HRU
                       deriv_data,                                & ! intent(in):    derivatives in model fluxes w.r.t. relevant state variables
                       indx_data,                                 & ! intent(in):    indices defining model states and layers
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
 implicit none 
 ! input
 logical(lgt),intent(in)         :: do_adjustTemp              ! flag to adjust temperature to account for the energy used in melt+freeze
 real(dp),intent(in),optional    :: fluxVec(:)                 ! model flux vector (mixed units)
 real(dp),intent(in)             :: stateVec(:)                ! model state vector (mixed units)
 type(var_d),      intent(in)    :: mpar_data                  ! model parameters for a local HRU
 type(var_dlength),intent(in)    :: diag_data                  ! diagnostic variables for a local HRU
 type(var_dlength),intent(in)    :: prog_data                  ! prognostic variables for a local HRU
 type(var_dlength),intent(in)    :: deriv_data                 ! derivatives in model fluxes w.r.t. relevant state variables 
 type(var_ilength),intent(in)    :: indx_data                  ! indices defining model states and layers                 
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
 integer(i4b)                    :: iState                     ! index of model state variable
 integer(i4b)                    :: ixFullVector               ! index within full state vector
 integer(i4b)                    :: ixDomainType               ! name of a given model domain
 integer(i4b)                    :: ixControlIndex             ! index within a given model domain
 integer(i4b)                    :: ixOther,ixOtherLocal       ! index of the coupled state variable within the (full, local) vector
 logical(lgt)                    :: isCoupled                  ! .true. if a given variable shared another state variable in the same control volume
 logical(lgt)                    :: isNrgState                 ! .true. if a given variable is an energy state 
 logical(lgt),allocatable        :: computedCoupling(:)        ! .true. if computed the coupling for a given state variable
 real(dp)                        :: scalarVolFracLiq           ! volumetric fraction of liquid water (-)
 real(dp)                        :: scalarVolFracIce           ! volumetric fraction of ice (-)
 real(dp)                        :: rhsNrgFlux                 ! the right-hand-side energy flux (J m-3), is zero for uncoupled and non-energy states
 character(len=256)              :: cMessage                   ! error message of downwind routine
 ! --------------------------------------------------------------------------------------------------------------------------------
 ! make association with variables in the data structures
 associate(&
 ! number of model layers, and layer type
 nSnow                   => indx_data%var(iLookINDEX%nSnow)%dat(1)                 ,& ! intent(in): [i4b]    total number of snow layers
 nSoil                   => indx_data%var(iLookINDEX%nSoil)%dat(1)                 ,& ! intent(in): [i4b]    total number of soil layers
 nLayers                 => indx_data%var(iLookINDEX%nLayers)%dat(1)               ,& ! intent(in): [i4b]    total number of snow and soil layers
 nVegState               => indx_data%var(iLookINDEX%nVegState)%dat(1)             ,& ! intent(in): [i4b]    number of vegetation state variables
 layerType               => indx_data%var(iLookINDEX%layerType)%dat                ,& ! intent(in): [i4b(:)] index defining type of layer (soil or snow)
 ! indices defining model states and layers
 ixCasNrg                => indx_data%var(iLookINDEX%ixCasNrg)%dat(1)              ,& ! intent(in): [i4b] index of canopy air space energy state variable
 ixVegNrg                => indx_data%var(iLookINDEX%ixVegNrg)%dat(1)              ,& ! intent(in): [i4b] index of canopy energy state variable
 ixVegWat                => indx_data%var(iLookINDEX%ixVegWat)%dat(1)              ,& ! intent(in): [i4b] index of canopy hydrology state variable (mass)
 ixSnowSoilNrg           => indx_data%var(iLookINDEX%ixSnowSoilNrg)%dat            ,& ! intent(in): [i4b(:)] indices IN THE FULL VECTOR for energy states in the snow-soil subdomain
 ixSnowSoilWat           => indx_data%var(iLookINDEX%ixSnowSoilWat)%dat            ,& ! intent(in): [i4b(:)] indices IN THE FULL VECTOR for total water states in the snow-soil subdomain
 ixSnowOnlyNrg           => indx_data%var(iLookINDEX%ixSnowOnlyNrg)%dat            ,& ! intent(in): [i4b(:)] indices IN THE FULL VECTOR for energy states in the snow subdomain
 ixSnowOnlyWat           => indx_data%var(iLookINDEX%ixSnowOnlyWat)%dat            ,& ! intent(in): [i4b(:)] indices IN THE FULL VECTOR for total water states in the snow subdomain
 ixSoilOnlyHyd           => indx_data%var(iLookINDEX%ixSoilOnlyHyd)%dat            ,& ! intent(in): [i4b(:)] indices IN THE FULL VECTOR for hydrology states in the soil subdomain
 ixVolFracWat            => indx_data%var(iLookINDEX%ixVolFracWat)%dat             ,& ! intent(in): [i4b(:)] indices IN THE SNOW+SOIL VECTOR for hydrology states in the soil subdomain
 ixMatricHead            => indx_data%var(iLookINDEX%ixMatricHead)%dat             ,& ! intent(in): [i4b(:)] indices IN THE SOIL VECTOR for hydrology states in the soil subdomain
 ! model diagnostic variables
 canopyDepth             => diag_data%var(iLookDIAG%scalarCanopyDepth)%dat(1)      ,& ! intent(in): [dp   ] canopy depth (m)
 scalarBulkVolHeatCapVeg => diag_data%var(iLookDIAG%scalarBulkVolHeatCapVeg)%dat(1),& ! intent(in): [dp   ] volumetric heat capacity of the vegetation (J m-3 K-1)
 mLayerVolHtCapBulk      => diag_data%var(iLookDIAG%mLayerVolHtCapBulk)%dat        ,& ! intent(in): [dp(:)] volumetric heat capacity in each layer (J m-3 K-1)
 ! model derivatives
 dTheta_dTkCanopy        => deriv_data%var(iLookDERIV%dTheta_dTkCanopy)%dat(1)     ,& ! intent(in): [dp   ] derivative of volumetric liquid water content w.r.t. temperature (K-1)
 mLayerdTheta_dTk        => deriv_data%var(iLookDERIV%mLayerdTheta_dTk)%dat        ,& ! intent(in): [dp(:)] derivative of volumetric liquid water content w.r.t. temperature (K-1)
 ! model states for the vegetation canopy
 scalarCanairTemp        => prog_data%var(iLookPROG%scalarCanairTemp)%dat(1)       ,& ! intent(in): [dp] temperature of the canopy air space (K)
 scalarCanopyTemp        => prog_data%var(iLookPROG%scalarCanopyTemp)%dat(1)       ,& ! intent(in): [dp] temperature of the vegetation canopy (K)
 scalarCanopyWat         => prog_data%var(iLookPROG%scalarCanopyWat)%dat(1)        ,& ! intent(in): [dp] mass of total water on the vegetation canopy (kg m-2)
 ! model state variable vectors for the snow-soil layers
 mLayerTemp              => prog_data%var(iLookPROG%mLayerTemp)%dat                ,& ! intent(in): [dp(:)] temperature of each snow/soil layer (K)
 mLayerVolFracWat        => prog_data%var(iLookPROG%mLayerVolFracWat)%dat          ,& ! intent(in): [dp(:)] volumetric fraction of total water (-)
 mLayerMatricHead        => prog_data%var(iLookPROG%mLayerMatricHead)%dat          ,& ! intent(in): [dp(:)] matric head (m)
 ! model diagnostic variables from a previous solution
 scalarCanopyIce         => prog_data%var(iLookPROG%scalarCanopyIce)%dat(1)        ,& ! intent(in): [dp(:)] mass of ice on the vegetation canopy (kg m-2)
 mLayerVolFracIce        => prog_data%var(iLookPROG%mLayerVolFracIce)%dat           & ! intent(in): [dp(:)] volumetric fraction of ice (-)
 ) ! association with variables in the data structures

 ! --------------------------------------------------------------------------------------------------------------------------------
 ! --------------------------------------------------------------------------------------------------------------------------------

 ! initialize error control
 err=0; message='varExtract/'

 ! check that the flux vector is present when the do_adjustTemp flag is true
 if(.not.present(fluxVec) .and. do_adjustTemp)then
  message=trim(message)//'flux vector must be present when do_adjustTemp=.true.'
  err=20; return
 endif

 if(do_adjustTemp) print*, 'fluxVec = ', fluxVec 

 ! *****
 ! * part 1: extract state variables...
 ! ************************************

 ! *** extract state variables for the vegetation canopy

 ! check if computing the vegetation flux
 if(ixCasNrg/=integerMissing .or. ixVegNrg/=integerMissing .or. ixVegWat/=integerMissing)then
 
  ! extract temperature of the canopy air space
  if(ixCasNrg/=integerMissing)then
   scalarCanairTempTrial = stateVec(ixCasNrg)
  else
   scalarCanairTempTrial = scalarCanairTemp     ! state variable from the last update
  endif
  
  ! extract canopy temperature
  if(ixVegNrg/=integerMissing)then
   scalarCanopyTempTrial = stateVec(ixVegNrg) 
  else
   scalarCanopyTempTrial = scalarCanopyTemp     ! state variable from the last update
  endif
  
  ! extract intercepted water
  if(ixVegWat/=integerMissing)then
   scalarCanopyWatTrial  = stateVec(ixVegWat)
  else
   scalarCanopyWatTrial  = scalarCanopyWat      ! state variable from the last update
  endif
  
 ! not computing the vegetation flux (veg buried with snow, or bare ground)
 else
  scalarCanairTempTrial = realMissing
  scalarCanopyTempTrial = realMissing
  scalarCanopyWatTrial  = realMissing
 endif  ! not computing the vegetation flux

 ! *** extract state variables from the snow+soil sub-domain

 ! initialize to the state variable from the last update
 mLayerTempTrial       = mLayerTemp
 mLayerVolFracWatTrial = mLayerVolFracWat
 mLayerMatricHeadTrial = mLayerMatricHead

 ! overwrite with the energy values from the state vector
 if(size(ixSnowSoilNrg)>0) mLayerTempTrial(1:nLayers)     = stateVec(ixSnowSoilNrg)

 ! overwrite the volumetric liquid water content values from the state vector
 ! NOTE: these could be in the snow and soil domains, implementing primary variable switching for Richards' equation
 if(size(ixSnowSoilWat)>0) mLayerVolFracWatTrial(ixVolFracWat) = stateVec(ixSnowSoilWat)

 ! overwrite the matric head values from the state vector
 if(size(ixSoilOnlyHyd)>0) mLayerMatricHeadTrial(ixMatricHead) = stateVec(ixSoilOnlyHyd) 

 ! --------------------------------------------------------------------------------------------------------------------------------
 ! --------------------------------------------------------------------------------------------------------------------------------

 ! *****
 ! * part 2: update volumetric fraction of liquid water and ice (and also update temperature if desired)... 
 ! ********************************************************************************************************

 ! make association with indices
 indices: associate(&
  ixNrgLayer          => indx_data%var(iLookINDEX%ixNrgLayer)%dat             ,& ! intent(in): [i4b(:)] indices IN THE FULL VECTOR for energy states in the snow+soil domain
  ixHydLayer          => indx_data%var(iLookINDEX%ixHydLayer)%dat             ,& ! intent(in): [i4b(:)] indices IN THE FULL VECTOR for hydrology states in the snow+soil domain
  ixMapFull2Subset    => indx_data%var(iLookINDEX%ixMapFull2Subset)%dat       ,& ! intent(in): [i4b(:)] list of indices in the state subset for each state in the full state vector
  ixDomainType_subset => indx_data%var(iLookINDEX%ixDomainType_subset)%dat    ,& ! intent(in): [i4b(:)] [state subset] id of domain for desired model state variables
  ixAllState_subset   => indx_data%var(iLookINDEX%ixAllState_subset)%dat      ,& ! intent(in): [i4b(:)] [state subset] list of indices in the full state vector
  ixControlVolume     => indx_data%var(iLookINDEX%ixControlVolume)%dat        ,& ! intent(in): [i4b(:)] index of the control volume for different domains (veg, snow, soil)
  ixStateType         => indx_data%var(iLookINDEX%ixStateType)%dat             & ! intent(in): [i4b(:)] indices defining the type of the state (iname_nrgLayer...)
 ) ! association with model indices

 ! make association with parameters in the data structures
 parameters: associate(&
  vGn_m               => diag_data%var(iLookDIAG%scalarVGn_m)%dat(1)          ,& ! intent(in): [dp] van Genutchen "m" parameter (-)
  vGn_n               => mpar_data%var(iLookPARAM%vGn_n)                      ,& ! intent(in): [dp] van Genutchen "n" parameter (-)
  vGn_alpha           => mpar_data%var(iLookPARAM%vGn_alpha)                  ,& ! intent(in): [dp] van Genutchen "alpha" parameter (m-1)
  theta_sat           => mpar_data%var(iLookPARAM%theta_sat)                  ,& ! intent(in): [dp] soil porosity (-)
  theta_res           => mpar_data%var(iLookPARAM%theta_res)                  ,& ! intent(in): [dp] soil residual volumetric water content (-)
  snowfrz_scale       => mpar_data%var(iLookPARAM%snowfrz_scale)               & ! intent(in): [dp] scaling parameter for the snow freezing curve (K-1)
 ) ! association with parameters in the data structures

 ! check that the dimensions are as expected
 if(size(ixDomainType_subset)/=size(stateVec) .or. size(ixAllState_subset)/=size(stateVec) )then
  message=trim(message)//'size mismatch for state subset'
  err=20; return
 endif

 ! allocate space and assign values to the flag vector
 allocate(computedCoupling(size(stateVec)),stat=err)        ! .true. if computed the coupling for a given state variable
 if(err/=0)then; message=trim(message)//'problem allocating computedCoupling'; return; endif
 computedCoupling(:)=.false.

 ! loop through model state variables
 do iState=1,size(stateVec)

  ! get domain type, and index of the control volume within the domain
  ixFullVector   = ixAllState_subset(iState)      ! index within full state vector
  ixDomainType   = ixDomainType_subset(iState)    ! named variables defining the domain (iname_cas, iname_veg, etc.)
  ixControlIndex = ixControlVolume(ixFullVector)  ! index within a given domain

  print*, 'ixFullVector   = ', ixFullVector
  print*, 'ixDomainType   = ', ixDomainType
  print*, 'ixControlIndex = ', ixControlIndex

  ! get the index of the other (energy or mass) state variable within the full state vector
  select case(ixDomainType)
   case(iname_cas);  cycle ! canopy air space: do nothing
   case(iname_veg);  ixOther = merge(ixVegWat,ixVegNrg,ixStateType(ixFullVector)==iname_nrgCanopy)  ! ixVegWat if stateType=iname_nrgCanopy; ixVegNrg otherwise
   case(iname_snow); ixOther = merge(ixHydLayer(ixControlIndex),ixNrgLayer(ixControlIndex),ixStateType(ixFullVector)==iname_nrgLayer)
   case(iname_soil); ixOther = merge(ixHydLayer(nSnow+ixControlIndex),ixNrgLayer(nSnow+ixControlIndex),ixStateType(ixFullVector)==iname_nrgLayer)
   case default; err=20; message=trim(message)//'expect case to be iname_cas, iname_veg, iname_snow, iname_soil'; return
  end select

  ! get the index in the local state vector
  ixOtherLocal = ixMapFull2Subset(ixOther)  ! ixOtherLocal could equal integerMissing

  print*, 'ixOther      = ', ixOther
  print*, 'ixOtherLocal = ', ixOtherLocal

  ! check if we have a coupled solution
  isCoupled    = (ixOtherLocal/=integerMissing)

  ! check if we are an energy state
  isNrgState   = (ixStateType(ixFullVector)==iname_nrgCanopy .or. ixStateType(ixFullVector)==iname_nrgLayer)

  print*, 'isCoupled  = ', isCoupled
  print*, 'isNrgState = ', isNrgState

  ! get the right-hand-side energy flux
  if(isCoupled .or. isNrgState)then
   rhsNrgFlux = merge(fluxVec(iState),fluxVec(ixOtherLocal),isNrgState)
  else
   rhsNrgFlux = 0._dp
  endif

  ! get desired variables for each control volume
  select case(ixDomainType)
  
   ! vegetation
   case(iname_veg)
    call updateVars(&
                    ! model control
                    ixDomainType,                                               & ! intent(in)   : named variable defining the domain type (iname_veg, iname_snow, iname_soil)
                    do_adjustTemp,                                              & ! intent(in)    : logical flag to adjust temperature to account for energy associated with melt+freeze
                    snowfrz_scale,                                              & ! intent(in)    : scaling parameter for the snow freezing curve (K-1)
                    vGn_alpha,vGn_n,theta_sat,theta_res,vGn_m,                  & ! intent(in)    : soil parameters
                    ! input
                    LH_fus*iden_water                                          ,& ! intent(in)    : energy for melt+freeze (J m-3)
                    dTheta_dTkCanopy                                           ,& ! intent(in)    : derivative of volumetric liquid water content w.r.t. temperature (K-1)
                    scalarBulkVolHeatCapVeg                                    ,& ! intent(in)    : volumetric heat capacity (J m-3 K-1)
                    scalarCanopyTemp                                           ,& ! intent(in)    : initial temperature (K) 
                    scalarCanopyIce/(iden_water*canopyDepth)                   ,& ! intent(in)    : initial volumetric fraction of ice (-)
                    scalarCanopyWatTrial/(iden_water*canopyDepth)              ,& ! intent(in)    : mass state variable = trial volumetric fraction of total water (-)
                    ! output 
                    scalarCanopyTempTrial                                      ,& ! intent(inout) : trial temperature (K)
                    scalarVolFracLiq                                           ,& ! intent(out)   : trial volumetric fraction of liquid water (-)
                    scalarVolFracIce                                           ,& ! intent(out)   : trial volumetric fraction if ice (-)
                    fracLiqTrial=fracLiqVeg                                    ,& ! intent(out)   : fraction of liquid water (-)
                    err=err,message=cmessage)                                     ! intent(out)   : error control
    if(err/=0)then; message=trim(message)//trim(cmessage); return; endif

    ! compute mass of water on the canopy
    scalarCanopyLiqTrial = scalarVolFracLiq*(iden_water*canopyDepth)
    scalarCanopyIceTrial = scalarVolFracIce*(iden_water*canopyDepth)

   ! snow layers
   case(iname_snow)
    call updateVars(&
                    ! model control
                    ixDomainType,                                               & ! intent(in)   : named variable defining the domain type (iname_veg, iname_snow, iname_soil)
                    do_adjustTemp,                                              & ! intent(in)    : logical flag to adjust temperature to account for energy associated with melt+freeze
                    snowfrz_scale,                                              & ! intent(in)    : scaling parameter for the snow freezing curve (K-1)
                    vGn_alpha,vGn_n,theta_sat,theta_res,vGn_m,                  & ! intent(in)    : soil parameters
                    ! input
                    LH_fus*iden_ice                                            ,& ! intent(in)    : energy for melt+freeze (J m-3)
                    mLayerdTheta_dTk(ixControlIndex)                           ,& ! intent(in)    : derivative of volumetric liquid water content w.r.t. temperature (K-1)
                    mLayerVolHtCapBulk(ixControlIndex)                         ,& ! intent(in)    : volumetric heat capacity (J m-3 K-1)
                    mLayerTemp(ixControlIndex)                                 ,& ! intent(in)    : initial temperature (K)
                    mLayerVolFracIce(ixControlIndex)                           ,& ! intent(in)    : initial volumetric fraction of ice (-)
                    mLayerVolFracWatTrial(ixControlIndex)                      ,& ! intent(in)    : mass state variable = trial volumetric fraction of water (-)
                    ! output 
                    mLayerTempTrial(ixControlIndex)                            ,& ! intent(inout) : trial temperature (K)
                    mLayerVolFracLiqTrial(ixControlIndex)                      ,& ! intent(out)   : trial volumetric fraction of liquid water (-)
                    mLayerVolFracIceTrial(ixControlIndex)                      ,& ! intent(out)   : trial volumetric fraction if ice (-)
                    fracLiqTrial=fracLiqSnow(ixControlIndex)                   ,& ! intent(out)   : fraction of liquid water (-)
                    err=err,message=cmessage)                                     ! intent(out)   : error control
    if(err/=0)then; message=trim(message)//trim(cmessage); return; endif

   ! soil layers
   case(iname_soil)
    call updateVars(&
                    ! model control
                    ixDomainType,                                               & ! intent(in)   : named variable defining the domain type (iname_veg, iname_snow, iname_soil)
                    do_adjustTemp,                                              & ! intent(in)    : logical flag to adjust temperature to account for energy associated with melt+freeze
                    snowfrz_scale,                                              & ! intent(in)    : scaling parameter for the snow freezing curve (K-1)
                    vGn_alpha,vGn_n,theta_sat,theta_res,vGn_m,                  & ! intent(in)    : soil parameters
                    ! input
                    LH_fus*iden_water                                          ,& ! intent(in)    : energy for melt+freeze (J m-3)
                    mLayerdTheta_dTk(nSnow+ixControlIndex)                     ,& ! intent(in)    : derivative of volumetric liquid water content w.r.t. temperature (K-1)
                    mLayerVolHtCapBulk(nSnow+ixControlIndex)                   ,& ! intent(in)    : volumetric heat capacity (J m-3 K-1)
                    mLayerTemp(nSnow+ixControlIndex)                           ,& ! intent(in)    : initial temperature (K)
                    mLayerVolFracIce(nSnow+ixControlIndex)                     ,& ! intent(in)    : initial volumetric fraction of ice (-)
                    mLayerMatricHeadTrial(ixControlIndex)                      ,& ! intent(in)    : mass state variable = trial matric head (m)
                    ! output 
                    mLayerTempTrial(nSnow+ixControlIndex)                      ,& ! intent(inout) : trial temperature (K)
                    mLayerVolFracLiqTrial(nSnow+ixControlIndex)                ,& ! intent(out)   : trial volumetric fraction of liquid water (-)
                    mLayerVolFracIceTrial(nSnow+ixControlIndex)                ,& ! intent(out)   : trial volumetric fraction if ice (-)
                    volFracWatTrial=mLayerVolFracWatTrial(nSnow+ixControlIndex),& ! intent(out)   : volumetric fraction of total water (-)
                    err=err,message=cmessage)                                     ! intent(out)   : error control
    if(err/=0)then; message=trim(message)//trim(cmessage); return; endif

   case default; err=20; message=trim(message)//'expect case to be iname_veg, iname_snow, iname_soil'; return

  endselect  ! domain type

 end do  ! looping through state variables

 ! deallocate space
 deallocate(computedCoupling,stat=err)        ! .true. if computed the coupling for a given state variable
 if(err/=0)then; message=trim(message)//'problem deallocating computedCoupling'; return; endif

 ! end association to the variables in the data structures
 end associate parameters
 end associate indices
 end associate

 end subroutine varExtract

 
 ! **********************************************************************************************************
 ! private subroutine updateProg: update prognostoc variables
 !  * compute volumetric fraction of liquid water and ice and update temperature
 ! **********************************************************************************************************
 subroutine updateVars(&
                       ! model control
                       ixDomainType,                                  & ! intent(in)    : named variable defining the domain type (iname_veg, iname_snow, iname_soil)
                       do_adjustTemp,                                 & ! intent(in)    : logical flag to adjust temperature to account for energy associated with melt+freeze
                       snowfrz_scale,                                 & ! intent(in)    : scaling parameter for the snow freezing curve (K-1)
                       vGn_alpha,vGn_n,theta_sat,theta_res,vGn_m,     & ! intent(in)    : soil parameters
                       ! input
                       meltNrg                                       ,& ! intent(in)    : energy for melt+freeze (J m-3)
                       dLiq_dT                                       ,& ! intent(in)    : derivative of volumetric liquid water content w.r.t. temperature (K-1)
                       heatCap                                       ,& ! intent(in)    : volumetric heat capacity (J m-3 K-1)
                       tempInit                                      ,& ! intent(in)    : initial temperature (K)
                       volFracIceInit                                ,& ! intent(in)    : initial volumetric fraction of ice (-)
                       massStateTrial                                ,& ! intent(in)    : trial mass state variable
                       ! output
                       tempTrial                                     ,& ! intent(inout) : trial temperature (K)
                       volFracLiqTrial                               ,& ! intent(out)  : volumetric fraction of liquid water (-)
                       volFracIceTrial                               ,& ! intent(out)  : volumetric fraction of ice (-)
                       volFracWatTrial ,&  ! OPTIONAL                   ! intent(out)  : volumetric fraction of total water (-)
                       fracLiqTrial    ,&  ! OPTIONAL                   ! intent(out)  : fraction liquid water (-)
                       ! error handling
                       err,message)                                     ! intent(out)  : error control
 ! --------------------------------------------------------------------------------------------------------------------------------
 USE updatState_module,only:updateSnow  ! update snow states
 USE updatState_module,only:updateSoil  ! update soil states
 ! --------------------------------------------------------------------------------------------------------------------------------
 implicit none
 ! model control
 integer(i4b),intent(in)        :: ixDomainType     ! named variable defining the domain type (iname_veg, iname_snow, iname_soil)
 logical(lgt),intent(in)        :: do_adjustTemp    ! logical flag to adjust temperature to account for energy associated with melt+freeze
 ! parameters
 real(dp),intent(in)            :: vGn_m            ! van Genutchen "m" parameter (-)
 real(dp),intent(in)            :: vGn_n            ! van Genutchen "n" parameter (-)
 real(dp),intent(in)            :: vGn_alpha        ! van Genutchen "alpha" parameter (m-1)
 real(dp),intent(in)            :: theta_sat        ! soil porosity (-)
 real(dp),intent(in)            :: theta_res        ! soil residual volumetric water content (-)
 real(dp),intent(in)            :: snowfrz_scale    ! scaling parameter for the snow freezing curve (K-1)
 ! input
 real(dp),intent(in)            :: meltNrg          ! energy for melt+freeze (J m-3)
 real(dp),intent(in)            :: dLiq_dT          ! derivative of volumetric liquid water content w.r.t. temperature (K-1)
 real(dp),intent(in)            :: heatCap          ! volumetric heat capacity (J m-3 K-1)
 real(dp),intent(in)            :: tempInit         ! initial temperature (K)
 real(dp),intent(in)            :: volFracIceInit   ! initial volumetric fraction of ice (-)
 real(dp),intent(in)            :: massStateTrial   ! trial mass state variable
 ! output
 real(dp),intent(inout)         :: tempTrial        ! trial temperature (K)
 real(dp),intent(out)           :: volFracLiqTrial  ! volumetric fraction of liquid water (-)
 real(dp),intent(out)           :: volFracIceTrial  ! volumetric fraction of ice (-)
 real(dp),intent(out),optional  :: volFracWatTrial  ! volumetric fraction of total water (-)
 real(dp),intent(out),optional  :: fracLiqTrial     ! fraction liquid water (-)
 ! error handling
 integer(i4b),intent(out)       :: err              ! error code
 character(*),intent(out)       :: message          ! error message
 ! -----------------------------------------------------------------------------------------------------------------
 ! local variables
 integer(i4b)                   :: iter             ! iteration index
 integer(i4b),parameter         :: maxiter=20       ! iteration index
 character(len=256)             :: cMessage         ! error message of downwind routine
 ! -----------------------------------------------------------------------------------------------------------------
 ! initialize subroutine
 err=0; message='updateVars/'

 ! check arguments for the snow domain
 if(ixDomainType==iname_snow)then
  if(.not.present(fracLiqTrial))then; err=20; message=trim(message)//'optional argument "fracLiqTrial" expected for the snow domain';        return; endif
  if(present(volFracWatTrial)  )then; err=20; message=trim(message)//'optional argument "volFracWatTrial" NOT expected for the snow domain'; return; endif
 endif

 ! check arguments for the soil domain
 if(ixDomainType==iname_soil)then
  if(.not.present(volFracWatTrial))then; err=20; message=trim(message)//'optional argument "volFracWatTrial" expected for the soil domain';  return; endif
  if(present(fracLiqTrial)        )then; err=20; message=trim(message)//'optional argument "fracLiqTrial" NOT expected for the soil domain'; return; endif
 endif

 ! iterate
 do iter=1,maxiter

  ! -----
  ! - diaggregate total water content into liquid water and ice...
  ! --------------------------------------------------------------

  ! identify domain type
  ! NOTE: cycle statement for iname_cas above, so that should not be reached here
  select case(ixDomainType)

   ! vegetation canopy and snow layers
   case(iname_veg, iname_snow)  
    call updateSnow(tempTrial,                                  & ! intent(in)   : temperature (K)
                    massStateTrial,                             & ! intent(in)   : volumetric fraction of total water (-)
                    snowfrz_scale,                              & ! intent(in)   : scaling parameter for the snow freezing curve (K-1)
                    volFracLiqTrial,                            & ! intent(out)  : volumetric fraction of liquid water (-)
                    volFracIceTrial,                            & ! intent(out)  : volumetric fraction of ice (-)
                    fracLiqTrial,                               & ! intent(out)  : fraction of liquid water (-)
                    err,cmessage)                                 ! intent(out)  : error control
    if(err/=0)then; message=trim(message)//trim(cmessage); return; endif

   ! soil layers
   case(iname_soil)
    call updateSoil(tempTrial,                                  & ! intent(in)   : temperature (K)
                    massStateTrial,                             & ! intent(in)   : matric head (m)
                    vGn_alpha,vGn_n,theta_sat,theta_res,vGn_m,  & ! intent(in)   : soil parameters
                    volFracWatTrial,                            & ! intent(out)  : volumetric fraction of total water (-)
                    volFracLiqTrial,                            & ! intent(out)  : volumetric fraction of liquid water (-)
                    volFracIceTrial,                            & ! intent(out)  : volumetric fraction of ice (-)
                    err,cmessage)                                 ! intent(out)  : error control
    if(err/=0)then; message=trim(message)//trim(cmessage); return; endif
  
   ! check
   case default; err=20; message=trim(message)//'expect case to be iname_veg, iname_snow, iname_soil'; return

  end select  ! domain type

  ! -----
  ! - compute new temperature...
  ! ----------------------------

  ! check if the need to update temperatures
  if(.not.do_adjustTemp) return

  !print*, 'dTheta_dTk = ', dTheta_dTk
  print*, 'PAUSE: adjusting temperature'; read(*,*)

 end do  ! iterating

 end subroutine updateVars


 ! **********************************************************************************************************
 ! private subroutine indxSubset: get a subset of indices for a given mask
 ! **********************************************************************************************************
 subroutine indxSubset(ixSubset,ixMaster,mask,err,message)
 implicit none
 ! input-output: subset of indices for allocation/population
 integer(i4b),intent(inout),allocatable :: ixSubset(:)           ! subset of indices
 ! input
 integer(i4b),intent(in)                :: ixMaster(:)           ! full list of indices
 logical(lgt),intent(in)                :: mask(:)               ! desired indices
 ! error control
 integer(i4b),intent(out)               :: err                   ! error code
 character(*),intent(out)               :: message               ! error message
 ! -----------------------------------------------------------------------------------------------------------------------------------
 ! local variables
 integer(i4b)                           :: nSubset               ! length of the subset
 ! -----------------------------------------------------------------------------------------------------------------------------------
 ! initialize errors
 err=0; message="indxSubset/"

 ! check size match
 if(size(ixMaster)/=size(mask))then
  message=trim(message)//'size mismatch'
  err=20; return
 endif

 ! get the number of variables
 nSubset = count(mask)

 ! check if we need to reallocate space
 if(size(ixSubset)/=nSubset) then

  ! deallocate space
  deallocate(ixSubset,stat=err)
  if(err/=0)then; message=trim(message)//'unable to deallocate space for variable'; err=20; return; endif

  ! allocate space
  allocate(ixSubset(nSubset),stat=err)
  if(err/=0)then; message=trim(message)//'unable to deallocate space for variable'; err=20; return; endif

 endif  ! allocating space

 ! define indices for variable types in specific sub-domains
 if(nSubset>0) ixSubset = pack(ixMaster, mask)

 end subroutine indxSubset

end module getVectorz_module

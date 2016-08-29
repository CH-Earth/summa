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
USE globalData,only:iname_liqLayer  ! named variable defining the liquid  water state variable for snow+soil layers
USE globalData,only:iname_matLayer  ! named variable defining the matric head state variable for soil layers
USE globalData,only:iname_lmpLayer  ! named variable defining the liquid matric potential state variable for soil layers

! metadata for information in the data structures
USE globalData,only:indx_meta       ! metadata for the variables in the index structure

! constants
USE multiconst,only:&
                    gravity,      & ! acceleration of gravity              (m s-2)
                    Tfreeze,      & ! temperature at freezing              (K)
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
USE var_lookup,only:iLookDIAG             ! named variables for structure elements
USE var_lookup,only:iLookPROG             ! named variables for structure elements
USE var_lookup,only:iLookDERIV            ! named variables for structure elements
USE var_lookup,only:iLookPARAM            ! named variables for structure elements
USE var_lookup,only:iLookINDEX            ! named variables for structure elements

! provide access to routines to update states
USE updatState_module,only:updateSnow     ! update snow states
USE updatState_module,only:updateSoil     ! update soil states

! provide access to functions for the constitutive functions and derivatives
USE snow_utils_module,only:dFracLiq_dTk   ! differentiate the freezing curve w.r.t. temperature (snow)
USE soil_utils_module,only:dTheta_dTk     ! differentiate the freezing curve w.r.t. temperature (soil)
USE soil_utils_module,only:dTheta_dPsi    ! derivative in the soil water characteristic (soil)
USE soil_utils_module,only:dPsi_dTheta    ! derivative in the soil water characteristic (soil)
USE soil_utils_module,only:matricHead     ! compute the matric head based on volumetric water content
USE soil_utils_module,only:volFracLiq     ! compute volumetric fraction of liquid water
USE soil_utils_module,only:crit_soilT     ! compute critical temperature below which ice exists
USE soil_utils_module,only:liquidHead     ! compute the liquid water matric potential 

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
                        ! input: data structures
                        nState,                  & ! intent(in):    number of desired state variables
                        prog_data,               & ! intent(in):    model prognostic variables for a local HRU
                        diag_data,               & ! intent(in):    model diagnostic variables for a local HRU
                        indx_data,               & ! intent(in):    indices defining model states and layers
                        ! output
                        stateVec,                & ! intent(out):   model state vector
                        fScale,                  & ! intent(out):   function scaling vector (mixed units)
                        xScale,                  & ! intent(out):   variable scaling vector (mixed units)
                        sMul,                    & ! intent(out):   multiplier for state vector (used in the residual calculations)
                        dMat,                    & ! intent(out):   diagonal of the Jacobian matrix (excludes fluxes) 
                        err,message)               ! intent(out):   error control
 ! --------------------------------------------------------------------------------------------------------------------------------
 USE nr_utility_module,only:arth                   ! get a sequence of numbers arth(start, incr, count)
 USE f2008funcs_module,only:findIndex              ! finds the index of the first value within a vector
 ! --------------------------------------------------------------------------------------------------------------------------------
 ! input: data structures
 integer(i4b),intent(in)         :: nState                 ! number of desired state variables
 type(var_dlength),intent(in)    :: prog_data              ! prognostic variables for a local HRU
 type(var_dlength),intent(in)    :: diag_data              ! diagnostic variables for a local HRU
 type(var_ilength),intent(in)    :: indx_data              ! indices defining model states and layers
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
 ! state subsets
 integer(i4b)                    :: iLayer                 ! index of layer within the snow+soil domain
 integer(i4b)                    :: ixStateSubset          ! index within the state subset
 logical(lgt),dimension(nState)  :: stateFlag              ! flag to denote that the state is populated
 ! --------------------------------------------------------------------------------------------------------------------------------
 ! --------------------------------------------------------------------------------------------------------------------------------
 ! make association with variables in the data structures
 fixedLength: associate(&
 ! model states for the vegetation canopy
 scalarCanairTemp    => prog_data%var(iLookPROG%scalarCanairTemp)%dat(1)       ,& ! intent(in) : [dp]     temperature of the canopy air space (K)
 scalarCanopyTemp    => prog_data%var(iLookPROG%scalarCanopyTemp)%dat(1)       ,& ! intent(in) : [dp]     temperature of the vegetation canopy (K)
 scalarCanopyWat     => prog_data%var(iLookPROG%scalarCanopyWat)%dat(1)        ,& ! intent(in) : [dp]     mass of total water on the vegetation canopy (kg m-2)
 ! model state variable vectors for the snow-soil layers
 mLayerTemp          => prog_data%var(iLookPROG%mLayerTemp)%dat                ,& ! intent(in) : [dp(:)]  temperature of each snow/soil layer (K)
 mLayerVolFracWat    => prog_data%var(iLookPROG%mLayerVolFracWat)%dat          ,& ! intent(in) : [dp(:)]  volumetric fraction of total water (-)
 mLayerVolFracLiq    => prog_data%var(iLookPROG%mLayerVolFracLiq)%dat          ,& ! intent(in) : [dp(:)]  volumetric fraction of liquid water (-)
 mLayerMatricHead    => prog_data%var(iLookPROG%mLayerMatricHead)%dat          ,& ! intent(in) : [dp(:)]  matric head (m)
 ! model diagnostic variables
 canopyDepth         => diag_data%var(iLookDIAG%scalarCanopyDepth)%dat(1)      ,& ! intent(in) : [dp]     canopy depth (m)
 volHeatCapVeg       => diag_data%var(iLookDIAG%scalarBulkVolHeatCapVeg)%dat(1),& ! intent(in) : [dp]     bulk volumetric heat capacity of vegetation (J m-3 K-1)
 mLayerVolHeatCap    => diag_data%var(iLookDIAG%mLayerVolHtCapBulk)%dat        ,& ! intent(in) : [dp(:)]  bulk volumetric heat capacity in each snow and soil layer (J m-3 K-1)
 mLayerMatricHeadLiq => diag_data%var(iLookDIAG%mLayerMatricHeadLiq)%dat       ,& ! intent(in) : [dp(:)]  matric potential of liquid water (m)
 ! indices defining specific model states
 ixCasNrg            => indx_data%var(iLookINDEX%ixCasNrg)%dat(1)              ,& ! intent(in) : [i4b]    index of canopy air space energy state variable
 ixVegNrg            => indx_data%var(iLookINDEX%ixVegNrg)%dat(1)              ,& ! intent(in) : [i4b]    index of canopy energy state variable
 ixVegWat            => indx_data%var(iLookINDEX%ixVegWat)%dat(1)              ,& ! intent(in) : [i4b]    index of canopy hydrology state variable (mass)
 ! vector of energy and hydrology indices for the snow and soil domains
 ixSnowSoilNrg       => indx_data%var(iLookINDEX%ixSnowSoilNrg)%dat            ,& ! intent(in) : [i4b(:)] index in the state subset for energy state variables in the snow+soil domain
 ixSnowSoilHyd       => indx_data%var(iLookINDEX%ixSnowSoilHyd)%dat            ,& ! intent(in) : [i4b(:)] index in the state subset for hydrology state variables in the snow+soil domain
 nSnowSoilNrg        => indx_data%var(iLookINDEX%nSnowSoilNrg )%dat(1)         ,& ! intent(in) : [i4b]    number of energy state variables in the snow+soil domain
 nSnowSoilHyd        => indx_data%var(iLookINDEX%nSnowSoilNrg )%dat(1)         ,& ! intent(in) : [i4b]    number of hydrology state variables in the snow+soil domain
 ! type of hydrology states in the snow+soil domain
 ixStateType_subset  => indx_data%var(iLookINDEX%ixStateType_subset)%dat       ,& ! intent(in) : [i4b(:)] [state subset] type of desired model state variables
 ixHydType           => indx_data%var(iLookINDEX%ixHydType)%dat                ,& ! intent(in) : [i4b(:)] index of the type of hydrology states in snow+soil domain
 ! number of layers
 nSnow               => indx_data%var(iLookINDEX%nSnow)%dat(1)                 ,& ! intent(in) : [i4b]    number of snow layers
 nSoil               => indx_data%var(iLookINDEX%nSoil)%dat(1)                 ,& ! intent(in) : [i4b]    number of soil layers
 nLayers             => indx_data%var(iLookINDEX%nLayers)%dat(1)                & ! intent(in) : [i4b]    total number of layers
 )  ! end association with variables in the data structures
 ! --------------------------------------------------------------------------------------------------------------------------------
 ! initialize error control
 err=0; message='popStateVec/'

 ! -----
 ! * initialize state vectors...
 ! -----------------------------

 ! initialize flags
 stateFlag(:) = .false.

 ! build elements of the state vector for the vegetation canopy
 if(ixCasNrg/=integerMissing) stateVec(ixCasNrg) = scalarCanairTemp
 if(ixVegNrg/=integerMissing) stateVec(ixVegNrg) = scalarCanopyTemp
 if(ixVegWat/=integerMissing) stateVec(ixVegWat) = scalarCanopyWat  ! kg m-2

 ! build the energy state vector for the snow and soil domain
 if(nSnowSoilNrg>0)then
  do concurrent (iLayer=1:nLayers,ixSnowSoilNrg(iLayer)/=integerMissing)   ! (loop through non-missing energy state variables in the snow+soil domain)
   ixStateSubset            = ixSnowSoilNrg(iLayer)  ! index within the state vector
   stateVec(ixStateSubset)  = mLayerTemp(iLayer)     ! transfer temperature from a layer to the state vector
   stateFlag(ixStateSubset) = .true.                 ! flag to denote that the state is populated
  end do  ! looping through non-missing energy state variables in the snow+soil domain
 endif

 ! build the hydrology state vector for the snow+soil domains
 ! NOTE: ixVolFracWat  and ixVolFracLiq can also include states in the soil domain, hence enable primary variable switching
 if(nSnowSoilHyd>0)then
  do concurrent (iLayer=1:nLayers,ixSnowSoilHyd(iLayer)/=integerMissing)   ! (loop through non-missing hydrology state variables in the snow+soil domain)
   ixStateSubset                                  = ixSnowSoilHyd(iLayer)   ! index within the state vector
   select case( ixHydType(iLayer) )
    case(iname_watLayer); stateVec(ixStateSubset) = mLayerVolFracWat(iLayer)           ! total water state variable for snow+soil layers
    case(iname_liqLayer); stateVec(ixStateSubset) = mLayerVolFracLiq(iLayer)           ! liquid water state variable for snow+soil layers
    case(iname_matLayer); stateVec(ixStateSubset) = mLayerMatricHead(iLayer-nSnow)     ! total water matric potential variable for soil layers
    case(iname_lmpLayer); stateVec(ixStateSubset) = mLayerMatricHeadLiq(iLayer-nSnow)  ! liquid matric potential state variable for soil layers
    case default; cycle
   end select
   stateFlag(ixStateSubset) = .true. ! flag to denote that the state is populated
  end do  ! looping through non-missing energy state variables in the snow+soil domain
 endif

 ! check that we populated all state variables
 if(count(stateFlag)/=nState)then
  message=trim(message)//'some state variables unpopulated'
  err=20; return
 endif

 ! -----
 ! * define scaling vectors...
 ! ---------------------------

 ! define the function and variable scaling factors for energy
 where(ixStateType_subset==iname_nrgCanair .or. ixStateType_subset==iname_nrgCanopy .or. ixStateType_subset==iname_nrgLayer)
  fScale = 1._dp / fScaleNrg  ! 1/(J m-3)
  xScale = 1._dp  ! K
 endwhere

 ! define the function and variable scaling factors for water on the vegetation canopy
 where(ixStateType_subset==iname_watCanopy)
  fScale = 1._dp / (fScaleLiq*canopyDepth*iden_water)  ! 1/(kg m-2)
  xScale = 1._dp  ! (kg m-2)
 endwhere

 ! define the function and variable scaling factors for water in the snow+soil domain
 where(ixStateType_subset==iname_watLayer .or. ixStateType_subset==iname_liqLayer)
  fScale = 1._dp / fScaleLiq  ! (-)
  xScale = 1._dp  ! (-)
 end where
 
 ! define the function and variable scaling factors for water in the snow+soil domain
 where(ixStateType_subset==iname_matLayer .or. ixStateType_subset==iname_lmpLayer)
  fScale = 1._dp / fScaleLiq  ! (-)
  xScale = 1._dp  ! (m)
 end where

 ! -----
 ! * define components of derivative matrices that are constant over a time step (substep)...
 ! ------------------------------------------------------------------------------------------

 ! define the multiplier for the state vector for residual calculations (vegetation canopy)
 if(ixCasNrg/=integerMissing) sMul(ixCasNrg) = Cp_air*iden_air   ! volumetric heat capacity of air (J m-3 K-1)
 if(ixVegNrg/=integerMissing) sMul(ixVegNrg) = volHeatCapVeg     ! volumetric heat capacity of the vegetation (J m-3 K-1)
 if(ixVegWat/=integerMissing) sMul(ixVegWat) = 1._dp             ! nothing else on the left hand side

 ! compute terms in the Jacobian for vegetation (excluding fluxes)
 ! NOTE: this is computed outside the iteration loop because it does not depend on state variables
 ! NOTE: energy for vegetation is computed *within* the iteration loop as it includes phase change
 if(ixCasNrg/=integerMissing) dMat(ixCasNrg) = Cp_air*iden_air   ! volumetric heat capacity of air (J m-3 K-1)
 if(ixVegNrg/=integerMissing) dMat(ixVegNrg) = realMissing       ! populated within the iteration loop 
 if(ixVegWat/=integerMissing) dMat(ixVegWat) = 1._dp             ! nothing else on the left hand side

 ! define the energy multiplier and diagonal elements for the state vector for residual calculations (snow-soil domain)
 if(nSnowSoilNrg>0)then
  do concurrent (iLayer=1:nLayers,ixSnowSoilNrg(iLayer)/=integerMissing)   ! (loop through non-missing energy state variables in the snow+soil domain)
   ixStateSubset        = ixSnowSoilNrg(iLayer)      ! index within the state vector
   sMul(ixStateSubset)  = mLayerVolHeatCap(iLayer)   ! transfer volumetric heat capacity to the state multiplier
   dMat(ixStateSubset)  = realMissing                ! diagonal element populated within the iteration loop
  end do  ! looping through non-missing energy state variables in the snow+soil domain
 endif

 ! define the hydrology multiplier and diagonal elements for the state vector for residual calculations (snow-soil domain)
 if(nSnowSoilHyd>0)then
  do concurrent (iLayer=1:nLayers,ixSnowSoilHyd(iLayer)/=integerMissing)   ! (loop through non-missing energy state variables in the snow+soil domain)
   ixStateSubset        = ixSnowSoilHyd(iLayer)      ! index within the state vector
   sMul(ixStateSubset)  = 1._dp                      ! state multiplier = 1 (nothing else on the left-hand-side) 
   dMat(ixStateSubset)  = 1._dp                      ! diagonal element = 1 (nothing else on the left-hand-side) 
  end do  ! looping through non-missing energy state variables in the snow+soil domain
 endif

 ! ------------------------------------------------------------------------------------------
 ! ------------------------------------------------------------------------------------------
 
 end associate fixedLength      ! end association to variables in the data structure where vector length does not change
 end subroutine popStateVec




 ! **********************************************************************************************************
 ! public subroutine varExtract: extract variables from the state vector and compute diagnostic variables
 ! **********************************************************************************************************
 subroutine varExtract(&
                       ! input
                       do_adjustTemp,                             & ! intent(in):    logical flag to adjust temperature to account for the energy used in melt+freeze
                       stateVec,                                  & ! intent(in):    model state vector (mixed units)
                       mpar_data,                                 & ! intent(in):    model parameters for a local HRU
                       diag_data,                                 & ! intent(in):    model diagnostic variables for a local HRU
                       prog_data,                                 & ! intent(in):    model prognostic variables for a local HRU
                       indx_data,                                 & ! intent(in):    indices defining model states and layers
                       deriv_data,                                & ! intent(inout): derivatives in model fluxes w.r.t. relevant state variables
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
                       mLayerMatricHeadTrial,                     & ! intent(out):   trial vector of total water matric potential (m)
                       mLayerMatricHeadLiqTrial,                  & ! intent(out):   trial vector of liquid water matric potential (m)
                       ! output: error control 
                       err,message)                                 ! intent(out):   error control
 ! --------------------------------------------------------------------------------------------------------------------------------
 ! --------------------------------------------------------------------------------------------------------------------------------
 implicit none 
 ! input
 logical(lgt),intent(in)         :: do_adjustTemp                   ! flag to adjust temperature to account for the energy used in melt+freeze
 real(dp),intent(in)             :: stateVec(:)                     ! model state vector (mixed units)
 type(var_d),      intent(in)    :: mpar_data                       ! model parameters for a local HRU
 type(var_dlength),intent(in)    :: diag_data                       ! diagnostic variables for a local HRU
 type(var_dlength),intent(in)    :: prog_data                       ! prognostic variables for a local HRU
 type(var_ilength),intent(in)    :: indx_data                       ! indices defining model states and layers                 
 type(var_dlength),intent(inout) :: deriv_data                      ! derivatives in model fluxes w.r.t. relevant state variables
 ! output: variables for the vegetation canopy
 real(dp),intent(out)            :: fracLiqVeg                      ! fraction of liquid water on the vegetation canopy (-)
 real(dp),intent(out)            :: scalarCanairTempTrial           ! trial value of canopy air temperature (K)
 real(dp),intent(out)            :: scalarCanopyTempTrial           ! trial value of canopy temperature (K)
 real(dp),intent(out)            :: scalarCanopyWatTrial            ! trial value of canopy total water (kg m-2)
 real(dp),intent(out)            :: scalarCanopyLiqTrial            ! trial value of canopy liquid water (kg m-2)
 real(dp),intent(out)            :: scalarCanopyIceTrial            ! trial value of canopy ice content (kg m-2)
 ! output: variables for the snow-soil domain
 real(dp),intent(out)            :: fracLiqSnow(:)                  ! volumetric fraction of water in each snow layer (-)
 real(dp),intent(out)            :: mLayerTempTrial(:)              ! trial vector of layer temperature (K)
 real(dp),intent(out)            :: mLayerVolFracWatTrial(:)        ! trial vector of volumetric total water content (-)
 real(dp),intent(out)            :: mLayerVolFracLiqTrial(:)        ! trial vector of volumetric liquid water content (-)
 real(dp),intent(out)            :: mLayerVolFracIceTrial(:)        ! trial vector of volumetric ice water content (-)
 real(dp),intent(out)            :: mLayerMatricHeadTrial(:)        ! trial vector of total water matric potential (m)
 real(dp),intent(out)            :: mLayerMatricHeadLiqTrial(:)     ! trial vector of liquid water matric potential (m)
 ! output: error control 
 integer(i4b),intent(out)        :: err                             ! error code
 character(*),intent(out)        :: message                         ! error message
 ! --------------------------------------------------------------------------------------------------------------------------------
 ! local variables
 integer(i4b)                    :: iState                          ! index of model state variable
 integer(i4b)                    :: iLayer                          ! index of layer within the snow+soil domain
 integer(i4b)                    :: ixFullVector                    ! index within full state vector
 integer(i4b)                    :: ixDomainType                    ! name of a given model domain
 integer(i4b)                    :: ixControlIndex                  ! index within a given model domain
 integer(i4b)                    :: ixOther,ixOtherLocal            ! index of the coupled state variable within the (full, local) vector
 logical(lgt)                    :: isCoupled                       ! .true. if a given variable shared another state variable in the same control volume
 logical(lgt)                    :: isNrgState                      ! .true. if a given variable is an energy state 
 logical(lgt),allocatable        :: computedCoupling(:)             ! .true. if computed the coupling for a given state variable
 real(dp)                        :: scalarVolFracLiq                ! volumetric fraction of liquid water (-)
 real(dp)                        :: scalarVolFracIce                ! volumetric fraction of ice (-)
 real(dp)                        :: Tcrit                           ! critical soil temperature below which ice exists (K)
 real(dp)                        :: xTemp                           ! temporary temperature (K)
 character(len=256)              :: cMessage                        ! error message of downwind routine
 logical(lgt),parameter          :: printFlag=.false.               ! flag to turn on printing
 ! --------------------------------------------------------------------------------------------------------------------------------
 ! make association with variables in the data structures
 associate(&
 ! number of model layers, and layer type
 nSnow                   => indx_data%var(iLookINDEX%nSnow)%dat(1)                 ,& ! intent(in):  [i4b]    total number of snow layers
 nSoil                   => indx_data%var(iLookINDEX%nSoil)%dat(1)                 ,& ! intent(in):  [i4b]    total number of soil layers
 nLayers                 => indx_data%var(iLookINDEX%nLayers)%dat(1)               ,& ! intent(in):  [i4b]    total number of snow and soil layers
 nVegState               => indx_data%var(iLookINDEX%nVegState)%dat(1)             ,& ! intent(in):  [i4b]    number of vegetation state variables
 layerType               => indx_data%var(iLookINDEX%layerType)%dat                ,& ! intent(in):  [i4b(:)] index defining type of layer (soil or snow)
 ! indices defining model states and layers
 ixCasNrg                => indx_data%var(iLookINDEX%ixCasNrg)%dat(1)              ,& ! intent(in):  [i4b]    index of canopy air space energy state variable
 ixVegNrg                => indx_data%var(iLookINDEX%ixVegNrg)%dat(1)              ,& ! intent(in):  [i4b]    index of canopy energy state variable
 ixVegWat                => indx_data%var(iLookINDEX%ixVegWat)%dat(1)              ,& ! intent(in):  [i4b]    index of canopy hydrology state variable (mass)
 ixSnowSoilNrg           => indx_data%var(iLookINDEX%ixSnowSoilNrg)%dat            ,& ! intent(in):  [i4b(:)] indices IN THE STATE SUBSET for energy states in the snow+soil subdomain
 ixSnowSoilHyd           => indx_data%var(iLookINDEX%ixSnowSoilNrg)%dat            ,& ! intent(in):  [i4b(:)] indices IN THE STATE SUBSET for hydrology states in the snow+soil subdomain
 nSnowSoilNrg            => indx_data%var(iLookINDEX%nSnowSoilNrg )%dat(1)         ,& ! intent(in):  [i4b]    number of energy state variables in the snow+soil domain
 nSnowSoilHyd            => indx_data%var(iLookINDEX%nSnowSoilHyd )%dat(1)         ,& ! intent(in):  [i4b]    number of hydrology variables in the snow+soil domain
 ! indices foe states within the snow+soil domain
 ixVolFracWat            => indx_data%var(iLookINDEX%ixVolFracWat)%dat             ,& ! intent(in):  [i4b(:)] indices IN THE SNOW+SOIL VECTOR for hydrology states in the snow+soil subdomain
 ixMatricHead            => indx_data%var(iLookINDEX%ixMatricHead)%dat             ,& ! intent(in):  [i4b(:)] indices IN THE SOIL VECTOR for hydrology states in the soil subdomain
 ixHydType               => indx_data%var(iLookINDEX%ixHydType)%dat                ,& ! intent(in) : [i4b(:)] index of the type of hydrology states in snow+soil domain
 ! model diagnostic variables
 canopyDepth             => diag_data%var(iLookDIAG%scalarCanopyDepth)%dat(1)      ,& ! intent(in):  [dp   ] canopy depth (m)
 scalarBulkVolHeatCapVeg => diag_data%var(iLookDIAG%scalarBulkVolHeatCapVeg)%dat(1),& ! intent(in):  [dp   ] volumetric heat capacity of the vegetation (J m-3 K-1)
 mLayerVolHtCapBulk      => diag_data%var(iLookDIAG%mLayerVolHtCapBulk)%dat        ,& ! intent(in):  [dp(:)] volumetric heat capacity in each layer (J m-3 K-1)
 ! model states for the vegetation canopy
 scalarCanairTemp        => prog_data%var(iLookPROG%scalarCanairTemp)%dat(1)       ,& ! intent(in):  [dp] temperature of the canopy air space (K)
 scalarCanopyTemp        => prog_data%var(iLookPROG%scalarCanopyTemp)%dat(1)       ,& ! intent(in):  [dp] temperature of the vegetation canopy (K)
 scalarCanopyWat         => prog_data%var(iLookPROG%scalarCanopyWat)%dat(1)        ,& ! intent(in):  [dp] mass of total water on the vegetation canopy (kg m-2)
 ! model state variable vectors for the snow-soil layers
 mLayerTemp              => prog_data%var(iLookPROG%mLayerTemp)%dat                ,& ! intent(in):  [dp(:)] temperature of each snow/soil layer (K)
 mLayerVolFracWat        => prog_data%var(iLookPROG%mLayerVolFracWat)%dat          ,& ! intent(in):  [dp(:)] volumetric fraction of total water (-)
 mLayerMatricHead        => prog_data%var(iLookPROG%mLayerMatricHead)%dat          ,& ! intent(in):  [dp(:)] total water matric potential (m)
 mLayerMatricHeadLiq     => diag_data%var(iLookDIAG%mLayerMatricHeadLiq)%dat       ,& ! intent(in):  [dp(:)] liquid water matric potential (m)
 ! model diagnostic variables from a previous solution
 scalarCanopyLiq         => prog_data%var(iLookPROG%scalarCanopyLiq)%dat(1)        ,& ! intent(in):  [dp(:)] mass of liquid water on the vegetation canopy (kg m-2)
 scalarCanopyIce         => prog_data%var(iLookPROG%scalarCanopyIce)%dat(1)        ,& ! intent(in):  [dp(:)] mass of ice on the vegetation canopy (kg m-2)
 mLayerVolFracLiq        => prog_data%var(iLookPROG%mLayerVolFracLiq)%dat          ,& ! intent(in):  [dp(:)] volumetric fraction of liquid water (-)
 mLayerVolFracIce        => prog_data%var(iLookPROG%mLayerVolFracIce)%dat          ,& ! intent(in):  [dp(:)] volumetric fraction of ice (-)
 ! derivatives
 dVolTot_dPsi0           => deriv_data%var(iLookDERIV%dVolTot_dPsi0   )%dat        ,& ! intent(out): [dp(:)] derivative in total water content w.r.t. total water matric potential
 dPsiLiq_dPsi0           => deriv_data%var(iLookDERIV%dPsiLiq_dPsi0   )%dat        ,& ! intent(out): [dp(:)] derivative in liquid water matric pot w.r.t. the total water matric pot (-)
 dPsiLiq_dTemp           => deriv_data%var(iLookDERIV%dPsiLiq_dTemp   )%dat        ,& ! intent(out): [dp(:)] derivative in the liquid water matric potential w.r.t. temperature
 mLayerdTheta_dTk        => deriv_data%var(iLookDERIV%mLayerdTheta_dTk)%dat        ,& ! intent(out): [dp(:)] derivative of volumetric liquid water content w.r.t. temperature
 dTheta_dTkCanopy        => deriv_data%var(iLookDERIV%dTheta_dTkCanopy)%dat(1)      & ! intent(out): [dp]    derivative of volumetric liquid water content w.r.t. temperature
 ) ! association with variables in the data structures

 ! --------------------------------------------------------------------------------------------------------------------------------
 ! --------------------------------------------------------------------------------------------------------------------------------

 ! initialize error control
 err=0; message='varExtract/'

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
 mLayerTempTrial          = mLayerTemp
 mLayerVolFracWatTrial    = mLayerVolFracWat
 mLayerVolFracLiqTrial    = mLayerVolFracLiq
 mLayerVolFracIceTrial    = mLayerVolFracIce
 mLayerMatricHeadTrial    = mLayerMatricHead      ! total water matric potential
 mLayerMatricHeadLiqTrial = mLayerMatricHeadLiq   ! liquid water matric potential

 ! overwrite with the energy values from the state vector
 if(nSnowSoilNrg>0)then
  do concurrent (iLayer=1:nLayers,ixSnowSoilNrg(iLayer)/=integerMissing)   ! (loop through non-missing energy state variables in the snow+soil domain)
   mLayerTempTrial(iLayer) = stateVec( ixSnowSoilNrg(iLayer) ) 
  end do  ! looping through non-missing energy state variables in the snow+soil domain
 endif

 ! overwrite with the energy values from the state vector
 ! NOTE: ixVolFracWat  and ixVolFracLiq can also include states in the soil domain, hence enable primary variable switching
 if(nSnowSoilHyd>0)then
  do concurrent (iLayer=1:nLayers,ixSnowSoilHyd(iLayer)/=integerMissing)   ! (loop through non-missing hydrology state variables in the snow+soil domain)
   select case( ixHydType(iLayer) )
    case(iname_watLayer); mLayerVolFracWatTrial(iLayer)          = stateVec( ixSnowSoilHyd(iLayer) ) ! total water state variable for snow+soil layers
    case(iname_liqLayer); mLayerVolFracLiqTrial(iLayer)          = stateVec( ixSnowSoilHyd(iLayer) ) ! liquid water state variable for snow+soil layers
    case(iname_matLayer); mLayerMatricHeadTrial(iLayer-nSnow)    = stateVec( ixSnowSoilHyd(iLayer) ) ! total water matric potential variable for soil layers
    case(iname_lmpLayer); mLayerMatricHeadLiqTrial(iLayer-nSnow) = stateVec( ixSnowSoilHyd(iLayer) ) ! liquid matric potential state variable for soil layers
    case default; cycle
   end select
  end do  ! looping through non-missing energy state variables in the snow+soil domain
 endif

 ! --------------------------------------------------------------------------------------------------------------------------------
 ! --------------------------------------------------------------------------------------------------------------------------------

 ! *****
 ! * part 2: update volumetric fraction of liquid water and ice (and also update temperature if desired)... 
 ! ********************************************************************************************************

 ! -----
 ! - preliminaries...
 ! ------------------

 ! make association with indices
 indices: associate(&
  ixNrgLayer          => indx_data%var(iLookINDEX%ixNrgLayer)%dat             ,& ! intent(in): [i4b(:)] indices IN THE FULL VECTOR for energy states in the snow+soil domain
  ixHydLayer          => indx_data%var(iLookINDEX%ixHydLayer)%dat             ,& ! intent(in): [i4b(:)] indices IN THE FULL VECTOR for hydrology states in the snow+soil domain
  ixMapFull2Subset    => indx_data%var(iLookINDEX%ixMapFull2Subset)%dat       ,& ! intent(in): [i4b(:)] list of indices in the state subset for each state in the full state vector
  ixMapSubset2Full    => indx_data%var(iLookINDEX%ixMapSubset2Full)%dat       ,& ! intent(in): [i4b(:)] [state subset] list of indices of the full state vector in the state subset
  ixDomainType_subset => indx_data%var(iLookINDEX%ixDomainType_subset)%dat    ,& ! intent(in): [i4b(:)] [state subset] id of domain for desired model state variables
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
 if( size(ixDomainType_subset)/=size(stateVec) )then
  message=trim(message)//'size mismatch for state subset'
  err=20; return
 endif

 ! allocate space and assign values to the flag vector
 allocate(computedCoupling(size(stateVec)),stat=err)        ! .true. if computed the coupling for a given state variable
 if(err/=0)then; message=trim(message)//'problem allocating computedCoupling'; return; endif
 computedCoupling(:)=.false.

 ! loop through model state variables
 do iState=1,size(stateVec)

  ! check the need for the computations
  if(computedCoupling(iState)) cycle

  ! -----
  ! - compute indices...
  ! --------------------

  ! get domain type, and index of the control volume within the domain
  ixFullVector   = ixMapSubset2Full(iState)       ! index within full state vector
  ixDomainType   = ixDomainType_subset(iState)    ! named variables defining the domain (iname_cas, iname_veg, etc.)
  ixControlIndex = ixControlVolume(ixFullVector)  ! index within a given domain

  ! get the layer index
  select case(ixDomainType)
   case(iname_cas);  cycle ! canopy air space: do nothing
   case(iname_veg);  iLayer = 0 
   case(iname_snow); iLayer = ixControlIndex
   case(iname_soil); iLayer = ixControlIndex + nSnow
   case default; err=20; message=trim(message)//'expect case to be iname_cas, iname_veg, iname_snow, iname_soil'; return
  end select
   
  ! get the index of the other (energy or mass) state variable within the full state vector
  select case(ixDomainType)
   case(iname_veg)             ; ixOther = merge(ixVegWat,ixVegNrg,ixStateType(ixFullVector)==iname_nrgCanopy)  ! ixVegWat if stateType=iname_nrgCanopy; ixVegNrg otherwise
   case(iname_snow, iname_soil); ixOther = merge(ixHydLayer(iLayer),ixNrgLayer(iLayer),ixStateType(ixFullVector)==iname_nrgLayer)
   case default; err=20; message=trim(message)//'expect case to be iname_veg, iname_snow, iname_soil'; return
  end select

  ! get the index in the local state vector
  ixOtherLocal = ixMapFull2Subset(ixOther)  ! ixOtherLocal could equal integerMissing
  if(ixOtherLocal/=integerMissing) computedCoupling(ixOtherLocal)=.true.

  ! check if we have a coupled solution
  isCoupled    = (ixOtherLocal/=integerMissing)

  ! check if we are an energy state
  isNrgState   = (ixStateType(ixFullVector)==iname_nrgCanopy .or. ixStateType(ixFullVector)==iname_nrgLayer)

  if(printFlag)then
   print*, 'ixFullVector   = ', ixFullVector
   print*, 'ixDomainType   = ', ixDomainType
   print*, 'ixControlIndex = ', ixControlIndex
   print*, 'ixOther        = ', ixOther
   print*, 'ixOtherLocal   = ', ixOtherLocal
   print*, 'do_adjustTemp  = ', do_adjustTemp
   print*, 'isCoupled      = ', isCoupled
   print*, 'isNrgState     = ', isNrgState
  endif

  ! -----
  ! - compute derivatives...
  ! ------------------------

  ! compute the derivative in total water content w.r.t. total water matric potential (m-1)
  ! NOTE: valid for frozen and unfrozen conditions
  if(ixDomainType==iname_soil) dVolTot_dPsi0(ixControlIndex) = dTheta_dPsi(mLayerMatricHeadTrial(ixControlIndex),vGn_alpha,theta_res,theta_sat,vGn_n,vGn_m)

  ! special case of hydrology state uncoupled with energy
  !  -- just compute the temperature derivatives and return
  if(.not.do_adjustTemp .and. .not.isNrgState .and. .not.isCoupled)then

   ! (compute the critical soil temperature below which ice exists)
   Tcrit = merge(crit_soilT(mLayerMatricHeadTrial(ixControlIndex)), Tfreeze, ixDomainType==iname_soil)
   xTemp = merge(scalarCanopyTemp, mLayerTempTrial(iLayer), ixDomainType==iname_veg)

   ! (compute the derivative in liquid water content w.r.t. temperature)
   select case(ixDomainType)
    case(iname_veg);  dTheta_dTkCanopy         = merge(dFracLiq_dTk(xTemp,snowfrz_scale)*scalarCanopyWat/(iden_water*canopyDepth),  0._dp, xTemp<Tcrit)
    case(iname_snow); mLayerdTheta_dTk(iLayer) = merge(dFracLiq_dTk(xTemp,snowfrz_scale)*mLayerVolFracWatTrial(iLayer),             0._dp, xTemp<Tcrit)
    case(iname_soil); mLayerdTheta_dTk(iLayer) = merge(dTheta_dTk(xTemp,theta_res,theta_sat,vGn_alpha,vGn_n,vGn_m),                 0._dp, xTemp<Tcrit) ! assume no volume expansion
    case default; err=20; message=trim(message)//'expect case to be iname_veg, iname_snow, iname_soil'; return
   end select  ! domain type

   ! (early return)
   return  ! do not update the volume fractions of liquid water and ice, and do not update liquid water matric potential
  endif

  ! -----
  ! - update volumetric fraction of liquid water and ice (and also update temperature if desired)...
  ! ------------------------------------------------------------------------------------------------

  ! get desired variables for each control volume
  select case(ixDomainType)
  
   ! vegetation
   case(iname_veg)
    call updateVars(&
                    ! model control
                    ixDomainType,                                               & ! intent(in)    : named variable defining the domain type (iname_veg, iname_snow, iname_soil)
                    (do_adjustTemp .and. .not.isCoupled .and. .not.isNrgState), & ! intent(in)    : logical flag to adjust temperature to account for energy associated with melt+freeze
                    snowfrz_scale,                                              & ! intent(in)    : scaling parameter for the snow freezing curve (K-1)
                    vGn_alpha,vGn_n,theta_sat,theta_res,vGn_m,                  & ! intent(in)    : soil parameters
                    ! input
                    LH_fus*iden_water                                          ,& ! intent(in)    : energy for melt+freeze (J m-3)
                    scalarBulkVolHeatCapVeg                                    ,& ! intent(in)    : volumetric heat capacity (J m-3 K-1)
                    scalarCanopyTemp                                           ,& ! intent(in)    : initial temperature (K) 
                    scalarCanopyIce/(iden_water*canopyDepth)                   ,& ! intent(in)    : initial volumetric fraction of ice (-)
                    scalarCanopyWatTrial/(iden_water*canopyDepth)              ,& ! intent(in)    : mass state variable = trial volumetric fraction of total water (-)
                    ! output 
                    scalarCanopyTempTrial                                      ,& ! intent(inout) : trial temperature (K)
                    scalarVolFracLiq                                           ,& ! intent(out)   : trial volumetric fraction of liquid water (-)
                    scalarVolFracIce                                           ,& ! intent(out)   : trial volumetric fraction if ice (-)
                    fracLiqTrial=fracLiqVeg                                    ,& ! intent(out)   : fraction of liquid water (-)
                    dLiq_dT=dTheta_dTkCanopy                                   ,& ! intent(out)   : derivative in liquid water w.r.t. temperature (K-1) 
                    err=err,message=cmessage)                                     ! intent(out)   : error control
    if(err/=0)then; message=trim(message)//trim(cmessage); return; endif

    ! compute mass of water on the canopy
    scalarCanopyLiqTrial = scalarVolFracLiq*(iden_water*canopyDepth)
    scalarCanopyIceTrial = scalarVolFracIce*(iden_water*canopyDepth)

   ! snow layers
   case(iname_snow)
    call updateVars(&
                    ! model control
                    ixDomainType,                                               & ! intent(in)    : named variable defining the domain type (iname_veg, iname_snow, iname_soil)
                    (do_adjustTemp .and. .not.isCoupled .and. .not.isNrgState), & ! intent(in)    : logical flag to adjust temperature to account for energy associated with melt+freeze
                    snowfrz_scale,                                              & ! intent(in)    : scaling parameter for the snow freezing curve (K-1)
                    vGn_alpha,vGn_n,theta_sat,theta_res,vGn_m,                  & ! intent(in)    : soil parameters
                    ! input
                    LH_fus*iden_ice                                            ,& ! intent(in)    : energy for melt+freeze (J m-3)
                    mLayerVolHtCapBulk(iLayer)                                 ,& ! intent(in)    : volumetric heat capacity (J m-3 K-1)
                    mLayerTemp(iLayer)                                         ,& ! intent(in)    : initial temperature (K)
                    mLayerVolFracIce(iLayer)                                   ,& ! intent(in)    : initial volumetric fraction of ice (-)
                    mLayerVolFracWatTrial(iLayer)                              ,& ! intent(in)    : mass state variable = trial volumetric fraction of water (-)
                    ! output 
                    mLayerTempTrial(iLayer)                                    ,& ! intent(inout) : trial temperature (K)
                    mLayerVolFracLiqTrial(iLayer)                              ,& ! intent(out)   : trial volumetric fraction of liquid water (-)
                    mLayerVolFracIceTrial(iLayer)                              ,& ! intent(out)   : trial volumetric fraction if ice (-)
                    fracLiqTrial=fracLiqSnow(iLayer)                           ,& ! intent(out)   : fraction of liquid water (-)
                    dLiq_dT=mLayerdTheta_dTk(iLayer)                           ,& ! intent(out)   : derivative in liquid water w.r.t. temperature (K-1) 
                    err=err,message=cmessage)                                     ! intent(out)   : error control
    if(err/=0)then; message=trim(message)//trim(cmessage); return; endif

   ! soil layers
   case(iname_soil)
    call updateVars(&
                    ! model control
                    ixDomainType,                                               & ! intent(in)    : named variable defining the domain type (iname_veg, iname_snow, iname_soil)
                    (do_adjustTemp .and. .not.isCoupled .and. .not.isNrgState), & ! intent(in)    : logical flag to adjust temperature to account for energy associated with melt+freeze
                    snowfrz_scale,                                              & ! intent(in)    : scaling parameter for the snow freezing curve (K-1)
                    vGn_alpha,vGn_n,theta_sat,theta_res,vGn_m,                  & ! intent(in)    : soil parameters
                    ! input
                    LH_fus*iden_water                                          ,& ! intent(in)    : energy for melt+freeze (J m-3)
                    mLayerVolHtCapBulk(iLayer)                                 ,& ! intent(in)    : volumetric heat capacity (J m-3 K-1)
                    mLayerTemp(iLayer)                                         ,& ! intent(in)    : initial temperature (K)
                    mLayerVolFracIce(iLayer)                                   ,& ! intent(in)    : initial volumetric fraction of ice (-)
                    mLayerMatricHeadTrial(ixControlIndex)                      ,& ! intent(in)    : mass state variable = trial matric head (m)
                    ! output 
                    mLayerTempTrial(iLayer)                                    ,& ! intent(inout) : trial temperature (K)
                    mLayerVolFracLiqTrial(iLayer)                              ,& ! intent(out)   : trial volumetric fraction of liquid water (-)
                    mLayerVolFracIceTrial(iLayer)                              ,& ! intent(out)   : trial volumetric fraction if ice (-)
                    volFracWatTrial=mLayerVolFracWatTrial(iLayer)              ,& ! intent(out)   : volumetric fraction of total water (-)
                    dLiq_dT=mLayerdTheta_dTk(iLayer)                           ,& ! intent(out)   : derivative in liquid water w.r.t. temperature (K-1) 
                    err=err,message=cmessage)                                     ! intent(out)   : error control
    if(err/=0)then; message=trim(message)//trim(cmessage); return; endif

    ! compute the liquid matric potential (and the derivatives w.r.t. total matric potential and temperature)
    call liquidHead(&
                    ! input
                    mLayerMatricHeadTrial(ixControlIndex)    ,& ! intent(in)    : total water matric potential (m)
                    mLayerVolFracLiqTrial(iLayer)            ,& ! intent(in)    : volumetric fraction of liquid water (-)
                    mLayerVolFracIceTrial(iLayer)            ,& ! intent(in)    : volumetric fraction of ice (-)
                    vGn_alpha,vGn_n,theta_sat,theta_res,vGn_m,& ! intent(in)    : soil parameters
                    dVolTot_dPsi0(ixControlIndex)            ,& ! intent(in)    : derivative in the soil water characteristic (m-1)
                    mLayerdTheta_dTk(iLayer)                 ,& ! intent(in)    : derivative in volumetric total water w.r.t. temperature (K-1)
                    ! output
                    mLayerMatricHeadLiqTrial(ixControlIndex) ,& ! intent(out)   : liquid water matric potential (m)
                    dPsiLiq_dPsi0(ixControlIndex)            ,& ! intent(out)   : derivative in the liquid water matric potential w.r.t. the total water matric potential (-)
                    dPsiLiq_dTemp(ixControlIndex)            ,& ! intent(out)   : derivative in the liquid water matric potential w.r.t. temperature (m K-1)
                    err,cmessage)                               ! intent(out)   : error control
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
 ! private subroutine updateVars: update prognostoc variables
 !  * compute volumetric fraction of liquid water and ice (and update temperature if desired)
 ! **********************************************************************************************************
 subroutine updateVars(&
                       ! model control
                       ixDomainType,                                  & ! intent(in)    : named variable defining the domain type (iname_veg, iname_snow, iname_soil)
                       do_adjustTemp,                                 & ! intent(in)    : logical flag to adjust temperature to account for energy associated with melt+freeze
                       snowfrz_scale,                                 & ! intent(in)    : scaling parameter for the snow freezing curve (K-1)
                       vGn_alpha,vGn_n,theta_sat,theta_res,vGn_m,     & ! intent(in)    : soil parameters
                       ! input
                       meltNrg                                       ,& ! intent(in)    : energy for melt+freeze (J m-3)
                       heatCap                                       ,& ! intent(in)    : volumetric heat capacity (J m-3 K-1)
                       tempInit                                      ,& ! intent(in)    : initial temperature (K)
                       volFracIceInit                                ,& ! intent(in)    : initial volumetric fraction of ice (-)
                       massStateTrial                                ,& ! intent(in)    : trial mass state variable
                       ! output
                       tempTrial                                     ,& ! intent(inout) : trial temperature (K)
                       volFracLiqTrial                               ,& ! intent(out)   : volumetric fraction of liquid water (-)
                       volFracIceTrial                               ,& ! intent(out)   : volumetric fraction of ice (-)
                       volFracWatTrial ,&  ! OPTIONAL                   ! intent(out)   : volumetric fraction of total water (-)
                       fracLiqTrial    ,&  ! OPTIONAL                   ! intent(out)   : fraction liquid water (-)
                       dLiq_dT                                       ,& ! intent(out)   : derivative in liquid water w.r.t. temperature (K-1)
                       ! error handling
                       err,message)                                     ! intent(out)   : error control
 ! --------------------------------------------------------------------------------------------------------------------------------
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
 real(dp),intent(out)           :: dLiq_dT          ! derivative in liquid water w.r.t. temperature (K-1)
 ! error handling
 integer(i4b),intent(out)       :: err              ! error code
 character(*),intent(out)       :: message          ! error message
 ! -----------------------------------------------------------------------------------------------------------------
 ! local variables
 integer(i4b)                   :: iter             ! iteration index
 integer(i4b),parameter         :: maxiter=20       ! iteration index
 real(dp)                       :: residual         ! residual in the energy equation (J m-3)
 real(dp)                       :: derivative       ! derivative in the energy equation (J m-3 K-1)
 real(dp)                       :: tempInc          ! iteration increment (K)
 real(dp)                       :: Tcrit            ! critical temperature above which all water is liquid (K)
 real(dp),parameter             :: nrgConvTol=1.e-4_dp ! convergence tolerance (J m-3)
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

  ! get the critical soil temperature where all water is frozen
  ! NOTE: the first argument is for soil, next argument for everything else
  Tcrit = merge(crit_soilT(massStateTrial), Tfreeze, ixDomainType==iname_soil) 

  ! get the derivative in total volume water w.r.t. temperature
  ! NOTE: derivative is zero if tempTrial>Tcrit
  select case(ixDomainType)
   case(iname_veg, iname_snow); dLiq_dT = merge(dFracLiq_dTk(tempTrial,snowfrz_scale)*massStateTrial, 0._dp, tempTrial<Tcrit)
   case(           iname_soil); dLiq_dT = merge(dTheta_dTk(tempTrial,theta_res,theta_sat,vGn_alpha,vGn_n,vGn_m), 0._dp, tempTrial<Tcrit) ! assume no volume expansion
   case default; err=20; message=trim(message)//'expect case to be iname_veg, iname_snow, iname_soil'; return
  end select  ! domain type

  ! -----
  ! - compute new temperature...
  ! ----------------------------

  ! check if the need to update temperatures
  if(.not.do_adjustTemp) return

  ! check failed convergence
  if(iter==maxiter)then
   message=trim(message)//'failed to converge'
   stop
   err=20; return
  endif

  ! compute the residual and check convergence criteria
  residual   = -heatCap*(tempTrial - tempInit) + meltNrg*(volFracIceTrial - volFracIceInit)  ! J m-3
  !if(iter>1) write(*,'(i4,1x,e20.10,1x,2(f20.10,1x))') iter, residual, tempTrial, tempInc
  if(abs(residual) < nrgConvTol) exit

  ! compute the derivative and the iteration increment
  derivative = heatCap + LH_fus*iden_water*dLiq_dT  ! J m-3 K-1
  tempInc    = residual/derivative        ! K

  ! add constraints for snow temperature
  if(ixDomainType==iname_snow)then
   if(tempInc > Tcrit - tempTrial) tempInc=(Tcrit - tempTrial)*0.5_dp  ! simple bi-section method
  endif

  ! update the temperature trial
  tempTrial = tempTrial + tempInc
  
 end do  ! iterating
 !print*, 'PAUSE: after iterations in '//trim(message); read(*,*)

 end subroutine updateVars

end module getVectorz_module

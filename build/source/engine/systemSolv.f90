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

module systemSolv_module

! data types
USE nrtype

! layer types
USE globalData,only:ix_soil,ix_snow ! named variables for snow and soil

! access the global print flag
USE globalData,only:globalPrintFlag

! access missing values
USE globalData,only:integerMissing  ! missing integer
USE globalData,only:realMissing     ! missing double precision number
USE globalData,only:quadMissing     ! missing quadruple precision number

! access matrix information
USE globalData,only: nBands         ! length of the leading dimension of the band diagonal matrix
USE globalData,only: ixFullMatrix   ! named variable for the full Jacobian matrix
USE globalData,only: ixBandMatrix   ! named variable for the band diagonal matrix

! named variables to describe the state variable type
USE globalData,only:iname_nrgCanair ! named variable defining the energy of the canopy air space
USE globalData,only:iname_nrgCanopy ! named variable defining the energy of the vegetation canopy
USE globalData,only:iname_watCanopy ! named variable defining the mass of water on the vegetation canopy
USE globalData,only:iname_nrgLayer  ! named variable defining the energy state variable for snow+soil layers
USE globalData,only:iname_watLayer  ! named variable defining the total water state variable for snow+soil layers
USE globalData,only:iname_matLayer  ! named variable defining the matric head state variable for soil layers

! constants
USE multiconst,only:&
                    gravity,      & ! acceleration of gravity              (m s-2)
                    Tfreeze,      & ! temperature at freezing              (K)
                    LH_fus,       & ! latent heat of fusion                (J kg-1)
                    LH_vap,       & ! latent heat of vaporization          (J kg-1)
                    LH_sub,       & ! latent heat of sublimation           (J kg-1)
                    Cp_air,       & ! specific heat of air                 (J kg-1 K-1)
                    iden_air,     & ! intrinsic density of air             (kg m-3)
                    iden_ice,     & ! intrinsic density of ice             (kg m-3)
                    iden_water      ! intrinsic density of liquid water    (kg m-3)

! provide access to indices that define elements of the data structures
USE var_lookup,only:iLookATTR       ! named variables for structure elements
USE var_lookup,only:iLookTYPE       ! named variables for structure elements
USE var_lookup,only:iLookPROG       ! named variables for structure elements
USE var_lookup,only:iLookDIAG       ! named variables for structure elements
USE var_lookup,only:iLookFLUX       ! named variables for structure elements
USE var_lookup,only:iLookFORCE      ! named variables for structure elements
USE var_lookup,only:iLookPARAM      ! named variables for structure elements
USE var_lookup,only:iLookINDEX      ! named variables for structure elements
USE var_lookup,only:iLookDECISIONS  ! named variables for elements of the decision structure

! provide access to the derived types to define the data structures
USE data_types,only:&
                    var_i,        & ! data vector (i4b)
                    var_d,        & ! data vector (dp)
                    var_ilength,  & ! data vector with variable length dimension (i4b)
                    var_dlength,  & ! data vector with variable length dimension (dp)
                    model_options   ! defines the model decisions

! look-up values for the choice of groundwater representation (local-column, or single-basin)
USE mDecisions_module,only:       &
 localColumn,                     & ! separate groundwater representation in each local soil column
 singleBasin                        ! single groundwater store over the entire basin

! look-up values for the choice of groundwater parameterization
USE mDecisions_module,only:      &
 qbaseTopmodel,                  & ! TOPMODEL-ish baseflow parameterization
 bigBucket,                      & ! a big bucket (lumped aquifer model)
 noExplicit                        ! no explicit groundwater parameterization

! safety: set private unless specified otherwise
implicit none
private
public::systemSolv

! control parameters
real(dp),parameter  :: valueMissing=-9999._dp     ! missing value
real(dp),parameter  :: verySmall=tiny(1.0_dp)     ! a very small number
real(dp),parameter  :: veryBig=1.e+20_dp          ! a very big number
real(dp),parameter  :: dx = 1.e-8_dp             ! finite difference increment

contains


 ! **********************************************************************************************************
 ! public subroutine systemSolv: run the coupled energy-mass model for one timestep
 ! **********************************************************************************************************
 subroutine systemSolv(&
                       ! input: model control
                       nSnow,          & ! intent(in): number of snow layers
                       nSoil,          & ! intent(in): number of soil layers
                       nLayers,        & ! intent(in): total number of layers
                       nState,         & ! intent(in): total number of state variables
                       dt,             & ! intent(in): time step (s)
                       maxiter,        & ! intent(in): maximum number of iterations
                       firstSubstep,   & ! intent(in): flag to denote first sub-step
                       computeVegFlux, & ! intent(in): flag to denote if computing energy flux over vegetation
                       ! input/output: data structures
                       type_data,      & ! intent(in):    type of vegetation and soil
                       attr_data,      & ! intent(in):    spatial attributes
                       forc_data,      & ! intent(in):    model forcing data
                       mpar_data,      & ! intent(in):    model parameters
                       indx_data,      & ! intent(inout): index data
                       prog_data,      & ! intent(inout): model prognostic variables for a local HRU
                       diag_data,      & ! intent(inout): model diagnostic variables for a local HRU
                       flux_data,      & ! intent(inout): model fluxes for a local HRU
                       bvar_data,      & ! intent(in):    model variables for the local basin
                       model_decisions,& ! intent(in):    model decisions
                       ! output: model control
                       niter,          & ! number of iterations taken
                       err,message)      ! error code and error message
 ! ---------------------------------------------------------------------------------------
 ! structure allocations
 USE globalData,only:deriv_meta                       ! metadata on the model derivatives
 USE allocspace_module,only:allocLocal                ! allocate local data structures
 ! simulation of fluxes and residuals given a trial state vector
 USE eval8summa_module,only:eval8summa                ! simulation of fluxes and residuals given a trial state vector
 USE summaSolve_module,only:summaSolve                ! calculate the iteration increment, evaluate the new state, and refine if necessary
 ! population/extracction of state vectors
 USE getVectorz_module,only:popStateVec               ! populate the state vector
 USE getVectorz_module,only:varExtract                ! extract variables from the state vector
 implicit none
 ! ---------------------------------------------------------------------------------------
 ! * dummy variables
 ! ---------------------------------------------------------------------------------------
 ! input: model control
 integer(i4b),intent(in)         :: nSnow                         ! number of snow layers
 integer(i4b),intent(in)         :: nSoil                         ! number of soil layers
 integer(i4b),intent(in)         :: nLayers                       ! total number of layers
 integer(i4b),intent(in)         :: nState                        ! total number of state variables
 real(dp),intent(in)             :: dt                            ! time step (seconds)
 integer(i4b),intent(in)         :: maxiter                       ! maximum number of iterations
 logical(lgt),intent(in)         :: firstSubStep                  ! flag to indicate if we are processing the first sub-step
 logical(lgt),intent(in)         :: computeVegFlux                ! flag to indicate if we are computing fluxes over vegetation (.false. means veg is buried with snow)
 ! input/output: data structures
 type(var_i),intent(in)          :: type_data                     ! type of vegetation and soil
 type(var_d),intent(in)          :: attr_data                     ! spatial attributes
 type(var_d),intent(in)          :: forc_data                     ! model forcing data
 type(var_d),intent(in)          :: mpar_data                     ! model parameters
 type(var_ilength),intent(inout) :: indx_data                     ! indices for a local HRU
 type(var_dlength),intent(inout) :: prog_data                     ! prognostic variables for a local HRU
 type(var_dlength),intent(inout) :: diag_data                     ! diagnostic variables for a local HRU
 type(var_dlength),intent(inout) :: flux_data                     ! model fluxes for a local HRU
 type(var_dlength),intent(in)    :: bvar_data                     ! model variables for the local basin
 type(model_options),intent(in)  :: model_decisions(:)            ! model decisions
 ! output: model control
 integer(i4b),intent(out)        :: niter                         ! number of iterations
 integer(i4b),intent(out)        :: err                           ! error code
 character(*),intent(out)        :: message                       ! error message
 ! *********************************************************************************************************************************************************
 ! *********************************************************************************************************************************************************
 ! ---------------------------------------------------------------------------------------
 ! * general local variables
 ! ---------------------------------------------------------------------------------------
 character(LEN=256)              :: cmessage                     ! error message of downwind routine
 integer(i4b)                    :: iter                         ! iteration index
 integer(i4b)                    :: nLeadDim                     ! length of the leading dimension of the Jacobian matrix (nBands or nState)
 integer(i4b)                    :: local_ixGroundwater          ! local index for groundwater representation
 real(dp),parameter              :: tempAccelerate=0.00_dp       ! factor to force initial canopy temperatures to be close to air temperature
 real(dp),parameter              :: xMinCanopyWater=0.0001_dp    ! minimum value to initialize canopy water (kg m-2)
 ! ------------------------------------------------------------------------------------------------------
 ! * state variables and diagniostic variables
 ! ------------------------------------------------------------------------------------------------------
 ! state variables
 real(dp)                        :: scalarCanairTempTrial ! trial value for temperature of the canopy air space (K)
 real(dp)                        :: scalarCanopyTempTrial ! trial value for temperature of the vegetation canopy (K)
 real(dp)                        :: scalarCanopyWatTrial  ! trial value for liquid water storage in the canopy (kg m-2)
 real(dp),dimension(nLayers)     :: mLayerTempTrial       ! trial value for temperature of layers in the snow and soil domains (K)
 real(dp),dimension(nLayers)     :: mLayerVolFracWatTrial ! trial value for volumetric fraction of total water (-)
 real(dp),dimension(nSoil)       :: mLayerMatricHeadTrial ! trial value for matric head (m)
 ! diagnostic variables
 real(dp)                        :: scalarCanopyLiqTrial  ! trial value for mass of liquid water on the vegetation canopy (kg m-2)
 real(dp)                        :: scalarCanopyIceTrial  ! trial value for mass of ice on the vegetation canopy (kg m-2)
 real(dp),dimension(nLayers)     :: mLayerVolFracLiqTrial ! trial value for volumetric fraction of liquid water (-)
 real(dp),dimension(nLayers)     :: mLayerVolFracIceTrial ! trial value for volumetric fraction of ice (-)
 ! ------------------------------------------------------------------------------------------------------
 ! * operator splitting
 ! ------------------------------------------------------------------------------------------------------
 integer(i4b),parameter          :: fullyImplicit=1001           ! named variable for the fully implicit solution
 integer(i4b),parameter          :: strangSplitting=1002         ! named variable for the Strang splitting solution
 integer(i4b)                    :: ixSplitOption=strangSplitting  ! selected operator splitting method
 integer(i4b)                    :: nOperSplit                   ! number of splitting operations
 integer(i4b)                    :: iSplit                       ! index of splitting operation
 integer(i4b),parameter          :: halfMass1=1                  ! order in sequence for the 1st mass operation
 integer(i4b),parameter          :: halfMass2=3                  ! order in sequence for the 2nd mass operation
 integer(i4b),parameter          :: fullNrg=2                    ! order in sequence for the energy operation
 logical(lgt),dimension(nState)  :: stateMask                    ! mask defining desired state variables
 integer(i4b)                    :: nSubset                      ! number of selected state variables for a given split
 ! ------------------------------------------------------------------------------------------------------
 ! * model solver
 ! ------------------------------------------------------------------------------------------------------
 logical(lgt),parameter          :: numericalJacobian=.false.    ! flag to compute the Jacobian matrix
 logical(lgt),parameter          :: testBandDiagonal=.false.     ! flag to test the band-diagonal matrix
 logical(lgt),parameter          :: forceFullMatrix=.false.      ! flag to force the use of the full Jacobian matrix
 logical(lgt)                    :: firstFluxCall                ! flag to define the first flux call
 integer(i4b)                    :: ixMatrix                     ! form of matrix (band diagonal or full matrix)
 type(var_dlength)               :: deriv_data                   ! derivatives in model fluxes w.r.t. relevant state variables 
 integer(i4b)                    :: ixSaturation                 ! index of the lowest saturated layer (NOTE: only computed on the first iteration)
 real(dp),allocatable            :: dBaseflow_dMatric(:,:)       ! derivative in baseflow w.r.t. matric head (s-1)  ! NOTE: allocatable, since not always needed
 real(dp),allocatable            :: stateVecInit(:)              ! initial state vector (mixed units)
 real(dp),allocatable            :: stateVecTrial(:)             ! trial state vector (mixed units)
 real(dp),allocatable            :: stateVecNew(:)               ! new state vector (mixed units)
 real(dp),allocatable            :: fluxVec0(:)                  ! flux vector (mixed units)
 real(dp),allocatable            :: fScale(:)                    ! characteristic scale of the function evaluations (mixed units)
 real(dp),allocatable            :: xScale(:)                    ! characteristic scale of the state vector (mixed units)
 real(dp),allocatable            :: dMat(:)                      ! diagonal matrix (excludes flux derivatives)
 real(qp),allocatable            :: sMul(:)       ! NOTE: qp     ! multiplier for state vector for the residual calculations
 real(qp),allocatable            :: rVec(:)       ! NOTE: qp     ! residual vector
 real(dp),allocatable            :: rAdd(:)                      ! additional terms in the residual vector
 real(dp)                        :: fOld,fNew                    ! function values (-); NOTE: dimensionless because scaled
 logical(lgt)                    :: feasible                     ! flag to define the feasibility of the solution
 logical(lgt)                    :: converged                    ! convergence flag
 real(dp),allocatable            :: resSinkNew(:)                ! additional terms in the residual vector
 real(dp),allocatable            :: fluxVecNew(:)                ! new flux vector
 real(qp),allocatable            :: resVecNew(:)  ! NOTE: qp     ! new residual vector
 ! ------------------------------------------------------------------------------------------------------
 ! * mass balance checks
 ! ------------------------------------------------------------------------------------------------------
 logical(lgt),parameter          :: checkMassBalance=.true.      ! flag to check the mass balance
 real(dp)                        :: balance0,balance1            ! storage at start and end of time step
 real(dp)                        :: vertFlux                     ! change in storage due to vertical fluxes
 real(dp)                        :: tranSink,baseSink,compSink   ! change in storage sue to sink terms
 real(dp)                        :: liqError                     ! water balance error
 ! ---------------------------------------------------------------------------------------
 ! point to variables in the data structures
 ! ---------------------------------------------------------------------------------------
 associate(&
 ! model decisions
 ixGroundwater           => model_decisions(iLookDECISIONS%groundwatr)%iDecision   ,& ! intent(in): [i4b] groundwater parameterization
 ixSpatialGroundwater    => model_decisions(iLookDECISIONS%spatial_gw)%iDecision   ,& ! intent(in): [i4b] spatial representation of groundwater (local-column or single-basin)
 ! domain boundary conditions
 airtemp                 => forc_data%var(iLookFORCE%airtemp)                      ,& ! intent(in): [dp] temperature of the upper boundary of the snow and soil domains (K)
 ! indices of model state variables
 ixCasNrg                => indx_data%var(iLookINDEX%ixCasNrg)%dat(1)             , & ! intent(in): [i4b]    index of canopy air space energy state variable
 ixVegNrg                => indx_data%var(iLookINDEX%ixVegNrg)%dat(1)             , & ! intent(in): [i4b]    index of canopy energy state variable
 ixVegWat                => indx_data%var(iLookINDEX%ixVegWat)%dat(1)             , & ! intent(in): [i4b]    index of canopy hydrology state variable (mass)
 ixSnowSoilNrg           => indx_data%var(iLookINDEX%ixSnowSoilNrg)%dat           , & ! intent(in): [i4b(:)] indices for energy states in the snow-soil subdomain
 ixSnowOnlyWat           => indx_data%var(iLookINDEX%ixSnowOnlyWat)%dat           , & ! intent(in): [i4b(:)] indices for total water states in the snow subdomain
 ixSoilOnlyHyd           => indx_data%var(iLookINDEX%ixSoilOnlyHyd)%dat           , & ! intent(in): [i4b(:)] indices for hydrology states in the soil subdomain
 ixStateType             => indx_data%var(iLookINDEX%ixStateType)%dat             , & ! intent(in): [i4b(:)] indices defining the type of the state (ixNrgState...)
 ! vegetation parameters
 canopyDepth             => diag_data%var(iLookDIAG%scalarCanopyDepth)%dat(1)      ,& ! intent(in): [dp] canopy depth (m)
 ! snow parameters
 snowfrz_scale           => mpar_data%var(iLookPARAM%snowfrz_scale)                ,& ! intent(in): [dp] scaling parameter for the snow freezing curve (K-1)
 ! soil parameters
 vGn_m                   => diag_data%var(iLookDIAG%scalarVGn_m)%dat(1)            ,& ! intent(in): [dp] van Genutchen "m" parameter (-)
 vGn_n                   => mpar_data%var(iLookPARAM%vGn_n)                        ,& ! intent(in): [dp] van Genutchen "n" parameter (-)
 vGn_alpha               => mpar_data%var(iLookPARAM%vGn_alpha)                    ,& ! intent(in): [dp] van Genutchen "alpha" parameter (m-1)
 theta_sat               => mpar_data%var(iLookPARAM%theta_sat)                    ,& ! intent(in): [dp] soil porosity (-)
 theta_res               => mpar_data%var(iLookPARAM%theta_res)                    ,& ! intent(in): [dp] soil residual volumetric water content (-)
 specificStorage         => mpar_data%var(iLookPARAM%specificStorage)              ,& ! intent(in): [dp] specific storage coefficient (m-1)
 ! convergence parameters
 absConvTol_liquid       => mpar_data%var(iLookPARAM%absConvTol_liquid)            ,& ! intent(in): [dp] absolute convergence tolerance for vol frac liq water (-)
 ! model diagnostic variables (fraction of liquid water)
 scalarFracLiqVeg        => diag_data%var(iLookDIAG%scalarFracLiqVeg)%dat(1)       ,& ! intent(out): [dp]    fraction of liquid water on vegetation (-)
 mLayerFracLiqSnow       => diag_data%var(iLookDIAG%mLayerFracLiqSnow)%dat         ,& ! intent(out): [dp(:)] fraction of liquid water in each snow layer (-)
 mLayerMeltFreeze        => diag_data%var(iLookDIAG%mLayerMeltFreeze)%dat          ,& ! intent(out): [dp(:)] melt-freeze in each snow and soil layer (kg m-3)
 ! model fluxes
 mLayerCompress          => diag_data%var(iLookDIAG%mLayerCompress)%dat            ,& ! intent(out): [dp(:)]  change in storage associated with compression of the soil matrix (-)
 iLayerLiqFluxSoil       => flux_data%var(iLookFLUX%iLayerLiqFluxSoil)%dat         ,& ! intent(out): [dp(0:)] vertical liquid water flux at soil layer interfaces (-)
 mLayerTranspire         => flux_data%var(iLookFLUX%mLayerTranspire)%dat           ,& ! intent(out): [dp(:)]  transpiration loss from each soil layer (m s-1)
 mLayerBaseflow          => flux_data%var(iLookFLUX%mLayerBaseflow)%dat            ,& ! intent(out): [dp(:)]  baseflow from each soil layer (m s-1)
 scalarExfiltration      => flux_data%var(iLookFLUX%scalarExfiltration)%dat(1)     ,& ! intent(out): [dp]     exfiltration from the soil profile (m s-1)
 ! sublimation (needed to check mass balance constraints)
 scalarCanopySublimation => flux_data%var(iLookFLUX%scalarCanopySublimation)%dat(1),& ! intent(out): [dp] sublimation of ice from the vegetation canopy (kg m-2 s-1)
 scalarSnowSublimation   => flux_data%var(iLookFLUX%scalarSnowSublimation)%dat(1)  ,& ! intent(out): [dp] sublimation of ice from the snow surface (kg m-2 s-1)
 ! model state variables (vegetation canopy)
 mLayerDepth             => prog_data%var(iLookPROG%mLayerDepth)%dat               ,& ! intent(in):    [dp(:)] depth of each layer in the snow-soil sub-domain (m)
 scalarCanairTemp        => prog_data%var(iLookPROG%scalarCanairTemp)%dat(1)       ,& ! intent(inout): [dp] temperature of the canopy air space (K)
 scalarCanopyTemp        => prog_data%var(iLookPROG%scalarCanopyTemp)%dat(1)       ,& ! intent(inout): [dp] temperature of the vegetation canopy (K)
 scalarCanopyIce         => prog_data%var(iLookPROG%scalarCanopyIce)%dat(1)        ,& ! intent(inout): [dp] mass of ice on the vegetation canopy (kg m-2)
 scalarCanopyLiq         => prog_data%var(iLookPROG%scalarCanopyLiq)%dat(1)        ,& ! intent(inout): [dp] mass of liquid water on the vegetation canopy (kg m-2)
 scalarCanopyWat         => prog_data%var(iLookPROG%scalarCanopyWat)%dat(1)        ,& ! intent(inout): [dp] mass of total water on the vegetation canopy (kg m-2)
 ! model state variables (snow and soil domains)
 mLayerTemp              => prog_data%var(iLookPROG%mLayerTemp)%dat                ,& ! intent(inout): [dp(:)] temperature of each snow/soil layer (K)
 mLayerVolFracIce        => prog_data%var(iLookPROG%mLayerVolFracIce)%dat          ,& ! intent(inout): [dp(:)] volumetric fraction of ice (-)
 mLayerVolFracLiq        => prog_data%var(iLookPROG%mLayerVolFracLiq)%dat          ,& ! intent(inout): [dp(:)] volumetric fraction of liquid water (-)
 mLayerVolFracWat        => prog_data%var(iLookPROG%mLayerVolFracWat)%dat          ,& ! intent(inout): [dp(:)] volumetric fraction of total water (-)
 mLayerMatricHead        => prog_data%var(iLookPROG%mLayerMatricHead)%dat          ,& ! intent(inout): [dp(:)] matric head (m)
 scalarAquiferStorage    => prog_data%var(iLookPROG%scalarAquiferStorage)%dat(1)    & ! intent(inout): [dp   ] aquifer storage (m)
 )
 ! ---------------------------------------------------------------------------------------
 ! initialize error control
 err=0; message="systemSolv/"

 ! *****
 ! (0) PRELIMINARIES...
 ! ********************

 ! -----
 ! * initialize...
 ! ---------------

 ! initialize the first flux call
 firstFluxCall=.true.

 ! compute the total water content in the vegetation canopy
 scalarCanopyWat = scalarCanopyLiq + scalarCanopyIce  ! kg m-2

 ! compute the total water content in snow and soil
 ! NOTE: no ice expansion allowed for soil
 if(nSnow>0)& 
 mLayerVolFracWat(1:      nSnow)   = mLayerVolFracLiq(1:      nSnow)   + mLayerVolFracIce(1:      nSnow)*(iden_ice/iden_water)
 mLayerVolFracWat(nSnow+1:nLayers) = mLayerVolFracLiq(nSnow+1:nLayers) + mLayerVolFracIce(nSnow+1:nLayers)

 ! define the number of operator splits
 select case(ixSplitOption)
  case(fullyImplicit);   nOperSplit=1
  case(strangSplitting); nOperSplit=3
  case default; err=20; message=trim(message)//'unable to identify operator splitting strategy'; return
 end select

 ! identify the matrix solution method
 ! (the type of matrix used to solve the linear system A.X=B)
 if(ixGroundwater==qbaseTopmodel .or. testBandDiagonal .or. forceFullMatrix)then
  nLeadDim=nState         ! length of the leading dimension
  ixMatrix=ixFullMatrix   ! named variable to denote the full Jacobian matrix
 else
  nLeadDim=nBands         ! length of the leading dimension
  ixMatrix=ixBandMatrix   ! named variable to denote the band-diagonal matrix
 endif

 ! modify the groundwater representation for this single-column implementation
 select case(ixSpatialGroundwater)
  case(singleBasin); local_ixGroundwater = noExplicit    ! force no explicit representation of groundwater at the local scale
  case(localColumn); local_ixGroundwater = ixGroundwater ! go with the specified decision
  case default; err=20; message=trim(message)//'unable to identify spatial representation of groundwater'; return
 end select ! (modify the groundwater representation for this single-column implementation)

 ! allocate space for the derivative structure
 call allocLocal(deriv_meta(:),deriv_data,nSnow,nSoil,err,cmessage)
 if(err/=0)then; err=20; message=trim(message)//trim(cmessage); return; endif

 ! allocate space for the baseflow derivatives
 if(ixGroundwater==qbaseTopmodel)then
  allocate(dBaseflow_dMatric(nSoil,nSoil),stat=err)  ! baseflow depends on total storage in the soil column, hence on matric head in every soil layer
 else
  allocate(dBaseflow_dMatric(0,0),stat=err)          ! allocate zero-length dimnensions to avoid passing around an unallocated matrix
 endif
 if(err/=0)then; err=20; message=trim(message)//'unable to allocate space for the baseflow derivatives'; return; endif

 ! ==========================================================================================================================================
 ! ==========================================================================================================================================
 ! ==========================================================================================================================================
 ! ==========================================================================================================================================

 ! operator splitting loop
 do iSplit=1,nOperSplit 

  ! define mask for the strang splitting
  if(ixSplitOption==strangSplitting)then
   select case(iSplit)
    case(fullNrg);             stateMask = (ixStateType==iname_nrgCanair .or. ixStateType==iname_nrgCanopy .or. ixStateType==iname_nrgLayer)
    case(halfMass1,halfMass2); stateMask = (ixStateType==iname_watCanopy .or. ixStateType==iname_watLayer .or. ixStateType==iname_matLayer)
    case default; err=20; message=trim(message)//'unable to identify splitting operation'; return
   end select
  else  ! (not strang splitting)
   stateMask(:) = .true.  ! use all state variables
  endif

  ! get the number of selected state variables
  nSubset = count(stateMask)

  ! allocate space for solution and scaling vectors
  allocate(stateVecInit(nSubset), stateVecTrial(nSubset), stateVecNew(nSubset), fluxVec0(nSubset), fluxVecNew(nSubset), fScale(nSubset), xScale(nSubset), stat=err)
  if(err/=0)then; err=20; message=trim(message)//'unable to allocate space for the solution and scaling vectors'; return; endif

  ! allocate space for the diagonal matrix, multipliers and residual vectors
  allocate(dMat(nSubset), sMul(nSubset), rVec(nSubset), rAdd(nSubset), resSinkNew(nSubset), resVecNew(nSubset), stat=err)
  if(err/=0)then; err=20; message=trim(message)//'unable to allocate space for the diagonal matrix, multipliers, and residual vectors'; return; endif

  ! initialize state vectors
  call popStateVec(&
                   ! input
                   computeVegFlux,              & ! intent(in):    flag to denote if computing energy flux over vegetation
                   pack(ixStateType,stateMask), & ! intent(in):    type of desired model state variables
                   nSubset,                     & ! intent(in):    number of desired state variables
                   ! input-output: data structures
                   prog_data,                   & ! intent(in):    model prognostic variables for a local HRU
                   diag_data,                   & ! intent(in):    model diagnostic variables for a local HRU
                   indx_data,                   & ! intent(in):    indices defining model states and layers
                   ! output
                   stateVecInit,                & ! intent(out):   initial model state vector (mixed units)
                   fScale,                      & ! intent(out):   function scaling vector (mixed units)
                   xScale,                      & ! intent(out):   variable scaling vector (mixed units)
                   sMul,                        & ! intent(out):   multiplier for state vector (used in the residual calculations)
                   dMat,                        & ! intent(out):   diagonal of the Jacobian matrix (excludes fluxes) 
                   err,cmessage)                  ! intent(out):   error control
  if(err/=0)then; message=trim(message)//trim(cmessage); return; endif  ! (check for errors)
  
  ! -----
  ! * compute the initial function evaluation...
  ! --------------------------------------------
  
  ! initialize the trial state vectors
  stateVecTrial = stateVecInit
  
  ! need to intialize canopy water at a positive value
  if(computeVegFlux)then
   if(scalarCanopyWat < xMinCanopyWater) stateVecTrial(ixVegWat) = scalarCanopyWat + xMinCanopyWater
  endif
  
  ! try to accelerate solution for energy
  if(computeVegFlux)then
   stateVecTrial(ixCasNrg) = stateVecInit(ixCasNrg) + (airtemp - stateVecInit(ixCasNrg))*tempAccelerate
   stateVecTrial(ixVegNrg) = stateVecInit(ixVegNrg) + (airtemp - stateVecInit(ixVegNrg))*tempAccelerate
  endif
  
  ! compute the flux and the residual vector for a given state vector
  ! NOTE 1: The derivatives computed in eval8summa are used to calculate the Jacobian matrix for the first iteration
  ! NOTE 2: The Jacobian matrix together with the residual vector is used to calculate the first iteration increment
  call eval8summa(&
                  ! input: model control
                  dt,                      & ! intent(in):    length of the time step (seconds)
                  nSnow,                   & ! intent(in):    number of snow layers
                  nSoil,                   & ! intent(in):    number of soil layers
                  nLayers,                 & ! intent(in):    total number of layers
                  nState,                  & ! intent(in):    total number of state variables
                  firstSubStep,            & ! intent(in):    flag to indicate if we are processing the first sub-step
                  firstFluxCall,           & ! intent(inout): flag to indicate if we are processing the first flux call
                  computeVegFlux,          & ! intent(in):    flag to indicate if we need to compute fluxes over vegetation
                  ! input: state vectors
                  stateVecTrial,           & ! intent(in):    model state vector
                  fScale,                  & ! intent(in):    function scaling vector
                  sMul,                    & ! intent(in):    state vector multiplier (used in the residual calculations)
                  ! input: data structures
                  model_decisions,         & ! intent(in):    model decisions
                  type_data,               & ! intent(in):    type of vegetation and soil
                  attr_data,               & ! intent(in):    spatial attributes
                  mpar_data,               & ! intent(in):    model parameters
                  forc_data,               & ! intent(in):    model forcing data
                  bvar_data,               & ! intent(in):    average model variables for the entire basin
                  prog_data,               & ! intent(in):    model prognostic variables for a local HRU
                  indx_data,               & ! intent(in):    index data
                  ! input-output: data structures
                  diag_data,               & ! intent(inout): model diagnostic variables for a local HRU
                  flux_data,               & ! intent(inout): model fluxes for a local HRU
                  deriv_data,              & ! intent(inout): derivatives in model fluxes w.r.t. relevant state variables
                  ! input-output: baseflow
                  ixSaturation,            & ! intent(inout): index of the lowest saturated layer (NOTE: only computed on the first iteration)
                  dBaseflow_dMatric,       & ! intent(out):   derivative in baseflow w.r.t. matric head (s-1)
                  ! output
                  feasible,                & ! intent(out):   flag to denote the feasibility of the solution
                  fluxVec0,                & ! intent(out):   flux vector
                  rAdd,                    & ! intent(out):   additional (sink) terms on the RHS of the state equation
                  rVec,                    & ! intent(out):   residual vector
                  fOld,                    & ! intent(out):   function evaluation
                  err,cmessage)              ! intent(out):   error control
  if(err/=0)then; message=trim(message)//trim(cmessage); return; endif  ! (check for errors)
  
  ! check feasibility (state vector SHOULD be feasible at this point)
  if(.not.feasible)then
   message=trim(message)//'unfeasible state vector'
   err=20; return
  endif
  
  ! ==========================================================================================================================================
  ! ==========================================================================================================================================
  ! ==========================================================================================================================================
  ! ==========================================================================================================================================
  
  ! (1) MAIN ITERATION LOOP...
  ! **************************
  
  ! iterate
  do iter=1,maxiter
  
   ! keep track of the number of iterations
   niter = iter+1  ! +1 because xFluxResid was moved outside the iteration loop (for backwards compatibility)
  
   ! compute the next trial state vector
   !  1) Computes the Jacobian matrix based on derivatives from the last flux evaluation
   !  2) Computes the iteration increment based on Jacobian and residuals from the last flux evaluation
   !  3) Computes new fluxes and derivatives, new residuals, and (if necessary) refines the state vector
   call summaSolve(&
                   ! input: model control
                   dt,                      & ! intent(in):    length of the time step (seconds)
                   iter,                    & ! intent(in):    iteration index
                   nSnow,                   & ! intent(in):    number of snow layers
                   nSoil,                   & ! intent(in):    number of soil layers
                   nLayers,                 & ! intent(in):    total number of layers
                   nLeadDim,                & ! intent(in):    length of the leading dimension of he Jacobian matrix (either nBands or nState) 
                   nState,                  & ! intent(in):    total number of state variables
                   ixMatrix,                & ! intent(in):    type of matrix (full or band diagonal)
                   firstSubStep,            & ! intent(in):    flag to indicate if we are processing the first sub-step
                   firstFluxCall,           & ! intent(inout): flag to indicate if we are processing the first flux call
                   computeVegFlux,          & ! intent(in):    flag to indicate if we need to compute fluxes over vegetation
                   ! input: state vectors
                   stateVecTrial,           & ! intent(in):    trial state vector
                   fScale,                  & ! intent(in):    function scaling vector
                   xScale,                  & ! intent(in):    "variable" scaling vector, i.e., for state variables
                   rVec,                    & ! intent(in):    residual vector
                   sMul,                    & ! intent(in):    state vector multiplier (used in the residual calculations)
                   dMat,                    & ! intent(inout): diagonal matrix (excludes flux derivatives)
                   fOld,                    & ! intent(in):    old function evaluation
                   ! input: data structures
                   model_decisions,         & ! intent(in):    model decisions
                   type_data,               & ! intent(in):    type of vegetation and soil
                   attr_data,               & ! intent(in):    spatial attributes
                   mpar_data,               & ! intent(in):    model parameters
                   forc_data,               & ! intent(in):    model forcing data
                   bvar_data,               & ! intent(in):    average model variables for the entire basin
                   prog_data,               & ! intent(in):    model prognostic variables for a local HRU
                   indx_data,               & ! intent(in):    index data
                   ! input-output: data structures
                   diag_data,               & ! intent(inout): model diagnostic variables for a local HRU
                   flux_data,               & ! intent(inout): model fluxes for a local HRU
                   deriv_data,              & ! intent(inout): derivatives in model fluxes w.r.t. relevant state variables
                   ! input-output: baseflow
                   ixSaturation,            & ! intent(inout): index of the lowest saturated layer (NOTE: only computed on the first iteration)
                   dBaseflow_dMatric,       & ! intent(inout): derivative in baseflow w.r.t. matric head (s-1)
                   ! output
                   stateVecNew,             & ! intent(out):   new state vector
                   fluxVecNew,              & ! intent(out):   new flux vector
                   resSinkNew,              & ! intent(out):   additional (sink) terms on the RHS of the state equation
                   resVecNew,               & ! intent(out):   new residual vector
                   fNew,                    & ! intent(out):   new function evaluation
                   converged,               & ! intent(out):   convergence flag
                   err,cmessage)              ! intent(out):   error control
   if(err/=0)then; message=trim(message)//trim(cmessage); return; endif  ! (check for errors)
  
   ! update function evaluation, residual vector, and states
   ! NOTE 1: The derivatives computed in summaSolve are used to calculate the Jacobian matrix at the next iteration
   ! NOTE 2: The Jacobian matrix together with the residual vector is used to calculate the new iteration increment
   fOld          = fNew
   rVec          = resVecNew 
   stateVecTrial = stateVecNew
  
   ! exit iteration loop if converged
   if(converged) exit
  
   ! check convergence
   if(niter==maxiter)then; err=-20; message=trim(message)//'failed to converge'; return; endif
   !print*, 'PAUSE: iterating'; read(*,*)
  
  end do  ! iterating
  !print*, 'PAUSE: after iterations'; read(*,*)

  ! deallocate space for temporary vectors
  deallocate(stateVecInit, stateVecTrial, stateVecNew, fluxVec0, fluxVecNew, fScale, xScale, dMat, sMul, rVec, rAdd, resSinkNew, resVecNew, stat=err)
  if(err/=0)then; err=20; message=trim(message)//'unable to deallocate space for temporary vectors'; return; endif
 
 end do  ! operator splitting loop 

 ! -----
 ! * update states and compute total volumetric melt...
 ! ----------------------------------------------------

 ! update temperatures (ensure new temperature is consistent with the fluxes)
 stateVecTrial(ixSnowSoilNrg) = stateVecInit(ixSnowSoilNrg) + (fluxVecNew(ixSnowSoilNrg)*dt + resSinkNew(ixSnowSoilNrg))/real(sMul(ixSnowSoilNrg), dp)

 ! update volumetric water content in the snow (ensure change in state is consistent with the fluxes)
 ! NOTE: for soil water balance is constrained within the iteration loop
 if(nSnow>0)&
 stateVecTrial(ixSnowOnlyWat) = stateVecInit(ixSnowOnlyWat) + (fluxVecNew(ixSnowOnlyWat)*dt + resSinkNew(ixSnowOnlyWat))

 ! update states: compute liquid water and ice content from total water content
 call varExtract(&
                 ! input
                 stateVecTrial,                             & ! intent(in):    model state vector (mixed units)
                 indx_data,                                 & ! intent(in):    indices defining model states and layers
                 snowfrz_scale,                             & ! intent(in):    scaling parameter for the snow freezing curve (K-1)
                 vGn_alpha,vGn_n,theta_sat,theta_res,vGn_m, & ! intent(in):    van Genutchen soil parameters
                 ! output: variables for the vegetation canopy
                 scalarFracLiqVeg,                          & ! intent(out):   fraction of liquid water on the vegetation canopy (-)
                 scalarCanairTempTrial,                     & ! intent(out):   trial value of canopy air temperature (K)
                 scalarCanopyTempTrial,                     & ! intent(out):   trial value of canopy temperature (K)
                 scalarCanopyWatTrial,                      & ! intent(out):   trial value of canopy total water (kg m-2)
                 scalarCanopyLiqTrial,                      & ! intent(out):   trial value of canopy liquid water (kg m-2)
                 scalarCanopyIceTrial,                      & ! intent(out):   trial value of canopy ice content (kg m-2)
                 ! output: variables for the snow-soil domain
                 mLayerFracLiqSnow,                         & ! intent(out):   volumetric fraction of water in each snow layer (-)
                 mLayerTempTrial,                           & ! intent(out):   trial vector of layer temperature (K)
                 mLayerVolFracWatTrial,                     & ! intent(out):   trial vector of volumetric total water content (-)
                 mLayerVolFracLiqTrial,                     & ! intent(out):   trial vector of volumetric liquid water content (-)
                 mLayerVolFracIceTrial,                     & ! intent(out):   trial vector of volumetric ice water content (-)
                 mLayerMatricHeadTrial,                     & ! intent(out):   trial vector of matric head (m)
                 ! output: error control
                 err,cmessage)                                ! intent(out):   error control
 if(err/=0)then; message=trim(message)//trim(cmessage); return; endif  ! (check for errors)

 ! check the mass balance for the soil domain
 ! NOTE: this should never fail since did not converge if water balance was not within tolerance=absConvTol_liquid
 if(checkMassBalance)then
  balance0 = sum( (mLayerVolFracLiq(nSnow+1:nLayers)      + mLayerVolFracIce(nSnow+1:nLayers)      )*mLayerDepth(nSnow+1:nLayers) )
  balance1 = sum( (mLayerVolFracLiqTrial(nSnow+1:nLayers) + mLayerVolFracIceTrial(nSnow+1:nLayers) )*mLayerDepth(nSnow+1:nLayers) )
  vertFlux = -(iLayerLiqFluxSoil(nSoil) - iLayerLiqFluxSoil(0))*dt  ! m s-1 --> m
  tranSink = sum(mLayerTranspire)*dt                                ! m s-1 --> m
  baseSink = sum(mLayerBaseflow)*dt                                 ! m s-1 --> m
  compSink = sum(mLayerCompress(1:nSoil) * mLayerDepth(nSnow+1:nLayers) ) ! dimensionless --> m
  liqError = balance1 - (balance0 + vertFlux + tranSink - baseSink - compSink)
  if(abs(liqError) > absConvTol_liquid*10._dp)then  ! *10 to avoid precision issues
   message=trim(message)//'water balance error in the soil domain'
   err=-20; return ! negative error code forces time step reduction and another trial
  endif  ! if there is a water balance error
 endif  ! checking mass balance

 ! compute the melt in each snow and soil layer
 if(nSnow>0) mLayerMeltFreeze(      1:nSnow  ) = -(mLayerVolFracIceTrial(      1:nSnow  ) - mLayerVolFracIce(      1:nSnow  ))*iden_ice
             mLayerMeltFreeze(nSnow+1:nLayers) = -(mLayerVolFracIceTrial(nSnow+1:nLayers) - mLayerVolFracIce(nSnow+1:nLayers))*iden_water

 ! -----
 ! * check that there is sufficient ice content to support the converged sublimation rate...
 ! -----------------------------------------------------------------------------------------

 ! check that sublimation does not exceed the available water on the canopy
 if(computeVegFlux)then
  if(-dt*scalarCanopySublimation > scalarCanopyLiqTrial + scalarCanopyIceTrial)then  ! try again
   message=trim(message)//'insufficient water to support converged canopy sublimation rate'
   err=-20; return  ! negative error code means "try again"
  endif  ! if insufficient water for sublimation
 endif  ! if computing the veg flux

 ! check that sublimation does not exceed the available ice in the top snow layer
 if(nSnow > 0)then ! snow layers exist
  if(-dt*(scalarSnowSublimation/mLayerDepth(1))/iden_ice > mLayerVolFracIceTrial(1))then  ! try again
   message=trim(message)//'insufficient water to support converged surface sublimation rate'
   err=-20; return  ! negative error code means "try again"
  endif  ! if insufficient water for sublimation
 endif  ! if computing the veg flux


 ! -----
 ! * extract state variables for the start of the next time step...
 ! ----------------------------------------------------------------

 ! extract the vegetation states from the state vector
 if(computeVegFlux)then
  scalarCanairTemp = stateVecTrial(ixCasNrg)
  scalarCanopyTemp = stateVecTrial(ixVegNrg)
  scalarCanopyLiq  = scalarCanopyLiqTrial
  scalarCanopyIce  = scalarCanopyIceTrial
 endif

 ! extract state variables for the snow and soil domain
 mLayerTemp(1:nLayers)     = stateVecTrial(ixSnowSoilNrg)
 mLayerMatricHead(1:nSoil) = stateVecTrial(ixSoilOnlyHyd)

 ! save the volumetric liquid water and ice content
 mLayerVolFracLiq = mLayerVolFracLiqTrial  ! computed in updatState
 mLayerVolFracIce = mLayerVolFracIceTrial  ! computed in updatState

 ! ==========================================================================================================================================
 ! ==========================================================================================================================================
 ! ==========================================================================================================================================
 ! ==========================================================================================================================================

 ! deallocate space for the baseflow derivative matrix
 deallocate(dBaseflow_dMatric,stat=err)
 if(err/=0)then; err=20; message=trim(message)//'unable to deallocate space for the baseflow derivatives'; return; endif

 ! end associate statement
 end associate

 end subroutine systemSolv

end module systemSolv_module

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

! access the global print flag
USE globalData,only:globalPrintFlag

! access missing values
USE multiconst,only:integerMissing  ! missing integer
USE multiconst,only:realMissing     ! missing double precision number
USE multiconst,only:quadMissing     ! missing quadruple precision number

! access matrix information
USE globalData,only: nBands         ! length of the leading dimension of the band diagonal matrix
USE globalData,only: ixFullMatrix   ! named variable for the full Jacobian matrix
USE globalData,only: ixBandMatrix   ! named variable for the band diagonal matrix
USE globalData,only: iJac1          ! first layer of the Jacobian to print
USE globalData,only: iJac2          ! last layer of the Jacobian to print

! domain types
USE globalData,only:iname_veg       ! named variables for vegetation
USE globalData,only:iname_snow      ! named variables for snow
USE globalData,only:iname_soil      ! named variables for soil

! state variable type
USE globalData,only:iname_nrgCanair ! named variable defining the energy of the canopy air space
USE globalData,only:iname_nrgCanopy ! named variable defining the energy of the vegetation canopy
USE globalData,only:iname_watCanopy ! named variable defining the mass of total water on the vegetation canopy
USE globalData,only:iname_liqCanopy ! named variable defining the mass of liquid water on the vegetation canopy
USE globalData,only:iname_nrgLayer  ! named variable defining the energy state variable for snow+soil layers
USE globalData,only:iname_watLayer  ! named variable defining the total water state variable for snow+soil layers
USE globalData,only:iname_liqLayer  ! named variable defining the liquid  water state variable for snow+soil layers
USE globalData,only:iname_matLayer  ! named variable defining the matric head state variable for soil layers
USE globalData,only:iname_lmpLayer  ! named variable defining the liquid matric potential state variable for soil layers

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

! provide access to the number of flux variables
USE var_lookup,only:nFlux=>maxvarFlux ! number of model flux variables

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
real(dp),parameter  :: verySmall=1.e-12_dp        ! a very small number (used to check consistency)
real(dp),parameter  :: veryBig=1.e+20_dp          ! a very big number
real(dp),parameter  :: dx = 1.e-8_dp              ! finite difference increment

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
                       firstSubStep,   & ! intent(in): flag to denote first sub-step
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
 USE globalData,only:flux_meta                        ! metadata on the model fluxes
 USE globalData,only:deriv_meta                       ! metadata on the model derivatives
 USE globalData,only:flux2state_meta                  ! metadata on flux-to-state mapping
 USE allocspace_module,only:allocLocal                ! allocate local data structures
 ! simulation of fluxes and residuals given a trial state vector
 USE soil_utils_module,only:matricHead                ! compute the matric head based on volumetric water content
 USE soil_utils_module,only:liquidHead                ! compute the liquid water matric potential
 USE eval8summa_module,only:eval8summa                ! simulation of fluxes and residuals given a trial state vector
 USE summaSolve_module,only:summaSolve                ! calculate the iteration increment, evaluate the new state, and refine if necessary
 ! population/extracction of state vectors
 USE indexState_module,only:indexSplit                ! get state indices
 USE getVectorz_module,only:popStateVec               ! populate the state vector
 USE getVectorz_module,only:varExtract                ! extract variables from the state vector
 USE updateVars_module,only:updateVars                ! update prognostic variables
 ! numerical recipes utility modules
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
 character(LEN=256)              :: cmessage                      ! error message of downwind routine
 integer(i4b)                    :: iter                          ! iteration index
 integer(i4b)                    :: iSoil                         ! index of soil layer
 integer(i4b)                    :: iLayer                        ! index of layer in the snow+soil domain
 integer(i4b)                    :: iState                        ! index of model state
 integer(i4b)                    :: nLeadDim                      ! length of the leading dimension of the Jacobian matrix (nBands or nState)
 integer(i4b)                    :: local_ixGroundwater           ! local index for groundwater representation
 real(dp),parameter              :: tempAccelerate=0.00_dp        ! factor to force initial canopy temperatures to be close to air temperature
 real(dp),parameter              :: xMinCanopyWater=0.0001_dp     ! minimum value to initialize canopy water (kg m-2)
 ! ------------------------------------------------------------------------------------------------------
 ! * state variables and diagniostic variables
 ! ------------------------------------------------------------------------------------------------------
 ! trial state variables
 real(dp)                        :: scalarCanairTempTrial         ! trial value for temperature of the canopy air space (K)
 real(dp)                        :: scalarCanopyTempTrial         ! trial value for temperature of the vegetation canopy (K)
 real(dp)                        :: scalarCanopyWatTrial          ! trial value for liquid water storage in the canopy (kg m-2)
 real(dp),dimension(nLayers)     :: mLayerTempTrial               ! trial vector for temperature of layers in the snow and soil domains (K)
 real(dp),dimension(nLayers)     :: mLayerVolFracWatTrial         ! trial vector for volumetric fraction of total water (-)
 real(dp),dimension(nSoil)       :: mLayerMatricHeadTrial         ! trial vector for total water matric potential (m)
 real(dp),dimension(nSoil)       :: mLayerMatricHeadLiqTrial      ! trial vector for liquid water matric potential (m)
 ! diagnostic variables
 real(dp)                        :: scalarCanopyLiqTrial          ! trial value for mass of liquid water on the vegetation canopy (kg m-2)
 real(dp)                        :: scalarCanopyIceTrial          ! trial value for mass of ice on the vegetation canopy (kg m-2)
 real(dp),dimension(nLayers)     :: mLayerVolFracLiqTrial         ! trial vector for volumetric fraction of liquid water (-)
 real(dp),dimension(nLayers)     :: mLayerVolFracIceTrial         ! trial vector for volumetric fraction of ice (-)
 real(dp),dimension(nLayers)     :: mLayerVolFracIceInit          ! initial vector for volumetric fraction of ice (-)
 ! ------------------------------------------------------------------------------------------------------
 ! * operator splitting
 ! ------------------------------------------------------------------------------------------------------
 real(dp)                        :: dtSplit                       ! time step for a given operator-splitting operation
 real(dp)                        :: dt_wght                       ! weight given to a given flux calculation
 integer(i4b),parameter          :: fullyCoupled=1001             ! named variable for the fully coupled solution
 integer(i4b),parameter          :: deCoupled_nrgMass=1002        ! named variable for the solution where energy and mass is decoupled
 integer(i4b)                    :: ixSplitOption=decoupled_nrgMass  ! selected operator splitting method
 !integer(i4b)                    :: ixSplitOption=fullyCoupled   ! selected operator splitting method
 integer(i4b),parameter          :: explicitEuler=2001            ! explicit Euler solution
 integer(i4b),parameter          :: implicitEuler=2002            ! implicit Euler solution
 integer(i4b)                    :: ixIterOption=explicitEuler    ! selected option for the iterations
 !integer(i4b)                    :: ixIterOption=implicitEuler    ! selected option for the iterations
 integer(i4b)                    :: nOperSplit                    ! number of splitting operations
 integer(i4b)                    :: iSplit                        ! index of splitting operation
 integer(i4b),parameter          :: nrgSplit=1                    ! order in sequence for the energy operation
 integer(i4b),parameter          :: massSplit=2                   ! order in sequence for the mass operation
 integer(i4b)                    :: iSubstep                      ! index of substep for a given variable
 integer(i4b)                    :: nSubstep                      ! number of substeps for a given variable
 logical(lgt),dimension(nState)  :: stateMask                     ! mask defining desired state variables
 logical(lgt),dimension(nFlux)   :: fluxMask                      ! mask defining desired flux variables
 integer(i4b)                    :: nSubset                       ! number of selected state variables for a given split
 integer(i4b)                    :: iVar                          ! index of variables in data structures
 ! ------------------------------------------------------------------------------------------------------
 ! * model solver
 ! ------------------------------------------------------------------------------------------------------
 logical(lgt),parameter          :: forceFullMatrix=.false.       ! flag to force the use of the full Jacobian matrix
 logical(lgt)                    :: firstFluxCall                 ! flag to define the first flux call
 integer(i4b)                    :: ixMatrix                      ! form of matrix (band diagonal or full matrix)
 type(var_dlength)               :: flux_init                     ! model fluxes at the start of the time step 
 type(var_dlength)               :: flux_temp                     ! temporary model fluxes 
 type(var_dlength)               :: deriv_data                    ! derivatives in model fluxes w.r.t. relevant state variables 
 integer(i4b)                    :: ixSaturation                  ! index of the lowest saturated layer (NOTE: only computed on the first iteration)
 real(dp),allocatable            :: dBaseflow_dMatric(:,:)        ! derivative in baseflow w.r.t. matric head (s-1)  ! NOTE: allocatable, since not always needed
 real(dp),allocatable            :: stateVecInit(:)               ! initial state vector (mixed units)
 real(dp),allocatable            :: stateVecTrial(:)              ! trial state vector (mixed units)
 real(dp),allocatable            :: stateVecNew(:)                ! new state vector (mixed units)
 real(dp),allocatable            :: fluxVec0(:)                   ! flux vector (mixed units)
 real(dp),allocatable            :: fScale(:)                     ! characteristic scale of the function evaluations (mixed units)
 real(dp),allocatable            :: xScale(:)                     ! characteristic scale of the state vector (mixed units)
 real(dp),allocatable            :: dMat(:)                       ! diagonal matrix (excludes flux derivatives)
 real(qp),allocatable            :: sMul(:)       ! NOTE: qp      ! multiplier for state vector for the residual calculations
 real(qp),allocatable            :: rVec(:)       ! NOTE: qp      ! residual vector
 real(dp),allocatable            :: rAdd(:)                       ! additional terms in the residual vector
 real(dp)                        :: fOld,fNew                     ! function values (-); NOTE: dimensionless because scaled
 logical(lgt)                    :: feasible                      ! flag to define the feasibility of the solution
 logical(lgt)                    :: converged                     ! convergence flag
 real(dp),allocatable            :: resSinkNew(:)                 ! additional terms in the residual vector
 real(dp),allocatable            :: fluxVecNew(:)                 ! new flux vector
 real(qp),allocatable            :: resVecNew(:)  ! NOTE: qp      ! new residual vector
 real(dp),allocatable            :: solutionError(:)              ! vector of errors in the model solution
 real(dp),parameter              :: safety=0.85_dp                ! safety factor in adaptive sub-stepping
 real(dp),parameter              :: reduceMin=0.1_dp              ! mimimum factor that time step is reduced
 real(dp),parameter              :: increaseMax=4.0_dp            ! maximum factor that time step is increased
 real(dp),parameter              :: errorTol=100._dp               ! error tolerance in the explicit solution
 real(dp),dimension(1)           :: errorMax                      ! maximum error in explicit solution
 real(dp)                        :: dtSubstep                     ! length of the substep (seconds)
 real(dp)                        :: dtSum                         ! keep track of the portion of the substep that is completed
 ! ------------------------------------------------------------------------------------------------------
 ! * mass balance checks
 ! ------------------------------------------------------------------------------------------------------
 logical(lgt),parameter          :: checkTranspire=.false.        ! flag to check transpiration
 logical(lgt)                    :: checkMassBalance              ! flag to check the mass balance
 real(dp)                        :: soilBalance0,soilBalance1     ! soil storage at start and end of time step
 real(dp)                        :: canopyBalance0,canopyBalance1 ! soil storage at start and end of time step
 real(dp)                        :: vertFlux                      ! change in storage due to vertical fluxes
 real(dp)                        :: tranSink,baseSink,compSink    ! change in storage sue to sink terms
 real(dp)                        :: liqError                      ! water balance error
 ! ---------------------------------------------------------------------------------------
 ! point to variables in the data structures
 ! ---------------------------------------------------------------------------------------
 globalVars: associate(&
 ! model decisions
 ixGroundwater           => model_decisions(iLookDECISIONS%groundwatr)%iDecision   ,& ! intent(in):    [i4b]    groundwater parameterization
 ixSpatialGroundwater    => model_decisions(iLookDECISIONS%spatial_gw)%iDecision   ,& ! intent(in):    [i4b]    spatial representation of groundwater (local-column or single-basin)
 ! domain boundary conditions
 airtemp                 => forc_data%var(iLookFORCE%airtemp)                      ,& ! intent(in):    [dp]     temperature of the upper boundary of the snow and soil domains (K)
 ! vector of energy and hydrology indices for the snow and soil domains
 ixSnowSoilNrg           => indx_data%var(iLookINDEX%ixSnowSoilNrg)%dat            ,& ! intent(in):    [i4b(:)] index in the state subset for energy state variables in the snow+soil domain
 ixSnowSoilHyd           => indx_data%var(iLookINDEX%ixSnowSoilHyd)%dat            ,& ! intent(in):    [i4b(:)] index in the state subset for hydrology state variables in the snow+soil domain
 nSnowSoilNrg            => indx_data%var(iLookINDEX%nSnowSoilNrg )%dat(1)         ,& ! intent(in):    [i4b]    number of energy state variables in the snow+soil domain
 nSnowSoilHyd            => indx_data%var(iLookINDEX%nSnowSoilHyd )%dat(1)         ,& ! intent(in):    [i4b]    number of hydrology state variables in the snow+soil domain
 ! indices of model state variables
 ixStateType             => indx_data%var(iLookINDEX%ixStateType)%dat              ,& ! intent(in):    [i4b(:)] indices defining the type of the state (ixNrgState...)
 ixHydCanopy             => indx_data%var(iLookINDEX%ixHydCanopy)%dat              ,& ! intent(in):    [i4b(:)] indices IN THE FULL VECTOR for hydrology states in the canopy domain
 ixHydLayer              => indx_data%var(iLookINDEX%ixHydLayer)%dat               ,& ! intent(in):    [i4b(:)] indices IN THE FULL VECTOR for hydrology states in the snow+soil domain
 ixCasNrg                => indx_data%var(iLookINDEX%ixCasNrg)%dat(1)              ,& ! intent(in):    [i4b]    index of canopy air space energy state variable
 ixVegNrg                => indx_data%var(iLookINDEX%ixVegNrg)%dat(1)              ,& ! intent(in):    [i4b]    index of canopy energy state variable
 ixVegHyd                => indx_data%var(iLookINDEX%ixVegHyd)%dat(1)              ,& ! intent(in):    [i4b]    index of canopy hydrology state variable (mass)
 ! vegetation parameters
 canopyDepth             => diag_data%var(iLookDIAG%scalarCanopyDepth)%dat(1)      ,& ! intent(in):    [dp]     canopy depth (m)
 ! snow parameters
 snowfrz_scale           => mpar_data%var(iLookPARAM%snowfrz_scale)                ,& ! intent(in):    [dp]     scaling parameter for the snow freezing curve (K-1)
 ! soil parameters
 vGn_m                   => diag_data%var(iLookDIAG%scalarVGn_m)%dat(1)            ,& ! intent(in):    [dp]     van Genutchen "m" parameter (-)
 vGn_n                   => mpar_data%var(iLookPARAM%vGn_n)                        ,& ! intent(in):    [dp]     van Genutchen "n" parameter (-)
 vGn_alpha               => mpar_data%var(iLookPARAM%vGn_alpha)                    ,& ! intent(in):    [dp]     van Genutchen "alpha" parameter (m-1)
 theta_sat               => mpar_data%var(iLookPARAM%theta_sat)                    ,& ! intent(in):    [dp]     soil porosity (-)
 theta_res               => mpar_data%var(iLookPARAM%theta_res)                    ,& ! intent(in):    [dp]     soil residual volumetric water content (-)
 specificStorage         => mpar_data%var(iLookPARAM%specificStorage)              ,& ! intent(in):    [dp]     specific storage coefficient (m-1)
 ! convergence parameters
 absConvTol_liquid       => mpar_data%var(iLookPARAM%absConvTol_liquid)            ,& ! intent(in):    [dp]     absolute convergence tolerance for vol frac liq water (-)
 ! model diagnostic variables (fraction of liquid water)
 scalarFracLiqVeg        => diag_data%var(iLookDIAG%scalarFracLiqVeg)%dat(1)       ,& ! intent(out):   [dp]     fraction of liquid water on vegetation (-)
 mLayerFracLiqSnow       => diag_data%var(iLookDIAG%mLayerFracLiqSnow)%dat         ,& ! intent(out):   [dp(:)]  fraction of liquid water in each snow layer (-)
 mLayerMeltFreeze        => diag_data%var(iLookDIAG%mLayerMeltFreeze)%dat          ,& ! intent(out):   [dp(:)]  melt-freeze in each snow and soil layer (kg m-3)
 ! model fluxes and derivatives
 mLayerCompress          => diag_data%var(iLookDIAG%mLayerCompress)%dat            ,& ! intent(out):   [dp(:)]  change in storage associated with compression of the soil matrix (-)
 ! model state variables (vegetation canopy)
 mLayerDepth             => prog_data%var(iLookPROG%mLayerDepth)%dat               ,& ! intent(in):    [dp(:)]  depth of each layer in the snow-soil sub-domain (m)
 scalarCanairTemp        => prog_data%var(iLookPROG%scalarCanairTemp)%dat(1)       ,& ! intent(inout): [dp]     temperature of the canopy air space (K)
 scalarCanopyTemp        => prog_data%var(iLookPROG%scalarCanopyTemp)%dat(1)       ,& ! intent(inout): [dp]     temperature of the vegetation canopy (K)
 scalarCanopyIce         => prog_data%var(iLookPROG%scalarCanopyIce)%dat(1)        ,& ! intent(inout): [dp]     mass of ice on the vegetation canopy (kg m-2)
 scalarCanopyLiq         => prog_data%var(iLookPROG%scalarCanopyLiq)%dat(1)        ,& ! intent(inout): [dp]     mass of liquid water on the vegetation canopy (kg m-2)
 scalarCanopyWat         => prog_data%var(iLookPROG%scalarCanopyWat)%dat(1)        ,& ! intent(inout): [dp]     mass of total water on the vegetation canopy (kg m-2)
 ! model state variables (snow and soil domains)
 mLayerTemp              => prog_data%var(iLookPROG%mLayerTemp)%dat                ,& ! intent(inout): [dp(:)]  temperature of each snow/soil layer (K)
 mLayerVolFracIce        => prog_data%var(iLookPROG%mLayerVolFracIce)%dat          ,& ! intent(inout): [dp(:)]  volumetric fraction of ice (-)
 mLayerVolFracLiq        => prog_data%var(iLookPROG%mLayerVolFracLiq)%dat          ,& ! intent(inout): [dp(:)]  volumetric fraction of liquid water (-)
 mLayerVolFracWat        => prog_data%var(iLookPROG%mLayerVolFracWat)%dat          ,& ! intent(inout): [dp(:)]  volumetric fraction of total water (-)
 mLayerMatricHead        => prog_data%var(iLookPROG%mLayerMatricHead)%dat          ,& ! intent(inout): [dp(:)]  matric head (m)
 mLayerMatricHeadLiq     => diag_data%var(iLookDIAG%mLayerMatricHeadLiq)%dat       ,& ! intent(inout): [dp(:)]  matric potential of liquid water (m)
 scalarAquiferStorage    => prog_data%var(iLookPROG%scalarAquiferStorage)%dat(1)    & ! intent(inout): [dp   ]  aquifer storage (m)
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

 ! save water storage at the start of the step
 canopyBalance0 = merge(scalarCanopyWat, realMissing, computeVegFlux)
 soilBalance0   = sum( (mLayerVolFracLiq(nSnow+1:nLayers)      + mLayerVolFracIce(nSnow+1:nLayers)      )*mLayerDepth(nSnow+1:nLayers) )
 
 ! save volumetric ice content at the start of the step
 ! NOTE: used for volumetric loss due to melt-freeze
 mLayerVolFracIceInit(:) = mLayerVolFracIce(:)

 ! compute the total water content in snow and soil
 ! NOTE: no ice expansion allowed for soil
 if(nSnow>0)& 
 mLayerVolFracWat(      1:nSnow  ) = mLayerVolFracLiq(      1:nSnow  ) + mLayerVolFracIce(      1:nSnow  )*(iden_ice/iden_water)
 mLayerVolFracWat(nSnow+1:nLayers) = mLayerVolFracLiq(nSnow+1:nLayers) + mLayerVolFracIce(nSnow+1:nLayers)

 ! compute the liquid water matric potential (m)
 ! NOTE: include ice content as part of the solid porosity - major effect of ice is to reduce the pore size; ensure that effSat=1 at saturation
 ! (from Zhao et al., J. Hydrol., 1997: Numerical analysis of simultaneous heat and mass transfer...)
 do iSoil=1,nSoil
  call liquidHead(mLayerMatricHead(iSoil),mLayerVolFracLiq(nSnow+iSoil),mLayerVolFracIce(nSnow+iSoil),  & ! input:  state variables
                  vGn_alpha,vGn_n,theta_sat,theta_res,vGn_m,                                            & ! input:  parameters
                  matricHeadLiq=mLayerMatricHeadLiq(iSoil),                                             & ! output: liquid water matric potential (m)
                  err=err,message=cmessage)                                                               ! output: error control
  if(err/=0)then; message=trim(message)//trim(cmessage); return; endif
 end do  ! looping through soil layers (computing liquid water matric potential)

 ! define the number of operator splits
 select case(ixSplitOption)
  case(fullyCoupled);      nOperSplit=1
  case(deCoupled_nrgMass); nOperSplit=2
  case default; err=20; message=trim(message)//'unable to identify operator splitting strategy'; return
 end select

 ! modify the groundwater representation for this single-column implementation
 select case(ixSpatialGroundwater)
  case(singleBasin); local_ixGroundwater = noExplicit    ! force no explicit representation of groundwater at the local scale
  case(localColumn); local_ixGroundwater = ixGroundwater ! go with the specified decision
  case default; err=20; message=trim(message)//'unable to identify spatial representation of groundwater'; return
 end select ! (modify the groundwater representation for this single-column implementation)

 ! allocate space for the model fluxes at the start of the time step
 call allocLocal(flux_meta(:),flux_init,nSnow,nSoil,err,cmessage)
 if(err/=0)then; err=20; message=trim(message)//trim(cmessage); return; endif

 ! allocate space for the temporary model flux structure
 call allocLocal(flux_meta(:),flux_temp,nSnow,nSoil,err,cmessage)
 if(err/=0)then; err=20; message=trim(message)//trim(cmessage); return; endif

 ! allocate space for the derivative structure
 call allocLocal(deriv_meta(:),deriv_data,nSnow,nSoil,err,cmessage)
 if(err/=0)then; err=20; message=trim(message)//trim(cmessage); return; end if

 ! allocate space for the baseflow derivatives
 if(ixGroundwater==qbaseTopmodel)then
  allocate(dBaseflow_dMatric(nSoil,nSoil),stat=err)  ! baseflow depends on total storage in the soil column, hence on matric head in every soil layer
 else
  allocate(dBaseflow_dMatric(0,0),stat=err)          ! allocate zero-length dimnensions to avoid passing around an unallocated matrix
 end if 
 if(err/=0)then; err=20; message=trim(message)//'unable to allocate space for the baseflow derivatives'; return; end if

 ! ==========================================================================================================================================
 ! ==========================================================================================================================================
 ! ==========================================================================================================================================
 ! ==========================================================================================================================================

 ! check
 if(nrgSplit > massSplit)then
  message=trim(message)//'currently energy split must be done first in order for canopy water balance to close: '//&
   'if implement Strang Splitting then need to include start-of-step energy fluxes in average flux calculations'
  err=20; return
 endif

 ! operator splitting loop
 do iSplit=1,nOperSplit 

  !print*, 'iSplit, nOperSplit = ', iSplit, nOperSplit

  ! -----
  ! * define subsets for a given split...
  ! -------------------------------------
 
  ! define need to check the mass balance
  checkMassBalance = ( (ixSplitOption==deCoupled_nrgMass .and. iSplit==massSplit) .or. ixSplitOption/=deCoupled_nrgMass )
 
  ! modify state variable names for the mass split
  if(ixSplitOption==deCoupled_nrgMass .and. iSplit==massSplit)then

   ! modify the state type names associated with the state vector
   if(computeVegFlux)then
    where(ixStateType(ixHydCanopy)==iname_watCanopy) ixStateType(ixHydCanopy)=iname_liqCanopy
   endif
   where(ixStateType(ixHydLayer) ==iname_watLayer)  ixStateType(ixHydLayer) =iname_liqLayer
   where(ixStateType(ixHydLayer) ==iname_matLayer)  ixStateType(ixHydLayer) =iname_lmpLayer

   ! modify the state type names associated with the flux mapping structure
   do iVar=1,size(flux_meta)
    ! (mass of total water on the vegetation canopy --> mass of liquid water)
    if(flux2state_meta(iVar)%state1==iname_watCanopy) flux2state_meta(iVar)%state1=iname_liqCanopy
    if(flux2state_meta(iVar)%state2==iname_watCanopy) flux2state_meta(iVar)%state2=iname_liqCanopy
    ! (volumetric total water in the snow+soil domain --> volumetric liquid water)
    if(flux2state_meta(iVar)%state1==iname_watLayer)  flux2state_meta(iVar)%state1=iname_liqLayer
    if(flux2state_meta(iVar)%state2==iname_watLayer)  flux2state_meta(iVar)%state2=iname_liqLayer
    ! (total water matric potential in the snow+soil domain --> liquid water matric potential)
    if(flux2state_meta(iVar)%state1==iname_matLayer)  flux2state_meta(iVar)%state1=iname_lmpLayer
    if(flux2state_meta(iVar)%state2==iname_matLayer)  flux2state_meta(iVar)%state2=iname_lmpLayer
   end do

  endif  ! if modifying state variables for the mass split

  ! define mask for the decoupled solutions
  if(ixSplitOption==deCoupled_nrgMass)then
   select case(iSplit)
    case(nrgSplit);  stateMask = (ixStateType==iname_nrgCanair .or. ixStateType==iname_nrgCanopy .or. ixStateType==iname_nrgLayer)
    case(massSplit); stateMask = (ixStateType==iname_liqCanopy .or. ixStateType==iname_liqLayer  .or. ixStateType==iname_lmpLayer)
    case default; err=20; message=trim(message)//'unable to identify splitting operation'; return
   end select
  else  ! (fully coupled solutions)
   stateMask(:) = .true.  ! use all state variables
  endif

  ! get the number of selected state variables
  nSubset = count(stateMask)

  ! check splitting operation
  if(nSubset==0)then
   select case(iSplit)
    case(nrgSplit);  err=20; message=trim(message)//'[nrgSplit] no state variables in the splitting operation';  return
    case(massSplit); err=20; message=trim(message)//'[massSplit] no state variables in the splitting operation'; return
    case default;    err=20; message=trim(message)//'unable to identify splitting operation'; return
   end select
  endif

  ! get indices for a given split
  call indexSplit(stateMask,                   & ! intent(in)    : logical vector (.true. if state is in the subset)
                  nSnow,nSoil,nLayers,nSubset, & ! intent(in)    : number of snow and soil layers, and total number of layers
                  indx_data,                   & ! intent(inout) : index data structure
                  err,cmessage)                  ! intent(out)   : error control
  if(err/=0)then; message=trim(message)//trim(cmessage); return; endif

  ! make association with model indices defined in indexSplit
  stateSubset: associate(&
  ixSnowSoilNrg       => indx_data%var(iLookINDEX%ixSnowSoilNrg)%dat,      & ! intent(in): [i4b(:)] indices for energy states in the snow+soil subdomain
  ixSnowOnlyHyd       => indx_data%var(iLookINDEX%ixSnowOnlyHyd)%dat,      & ! intent(in): [i4b(:)] indices for hydrology states in the snow subdomain
  ixStateType_subset  => indx_data%var(iLookINDEX%ixStateType_subset)%dat, & ! intent(in): [i4b(:)] [state subset] type of desired model state variables
  ixDomainType_subset => indx_data%var(iLookINDEX%ixDomainType_subset)%dat & ! intent(in): [i4b(:)] [state subset] type of desired model state variables
  ) ! associations
  
  ! allocate space for solution and scaling vectors
  allocate(stateVecInit(nSubset), stateVecTrial(nSubset), stateVecNew(nSubset), fluxVec0(nSubset), fluxVecNew(nSubset), fScale(nSubset), xScale(nSubset), stat=err)
  if(err/=0)then; err=20; message=trim(message)//'unable to allocate space for the solution and scaling vectors'; return; endif

  ! allocate space for the diagonal matrix, multipliers and residual vectors
  allocate(dMat(nSubset), sMul(nSubset), rVec(nSubset), rAdd(nSubset), resSinkNew(nSubset), resVecNew(nSubset), stat=err)
  if(err/=0)then; err=20; message=trim(message)//'unable to allocate space for the diagonal matrix, multipliers, and residual vectors'; return; endif

  ! if explicit Euler, then allocate space for the solution error
  if(ixIterOption==explicitEuler)then
   allocate(solutionError(nSubset), stat=err)
   if(err/=0)then; err=20; message=trim(message)//'unable to allocate space for the solution error'; return; endif
  endif

  ! identify the matrix solution method
  ! NOTE: uses nSubset so must be defined here
  ! (the type of matrix used to solve the linear system A.X=B)
  if(local_ixGroundwater==qbaseTopmodel .or. forceFullMatrix)then
   nLeadDim=nSubset        ! length of the leading dimension
   ixMatrix=ixFullMatrix   ! named variable to denote the full Jacobian matrix
  else
   nLeadDim=nBands         ! length of the leading dimension
   ixMatrix=ixBandMatrix   ! named variable to denote the band-diagonal matrix
  endif

  ! define the mask of the fluxes used
  do iVar=1,size(flux_meta)
   fluxMask(iVar) = any(ixStateType_subset==flux2state_meta(iVar)%state1) .or. any(ixStateType_subset==flux2state_meta(iVar)%state2)
  end do

  ! initialize the model fluxes (some model fluxes are not computed in the iterations)
  do iVar=1,size(flux_meta)
   ! copy over the master flux structure since some model fluxes are not computed in the iterations
   flux_init%var(iVar)%dat(:) = flux_data%var(iVar)%dat(:)
   flux_temp%var(iVar)%dat(:) = flux_data%var(iVar)%dat(:)
   ! set master flux vector to zero for the fluxes computed here
   if(fluxMask(iVar))then
    flux_data%var(iVar)%dat(:) = 0._dp
   endif
  end do

  ! -----
  ! * variable-dependent sub-stepping...
  ! ------------------------------------
  
  ! define the number of substeps for each solution
  if(ixSplitOption==deCoupled_nrgMass)then
   select case(iSplit)
    case(nrgSplit);  nSubstep=1
    case(massSplit); nSubstep=1
    case default; err=20; message=trim(message)//'unable to identify splitting operation'; return
   end select

  ! fully implicit, then nSubstep=1
  else
   nSubstep=1
  endif

  ! define the length of the substep
  dtSplit   = dt/real(nSubstep, kind(dp))
  dtSubstep = dtSplit

  ! initialize subStep
  dtSum    = 0._dp  ! keep track of the portion of the time step that is completed
  iSubstep = 0

  ! loop through substeps
  ! NOTE: continuous do statement with exit clause
  substeps: do 

   ! increment substep
   iSubstep = iSubstep+1

   ! test
   !write(*,'(a,1x,3(i5,1x),a)') 'iSubstep, iSplit, nSnow = ', iSubstep, iSplit, nSnow, merge('operSplitting','fully_coupled',ixSplitOption==deCoupled_nrgMass)

   ! **************************************************************************************************************************
   ! **************************************************************************************************************************
   ! **************************************************************************************************************************
   ! *** NUMERICAL SOLUTION FOR A GIVEN SUBSTEP AND SPLIT *********************************************************************
   ! **************************************************************************************************************************
   ! **************************************************************************************************************************
   ! **************************************************************************************************************************

   ! -----
   ! * populate state vectors...
   ! ---------------------------

   ! initialize the global print flag
   globalPrintFlag=.false.

   ! initialize state vectors
   call popStateVec(&
                    ! input
                    nSubset,                          & ! intent(in):    number of desired state variables
                    prog_data,                        & ! intent(in):    model prognostic variables for a local HRU
                    diag_data,                        & ! intent(in):    model diagnostic variables for a local HRU
                    indx_data,                        & ! intent(in):    indices defining model states and layers
                    ! output
                    stateVecInit,                     & ! intent(out):   initial model state vector (mixed units)
                    fScale,                           & ! intent(out):   function scaling vector (mixed units)
                    xScale,                           & ! intent(out):   variable scaling vector (mixed units)
                    sMul,                             & ! intent(out):   multiplier for state vector (used in the residual calculations)
                    dMat,                             & ! intent(out):   diagonal of the Jacobian matrix (excludes fluxes) 
                    err,cmessage)                       ! intent(out):   error control
   if(err/=0)then; message=trim(message)//trim(cmessage); return; endif  ! (check for errors)
 
   ! -----
   ! * compute the initial function evaluation...
   ! --------------------------------------------
   
   ! initialize the trial state vectors
   stateVecTrial = stateVecInit
  
   ! need to intialize canopy water at a positive value
   if(ixVegHyd/=integerMissing)then
    if(scalarCanopyWat < xMinCanopyWater) stateVecTrial(ixVegHyd) = scalarCanopyWat + xMinCanopyWater
   endif
   
   ! try to accelerate solution for energy
   if(iSplit==1 .and. iSubStep==1)then
    if(ixCasNrg/=integerMissing) stateVecTrial(ixCasNrg) = stateVecInit(ixCasNrg) + (airtemp - stateVecInit(ixCasNrg))*tempAccelerate
    if(ixVegNrg/=integerMissing) stateVecTrial(ixVegNrg) = stateVecInit(ixVegNrg) + (airtemp - stateVecInit(ixVegNrg))*tempAccelerate
   endif

   ! compute the flux and the residual vector for a given state vector
   ! NOTE 1: The derivatives computed in eval8summa are used to calculate the Jacobian matrix for the first iteration
   ! NOTE 2: The Jacobian matrix together with the residual vector is used to calculate the first iteration increment
   call eval8summa(&
                   ! input: model control
                   dtSubstep,               & ! intent(in):    length of the time step (seconds)
                   nSnow,                   & ! intent(in):    number of snow layers
                   nSoil,                   & ! intent(in):    number of soil layers
                   nLayers,                 & ! intent(in):    number of layers
                   nSubset,                 & ! intent(in):    number of state variables in the current subset
                   firstSubStep,            & ! intent(in):    flag to indicate if we are processing the first sub-step
                   firstFluxCall,           & ! intent(inout): flag to indicate if we are processing the first flux call
                   .true.,                  & ! intent(in):    flag to indicate if we are processing the first iteration in a splitting operation
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
                   flux_init,               & ! intent(inout): model fluxes for a local HRU (initial flux structure)
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
    message=trim(message)//'unfeasible state vector for the initial function evaluation'
    err=20; return
   endif
 
   ! copy over the initial flux structure since some model fluxes are not computed in the iterations
   do concurrent ( iVar=1:size(flux_meta) )
    flux_temp%var(iVar)%dat(:) = flux_init%var(iVar)%dat(:)
   end do

   ! ** if explicit Euler, then estimate state vector at the end of the time step
   if(ixIterOption==explicitEuler)then
    ! --> update state vector
    stateVecTrial(:) = stateVecInit(:) + (fluxVec0(:)*dtSubstep + rAdd(:))/real(sMul(:), dp)
    ! --> impose constraints
    do concurrent (iState=1:nSubset)
     select case( ixStateType_subset(iState) )
      ! impose non-negativity constraints for mass states
      case(iname_watCanopy,iname_liqCanopy,iname_watLayer,iname_liqLayer)
       if(stateVecTrial(iState) < 0._dp) stateVecTrial(iState)=0._dp
      ! impose below-freezing constraints for snow temperature
      case(iname_nrgLayer)
       if(ixDomainType_subset(iState)==iname_snow .and. stateVecTrial(iState) > Tfreeze) stateVecTrial(iState)=Tfreeze
      ! skip unconstrained states
      case default; cycle
     end select  ! selecting state variables
    end do ! looping through states
   endif  ! if explicit Euler

   ! ==========================================================================================================================================
   ! ==========================================================================================================================================
   ! ==========================================================================================================================================
   ! ==========================================================================================================================================
   
   ! (1) MAIN ITERATION LOOP...
   ! **************************

   ! iterate
   ! NOTE: this do loop is skipped in the explicitEuler solution (localMaxIter=0)
   do iter=1,maxIter

    ! print iteration count
    !print*, '*** iter, dt = ', iter, dtSubstep, merge('nrg','wat',iSplit==nrgSplit)
 
    ! keep track of the number of iterations
    niter = iter+1  ! +1 because xFluxResid was moved outside the iteration loop (for backwards compatibility)
  
    ! compute the next trial state vector
    !  1) Computes the Jacobian matrix based on derivatives from the last flux evaluation
    !  2) Computes the iteration increment based on Jacobian and residuals from the last flux evaluation
    !  3) Computes new fluxes and derivatives, new residuals, and (if necessary) refines the state vector
    ! NOTE: only returns the flux vector and function evaluation when ixIterOption==explicitEuler
    call summaSolve(&
                    ! input: model control
                    dtSubstep,                     & ! intent(in):    length of the time step (seconds)
                    (ixIterOption==explicitEuler), & ! intent(in):    logical flag to only return the flux and function evaluation
                    iter,                          & ! intent(in):    iteration index
                    nSnow,                         & ! intent(in):    number of snow layers
                    nSoil,                         & ! intent(in):    number of soil layers
                    nLayers,                       & ! intent(in):    total number of layers
                    nLeadDim,                      & ! intent(in):    length of the leading dimension of the Jacobian matrix (either nBands or nState) 
                    nSubset,                       & ! intent(in):    total number of state variables
                    ixMatrix,                      & ! intent(in):    type of matrix (full or band diagonal)
                    firstSubStep,                  & ! intent(in):    flag to indicate if we are processing the first sub-step
                    firstFluxCall,                 & ! intent(inout): flag to indicate if we are processing the first flux call
                    computeVegFlux,                & ! intent(in):    flag to indicate if we need to compute fluxes over vegetation
                    ! input: state vectors
                    stateVecTrial,                 & ! intent(in):    trial state vector
                    fScale,                        & ! intent(in):    function scaling vector
                    xScale,                        & ! intent(in):    "variable" scaling vector, i.e., for state variables
                    rVec,                          & ! intent(in):    residual vector
                    sMul,                          & ! intent(in):    state vector multiplier (used in the residual calculations)
                    dMat,                          & ! intent(inout): diagonal matrix (excludes flux derivatives)
                    fOld,                          & ! intent(in):    old function evaluation
                    ! input: data structures       
                    model_decisions,               & ! intent(in):    model decisions
                    type_data,                     & ! intent(in):    type of vegetation and soil
                    attr_data,                     & ! intent(in):    spatial attributes
                    mpar_data,                     & ! intent(in):    model parameters
                    forc_data,                     & ! intent(in):    model forcing data
                    bvar_data,                     & ! intent(in):    average model variables for the entire basin
                    prog_data,                     & ! intent(in):    model prognostic variables for a local HRU
                    indx_data,                     & ! intent(in):    index data
                    ! input-output: data structures
                    diag_data,                     & ! intent(inout): model diagnostic variables for a local HRU
                    flux_temp,                     & ! intent(inout): model fluxes for a local HRU (temporary structure)
                    deriv_data,                    & ! intent(inout): derivatives in model fluxes w.r.t. relevant state variables
                    ! input-output: baseflow       
                    ixSaturation,                  & ! intent(inout): index of the lowest saturated layer (NOTE: only computed on the first iteration)
                    dBaseflow_dMatric,             & ! intent(inout): derivative in baseflow w.r.t. matric head (s-1)
                    ! output
                    stateVecNew,                   & ! intent(out):   new state vector
                    fluxVecNew,                    & ! intent(out):   new flux vector
                    resSinkNew,                    & ! intent(out):   additional (sink) terms on the RHS of the state equa
                    resVecNew,                     & ! intent(out):   new residual vector
                    fNew,                          & ! intent(out):   new function evaluation
                    converged,                     & ! intent(out):   convergence flag
                    err,cmessage)                    ! intent(out):   error control
    if(err/=0)then; message=trim(message)//trim(cmessage); return; endif  ! (check for errors)
    
    ! update function evaluation, residual vector, and states
    ! NOTE 1: The derivatives computed in summaSolve are used to calculate the Jacobian matrix at the next iteration
    ! NOTE 2: The Jacobian matrix together with the residual vector is used to calculate the new iteration increment
    if(ixIterOption==implicitEuler)then
     fOld          = fNew
     rVec          = resVecNew 
     stateVecTrial = stateVecNew
    endif

    ! print progress
    !write(*,'(a,10(f16.14,1x))') 'rVec                  = ', rVec(iJac1:iJac2)
    !write(*,'(a,10(f16.10,1x))') 'fluxVecNew            = ', fluxVecNew(iJac1:iJac2)*dtSubstep
    !write(*,'(a,10(f16.10,1x))') 'stateVecTrial         = ', stateVecTrial(iJac1:iJac2)
    !print*, 'PAUSE: check states and fluxes'; read(*,*) 
  
    ! exit iteration loop if converged
    if(converged .or. ixIterOption==explicitEuler) exit
   
    ! check convergence
    if(niter==maxiter)then
     message=trim(message)//'failed to converge'
     err=-20; return
    endif
    !print*, 'PAUSE: iterating'; read(*,*)
   
   end do  ! iterating
   !print*, 'PAUSE: after iterations'; read(*,*)
  
   ! -----
   ! * update states...
   ! ------------------

   ! identify solution method
   select case(ixIterOption)

    ! * explicit Euler: update state vector based on the average of the start-of-step and end-of-step fluxes
    case(explicitEuler)
     stateVecTrial(:) = stateVecInit(:) + (0.5_dp*(fluxVec0(:) + fluxVecNew(:))*dtSubstep + 0.5_dp*(rAdd(:) + resSinkNew(:)) ) / real(sMul(:), dp)
     solutionError(:) = (fluxVec0(:)*dtSubstep + rAdd(:)) - (fluxVecNew(:)*dtSubstep + resSinkNew(:))
     write(*,'(a,1x,10(f20.10,1x))') 'fluxVec0      = ', fluxVec0(iJac1:iJac2)
     write(*,'(a,1x,10(f20.10,1x))') 'fluxVecNew    = ', fluxVecNew(iJac1:iJac2)
     write(*,'(a,1x,10(f20.10,1x))') 'rAdd          = ', rAdd(iJac1:iJac2)
     write(*,'(a,1x,10(f20.10,1x))') 'resSinkNew    = ', resSinkNew(iJac1:iJac2)
     write(*,'(a,1x,10(f20.10,1x))') 'stateVecInit  = ', stateVecInit(iJac1:iJac2)
     write(*,'(a,1x,10(f20.10,1x))') 'stateVecTrial = ', stateVecTrial(iJac1:iJac2)
     write(*,'(a,1x,10(f20.10,1x))') 'stateVecNew   = ', stateVecNew(iJac1:iJac2)
     write(*,'(a,1x,10(f20.10,1x))') 'solutionError = ', solutionError(iJac1:iJac2)
     print*, 'dt = ', dtSplit, dtSubstep, dtSum
     if(iSplit==2)then
      print*, 'PAUSE: checking state vector for the explicit Euler solution'; read(*,*)
     endif

     ! if the time step is rejected, then reduce the length of the substep and cycle
     errorMax = maxval( abs(solutionError) )
     if(errorMax(1) > errorTol)then
      dtSubstep = dtSubstep*max(safety*sqrt(errorTol/errorMax(1)), reduceMin)
      cycle substeps
     endif

     ! if the time step is accepted then check the need to increase the time step length
     dtSubstep = dtSubstep*min(safety*sqrt(errorTol/errorMax(1)), increaseMax)

    ! * implicit Euler: update state vector to be consistent with the fluxes 
    case(implicitEuler)

     ! update temperatures (ensure new temperature is consistent with the fluxes)
     if(nSnowSoilNrg>0)then
      do concurrent (iLayer=1:nLayers,ixSnowSoilNrg(iLayer)/=integerMissing)   ! (loop through non-missing energy state variables in the snow+soil domain)
       iState = ixSnowSoilNrg(iLayer)
       stateVecTrial(iState) = stateVecInit(iState) + (fluxVecNew(iState)*dtSubstep + resSinkNew(iState))/real(sMul(iState), dp)
      end do  ! looping through non-missing energy state variables in the snow+soil domain
     endif
   
     ! update volumetric water content in the snow (ensure change in state is consistent with the fluxes)
     ! NOTE: for soil water balance is constrained within the iteration loop
     if(nSnowSoilHyd>0)then
      do concurrent (iLayer=1:nSnow,ixSnowSoilHyd(iLayer)/=integerMissing)   ! (loop through non-missing energy state variables in the snow domain)
       iState = ixSnowSoilHyd(iLayer)
       stateVecTrial(iState) = stateVecInit(iState) + (fluxVecNew(iState)*dtSubstep + resSinkNew(iState))
      end do  ! looping through non-missing energy state variables in the snow+soil domain
     endif

    case default; err=20; message=trim(message)//'cannot find iterative option (expect explicitEuler or explicitEuler)'; return

   end select  ! iterative option
    
   ! -----
   ! * update model fluxes...
   ! ------------------------
 
   ! average start-of-step and end-of-step fluxes for explicit Euler
   if(ixIterOption==explicitEuler)then
    do iVar=1,size(flux_meta)
     if(fluxMask(iVar)) flux_temp%var(iVar)%dat(:) = 0.5_dp*(flux_init%var(iVar)%dat(:) + flux_temp%var(iVar)%dat(:) )
    end do
   endif

   ! increment fluxes
   dt_wght = dtSubstep/dt ! (define weight applied to each splitting operation) 
   do iVar=1,size(flux_meta)
    if(fluxMask(iVar)) flux_data%var(iVar)%dat(:) = flux_data%var(iVar)%dat(:) + flux_temp%var(iVar)%dat(:)*dt_wght
   end do

   ! -----
   ! * extract variables and update states...
   ! ----------------------------------------

   ! extract states from the state vector 
   call varExtract(&
                   ! input
                   stateVecTrial,            & ! intent(in):    model state vector (mixed units)
                   diag_data,                & ! intent(in):    model diagnostic variables for a local HRU
                   prog_data,                & ! intent(in):    model prognostic variables for a local HRU
                   indx_data,                & ! intent(in):    indices defining model states and layers
                   ! output: variables for the vegetation canopy
                   scalarCanairTempTrial,    & ! intent(out):   trial value of canopy air temperature (K)
                   scalarCanopyTempTrial,    & ! intent(out):   trial value of canopy temperature (K)
                   scalarCanopyWatTrial,     & ! intent(out):   trial value of canopy total water (kg m-2)
                   scalarCanopyLiqTrial,     & ! intent(out):   trial value of canopy liquid water (kg m-2)
                   scalarCanopyIceTrial,     & ! intent(out):   trial value of canopy ice content (kg m-2)
                   ! output: variables for the snow-soil domain
                   mLayerTempTrial,          & ! intent(inout): trial vector of layer temperature (K)
                   mLayerVolFracWatTrial,    & ! intent(out):   trial vector of volumetric total water content (-)
                   mLayerVolFracLiqTrial,    & ! intent(out):   trial vector of volumetric liquid water content (-)
                   mLayerVolFracIceTrial,    & ! intent(out):   trial vector of volumetric ice water content (-)
                   mLayerMatricHeadTrial,    & ! intent(out):   trial vector of total water matric potential (m)
                   mLayerMatricHeadLiqTrial, & ! intent(out):   trial vector of liquid water matric potential (m)
                   ! output: error control
                   err,cmessage)               ! intent(out):   error control
   if(err/=0)then; message=trim(message)//trim(cmessage); return; end if  ! (check for errors)
   !print*, 'PAUSE: after varExtract'; read(*,*)

   ! update diagnostic variables
   call updateVars(&
                   ! input
                   (iSplit==massSplit),                       & ! intent(in):    logical flag to adjust temperature to account for the energy used in melt+freeze
                   mpar_data,                                 & ! intent(in):    model parameters for a local HRU
                   indx_data,                                 & ! intent(in):    indices defining model states and layers
                   prog_data,                                 & ! intent(in):    model prognostic variables for a local HRU
                   diag_data,                                 & ! intent(inout): model diagnostic variables for a local HRU
                   deriv_data,                                & ! intent(inout): derivatives in model fluxes w.r.t. relevant state variables
                   ! output: variables for the vegetation canopy
                   scalarCanopyTempTrial,                     & ! intent(inout): trial value of canopy temperature (K)
                   scalarCanopyWatTrial,                      & ! intent(inout): trial value of canopy total water (kg m-2)
                   scalarCanopyLiqTrial,                      & ! intent(inout): trial value of canopy liquid water (kg m-2)
                   scalarCanopyIceTrial,                      & ! intent(inout): trial value of canopy ice content (kg m-2)
                   ! output: variables for the snow-soil domain
                   mLayerTempTrial,                           & ! intent(inout): trial vector of layer temperature (K)
                   mLayerVolFracWatTrial,                     & ! intent(inout): trial vector of volumetric total water content (-)
                   mLayerVolFracLiqTrial,                     & ! intent(inout): trial vector of volumetric liquid water content (-)
                   mLayerVolFracIceTrial,                     & ! intent(inout): trial vector of volumetric ice water content (-)
                   mLayerMatricHeadTrial,                     & ! intent(inout): trial vector of total water matric potential (m)
                   mLayerMatricHeadLiqTrial,                  & ! intent(inout): trial vector of liquid water matric potential (m)
                   ! output: error control
                   err,cmessage)                                ! intent(out):   error control
   if(err/=0)then; message=trim(message)//trim(cmessage); return; end if  ! (check for errors)
   !print*, 'PAUSE: after updateVars'; read(*,*)

   ! -----
   ! * extract state variables for the start of the next time step...
   ! ----------------------------------------------------------------
 
   ! build elements of the state vector for the vegetation canopy
   scalarCanairTemp    = scalarCanairTempTrial    ! trial value of canopy air temperature (K)
   scalarCanopyTemp    = scalarCanopyTempTrial    ! trial value of canopy temperature (K)
   scalarCanopyWat     = scalarCanopyWatTrial     ! trial value of canopy total water (kg m-2) kg m-2
   scalarCanopyLiq     = scalarCanopyLiqTrial     ! trial value of canopy liquid water (kg m-2)kg m-2
   scalarCanopyIce     = scalarCanopyIceTrial     ! trial value of canopy ice content (kg m-2) kg m-2

   ! build elements of the state vector for the snow+soil domain
   mLayerTemp          = mLayerTempTrial          ! trial vector of layer temperature (K)
   mLayerVolFracWat    = mLayerVolFracWatTrial    ! trial vector of volumetric total water content (-)
   mLayerVolFracLiq    = mLayerVolFracLiqTrial    ! trial vector of volumetric liquid water content (-)
   mLayerVolFracIce    = mLayerVolFracIceTrial    ! trial vector of volumetric ice water content (-)
   mLayerMatricHead    = mLayerMatricHeadTrial    ! trial vector of matric head (m)
   mLayerMatricHeadLiq = mLayerMatricHeadLiqTrial ! trial vector of matric head (m)

   ! ------------------------------------------------------
   ! ------------------------------------------------------

   ! increment sub-step
   dtSum = dtSum + dtSubstep

   ! check that we have completed the sub-step
   if(dtSum >= dt) exit substeps

   ! adjust length of the sub-step (make sure that we don't exceed the step)
   dtSubstep = min(dt - dtSum, dtSubstep)

  end do substeps  ! time steps for variable-dependent sub-stepping

  ! end associations with variables for the state update
  end associate stateSubset

  ! deallocate space for temporary vectors
  deallocate(stateVecInit, stateVecTrial, stateVecNew, fluxVec0, fluxVecNew, fScale, xScale, dMat, sMul, rVec, rAdd, resSinkNew, resVecNew, stat=err)
  if(err/=0)then; err=20; message=trim(message)//'unable to deallocate space for temporary vectors'; return; endif

  ! if explicit Euler, then deallocate space for the solution error
  if(ixIterOption==explicitEuler)then
   deallocate(solutionError, stat=err)
   if(err/=0)then; err=20; message=trim(message)//'unable to allocate space for the solution error'; return; endif
  endif

 end do  ! operator splitting loop 
 !print*, 'PAUSE: end of splitting loop'; read(*,*)

 ! ==========================================================================================================================================
 ! ==========================================================================================================================================
 ! ==========================================================================================================================================
 ! ==========================================================================================================================================

 ! -----
 ! * check mass balance...
 ! -----------------------

 ! point to flux variables in the data structure
 fluxVars: associate(&
 scalarRainfall            => flux_data%var(iLookFLUX%scalarRainfall)%dat(1)             ,& ! intent(in):  [dp]     rainfall rate (kg m-2 s-1)
 scalarThroughfallRain     => flux_data%var(iLookFLUX%scalarThroughfallRain)%dat(1)      ,& ! intent(out): [dp]     rain reaches ground without touching the canopy (kg m-2 s-1)
 scalarCanopyEvaporation   => flux_data%var(iLookFLUX%scalarCanopyEvaporation)%dat(1)    ,& ! intent(out): [dp]     canopy evaporation/condensation (kg m-2 s-1)
 scalarCanopyTranspiration => flux_data%var(iLookFLUX%scalarCanopyTranspiration)%dat(1)  ,& ! intent(out): [dp]     canopy transpiration (kg m-2 s-1)
 scalarCanopyLiqDrainage   => flux_data%var(iLookFLUX%scalarCanopyLiqDrainage)%dat(1)    ,& ! intent(out): [dp]     drainage liquid water from vegetation canopy (kg m-2 s-1)
 iLayerLiqFluxSoil         => flux_data%var(iLookFLUX%iLayerLiqFluxSoil)%dat             ,& ! intent(out): [dp(0:)] vertical liquid water flux at soil layer interfaces (-)
 mLayerTranspire           => flux_data%var(iLookFLUX%mLayerTranspire)%dat               ,& ! intent(out): [dp(:)]  transpiration loss from each soil layer (m s-1)
 mLayerBaseflow            => flux_data%var(iLookFLUX%mLayerBaseflow)%dat                ,& ! intent(out): [dp(:)]  baseflow from each soil layer (m s-1)
 scalarExfiltration        => flux_data%var(iLookFLUX%scalarExfiltration)%dat(1)         ,& ! intent(out): [dp]     exfiltration from the soil profile (m s-1)
 scalarCanopySublimation   => flux_data%var(iLookFLUX%scalarCanopySublimation)%dat(1)    ,& ! intent(out): [dp]     sublimation of ice from the vegetation canopy (kg m-2 s-1)
 scalarSnowSublimation     => flux_data%var(iLookFLUX%scalarSnowSublimation)%dat(1)       & ! intent(out): [dp]     sublimation of ice from the snow surface (kg m-2 s-1)
 ) ! associating flux variables in the data structure

 ! NOTE: This could be moved to the updateVars subroutine to avoid need to output Trial values

 ! check the mass balance
 ! NOTE: this should never fail since did not converge if water balance was not within tolerance=absConvTol_liquid
 if(checkMassBalance)then

  ! check mass balance for the canopy
  if(computeVegFlux)then
   canopyBalance1 = canopyBalance0 + (scalarRainfall + scalarCanopyEvaporation - scalarThroughfallRain - scalarCanopyLiqDrainage)*dtSplit
   liqError       = canopyBalance1 - scalarCanopyWatTrial
   if(abs(liqError) > absConvTol_liquid*10._dp)then  ! *10 to avoid precision issues
    write(*,'(a,1x,f30.20)')  'canopyBalance0          = ', canopyBalance0
    write(*,'(a,1x,f30.20)')  'canopyBalance1          = ', canopyBalance1
    write(*,'(a,1x,f30.20)')  'scalarCanopyWatTrial    = ', scalarCanopyWatTrial
    write(*,'(a,1x,f30.20)')  'scalarCanopyLiqTrial    = ', scalarCanopyLiqTrial
    write(*,'(a,1x,f30.20)')  'scalarCanopyIceTrial    = ', scalarCanopyIceTrial
    write(*,'(a,1x,f30.20)')  'scalarRainfall          = ', scalarRainfall*dtSplit 
    write(*,'(a,1x,f30.20)')  'scalarCanopyEvaporation = ', scalarCanopyEvaporation*dtSplit   
    write(*,'(a,1x,f30.20)')  'scalarThroughfallRain   = ', scalarThroughfallRain*dtSplit
    write(*,'(a,1x,f30.20)')  'scalarCanopyLiqDrainage = ', scalarCanopyLiqDrainage*dtSplit
    write(*,'(a,1x,f30.20)')  'liqError                = ', liqError
    message=trim(message)//'water balance error in the canopy domain'
    print*, trim(message)
    err=-20; return ! negative error code forces time step reduction and another trial
   endif  ! if there is a water balance error
  endif  ! if computing the water balance error

  ! check mass balance for soil
  soilBalance1 = sum( (mLayerVolFracLiqTrial(nSnow+1:nLayers) + mLayerVolFracIceTrial(nSnow+1:nLayers) )*mLayerDepth(nSnow+1:nLayers) )
  vertFlux     = -(iLayerLiqFluxSoil(nSoil) - iLayerLiqFluxSoil(0))*dt  ! m s-1 --> m
  tranSink     = sum(mLayerTranspire)*dt                                ! m s-1 --> m
  baseSink     = sum(mLayerBaseflow)*dt                                 ! m s-1 --> m
  compSink     = sum(mLayerCompress(1:nSoil) * mLayerDepth(nSnow+1:nLayers) ) ! dimensionless --> m
  liqError     = soilBalance1 - (soilBalance0 + vertFlux + tranSink - baseSink - compSink)
  if(abs(liqError) > absConvTol_liquid*10._dp)then  ! *10 to avoid precision issues
   write(*,'(a,1x,f30.20)')  'soilBalance0 = ', soilBalance0
   write(*,'(a,1x,f30.20)')  'soilBalance1 = ', soilBalance1
   write(*,'(a,1x,f30.20)')  'liqWaterSoil = ', sum( mLayerVolFracLiqTrial(nSnow+1:nLayers)*mLayerDepth(nSnow+1:nLayers) ) 
   write(*,'(a,1x,f30.20)')  'iceWaterSoil = ', sum( mLayerVolFracIceTrial(nSnow+1:nLayers)*mLayerDepth(nSnow+1:nLayers) ) 
   write(*,'(a,1x,f30.20)')  'vertFlux     = ', vertFlux
   write(*,'(a,1x,f30.20)')  'tranSink     = ', tranSink
   write(*,'(a,1x,f30.20)')  'baseSink     = ', baseSink
   write(*,'(a,1x,f30.20)')  'compSink     = ', compSink
   write(*,'(a,1x,f30.20)')  'liqError     = ', liqError
   message=trim(message)//'water balance error in the soil domain'
   print*, trim(message)
   stop
   err=-20; return ! negative error code forces time step reduction and another trial
  endif  ! if there is a water balance error

 endif  ! checking mass balance
 
 ! -----
 ! * check that there is sufficient ice content to support the converged sublimation rate...
 ! -----------------------------------------------------------------------------------------

 ! NOTE: This could be moved to the varExtract subroutine to avoid need to output Trial values

 ! check that sublimation does not exceed the available water on the canopy
 if(computeVegFlux)then
  if(-dt*scalarCanopySublimation > scalarCanopyLiqTrial + scalarCanopyIceTrial)then  ! try again
   print*, 'scalarCanopySublimation = ', scalarCanopySublimation
   print*, 'scalarCanopyLiqTrial    = ', scalarCanopyLiqTrial
   print*, 'scalarCanopyIceTrial    = ', scalarCanopyIceTrial
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
 
 ! end association to flux variables
 end associate fluxVars

 ! compute the melt in each snow and soil layer
 if(nSnow>0) mLayerMeltFreeze(      1:nSnow  ) = -(mLayerVolFracIce(      1:nSnow  ) - mLayerVolFracIceInit(      1:nSnow  ))*iden_ice
             mLayerMeltFreeze(nSnow+1:nLayers) = -(mLayerVolFracIce(nSnow+1:nLayers) - mLayerVolFracIceInit(nSnow+1:nLayers))*iden_water
   
 ! deallocate space for the baseflow derivative matrix
 deallocate(dBaseflow_dMatric,stat=err)
 if(err/=0)then; err=20; message=trim(message)//'unable to deallocate space for the baseflow derivatives'; return; end if

 ! end associate statements
 end associate globalVars

 end subroutine systemSolv

end module systemSolv_module

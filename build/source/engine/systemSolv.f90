module systemSolv_module
! data types
USE nrtype
! layer types
USE data_struc,only:ix_soil,ix_snow ! named variables for snow and soil
! access the number of snow and soil layers
USE data_struc,only:&
                    nSnow,        & ! number of snow layers  
                    nSoil,        & ! number of soil layers  
                    nLayers         ! total number of layers

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
! look-up values for the choice of groundwater representation (local-column, or single-basin)
USE mDecisions_module,only:  &
 localColumn,                & ! separate groundwater representation in each local soil column
 singleBasin                   ! single groundwater store over the entire basin
! look-up values for the choice of groundwater parameterization
USE mDecisions_module,only:  &
 qbaseTopmodel,              & ! TOPMODEL-ish baseflow parameterization
 bigBucket,                  & ! a big bucket (lumped aquifer model)
 noExplicit                    ! no explicit groundwater parameterization
! look-up values for the form of Richards' equation
USE mDecisions_module,only:  &
 moisture,                   & ! moisture-based form of Richards' equation
 mixdform                      ! mixed form of Richards' equation
! look-up values for the choice of boundary conditions for hydrology
USE mDecisions_module,only:  &
 prescribedHead,             & ! prescribed head (volumetric liquid water content for mixed form of Richards' eqn)
 funcBottomHead,             & ! function of matric head in the lower-most layer
 freeDrainage,               & ! free drainage
 liquidFlux,                 & ! liquid water flux
 zeroFlux                      ! zero flux
implicit none
private
public::systemSolv
! control parameters
real(dp),parameter  :: valueMissing=-9999._dp     ! missing value
real(dp),parameter  :: verySmall=1.e-6_dp         ! a very small number
real(dp),parameter  :: dx = 1.e-8_dp              ! finite difference increment
contains

 ! ************************************************************************************************
 ! new subroutine: run the coupled energy-mass model for one timestep
 ! ************************************************************************************************
 subroutine systemSolv(&
                       ! input: model control
                       dt,             & ! time step (s)
                       maxiter,        & ! maximum number of iterations
                       firstSubstep,   & ! flag to denote first sub-step
                       computeVegFlux, & ! flag to denote if computing energy flux over vegetation
                       ! input/output: data structures
                       type_data,      & ! intent(in):    type of vegetation and soil
                       attr_data,      & ! intent(in):    spatial attributes
                       forc_data,      & ! intent(in):    model forcing data
                       mpar_data,      & ! intent(in):    model parameters
                       mvar_data,      & ! intent(inout): model variables for a local HRU
                       bvar_data,      & ! intent(in):    model variables for the local basin
                       model_decisions,& ! intent(in):    model decisions
                       ! output: model control
                       niter,          & ! number of iterations taken
                       err,message)      ! error code and error message
 ! ---------------------------------------------------------------------------------------
 ! provide access to the derived types to define the data structures
 USE data_struc,only:&
                     var_i,            & ! data vector (i4b)
                     var_d,            & ! data vector (dp)
                     var_dlength,      & ! data vector with variable length dimension (dp)
                     model_options       ! defines the model decisions
 ! provide access to indices that define elements of the data structures
 USE var_lookup,only:iLookATTR,iLookTYPE,iLookPARAM,iLookFORCE,iLookMVAR,iLookINDEX ! named variables for structure elements
 USE var_lookup,only:iLookDECISIONS                                                 ! named variables for elements of the decision structure
 ! provide access to the numerical recipes modules
 USE nr_utility_module,only:arth                      ! creates a sequence of numbers (start, incr, n)
 ! provide access to utility modules
 USE soil_utils_module,only:volFracLiq                ! compute volumetric fraction of liquid water
 ! provide access to the flux modules
 USE vegnrgflux_module,only:vegnrgflux                ! compute energy fluxes over vegetation
 USE ssdnrgflux_module,only:ssdnrgflux                ! compute energy fluxes throughout the snow and soil subdomains
 USE vegliqflux_module,only:vegliqflux                ! compute liquid water fluxes through vegetation
 USE snowliqflx_module,only:snowliqflx                ! compute liquid water fluxes through snow 
 USE soilliqflx_module,only:soilliqflx                ! compute liquid water fluxes through soil
 USE groundwatr_module,only:satstorage                ! compute saturated storage in the soil profile
 USE groundwatr_module,only:soilbsflow                ! compute total baseflow from the soil profile
 USE groundwatr_module,only:disaggflow                ! diaggregate total baseflow to individual soil layers
 USE groundwatr_module,only:aquifrflux                ! compute liquid water fluxes in the aquifer
 USE snwDensify_module,only:snwDensify                ! snow densification
 ! provide access to the solver modules
 USE matrixSolv_module,only:matrixSolv                ! solve full matrix
 implicit none
 ! ---------------------------------------------------------------------------------------
 ! * dummy variables
 ! ---------------------------------------------------------------------------------------
 ! input: model control
 real(dp),intent(in)             :: dt                            ! time step (seconds)
 integer(i4b),intent(in)         :: maxiter                       ! maximum number of iterations
 logical(lgt),intent(in)         :: firstSubStep                  ! flag to indicate if we are processing the first sub-step
 logical(lgt),intent(in)         :: computeVegFlux                ! flag to indicate if we are computing fluxes over vegetation (.false. means veg is buried with snow)
 ! input/output: data structures
 type(var_i),intent(in)          :: type_data                     ! type of vegetation and soil
 type(var_d),intent(in)          :: attr_data                     ! spatial attributes
 type(var_d),intent(in)          :: forc_data                     ! model forcing data
 type(var_d),intent(in)          :: mpar_data                     ! model parameters
 type(var_dlength),intent(inout) :: mvar_data                     ! model variables for a local HRU
 type(var_dlength),intent(in)    :: bvar_data                     ! model variables for the local basin
 type(model_options),intent(in)  :: model_decisions(:)            ! model decisions
 ! output: model control
 integer(i4b),intent(out)        :: niter                         ! number of iterations
 integer(i4b),intent(out)        :: err                           ! error code
 character(*),intent(out)        :: message                       ! error message
 ! ---------------------------------------------------------------------------------------
 ! * variables in the data structures
 ! ---------------------------------------------------------------------------------------
 ! model decisions structure
 integer(i4b)                    :: ixRichards                   ! intent(in): choice of option for Richards eqn
 integer(i4b)                    :: ixGroundwater                ! intent(in): choice of groundwater parameterization
 integer(i4b)                    :: ixSpatialGroundwater         ! intent(in): spatial representation of groundwater (local-column or single-basin)
 ! domain boundary conditions
 real(dp)                        :: upperBoundTemp               ! intent(in): temperature of the upper boundary of the snow and soil domains (K)
 real(dp)                        :: scalarRainfall               ! intent(in): rainfall (kg m-2 s-1)
 real(dp)                        :: scalarSfcMeltPond            ! intent(in): ponded water caused by melt of the "snow without a layer" (kg m-2)
 ! diagnostic variables
 real(dp),dimension(nLayers)     :: mLayerDepth                  ! intent(in): depth of each layer in the snow-soil sub-domain (m)
 real(dp)                        :: scalarBulkVolHeatCapVeg      ! intent(in): bulk volumetric heat capacity of vegetation (J m-3 K-1)
 real(dp),dimension(nLayers)     :: mLayerVolHtCapBulk           ! intent(in): bulk volumetric heat capacity in each snow and soil layer (J m-3 K-1)
 ! vegetation parameters
 real(dp)                        :: heightCanopyTop              ! intent(in): height of the top of the vegetation canopy (m)
 real(dp)                        :: heightCanopyBottom           ! intent(in): height of the bottom of the vegetation canopy (m)
 ! soil parameters
 real(dp)                        :: vGn_alpha                    ! intent(in): van Genutchen "alpha" parameter (m-1)
 real(dp)                        :: vGn_n                        ! intent(in): van Genutchen "n" parameter (-)
 real(dp)                        :: vGn_m                        ! intent(in): van Genutchen "m" parameter (-)
 real(dp)                        :: theta_sat                    ! intent(in): soil porosity (-)
 real(dp)                        :: theta_res                    ! intent(in): soil residual volumetric water content (-)
 real(dp)                        :: specificStorage              ! intent(in): specific storage coefficient (m-1)
 ! model state variables (vegetation canopy)
 real(dp)                        :: scalarCanairTemp             ! intent(inout): temperature of the canopy air space (K)
 real(dp)                        :: scalarCanopyTemp             ! intent(inout): temperature of the vegetation canopy (K)
 real(dp)                        :: scalarCanopyIce              ! intent(inout): mass of ice on the vegetation canopy (kg m-2)
 real(dp)                        :: scalarCanopyLiq              ! intent(inout): mass of liquid water on the vegetation canopy (kg m-2)
 ! model state variables (snow and soil domains)
 real(dp),dimension(nLayers)     :: mLayerTemp                   ! intent(inout): temperature of each snow/soil layer (K)
 real(dp),dimension(nLayers)     :: mLayerVolFracIce             ! intent(inout): volumetric fraction of ice (-)
 real(dp),dimension(nLayers)     :: mLayerVolFracLiq             ! intent(inout): volumetric fraction of liquid water (-)
 real(dp),dimension(nLayers)     :: mLayerMatricHead             ! intent(inout): matric head (m)
 real(dp)                        :: scalarAquiferStorage         ! intent(inout): aquifer storage (m)
 ! *********************************************************************************************************************************************************
 ! *********************************************************************************************************************************************************
 ! ---------------------------------------------------------------------------------------
 ! * general local variables
 ! ---------------------------------------------------------------------------------------
 character(LEN=256)              :: cmessage                     ! error message of downwind routine
 real(dp)                        :: canopyDepth                  ! depth of the vegetation canopy (m)
 real(dp)                        :: scalarCanairTempNew          ! temperature of the canopy air space at the end of the sub-step (K)
 real(dp)                        :: scalarCanopyTempNew          ! temperature of the vegetation canopy at the end of the sub-step (K)
 real(dp)                        :: scalarCanopyIceNew           ! mass of ice on the vegetation canopy at the end of the sub-step (kg m-2)
 real(dp)                        :: scalarCanopyLiqNew           ! mass of liquid water on the vegetation canopy at the end of the sub-step (kg m-2)
 real(dp),allocatable            :: mLayerTempNew(:)             ! temperature of each snow/soil layer at the end of the sub-step (K)
 real(dp),allocatable            :: mLayerVolFracIceNew(:)       ! volumetric fraction of ice at the end of the sub-step (-)
 real(dp),allocatable            :: mLayerVolFracLiqNew(:)       ! volumetric fraction of liquid water at the end of the sub-step (-)
 real(dp),allocatable            :: mLayerMatricHeadNew(:)       ! matric head at the end of the sub-step (m)
 real(dp)                        :: scalarAquiferStorageNew      ! aquifer storage at the end of the sub-step (m)
 real(dp)                        :: scalarCanopyWater            ! total storage of wate rin the canopy; liquidwater plus ice (kg m-2)
 real(dp)                        :: volSub                       ! volumetric sublimation (kg m-3)
 real(dp)                        :: testSWE                      ! test SWE (kg m-2) 
 integer(i4b)                    :: iter                         ! iteration index
 integer(i4b)                    :: iLayer                       ! index of model layer
 integer(i4b)                    :: jLayer                       ! index of model layer within the full state vector
 integer(i4b)                    :: kLayer                       ! index of model layer within the snow-soil domain
 integer(i4b)                    :: local_ixGroundwater          ! local index for groundwater representation
 ! ------------------------------------------------------------------------------------------------------
 ! * trial state variables
 ! ------------------------------------------------------------------------------------------------------
 ! trial state variables (vegetation canopy)
 real(dp)                       :: scalarCanairTempTrial         ! trial value for temperature of the canopy air space (K)
 real(dp)                       :: scalarCanopyTempTrial         ! trial value for temperature of the vegetation canopy (K)
 real(dp)                       :: scalarCanopyIceTrial          ! trial value for mass of ice on the vegetation canopy (kg m-2)
 real(dp)                       :: scalarCanopyLiqTrial          ! trial value for mass of liquid water on the vegetation canopy (kg m-2)
 ! trial state variables (snow and soil domains)
 real(dp),dimension(nLayers)    :: mLayerTempTrial               ! trial value for temperature of each snow/soil layer (K)
 real(dp),dimension(nLayers)    :: mLayerVolFracIceTrial         ! trial value for volumetric fraction of ice (-)
 real(dp),dimension(nLayers)    :: mLayerVolFracLiqTrial         ! trial value for volumetric fraction of liquid water (-)
 real(dp),dimension(nSoil)      :: mLayerMatricHeadTrial         ! trial value for matric head (m)
 real(dp)                       :: scalarAquiferStorageTrial     ! trial value for aquifer storage (m)
 ! ------------------------------------------------------------------------------------------------------
 ! * model fluxes and derivatives
 ! ------------------------------------------------------------------------------------------------------
 ! energy fluxes and derivatives for the vegetation domain
 real(dp)                        :: canairNetNrgFlux             ! net energy flux for the canopy air space (W m-2)
 real(dp)                        :: canopyNetNrgFlux             ! net energy flux for the vegetation canopy (W m-2)
 real(dp)                        :: groundNetNrgFlux             ! net energy flux for the ground surface (W m-2)
 real(dp)                        :: dCanairNetFlux_dCanairTemp   ! derivative in net canopy air space flux w.r.t. canopy air temperature (W m-2 K-1)
 real(dp)                        :: dCanairNetFlux_dCanopyTemp   ! derivative in net canopy air space flux w.r.t. canopy temperature (W m-2 K-1)
 real(dp)                        :: dCanairNetFlux_dGroundTemp   ! derivative in net canopy air space flux w.r.t. ground temperature (W m-2 K-1)
 real(dp)                        :: dCanopyNetFlux_dCanairTemp   ! derivative in net canopy flux w.r.t. canopy air temperature (W m-2 K-1)
 real(dp)                        :: dCanopyNetFlux_dCanopyTemp   ! derivative in net canopy flux w.r.t. canopy temperature (W m-2 K-1)
 real(dp)                        :: dCanopyNetFlux_dGroundTemp   ! derivative in net canopy flux w.r.t. ground temperature (W m-2 K-1)
 real(dp)                        :: dGroundNetFlux_dCanairTemp   ! derivative in net ground flux w.r.t. canopy air temperature (W m-2 K-1)
 real(dp)                        :: dGroundNetFlux_dCanopyTemp   ! derivative in net ground flux w.r.t. canopy temperature (W m-2 K-1)
 real(dp)                        :: dGroundNetFlux_dGroundTemp   ! derivative in net ground flux w.r.t. ground temperature (W m-2 K-1)
 real(dp)                        :: dCanopyNetFlux_dCanLiq       ! derivative in net canopy fluxes w.r.t. canopy liquid water content (J kg-1 s-1)
 real(dp)                        :: dGroundNetFlux_dCanLiq       ! derivative in net ground fluxes w.r.t. canopy liquid water content (J kg-1 s-1)
 ! liquid water fluxes and derivatives associated with transpiration
 real(dp)                        :: scalarCanopyTranspiration    ! canopy transpiration (kg m-2 s-1)
 real(dp)                        :: scalarCanopyEvaporation      ! canopy evaporation/condensation (kg m-2 s-1)
 real(dp)                        :: scalarGroundEvaporation      ! ground evaporation/condensation -- below canopy or non-vegetated (kg m-2 s-1)
 real(dp)                        :: dCanopyEvaporation_dCanLiq   ! derivative in canopy evaporation w.r.t. canopy liquid water content (s-1)
 ! energy fluxes and derivatives for the snow and soil domains
 real(dp),dimension(nLayers)     :: ssdNetNrgFlux                ! net energy flux for each layer (J m-3 s-1)
 real(dp),dimension(0:nLayers)   :: iLayerNrgFlux                ! energy flux at the layer interfaces (W m-2)
 real(dp),dimension(0:nLayers)   :: dNrgFlux_dTempAbove          ! derivatives in the flux w.r.t. temperature in the layer above (J m-2 s-1 K-1)
 real(dp),dimension(0:nLayers)   :: dNrgFlux_dTempBelow          ! derivatives in the flux w.r.t. temperature in the layer below (J m-2 s-1 K-1)
 ! liquid water fluxes and derivatives for the vegetation domain
 real(dp)                        :: canopyNetLiqFlux             ! net liquid water flux for the vegetation canopy (kg m-2 s-1)
 real(dp)                        :: scalarThroughfallRain        ! rain that reaches the ground without ever touching the canopy (kg m-2 s-1)
 real(dp)                        :: scalarCanopyLiqDrainage      ! drainage of liquid water from the vegetation canopy (kg m-2 s-1)
 real(dp)                        :: scalarCanopyLiqDrainageDeriv ! derivative in canopy drainage w.r.t. canopy liquid water (s-1)
 real(dp)                        :: dCanopyEvaporation_dTCanair  ! derivative in canopy evaporation w.r.t. canopy air temperature (kg m-2 s-1 K-1)
 real(dp)                        :: dCanopyEvaporation_dTCanopy  ! derivative in canopy evaporation w.r.t. canopy temperature (kg m-2 s-1 K-1)
 real(dp)                        :: dCanopyEvaporation_dTGround  ! derivative in canopy evaporation w.r.t. ground temperature (kg m-2 s-1 K-1)
 ! liquid water fluxes and derivatives for the snow domain 
 real(dp),dimension(0:nSnow)     :: iLayerLiqFluxSnow            ! vertical liquid water flux at layer interfaces (m s-1)
 real(dp),dimension(0:nSnow)     :: iLayerLiqFluxSnowDeriv       ! derivative in vertical liquid water flux at layer interfaces (m s-1)
 real(dp)                        :: scalarRainPlusMelt           ! surface water input to the soil zone (m s-1)
 ! liquid water fluxes and derivatives for the soil domain
 real(dp),dimension(0:nSoil)     :: iLayerLiqFluxSoil            ! liquid flux at soil layer interfaces (m s-1)
 real(dp),dimension(nSoil)       :: mLayerTranspire              ! transpiration loss from each soil layer (m s-1)
 real(dp),dimension(nSoil)       :: mLayerBaseflow               ! baseflow from each soil layer -- only compute at the start of the step (m s-1)
 real(dp),dimension(0:nSoil)     :: dq_dStateAbove               ! change in the flux in layer interfaces w.r.t. state variables in the layer above
 real(dp),dimension(0:nSoil)     :: dq_dStateBelow               ! change in the flux in layer interfaces w.r.t. state variables in the layer below
 real(dp),dimension(nSoil)       :: mLayerHydCond                ! hydraulic conductivity in each soil layer (m s-1)
 real(dp),dimension(nSoil)       :: mLayerdTheta_dPsi            ! derivative in the soil water characteristic w.r.t. psi (m-1)
 real(dp),dimension(nSoil)       :: mLayerdPsi_dTheta            ! derivative in the soil water characteristic w.r.t. theta (m)
 real(dp),dimension(nSoil)       :: compress                     ! compressibility of soil (-)
 real(dp),dimension(nSoil)       :: dCompress_dPsi               ! derivative in compressibility w.r.t matric head (m-1)
 real(dp),dimension(nSnow)       :: snowNetLiqFlux               ! net liquid water flux for each snow layer (s-1)
 real(dp),dimension(nSoil)       :: soilNetLiqFlux               ! net liquid water flux for each soil layer (s-1)
 ! liquid water fluxes and derivatives for the aquifer
 real(dp)                        :: scalarAquiferTranspire       ! transpiration loss from the aquifer at the start-of-step (m s-1)
 real(dp)                        :: scalarAquiferRecharge        ! recharge to the aquifer (m s-1)
 real(dp)                        :: scalarAquiferBaseflow        ! total baseflow from the aquifer (m s-1)
 real(dp)                        :: scalarAquiferBaseflowDeriv   ! derivative in baseflow from the aquifer w.r.t. aquifer storage (s-1)
 ! ------------------------------------------------------------------------------------------------------
 ! * model indices
 ! ------------------------------------------------------------------------------------------------------
 integer(i4b),parameter          :: iJac1=1                      ! first layer of the Jacobian to print
 integer(i4b),parameter          :: iJac2=5                      ! last layer of the Jacobian to print
 integer(i4b),parameter          :: nVegNrg=2                    ! number of energy state variables for vegetation
 integer(i4b),parameter          :: nVegLiq=1                    ! number of hydrology state variables for vegetation
 integer(i4b)                    :: nVegState                    ! number of vegetation state variables (defines position of snow-soil states in the state vector)
 integer(i4b)                    :: nState                       ! total number of model state variables
 integer(i4b),parameter          :: ixCasNrg=1                   ! index of the canopy air space state variable
 integer(i4b),parameter          :: ixVegNrg=2                   ! index of the canopy energy state variable
 integer(i4b),parameter          :: ixVegLiq=3                   ! index of the canopy liquid water state variable
 integer(i4b)                    :: ixTopNrg                     ! index of the upper-most energy state variable in the snow-soil subdomain
 integer(i4b)                    :: ixTopLiq                     ! index of the upper-most liquid water state variable in the snow subdomain
 integer(i4b)                    :: ixTopMat                     ! index of the upper-most matric head state variable in the soil subdomain
 integer(i4b),dimension(nLayers) :: ixSnowSoilNrg                ! indices for energy state variables in the snow-soil subdomain
 integer(i4b),dimension(nLayers) :: ixSnowSoilLiq                ! indices for liquid warer state variables in the snow-soil subdomain
 integer(i4b),dimension(nSnow)   :: ixSnowOnlyLiq                ! indices for liquid water state variables in the snow subdomain
 integer(i4b),dimension(nSoil)   :: ixSoilOnlyMat                ! indices for matric head state variables in the soil subdomain
 integer(i4b),parameter          :: nVarSnowSoil=2               ! number of state variables in the snow and soil domain (energy and liquid water/matric head)
 integer(i4b),parameter          :: nRHS=1                       ! number of unknown variables on the RHS of the linear system A.X=B
 integer(i4b),parameter          :: ku=2                         ! number of super-diagonal bands
 integer(i4b),parameter          :: kl=2                         ! number of sub-diagonal bands
 integer(i4b),parameter          :: ixSup2=kl+ku-1               ! index for the 2nd super-diagonal band
 integer(i4b),parameter          :: ixSup1=kl+ku-0               ! index for the 1st super-diagonal band
 integer(i4b),parameter          :: ixDiag=kl+ku+1               ! index for the diagonal band
 integer(i4b),parameter          :: ixSub1=kl+ku+2               ! index for the 1st sub-diagonal band
 integer(i4b),parameter          :: ixSub2=kl+ku+3               ! index for the 2nd sub-diagonal band
 integer(i4b),parameter          :: nBands=2*kl+ku+1             ! length of tyhe leading dimension of the band diagonal matrix
 integer(i4b),parameter          :: ixFullMatrix=1001            ! named variable for the full Jacobian matrix
 integer(i4b),parameter          :: ixBandMatrix=1002            ! named variable for the band diagonal matrix
 integer(i4b)                    :: ixSolve=ixFullMatrix         ! option selected for the type of matrix used to solve the linear system A.X=B
 !integer(i4b)                    :: ixSolve=ixBandMatrix         ! option selected for the type of matrix used to solve the linear system A.X=B
 ! ------------------------------------------------------------------------------------------------------
 ! * model solver
 ! ------------------------------------------------------------------------------------------------------
 logical(lgt),parameter          :: numericalJacobian=.false.    ! flag to compute the Jacobian matrix
 real(dp),allocatable            :: stateVecInit(:)              ! initial state vector for the of the state equations (mixed units)
 real(dp),allocatable            :: stateVecTrial(:)             ! trial state vector for the of the state equations (mixed units)
 real(dp),allocatable            :: fluxVec0(:)                  ! flux vector (mixed units)
 real(dp),allocatable            :: fluxVec1(:)                  ! flux vector used in the numerical Jacobian calculations (mixed units)
 real(dp),allocatable            :: aJac(:,:)                    ! analytical Jacobian matrix
 real(dp),allocatable            :: nJac(:,:)                    ! numerical Jacobian matrix
 real(dp),allocatable            :: dMat(:)                      ! diagonal matrix (excludes flux derivatives)
 real(dp),allocatable            :: sMul(:)                      ! multiplier for state vector for the residual calculations
 real(dp),allocatable            :: rAdd(:)                      ! additional terms in the residual vector
 real(dp),allocatable            :: rVec(:)                      ! residual vector
 real(dp),allocatable            :: xInc(:)                      ! iteration increment
 real(dp),allocatable            :: rhs(:,:)                     ! the nState-by-nRHS matrix of matrix B, for the linear system A.X=B
 integer(i4b),allocatable        :: iPiv(:)                      ! defines if row i of the matrix was interchanged with row iPiv(i)
 real(dp),dimension(1)           :: energy_max                   ! maximum absolute value of the energy residual (J m-3)
 real(dp),dimension(1)           :: liquid_max                   ! maximum absolute value of the volumetric liquid water content residual (-)
 real(dp),parameter              :: absConvTol_energy=1.e-1_dp   ! convergence tolerance for energy (J m-3)
 real(dp),parameter              :: absConvTol_liquid=1.e-8_dp   ! convergence tolerance for volumetric liquid water content (-)
 ! ---------------------------------------------------------------------------------------
 ! point to variables in the data structures
 ! ---------------------------------------------------------------------------------------
 associate(&

 ! model decisions
 ixRichards              => model_decisions(iLookDECISIONS%f_Richards)%iDecision   ,&  ! intent(in): [i4b] index of the form of Richards' equation
 ixGroundwater           => model_decisions(iLookDECISIONS%groundwatr)%iDecision   ,&  ! intent(in): [i4b] groundwater parameterization
 ixSpatialGroundwater    => model_decisions(iLookDECISIONS%spatial_gw)%iDecision   ,&  ! intent(in): [i4b] spatial representation of groundwater (local-column or single-basin)

 ! domain boundary conditions
 upperBoundTemp          => forc_data%var(iLookFORCE%airtemp)                      ,&  ! intent(in): [dp] temperature of the upper boundary of the snow and soil domains (K)
 scalarRainfall          => mvar_data%var(iLookMVAR%scalarRainfall)%dat(1)         ,&  ! intent(in): [dp] rainfall rate (kg m-2 s-1)
 scalarSfcMeltPond       => mvar_data%var(iLookMVAR%scalarSfcMeltPond)%dat(1)      ,&  ! intent(in): [dp] ponded water caused by melt of the "snow without a layer" (kg m-2)

 ! diagnostic variables
 mLayerDepth             => mvar_data%var(iLookMVAR%mLayerDepth)%dat               ,&  ! intent(in): [dp(:)] depth of each layer in the snow-soil sub-domain (m)
 scalarBulkVolHeatCapVeg => mvar_data%var(iLookMVAR%scalarBulkVolHeatCapVeg)%dat(1),&  ! intent(in): [dp   ] bulk volumetric heat capacity of vegetation (J m-3 K-1)
 mLayerVolHtCapBulk      => mvar_data%var(iLookMVAR%mLayerVolHtCapBulk)%dat        ,&  ! intent(in): [dp(:)] bulk volumetric heat capacity in each snow and soil layer (J m-3 K-1)

 ! vegetation parameters
 heightCanopyTop         => mpar_data%var(iLookPARAM%heightCanopyTop)              ,&  ! intent(in): [dp] height of the top of the vegetation canopy (m)
 heightCanopyBottom      => mpar_data%var(iLookPARAM%heightCanopyBottom)           ,&  ! intent(in): [dp] height of the bottom of the vegetation canopy (m)

 ! soil parameters
 vGn_alpha               => mpar_data%var(iLookPARAM%vGn_alpha)                    ,&  ! intent(in): [dp] van Genutchen "alpha" parameter (m-1)
 vGn_n                   => mpar_data%var(iLookPARAM%vGn_n)                        ,&  ! intent(in): [dp] van Genutchen "n" parameter (-)
 vGn_m                   => mvar_data%var(iLookMVAR%scalarVGn_m)%dat(1)            ,&  ! intent(in): [dp] van Genutchen "m" parameter (-)
 theta_sat               => mpar_data%var(iLookPARAM%theta_sat)                    ,&  ! intent(in): [dp] soil porosity (-)
 theta_res               => mpar_data%var(iLookPARAM%theta_res)                    ,&  ! intent(in): [dp] soil residual volumetric water content (-)
 specificStorage         => mpar_data%var(iLookPARAM%specificStorage)              ,&  ! intent(in): [dp] specific storage coefficient (m-1)

 ! model state variables (vegetation canopy)
 scalarCanairTemp        => mvar_data%var(iLookMVAR%scalarCanairTemp)%dat(1)       ,&  ! intent(inout): [dp] temperature of the canopy air space (K)
 scalarCanopyTemp        => mvar_data%var(iLookMVAR%scalarCanopyTemp)%dat(1)       ,&  ! intent(inout): [dp] temperature of the vegetation canopy (K)
 scalarCanopyIce         => mvar_data%var(iLookMVAR%scalarCanopyIce)%dat(1)        ,&  ! intent(inout): [dp] mass of ice on the vegetation canopy (kg m-2)
 scalarCanopyLiq         => mvar_data%var(iLookMVAR%scalarCanopyLiq)%dat(1)        ,&  ! intent(inout): [dp] mass of liquid water on the vegetation canopy (kg m-2)

    ! model state variables (snow and soil domains)
 mLayerTemp              => mvar_data%var(iLookMVAR%mLayerTemp)%dat                ,&  ! intent(inout): [dp(:)] temperature of each snow/soil layer (K)
 mLayerVolFracIce        => mvar_data%var(iLookMVAR%mLayerVolFracIce)%dat          ,&  ! intent(inout): [dp(:)] volumetric fraction of ice (-)
 mLayerVolFracLiq        => mvar_data%var(iLookMVAR%mLayerVolFracLiq)%dat          ,&  ! intent(inout): [dp(:)] volumetric fraction of liquid water (-)
 mLayerMatricHead        => mvar_data%var(iLookMVAR%mLayerMatricHead)%dat          ,&  ! intent(inout): [dp(:)] matric head (m)
 scalarAquiferStorage    => mvar_data%var(iLookMVAR%scalarAquiferStorage)%dat(1)    &  ! intent(inout): [dp   ] aquifer storage (m)
 )
 ! ---------------------------------------------------------------------------------------
 ! initialize error control
 err=0; message="systemSolv/"

 ! *****
 ! (0) PRELIMINARIES...
 ! ********************

 ! modify the groundwater representation for this single-column implementation
 select case(ixSpatialGroundwater)
  case(singleBasin); local_ixGroundwater = noExplicit    ! force no explicit representation of groundwater at the local scale
  case(localColumn); local_ixGroundwater = ixGroundwater ! go with the specified decision
  case default; err=20; message=trim(message)//'unable to identify spatial representation of groundwater'; return
 end select ! (modify the groundwater representation for this single-column implementation) 

 ! define canopy depth (m)
 canopyDepth = heightCanopyTop - heightCanopyBottom

 ! get an initial canopy temperature if veg just starts protruding through snow on the ground
 if(computeVegFlux)then
  ! (NOTE: if canopy temperature is below absolute zero then canopy was previously buried by snow)
  if(scalarCanopyTemp < 0._dp .or. scalarCanairTemp < 0._dp)then
   ! check there is snow (there really has to be)
   if(nSnow == 0)then
    message=trim(message)//'no snow when canopy temperature or canopy air temperature is undefined -- canopy temps can only be undefined when buried with snow'
    err=20; return
   endif
   ! set canopy temperature to the temperature of the top snow layer + small offset to check derivative calculations
   scalarCanairTemp = mLayerTemp(1) + 0.1_dp
   scalarCanopyTemp = mLayerTemp(1) + 0.1_dp
  endif  ! (if canopy temperature undefined -- means canopy previously buried with snow)
 endif  ! (if computing vegetation fluxes -- canopy exposed)

 ! define the number of vegetation state variables (defines position of snow-soil states in the state vector)
 if(computeVegFlux)then
  nVegState = nVegNrg + nVegLiq
 else
  nVegState = 0
 endif

 ! define the number of model state variables
 nState = nVegState + nLayers*nVarSnowSoil   ! *nVarSnowSoil (both energy and liquid water)

 ! allocate space for the state vectors
 allocate(stateVecInit(nState),stateVecTrial(nState),stat=err)
 if(err/=0)then; err=20; message=trim(message)//'unable to allocate space for the state vector'; return; endif

 ! allocate space for the Jacobian matrix
 select case(ixSolve)
  case(ixFullMatrix); allocate(aJac(nState,nState),stat=err)
  case(ixBandMatrix); allocate(aJac(nBands,nState),stat=err)
  case default; err=20; message=trim(message)//'unable to identify option for the type of matrix'
 end select
 if(err/=0)then; err=20; message=trim(message)//'unable to allocate space for the Jacobian matrix'; return; endif

 ! allocate space for the flux vectors and Jacobian matrix
 allocate(dMat(nState),sMul(nState),rAdd(nState),fluxVec0(nState),rVec(nState),rhs(1,nState),iPiv(nState),xInc(nState),stat=err)
 if(err/=0)then; err=20; message=trim(message)//'unable to allocate space for the solution vectors'; return; endif

 ! define variables to calculate the numerical Jacobian matrix
 if(numericalJacobian)then
  ! (allocate space for the flux vector and Jacobian matrix
  allocate(fluxVec1(nState),nJac(nState,nState),stat=err)
  if(err/=0)then; err=20; message=trim(message)//'unable to allocate space for the flux vector and numerical Jacobian matrix'; return; endif
 endif  ! if calculating the numerical approximation of the Jacobian matrix

 ! define the index of the top layer
 ixTopNrg = nVegState + 1                       ! energy
 ixTopLiq = nVegState + 2                       ! liquid water (only snow)
 ixTopMat = nVegState + nSnow*nVarSnowSoil + 2  ! matric head (only soil)

 ! define the indices within the snow-soil domain
 ixSnowSoilNrg = arth(ixTopNrg,nVarSnowSoil,nLayers)  ! energy
 ixSnowSoilLiq = arth(ixTopLiq,nVarSnowSoil,nLayers)  ! liquid water

 ! define indices just for the snow and soil domains
 ixSoilOnlyMat = arth(ixTopMat,nVarSnowSoil,nSoil)    ! matric head
 if(nSnow>1)&  ! (liquid water in snow only defined if snow layers exist)
 ixSnowOnlyLiq = arth(ixTopLiq,nVarSnowSoil,nSnow)    ! liquid water
 print*, 'nVegState     = ', nVegState
 print*, 'ixSnowSoilNrg = ', ixSnowSoilNrg
 print*, 'ixSoilOnlyMat = ', ixSoilOnlyMat
 print*, 'ixSnowOnlyLiq = ', ixSnowOnlyLiq

 ! define additional vectors used in the residual calculations
 sMul(:) = 1._dp  ! multiplier for the state vector
 rAdd(:) = 0._dp  ! additional terms in the residual calculations (phase change, compressibility, etc.)

 ! define the multiplier for the state vector for residual calculations (vegetation canopy)
 if(computeVegFlux)then
  sMul(ixCasNrg) = Cp_air*iden_air          ! volumetric heat capacity of air (J m-3 K-1)
  sMul(ixVegNrg) = scalarBulkVolHeatCapVeg  ! volumetric heat capacity of the vegetation (J m-3 K-1)
  sMul(ixVegLiq) = 1._dp                    ! nothing else on the left hand side
 endif

 ! define the multiplier for the state vector for residual calculations (snow-soil domain)
 sMul(ixSnowSoilNrg) = mLayerVolHtCapBulk(1:nLayers)
 sMul(ixSnowSoilLiq) = 1._dp

 ! compute terms in the Jacobian for vegetation (excluding fluxes)
 ! NOTE: this is computed outside the iteration loop because it does not depend on state variables
 ! NOTE: energy for vegetation is computed *within* the iteration loop as it includes phase change
 if(computeVegFlux)then
  dMat(ixCasNrg) = Cp_air*iden_air          ! volumetric heat capacity of air (J m-3 K-1)
  dMat(ixVegLiq) = 1._dp                    ! nothing else on the left hand side
 endif

 ! compute terms in the Jacobian for the snow domain (excluding fluxes)
 ! NOTE: this is computed outside the iteration loop because it does not depend on state variables
 if(nSnow>0)&  ! (liquid water in snow only defined if snow layers exist)
 dMat(ixSnowOnlyLiq) = 1._dp

 ! initialize the analytical Jacobian matrix
 aJac(:,:) = 0._dp

 ! build the state vector for the vegetation canopy
 if(computeVegFlux)then
  stateVecInit(ixCasNrg) = scalarCanairTemp
  stateVecInit(ixVegNrg) = scalarCanopyTemp
  stateVecInit(ixVegLiq) = scalarCanopyLiq
 endif

 ! build the state vector for the snow and soil domain
 stateVecInit(ixSnowSoilNrg) = mLayerTemp(1:nLayers)
 stateVecInit(ixSoilOnlyMat) = mLayerMatricHead(1:nSoil)
 if(nSnow>0)&
 stateVecInit(ixSnowOnlyLiq) = mLayerVolFracLiq(1:nSnow)

 ! initialize the trial state vectors
 stateVecTrial = stateVecInit 

 ! ==========================================================================================================================================
 ! ==========================================================================================================================================
 ! ==========================================================================================================================================
 ! ==========================================================================================================================================

 ! (1) MAIN ITERATION LOOP...
 ! **************************

 ! iterate
 do iter=1,maxiter

  ! keep track of the number of iterations
  niter = iter

  ! test
  !write(*,'(a,1x,10(e15.5,1x))') 'stateVecInit(1:10)  = ', stateVecInit(1:10)
  !write(*,'(a,1x,10(e15.5,1x))') 'stateVecTrial(1:10) = ', stateVecTrial(1:10)

  ! -----
  ! * compute model fluxes...
  ! -------------------------

  ! compute model flux for a given state vector
  call computFlux(&
                  ! input:  full state vector
                  stateVecTrial,                      & ! intent(in): full state vector (mixed units)
                  ! output: state variables for the vegetation
                  scalarCanairTempTrial,              & ! intent(out): trial value for the temperature of the canopy air space (K)
                  scalarCanopyTempTrial,              & ! intent(out): trial value for the temperature of the vegetation canopy (K)
                  scalarCanopyLiqTrial,               & ! intent(out): trial value for the liquid water on the vegetation canopy (kg m-2)
                  ! output: state variables for the snow and soil domains
                  mLayerTempTrial,                    & ! intent(out): trial value for the temperature of each snow and soil layer (K)
                  mLayerVolFracLiqTrial,              & ! intent(out): trial value for the volumetric liquid water content in each snow and soil layer (-)
                  mLayerMatricHeadTrial,              & ! intent(out): trial value for the matyric head in each soil layer (m)
                  ! output: flux vector
                  fluxVec0,                           & ! intent(out): flux vector (mixed units)
                  ! output: error control
                  err,cmessage)                         ! intent(out): error code and error message
  if(err/=0)then; message=trim(message)//trim(cmessage); return; endif  ! (check for errors)
  !write(*,'(a,1x,10(e15.5,1x))') 'fluxVec0(1:10) = ', fluxVec0(1:10)

  ! compute soil compressibility (-) and its derivative w.r.t. matric head (m)
  ! NOTE: we already extracted trial matrix head and volumetric liquid water as part of the flux calculations
  call soilCmpres(&
                  ! input:
                  ixRichards,                             & ! intent(in): choice of option for Richards' equation
                  mLayerMatricHead(1:nSoil),              & ! intent(in): matric head at the start of the time step (m)
                  mLayerMatricHeadTrial(1:nSoil),         & ! intent(in): trial value of matric head (m)
                  mLayerVolFracLiqTrial(nSnow+1:nLayers), & ! intent(in): trial value for the volumetric liquid water content in each soil layer (-)
                  mLayerdTheta_dPsi,                      & ! intent(in): derivative in the soil water characteristic (m-1)
                  specificStorage,                        & ! intent(in): specific storage coefficient (m-1)
                  theta_sat,                              & ! intent(in): soil porosity (-)
                  ! output:
                  compress,                               & ! intent(out): compressibility of the soil matrix (-)
                  dCompress_dPsi,                         & ! intent(out): derivative in compressibility w.r.t. matric head (m-1)
                  err,cmessage)                             ! intent(out): error code and error message
  if(err/=0)then; message=trim(message)//trim(cmessage); return; endif  ! (check for errors)


  ! -----
  ! * compute Jacobian...
  ! ---------------------

  ! compute terms in the Jacobian for vegetation (excluding fluxes)
  ! NOTE: energy for vegetation is computed *within* the iteration loop as it includes phase change
  if(computeVegFlux)then
   dMat(ixVegNrg) = scalarBulkVolHeatCapVeg  ! volumetric heat capacity of the vegetation (J m-3 K-1)
  endif

  ! compute additional terms for the Jacobian for the snow-soil domain (excluding fluxes)
  ! NOTE: energy for vegetation is computed *within* the iteration loop as it includes phase change
  dMat(ixSnowSoilNrg) = mLayerVolHtCapBulk(1:nLayers)  

  ! compute additional terms for the Jacobian for the soil domain (excluding fluxes)
  if(ixRichards==moisture)then; err=20; message=trim(message)//'have not implemented the moisture-based form of RE yet'; return; endif
  dMat(ixSoilOnlyMat) = mLayerdTheta_dPsi(1:nSoil) + dCompress_dPsi(1:nSoil)

  ! compute the analytical Jacobian matrix
  select case(ixSolve)
   case(ixFullMatrix); call analJacob(err,cmessage)
   case(ixBandMatrix); call cpactBand(err,cmessage)
   case default; err=20; message=trim(message)//'unable to identify option for the type of matrix'
  end select
  if(err/=0)then; message=trim(message)//trim(cmessage); return; endif  ! (check for errors)
  !pause ' after analytical jacobian'

  ! *** testing: compute the numerical approximation of the Jacobian matrix
  if(numericalJacobian)then
   call numlJacob(stateVecTrial,fluxVec0,err,cmessage)
   if(err/=0)then; message=trim(message)//trim(cmessage); return; endif  ! (check for errors)
  endif  ! if computing the numerical Jacobian matrix

  ! -----
  ! * compute residual vector...
  ! ----------------------------

  ! compute additional terms on the RHS of the residual vector
  rAdd(ixSoilOnlyMat) = -compress(1:nSoil)

  ! compute the residual vector for the vegetation canopy
  ! NOTE: ignore phase change for now
  if(computeVegFlux)then
   rVec(ixCasNrg) = sMul(ixCasNrg)*scalarCanairTempTrial - ( (sMul(ixCasNrg)*scalarCanairTemp + fluxVec0(ixCasNrg)*dt) + rAdd(ixCasNrg) )
   rVec(ixVegNrg) = sMul(ixVegNrg)*scalarCanopyTempTrial - ( (sMul(ixVegNrg)*scalarCanopyTemp + fluxVec0(ixVegNrg)*dt) + rAdd(ixVegNrg) )
   rVec(ixVegLiq) =                scalarCanopyLiqTrial  - ( (               scalarCanopyLiq  + fluxVec0(ixVegLiq)*dt) + rAdd(ixVegLiq) )
  endif

  ! compute the residual vector for the snow and soil sub-domains for energy
  rVec(ixSnowSoilNrg) = sMul(ixSnowSoilNrg)*mLayerTempTrial(1:nLayers) - ( (sMul(ixSnowSoilNrg)*mLayerTemp(1:nLayers)  + fluxVec0(ixSnowSoilNrg)*dt) + rAdd(ixSnowSoilNrg) )

  ! compute the residual vector for the snow and soil sub-domains for liquid water
  ! NOTE: this is done regardless of the form of Richards' equation
  rVec(ixSnowSoilLiq) = mLayerVolFracLiqTrial(1:nLayers) - ( (mLayerVolFracLiq(1:nLayers)  + fluxVec0(ixSnowSoilLiq)*dt) + rAdd(ixSnowSoilLiq) )

  ! test
  !write(*,'(a,1x,10(e15.5,1x))') 'rVec(1:10) = ',     rVec(1:10)
  !write(*,'(a,1x,10(e15.5,1x))') 'rAdd(1:10) = ',     rAdd(1:10)
  !write(*,'(a,1x,10(e15.5,1x))') 'fluxVec0(1:10) = ', fluxVec0(1:10)

  !print*, '***'
  !write(*,'(a,1x,10(e15.5,1x))') 'mLayerVolFracLiqTrial(1:10) = ', mLayerVolFracLiqTrial(1:10)
  
  ! -----
  ! * solve linear system and update states...
  ! ------------------------------------------

  ! use the lapack routines to solve the linear system A.X=B
  call lapackSolv(aJac,rVec,xInc,err,cmessage)
  if(err/=0)then; message=trim(message)//trim(cmessage); return; endif  ! (check for errors)

  ! update the state vector (mixed units)
  stateVecTrial(:) = stateVecTrial(:) + xInc(:)

  ! -----
  ! * check convergence...
  ! ----------------------

  ! check convergence based on the residuals for energy (J m-3)
  if(computeVegFlux)then
   energy_max = maxval(abs( (/rVec(ixCasNrg), rVec(ixVegNrg), rVec(ixSnowSoilNrg)/) ) )
  else
   energy_max = maxval(abs( rVec(ixSnowSoilNrg) ) )
  endif

  ! check convergence based on the residuals for volumetric liquid water content (-)
  liquid_max = maxval(abs( rVec(ixSnowSoilLiq) ) )

  ! convergence check: 
  if( liquid_max(1) < absConvTol_liquid .and. energy_max(1) < absConvTol_energy) exit

  ! print progress towards solution
  print*, 'dt = ', dt
  write(*,'(a,1x,10(e15.5,1x))') 'liquid_max(1), energy_max(1) = ', liquid_max(1), energy_max(1)
  pause 'iterating'

  ! check convergence
  if(niter==maxiter)then; err=20; message=trim(message)//'failed to converge'; return; endif

 end do  ! iterating

 ! -----
 ! * extract state variables for the start of the next time step...
 ! ----------------------------------------------------------------

 ! extract the vegetation states from the state vector
 if(computeVegFlux)then
  scalarCanairTemp = stateVecTrial(ixCasNrg)
  scalarCanopyTemp = stateVecTrial(ixVegNrg)
  scalarCanopyLiq  = stateVecTrial(ixVegLiq)
 endif

 ! extract state variables for the snow and soil domain
 mLayerTemp(1:nLayers)     = stateVecTrial(ixSnowSoilNrg)
 mLayerMatricHead(1:nSoil) = stateVecTrial(ixSoilOnlyMat)
 if(nSnow>0)&  ! (liquid water in snow only defined if snow layers exist)
 mLayerVolFracLiq(1:nSnow) = stateVecTrial(ixSnowOnlyLiq)

 ! compute the volumetric liquid water content
 do iLayer=1,nSoil
  mLayerVolFracLiq(iLayer+nSnow) = volFracLiq(mLayerMatricHead(iLayer),vGn_alpha,theta_res,theta_sat,vGn_n,vGn_m)
 end do

 ! NOTE: still need to compute the volumetric fraction of ice

 ! ==========================================================================================================================================
 ! ==========================================================================================================================================
 ! ==========================================================================================================================================
 ! ==========================================================================================================================================

 ! deallocate space for the state vectors etc.
 deallocate(stateVecInit,stateVecTrial,dMat,sMul,rAdd,fluxVec0,aJac,rVec,rhs,iPiv,xInc,stat=err)
 if(err/=0)then; err=20; message=trim(message)//'unable to deallocate space for the state/flux vectors and analytical Jacobian matrix'; return; endif

 ! deallocate space for the variables used to create the numerical Jacobian matrix
 if(numericalJacobian)then
  deallocate(fluxVec1,nJac,stat=err)
  if(err/=0)then; err=20; message=trim(message)//'unable to allocate space for the flux vector and numerical Jacobian matrix'; return; endif
 endif

 ! end associate statement
 end associate

 contains

  ! ************************************************************************************************
  ! ************************************************************************************************
  ! ************************************************************************************************
  ! *** INTERNAL SUBROUTINES ***********************************************************************
  ! ************************************************************************************************
  ! ************************************************************************************************
  ! ************************************************************************************************

  ! ************************************************************************************************
  ! internal subroutine: compute model fluxes
  ! ************************************************************************************************
  subroutine computFlux(&
                        ! input:  full state vector
                        stateVec,                           & ! intent(in): full state vector (mixed units)
                        ! output: state variables for the vegetation
                        scalarCanairTempTrial,              & ! intent(out): trial value for the temperature of the canopy air space (K)
                        scalarCanopyTempTrial,              & ! intent(out): trial value for the temperature of the vegetation canopy (K)
                        scalarCanopyLiqTrial,               & ! intent(out): trial value for the liquid water on the vegetation canopy (kg m-2)
                        ! output: state variables for the snow and soil domains
                        mLayerTempTrial,                    & ! intent(out): trial value for the temperature of each snow and soil layer (K)
                        mLayerVolFracLiqTrial,              & ! intent(out): trial value for the volumetric liquid water content in each snow and soil layer (-)
                        mLayerMatricHeadTrial,              & ! intent(out): trial value for the matyric head in each soil layer (m)
                        ! output: flux vector
                        fluxVec,                            & ! intent(out): flux vector (mixed units)
                        ! output: error control
                        err,message)                          ! intent(out): error code and error message
  ! --------------------------------------------------------------
  implicit none
  ! input:  state vector
  real(dp),intent(in)            :: stateVec(:)               ! model state vector (mixed units)
  ! output: trial state variables (vegetation canopy)
  real(dp),intent(out)           :: scalarCanairTempTrial     ! trial value for temperature of the canopy air space (K)
  real(dp),intent(out)           :: scalarCanopyTempTrial     ! trial value for temperature of the vegetation canopy (K)
  real(dp),intent(out)           :: scalarCanopyLiqTrial      ! trial value for mass of liquid water on the vegetation canopy (kg m-2)
  ! output: trial state variables (snow and soil domains)
  real(dp),intent(out)           :: mLayerTempTrial(:)        ! trial value for temperature of each snow/soil layer (K)
  real(dp),intent(out)           :: mLayerVolFracLiqTrial(:)  ! trial value for volumetric fraction of liquid water (-)
  real(dp),intent(out)           :: mLayerMatricHeadTrial(:)  ! trial value for matric head (m)
  ! output: flux vector
  real(dp),intent(out)           :: fluxVec(:)                ! model flux vector (mixed units)
  ! output: error control
  integer(i4b),intent(out)       :: err                       ! error code
  character(*),intent(out)       :: message                   ! error message
  ! general local variables
  character(LEN=256)             :: cmessage                       ! error message of downwind routine
  ! --------------------------------------------------------------
  ! initialize error control
  err=0; message='computFlux/'

  ! *****
  ! (0) PRELIMINARIES...
  ! ********************

  ! get the necessary variables for the flux computations
  associate(&

  ! domain boundary conditions
  upperBoundTemp          => forc_data%var(iLookFORCE%airtemp)                      ,&  ! intent(in): [dp] temperature of the upper boundary of the snow and soil domains (K)
  scalarRainfall          => mvar_data%var(iLookMVAR%scalarRainfall)%dat(1)         ,&  ! intent(in): [dp] rainfall rate (kg m-2 s-1)
  scalarSfcMeltPond       => mvar_data%var(iLookMVAR%scalarSfcMeltPond)%dat(1)      ,&  ! intent(in): [dp] ponded water caused by melt of the "snow without a layer" (kg m-2)

  ! layer depth
  mLayerDepth             => mvar_data%var(iLookMVAR%mLayerDepth)%dat               ,&  ! intent(in): [dp(:)] depth of each layer in the snow-soil sub-domain (m)

  ! soil parameters
  vGn_alpha               => mpar_data%var(iLookPARAM%vGn_alpha)                    ,&  ! intent(in): [dp] van Genutchen "alpha" parameter (m-1)
  vGn_n                   => mpar_data%var(iLookPARAM%vGn_n)                        ,&  ! intent(in): [dp] van Genutchen "n" parameter (-)
  vGn_m                   => mvar_data%var(iLookMVAR%scalarVGn_m)%dat(1)            ,&  ! intent(in): [dp] van Genutchen "m" parameter (-)
  theta_sat               => mpar_data%var(iLookPARAM%theta_sat)                    ,&  ! intent(in): [dp] soil porosity (-)
  theta_res               => mpar_data%var(iLookPARAM%theta_res)                    ,&  ! intent(in): [dp] soil residual volumetric water content (-)
  specificStorage         => mpar_data%var(iLookPARAM%specificStorage)               &  ! intent(in): [dp] specific storage coefficient (m-1)
  )

  ! NOTE: Need to increase cleverness and avoid copying vectors

  ! extract the vegetation states from the state vector
  if(computeVegFlux)then
   scalarCanairTempTrial = stateVec(ixCasNrg)
   scalarCanopyTempTrial = stateVec(ixVegNrg)
   scalarCanopyLiqTrial  = stateVec(ixVegLiq)
  endif

  ! extract state variables for the snow and soil domain
  mLayerTempTrial(1:nLayers)     = stateVec(ixSnowSoilNrg)
  mLayerMatricHeadTrial(1:nSoil) = stateVec(ixSoilOnlyMat)
  if(nSnow>0)&  ! (liquid water in snow only defined if snow layers exist)
  mLayerVolFracLiqTrial(1:nSnow) = stateVec(ixSnowOnlyLiq)

  ! compute the volumetric liquid water content
  do iLayer=1,nSoil
   mLayerVolFracLiqTrial(iLayer+nSnow) = volFracLiq(mLayerMatricHeadTrial(iLayer),vGn_alpha,theta_res,theta_sat,vGn_n,vGn_m)
  end do

  ! compute the volumetric ice content
  ! *** NOTE: we still need to do this
  mLayerVolFracIceTrial(:) = 0._dp

  ! initialize liquid water fluxes throughout the snow and soil domains
  ! NOTE: used in the energy routines, which is called before the hydrology routines
  if(iter==1)then
   if(nSnow > 0)&
   iLayerLiqFluxSnow(0:nSnow) = 0._dp
   iLayerLiqFluxSoil(0:nSoil) = 0._dp
  endif

  ! *****
  ! (1) CALCULATE ENERGY FLUXES OVER VEGETATION...
  ! **********************************************
  call vegNrgFlux(&
                  ! input: model control
                  iter,                                   & ! intent(in): iteration index
                  firstSubStep,                           & ! intent(in): flag to indicate if we are processing the first sub-step
                  computeVegFlux,                         & ! intent(in): flag to indicate if we need to compute fluxes over vegetation
                  ! input: model state variables
                  upperBoundTemp,                         & ! intent(in): temperature of the upper boundary (K) --> NOTE: use air temperature
                  scalarCanairTempTrial,                  & ! intent(in): trial value of the canopy air space temperature (K)
                  scalarCanopyTempTrial,                  & ! intent(in): trial value of canopy temperature (K)
                  mLayerTempTrial(1),                     & ! intent(in): trial value of ground temperature (K)
                  scalarCanopyIceTrial,                   & ! intent(in): trial value of mass of ice on the vegetation canopy (kg m-2)
                  scalarCanopyLiqTrial,                   & ! intent(in): trial value of mass of liquid water on the vegetation canopy (kg m-2)
                  ! input/output: data structures
                  type_data,                              & ! intent(in):    type of vegetation and soil
                  attr_data,                              & ! intent(in):    spatial attributes
                  forc_data,                              & ! intent(in):    model forcing data
                  mpar_data,                              & ! intent(in):    model parameters
                  mvar_data,                              & ! intent(inout): model variables for a local HRU
                  bvar_data,                              & ! intent(in):    model variables for the local basin
                  model_decisions,                        & ! intent(in):    model decisions
                  ! output: liquid water fluxes associated with evaporation/transpiration
                  scalarCanopyTranspiration,              & ! intent(out): canopy transpiration (kg m-2 s-1)
                  scalarCanopyEvaporation,                & ! intent(out): canopy evaporation/condensation (kg m-2 s-1)
                  scalarGroundEvaporation,                & ! intent(out): ground evaporation/condensation -- below canopy or non-vegetated (kg m-2 s-1)
                  ! output: fluxes
                  canairNetNrgFlux,                       & ! intent(out): net energy flux for the canopy air space (W m-2)
                  canopyNetNrgFlux,                       & ! intent(out): net energy flux for the vegetation canopy (W m-2)
                  groundNetNrgFlux,                       & ! intent(out): net energy flux for the ground surface (W m-2)
                  ! output: flux derivatives
                  dCanairNetFlux_dCanairTemp,             & ! intent(out): derivative in net canopy air space flux w.r.t. canopy air temperature (W m-2 K-1)
                  dCanairNetFlux_dCanopyTemp,             & ! intent(out): derivative in net canopy air space flux w.r.t. canopy temperature (W m-2 K-1)
                  dCanairNetFlux_dGroundTemp,             & ! intent(out): derivative in net canopy air space flux w.r.t. ground temperature (W m-2 K-1)
                  dCanopyNetFlux_dCanairTemp,             & ! intent(out): derivative in net canopy flux w.r.t. canopy air temperature (W m-2 K-1)
                  dCanopyNetFlux_dCanopyTemp,             & ! intent(out): derivative in net canopy flux w.r.t. canopy temperature (W m-2 K-1)
                  dCanopyNetFlux_dGroundTemp,             & ! intent(out): derivative in net canopy flux w.r.t. ground temperature (W m-2 K-1)
                  dGroundNetFlux_dCanairTemp,             & ! intent(out): derivative in net ground flux w.r.t. canopy air temperature (W m-2 K-1)
                  dGroundNetFlux_dCanopyTemp,             & ! intent(out): derivative in net ground flux w.r.t. canopy temperature (W m-2 K-1)
                  dGroundNetFlux_dGroundTemp,             & ! intent(out): derivative in net ground flux w.r.t. ground temperature (W m-2 K-1)
                  ! output: liquid water flux derivarives
                  dCanopyEvaporation_dCanLiq,             & ! intent(out): derivative in canopy evaporation w.r.t. canopy liquid water content (s-1)
                  dCanopyEvaporation_dTCanair,            & ! intent(out): derivative in canopy evaporation w.r.t. canopy air temperature (kg m-2 s-1 K-1)
                  dCanopyEvaporation_dTCanopy,            & ! intent(out): derivative in canopy evaporation w.r.t. canopy temperature (kg m-2 s-1 K-1)
                  dCanopyEvaporation_dTGround,            & ! intent(out): derivative in canopy evaporation w.r.t. ground temperature (kg m-2 s-1 K-1)
                  ! output: cross derivative terms
                  dCanopyNetFlux_dCanLiq,                 & ! intent(out): derivative in net canopy fluxes w.r.t. canopy liquid water content (J kg-1 s-1)
                  dGroundNetFlux_dCanLiq,                 & ! intent(out): derivative in net ground fluxes w.r.t. canopy liquid water content (J kg-1 s-1)
                  ! output: error control
                  err,cmessage)                             ! intent(out): error control
  if(err/=0)then; message=trim(message)//trim(cmessage); return; endif  ! (check for errors)
  !print*, 'canairNetNrgFlux = ', canairNetNrgFlux
  !print*, 'canopyNetNrgFlux = ', canopyNetNrgFlux
  !print*, 'groundNetNrgFlux = ', groundNetNrgFlux
 
  !print*, 'in systemSolv: scalarGroundEvaporation = ', scalarGroundEvaporation
  !print*, 'in systemSolv: scalarCanopyEvaporation = ', scalarCanopyEvaporation
  !print*, 'in systemSolv: dCanopyEvaporation_dCanLiq = ', dCanopyEvaporation_dCanLiq


  ! *****
  ! (2) CALCULATE ENERGY FLUXES THROUGH THE SNOW-SOIL DOMAIN...
  ! ***********************************************************
  ! calculate energy fluxes at layer interfaces through the snow and soil domain
  call ssdNrgFlux(&
                  ! input: fluxes and derivatives at the upper boundary
                  groundNetNrgFlux,                       & ! intent(in): total flux at the ground surface (W m-2)
                  dGroundNetFlux_dGroundTemp,             & ! intent(in): derivative in total ground surface flux w.r.t. ground temperature (W m-2 K-1)
                  ! input: liquid water fluxes throughout the snow and soil domains
                  iLayerLiqFluxSnow,                      & ! intent(in): liquid flux at the interface of each snow layer (m s-1)
                  iLayerLiqFluxSoil,                      & ! intent(in): liquid flux at the interface of each soil layer (m s-1)
                  ! input: trial value of model state variabes
                  mLayerTempTrial,                        & ! intent(in): trial temperature at the current iteration (K)
                  ! output: fluxes and derivatives at all layer interfaces
                  iLayerNrgFlux,                          & ! intent(out): energy flux at the layer interfaces (W m-2)
                  dNrgFlux_dTempAbove,                    & ! intent(out): derivatives in the flux w.r.t. temperature in the layer above (W m-2 K-1)
                  dNrgFlux_dTempBelow,                    & ! intent(out): derivatives in the flux w.r.t. temperature in the layer below (W m-2 K-1)
                  ! output: error control
                  err,cmessage)                             ! intent(out): error control
  if(err/=0)then; message=trim(message)//trim(cmessage); return; endif
  ! calculate net energy fluxes for each snow and soil layer (J m-3 s-1)
  do iLayer=1,nLayers
   ssdNetNrgFlux(iLayer) = -(iLayerNrgFlux(iLayer) - iLayerNrgFlux(iLayer-1))/mLayerDepth(iLayer)
  end do
  print*, 'iLayerNrgFlux = ', iLayerNrgFlux
  print*, 'ssdNetNrgFlux = ', ssdNetNrgFlux

  ! *****
  ! (3) CALCULATE THE LIQUID FLUX THROUGH VEGETATION...
  ! ***************************************************
  call vegLiqFlux(&
                  ! input
                  computeVegFlux,                         & ! intent(in): flag to denote if computing energy flux over vegetation
                  scalarCanopyLiqTrial,                   & ! intent(in): trial mass of liquid water on the vegetation canopy at the current iteration (kg m-2)
                  scalarRainfall,                         & ! intent(in): rainfall rate (kg m-2 s-1)
                  ! output
                  scalarThroughfallRain,                  & ! intent(out): rain that reaches the ground without ever touching the canopy (kg m-2 s-1)
                  scalarCanopyLiqDrainage,                & ! intent(out): drainage of liquid water from the vegetation canopy (kg m-2 s-1)
                  scalarCanopyLiqDrainageDeriv,           & ! intent(out): derivative in canopy drainage w.r.t. canopy liquid water (s-1)
                  err,cmessage)                             ! intent(out): error control
  if(err/=0)then; message=trim(message)//trim(cmessage); return; endif
  ! calculate the net liquid water flux for the vegetation canopy
  canopyNetLiqFlux = scalarRainfall + scalarCanopyEvaporation - scalarThroughfallRain - scalarCanopyLiqDrainage
  !print*, 'scalarRainfall = ', scalarRainfall
  !print*, 'scalarCanopyEvaporation = ', scalarCanopyEvaporation
  !print*, 'scalarThroughfallRain = ', scalarThroughfallRain
  !print*, 'scalarCanopyLiqDrainage = ', scalarCanopyLiqDrainage

  ! *****
  ! (4) CALCULATE THE LIQUID FLUX THROUGH SNOW...
  ! *********************************************

  if(nSnow > 0)then
   ! compute liquid fluxes
   call snowLiqFlx(&
                   ! input: model control
                   iter,                                  & ! intent(in): iteration index
                   ! input: forcing for the snow domain
                   scalarThroughfallRain,                 & ! intent(in): rain that reaches the snow surface without ever touching vegetation (kg m-2 s-1)
                   scalarCanopyLiqDrainage,               & ! intent(in): liquid drainage from the vegetation canopy (kg m-2 s-1)
                   ! input: model state vector
                   mLayerVolFracLiqTrial(1:nSnow),        & ! intent(in): trial value of volumetric fraction of liquid water at the current iteration (-)
                   ! output: fluxes and derivatives
                   iLayerLiqFluxSnow(0:nSnow),            & ! intent(out): vertical liquid water flux at layer interfaces (m s-1)
                   iLayerLiqFluxSnowDeriv(0:nSnow),       & ! intent(out): derivative in vertical liquid water flux at layer interfaces (m s-1)
                   ! output: error control
                   err,cmessage)                            ! intent(out): error control
   if(err/=0)then; message=trim(message)//trim(cmessage); return; endif
   ! define forcing for the soil domain
   scalarRainPlusMelt = iLayerLiqFluxSnow(nSnow)    ! drainage from the base of the snowpack
   ! calculate net liquid water fluxes for each soil layer (s-1)
   do iLayer=1,nSnow
    snowNetLiqFlux(iLayer) = -(iLayerLiqFluxSnow(iLayer) - iLayerLiqFluxSnow(iLayer-1))/mLayerDepth(iLayer)
   end do
  else
   ! define forcing for the soil domain
   scalarRainPlusMelt = (scalarThroughfallRain + scalarCanopyLiqDrainage)/iden_water + &  ! liquid flux from the canopy (m s-1)
                         + (scalarSfcMeltPond/dt)/iden_water  ! melt of the snow without a layer (m s-1)
  endif

  ! *****
  ! (5) CALCULATE THE LIQUID FLUX THROUGH SOIL...
  ! *********************************************
  call soilLiqFlx(&
                  ! input: model control
                  iter,                                   & ! intent(in): iteration index
                  .true.,                                 & ! intent(in): flag indicating if derivatives are desired
                  ! input: trial state variables
                  mLayerMatricHeadTrial,                  & ! intent(in): matric head (m)
                  mLayerVolFracLiqTrial(nSnow+1:nLayers), & ! intent(in): volumetric fraction of liquid water (-)
                  mLayerVolFracIceTrial(nSnow+1:nLayers), & ! intent(in): volumetric fraction of ice (-)
                  ! input: fluxes
                  scalarCanopyTranspiration,              & ! intent(in): canopy transpiration (kg m-2 s-1)
                  scalarGroundEvaporation,                & ! intent(in): ground evaporation (kg m-2 s-1)
                  scalarRainPlusMelt,                     & ! intent(in): rain plus melt (m s-1)
                  ! output: diagnostic variables for model layers
                  mLayerdTheta_dPsi,                      & ! intent(out): derivative in the soil water characteristic w.r.t. psi (m-1)
                  mLayerdPsi_dTheta,                      & ! intent(out): derivative in the soil water characteristic w.r.t. theta (m)
                  ! output: fluxes
                  iLayerLiqFluxSoil,                      & ! intent(out): liquid fluxes at layer interfaces (m s-1)
                  mLayerTranspire,                        & ! intent(out): transpiration loss from each soil layer (m s-1)
                  mLayerHydCond,                          & ! intent(out): hydraulic conductivity in each layer (m s-1)
                  ! output: derivatives in fluxes w.r.t. state variables -- matric head or volumetric lquid water -- in the layer above and layer below (m s-1 or s-1)
                  dq_dStateAbove,                         & ! intent(out): derivatives in the flux w.r.t. volumetric liquid water content in the layer above (m s-1)
                  dq_dStateBelow,                         & ! intent(out): derivatives in the flux w.r.t. volumetric liquid water content in the layer below (m s-1)
                  ! output: error control
                  err,cmessage)                             ! intent(out): error control
  if(err/=0)then; message=trim(message)//trim(cmessage); return; endif
  ! calculate net liquid water fluxes for each soil layer (s-1)
  do iLayer=1,nSoil
   soilNetLiqFlux(iLayer) = -(iLayerLiqFluxSoil(iLayer) - iLayerLiqFluxSoil(iLayer-1))/mLayerDepth(iLayer+nSnow)
  end do
  !print*, 'iLayerLiqFluxSoil = ', iLayerLiqFluxSoil
  !print*, 'soilNetLiqFlux    = ', soilNetLiqFlux

  ! *****
  ! (X) WRAP UP...
  ! **************

  ! define model flux vector for the vegetation sub-domain
  if(computeVegFlux)then
   fluxVec(ixCasNrg) = canairNetNrgFlux/canopyDepth
   fluxVec(ixVegNrg) = canopyNetNrgFlux/canopyDepth
   fluxVec(ixVegLiq) = canopyNetLiqFlux
  endif

  ! define the model flux vector for the snow and soil sub-domains
  fluxVec(ixSnowSoilNrg) = ssdNetNrgFlux(1:nLayers)
  fluxVec(ixSoilOnlyMat) = soilNetLiqFlux(1:nSoil)
  if(nSnow>0)&
  fluxVec(ixSnowOnlyLiq) = snowNetLiqFlux(1:nSnow) 

  ! print progress
  !print*, '**'
  !write(*,'(a,1x,100(f15.9,1x))') 'stateVec(:) = ', stateVec(:)
  !write(*,'(a,1x,100(f15.9,1x))') 'fluxVec(:)  = ', fluxVec(:)

  ! end associate statement
  end associate

  end subroutine computFlux

  ! ************************************************************************************************
  ! internal subroutine: compute the compact band-diagonal matric
  ! ************************************************************************************************
  subroutine cpactBand(err,message)
  ! dummy variables
  integer(i4b),intent(out)       :: err                     ! error code
  character(*),intent(out)       :: message                 ! error message
  ! --------------------------------------------------------------
  ! initialize error control
  err=0; message='cpactBand/'

  ! associate variables from data structures
  associate(mLayerDepth => mvar_data%var(iLookMVAR%mLayerDepth)%dat) ! intent(in): [dp(:)] depth of each layer in the snow-soil sub-domain (m)

  ! compute elements that are only defined when vegetation protrudes over the snow surface
  if(computeVegFlux)then
   message=trim(message)//'band diagonal matrix not yet defined for vegetation'
   err=20; return
  endif

  ! energy fluxes for the snow-soil domain
  do iLayer=1,nLayers  ! loop through layers in the snow-soil domain
   ! (define layer indices)
   jLayer = ixSnowSoilNrg(iLayer)   ! layer index within the full state vector
   ! (define the compact band-diagonal matrix)
   if(iLayer > 1)       aJac(ixSup2,jLayer) = (dt/mLayerDepth(iLayer-1))*( dNrgFlux_dTempBelow(iLayer-1) )
                        aJac(ixSup1,jLayer) = 0._dp   ! cross derivative term: not yet computed
                        aJac(ixDiag,jLayer) = (dt/mLayerDepth(iLayer))*(-dNrgFlux_dTempBelow(iLayer-1) + dNrgFlux_dTempAbove(iLayer)) + dMat(jLayer)
                        aJac(ixSub1,jLayer) = 0._dp   ! cross derivative term: not yet computed
   if(iLayer < nLayers) aJac(ixSub2,jLayer) = (dt/mLayerDepth(iLayer+1))*(-dNrgFlux_dTempAbove(iLayer  ) )
  end do  ! (looping through layers in the snow-soil system)

  ! liquid water fluxes for the snow domain
  message=trim(message)//'have not yet defined liquid fluxes for the snow domain'
  err=20; return
 
   



  ! liquid water fluxes for the soil domain
  do iLayer=1,nSoil    ! loop through layers in the soil domain
   ! (define layer indices)
   jLayer = ixSoilOnlyMat(iLayer)  ! layer index within the full state vector
   kLayer = iLayer+nSnow           ! layer index within the full snow-soil vector
   ! (define the compact band-diagonal matrix)
   if(kLayer > nSnow+1) aJac(ixSup2,jLayer) = (dt/mLayerDepth(kLayer-1))*( dq_dStateBelow(iLayer-1))
                        aJac(ixSup1,jLayer) = 0._dp   ! cross derivative term: not yet computed
                        aJac(ixDiag,jLayer) = (dt/mLayerDepth(kLayer))*(-dq_dStateBelow(iLayer-1) + dq_dStateAbove(iLayer)) + dMat(jLayer)
                        aJac(ixSub1,jLayer) = 0._dp   ! cross derivative term: not yet computed
   if(kLayer < nLayers) aJac(ixSub2,jLayer) = (dt/mLayerDepth(kLayer+1))*(-dq_dStateAbove(iLayer))
  end do  ! (looping through soil layers)

  ! end association ro variables in the data structures
  end associate

  end subroutine cpactBand


  ! ************************************************************************************************
  ! internal subroutine: compute the Jacobian matrix (analytical)
  ! ************************************************************************************************
  subroutine analJacob(err,message)
  implicit none
  ! dummy variables
  integer(i4b),intent(out)       :: err                     ! error code
  character(*),intent(out)       :: message                 ! error message
  ! --------------------------------------------------------------
  ! initialize error control
  err=0; message='analJacob/'

  ! associate variables from data structures
  associate(mLayerDepth => mvar_data%var(iLookMVAR%mLayerDepth)%dat) ! intent(in): [dp(:)] depth of each layer in the snow-soil sub-domain (m)

  ! compute elements that are only defined when vegetation protrudes over the snow surface
  if(computeVegFlux)then

   ! liquid water fluxes for vegetation canopy (-)
   aJac(ixVegLiq,ixVegLiq) = -(dCanopyEvaporation_dCanLiq - scalarCanopyLiqDrainageDeriv)*dt + 1._dp

   ! cross-derivative terms w.r.t. system temperatures (kg m-2 K-1)
   aJac(ixVegLiq,ixCasNrg) = -dCanopyEvaporation_dTCanair*dt
   aJac(ixVegLiq,ixVegNrg) = -dCanopyEvaporation_dTCanopy*dt
   aJac(ixVegLiq,ixTopNrg) = -dCanopyEvaporation_dTGround*dt

   ! cross-derivative terms w.r.t. canopy liquid water (J m-1 kg-1)
   aJac(ixVegNrg,ixVegLiq) = (dt/canopyDepth)   *(-dCanopyNetFlux_dCanLiq)
   aJac(ixTopNrg,ixVegLiq) = (dt/mLayerDepth(1))*(-dGroundNetFlux_dCanLiq)

   ! energy fluxes with the canopy air space (J m-3 K-1)
   aJac(ixCasNrg,ixCasNrg) = (dt/canopyDepth)*(-dCanairNetFlux_dCanairTemp) + dMat(ixCasNrg)
   aJac(ixCasNrg,ixVegNrg) = (dt/canopyDepth)*(-dCanairNetFlux_dCanopyTemp)
   aJac(ixCasNrg,ixTopNrg) = (dt/canopyDepth)*(-dCanairNetFlux_dGroundTemp)
 
   ! energy fluxes with the vegetation canopy (J m-3 K-1)
   aJac(ixVegNrg,ixCasNrg) = (dt/canopyDepth)*(-dCanopyNetFlux_dCanairTemp)
   aJac(ixVegNrg,ixVegNrg) = (dt/canopyDepth)*(-dCanopyNetFlux_dCanopyTemp) + dMat(ixVegNrg)
   aJac(ixVegNrg,ixTopNrg) = (dt/canopyDepth)*(-dCanopyNetFlux_dGroundTemp)

   ! energy fluxes with the surface (J m-3 K-1)
   aJac(ixTopNrg,ixCasNrg) = (dt/mLayerDepth(1))*(-dGroundNetFlux_dCanairTemp)
   aJac(ixTopNrg,ixVegNrg) = (dt/mLayerDepth(1))*(-dGroundNetFlux_dCanopyTemp)

  endif  ! if there is a need to compute energy fluxes within vegetation

  ! energy fluxes for the snow-soil domain
  do iLayer=1,nLayers  ! loop through layers in the snow-soil domain
   ! - define layer indices
   jLayer = ixSnowSoilNrg(iLayer)
   ! - compute the Jacobian
   aJac(jLayer,jLayer)   = (dt/mLayerDepth(iLayer))*(-dNrgFlux_dTempBelow(iLayer-1) + dNrgFlux_dTempAbove(iLayer)) + dMat(jLayer)
   if(iLayer > 1)       aJac(jLayer-nVarSnowSoil,jLayer) = (dt/mLayerDepth(iLayer-1))*( dNrgFlux_dTempBelow(iLayer-1) )
   if(iLayer < nLayers) aJac(jLayer+nVarSnowSoil,jLayer) = (dt/mLayerDepth(iLayer+1))*(-dNrgFlux_dTempAbove(iLayer  ) )
  end do  ! (looping through layers in the snow-soil system)

  ! liquid water fluxes for the snow domain
  do iLayer=1,nSnow
   ! - define layer indices
   jLayer = ixSnowOnlyLiq(iLayer)   ! layer index within the full state vector
   ! - compute the Jacobian
   aJac(jLayer,jLayer) = (dt/mLayerDepth(iLayer))*iLayerLiqFluxSnowDeriv(iLayer) + dMat(jLayer)
   if(iLayer > 1)       aJac(jLayer-nVarSnowSoil,jLayer) = 0._dp  ! sub-diagonal: no dependence on other layers
   if(iLayer < nLayers) aJac(jLayer+nVarSnowSoil,jLayer) = 0._dp  ! super-diagonal: no dependence on other layers
  end do

  ! liquid water fluxes for the soil domain
  do iLayer=1,nSoil    ! loop through layers in the soil domain
   ! - define layer indices
   jLayer = ixSoilOnlyMat(iLayer)  ! layer index within the full state vector
   kLayer = iLayer+nSnow           ! layer index within the full snow-soil vector
   ! - compute the Jacobian
   aJac(jLayer,jLayer) = (dt/mLayerDepth(kLayer))*(-dq_dStateBelow(iLayer-1) + dq_dStateAbove(iLayer)) + dMat(jLayer)
   if(kLayer > nSnow+1) aJac(jLayer-nVarSnowSoil,jLayer) = (dt/mLayerDepth(kLayer-1))*( dq_dStateBelow(iLayer-1))
   if(kLayer < nLayers) aJac(jLayer+nVarSnowSoil,jLayer) = (dt/mLayerDepth(kLayer+1))*(-dq_dStateAbove(iLayer))
  end do  ! (looping through soil layers)

  ! print the Jacobian
  !print*, '** analytical Jacobian:'; do iLayer=iJac1,iJac2; write(*,'(i4,1x,100(e15.5,1x))') iLayer, aJac(iJac1:iJac2,iLayer); end do
  !pause

  ! end the association to data structures
  end associate

  end subroutine analJacob



  ! ************************************************************************************************
  ! internal subroutine: compute the Jacobian matrix (numerical)
  ! ************************************************************************************************
  subroutine numlJacob(stateVec,fluxVec,err,message)
  implicit none
  ! dummy
  real(dp),intent(in)            :: stateVec(:)             ! model state vector (mixed units)
  real(dp),intent(in)            :: fluxVec(:)              ! model flux vector (mixed units)
  integer(i4b),intent(out)       :: err                     ! error code
  character(*),intent(out)       :: message                 ! error message
  ! local
  real(dp),dimension(nState)     :: stateVecPerturbed       ! perturbed state vector
  real(dp),dimension(nState)     :: fluxVecJac              ! flux vector
  integer(i4b)                   :: iJac                    ! index of row of the Jacobian matrix
  integer(i4b),parameter         :: iTry=2                  ! index of trial model state variable (used for testing)
  ! trial state variables (vegetation canopy)
  real(dp)                       :: scalarCanairTempLocal   ! trial value for temperature of the canopy air space (K)
  real(dp)                       :: scalarCanopyTempLocal   ! trial value for temperature of the vegetation canopy (K)
  real(dp)                       :: scalarCanopyLiqLocal    ! trial value for mass of liquid water on the vegetation canopy (kg m-2)
  ! trial state variables (snow and soil domains)
  real(dp),dimension(nLayers)    :: mLayerTempLocal         ! trial value for temperature of each snow/soil layer (K)
  real(dp),dimension(nLayers)    :: mLayerVolFracLiqLocal   ! trial value for volumetric fraction of liquid water (-)
  real(dp),dimension(nSoil)      :: mLayerMatricHeadLocal   ! trial value for matric head (m)
  ! --------------------------------------------------------------
  ! initialize error control
  err=0; message='numlJacob/'

  ! get a copy of the state vector to perturb
  stateVecPerturbed(:) = stateVec(:)
 
  ! loop through state variables
  do iJac=1,nState

   ! (perturb state vector)
   stateVecPerturbed(iJac) = stateVec(iJac) + dx

   ! (compute fluxes)
   call computFlux(&
                   ! input:  full state vector
                   stateVecPerturbed,                  & ! intent(in): full state vector (mixed units)
                   ! output: state variables for the vegetation
                   scalarCanairTempLocal,              & ! intent(out): trial value for the temperature of the canopy air space (K)
                   scalarCanopyTempLocal,              & ! intent(out): trial value for the temperature of the vegetation canopy (K)
                   scalarCanopyLiqLocal,               & ! intent(out): trial value for the liquid water on the vegetation canopy (kg m-2)
                   ! output: state variables for the snow and soil domains
                   mLayerTempLocal,                    & ! intent(out): trial value for the temperature of each snow and soil layer (K)
                   mLayerVolFracLiqLocal,              & ! intent(out): trial value for the volumetric liquid water content in each snow and soil layer (-)
                   mLayerMatricHeadLocal,              & ! intent(out): trial value for the matyric head in each soil layer (m)
                   ! output: flux vector
                   fluxVecJac,                         & ! intent(out): flux vector (mixed units)
                   ! output: error control
                   err,cmessage)                         ! intent(out): error code and error message
   if(err/=0)then; message=trim(message)//trim(cmessage); return; endif  ! (check for errors)

   ! (compute the row of the Jacobian matrix)
   nJac(:,iJac) = -dt*(fluxVecJac(:) - fluxVec(:))/dx

   ! (add in the diagonal matrix)
   nJac(iJac,iJac) = nJac(iJac,iJac) + dMat(iJac)

   ! (test)
   !if(iJac<10) write(*,'(a,1x,10(e15.5,1x))') 'fluxVecJac(1:10) = ', fluxVecJac(1:10)
   !if(iJac==iTry) print*, 'iTry, stateVec(iTry), stateVecPerturbed(iTry) = ', iTry, stateVec(iTry), stateVecPerturbed(iTry) 
   !if(iJac==iTry) write(*,'(a,1x,i4,1x,10(f20.8,1x))'), 'iTry, -dt*(fluxVecJac(iTry) - fluxVec(iTry))/dx = ', iTry, -dt*(fluxVecJac(iTry) - fluxVec(iTry))/dx
   !if(iJac==iTry) pause ' in numerical Jacobian calculations'

   ! (set the state back to the input value)
   stateVecPerturbed(iJac) = stateVec(iJac)

  end do  ! (looping through state variables)

  ! print the Jacobian
  print*, '** numerical Jacobian:'; do iJac=iJac1,iJac2; write(*,'(i4,1x,100(e15.5,1x))') iJac, nJac(iJac1:iJac2,iJac); end do
  pause 'testing Jacobian'

  end subroutine numlJacob


  ! ************************************************************************************************
  ! internal subroutine: use the lapack routines to solve the linear system A.X=B
  ! ************************************************************************************************
  subroutine lapackSolv(aJac,rVec,xInc,err,message)
  implicit none
  ! dummy
  real(dp),intent(inout)         :: aJac(:,:)     ! input = the Jacobian matrix A; output = decomposed matrix
  real(dp),intent(in)            :: rVec(:)       ! the residual vector B
  real(dp),intent(out)           :: xInc(:)       ! the solution vector X
  integer(i4b),intent(out)       :: err           ! error code
  character(*),intent(out)       :: message       ! error message

  ! initialize error control
  select case(ixSolve)
   case(ixFullMatrix); message='lapackSolv/dgesv/'
   case(ixBandMatrix); message='lapackSolv/dgbsv/'
   case default; err=20; message=trim(message)//'unable to identify option for the type of matrix'
  end select

  ! form the rhs matrix
  rhs(1,1:nState) = -rVec(1:nState)

  ! identify option to solve the linear system A.X=B
  select case(ixSolve)

   ! lapack: use the full Jacobian matrix to solve the linear system A.X=B
   case(ixFullMatrix)
    call dgesv(nState, &  ! intent(in):    [i4b]               number of state variables
               nRHS,   &  ! intent(in):    [i4b]               number of columns of the matrix B
               aJac,   &  ! intent(inout): [dp(nState,nState)] input = the nState-by-nState Jacobian matrix A; output = decomposed matrix
               nState, &  ! intent(in):    [i4b]               the leading dimension of aJac
               iPiv,   &  ! intent(out):   [i4b(nState)]       defines if row i of the matrix was interchanged with row iPiv(i)
               rhs,    &  ! intent(inout): [dp(nState,nRHS)]   input = the nState-by-nRHS matrix of matrix B; output: the solution matrix X
               nState, &  ! intent(in):    [i4b]               the leading dimension of matrix rhs
               err)       ! intent(out)    [i4b]               error code

   ! lapack: use the band diagonal matrix to solve the linear system A.X=B
   case(ixBandMatrix)
    call dgbsv(nState, &  ! intent(in):    [i4b]               number of state variables
               kl,     &  ! intent(in):    [i4b]               number of subdiagonals within the band of A
               ku,     &  ! intent(in):    [i4b]               number of superdiagonals within the band of A
               nRHS,   &  ! intent(in):    [i4b]               number of columns of the matrix B
               aJac,   &  ! intent(inout): [dp(nBands,nState)] input = the nBands-by-nState Jacobian matrix A; output = decomposed matrix
               nBands, &  ! intent(in):    [i4b]               the leading dimension of aJac
               iPiv,   &  ! intent(out):   [i4b(nState)]       defines if row i of the matrix was interchanged with row iPiv(i)
               rhs,    &  ! intent(inout): [dp(nState,nRHS)]   input = the nState-by-nRHS matrix of matrix B; output: the solution matrix X
               nState, &  ! intent(in):    [i4b]               the leading dimension of matrix rhs
               err)       ! intent(out)    [i4b]               error code

   ! check that we found a valid option (should not get here because of the check above; included for completeness)
   case default; err=20; message=trim(message)//'unable to identify option for the type of matrix'

  end select  ! (option to solve the linear system A.X=B)

  ! identify any errors
  if(err/=0)then
   if(err<0)then
    write(message,'(a,i0,a)') trim(message)//'the ',err,'-th argument had an illegal value'
    err=abs(err); return
   else
    write(message,'(a,i0,a,i0,a)') trim(message)//'U(',err,',',err,') is exactly zero - factorization complete, but U is singular so the solution could not be completed'
    return
   endif
  endif

  ! extract the iteration increment
  xInc(1:nState) = rhs(1,1:nState)

  end subroutine lapackSolv


  ! ================================================================================================
  ! ================================================================================================


 end subroutine systemSolv


 ! ************************************************************************************************
 ! private subroutine: compute soil compressibility (-) and its derivative w.r.t matric head (m-1)
 ! ************************************************************************************************
 subroutine soilCmpres(&
                       ! input:
                       ixRichards,                         & ! intent(in): choice of option for Richards' equation
                       mLayerMatricHead,                   & ! intent(in): matric head at the start of the time step (m)
                       mLayerMatricHeadTrial,              & ! intent(in): trial value of matric head (m)
                       mLayerVolFracLiqTrial,              & ! intent(in): trial value for the volumetric liquid water content in each soil layer (-)
                       mLayerdTheta_dPsi,                  & ! intent(in): derivative in the soil water characteristic (m-1)
                       specificStorage,                    & ! intent(in): specific storage coefficient (m-1)
                       theta_sat,                          & ! intent(in): soil porosity (-)
                       ! output:
                       compress,                           & ! intent(out): compressibility of the soil matrix (-)
                       dCompress_dPsi,                     & ! intent(out): derivative in compressibility w.r.t. matric head (m-1)
                       err,message)                          ! intent(out): error code and error message
 implicit none
 ! input:
 integer(i4b),intent(in)        :: ixRichards                ! choice of option for Richards' equation
 real(dp),intent(in)            :: mLayerMatricHead(:)       ! matric head at the start of the time step (m)
 real(dp),intent(in)            :: mLayerMatricHeadTrial(:)  ! trial value for matric head (m)
 real(dp),intent(in)            :: mLayerVolFracLiqTrial(:)  ! trial value for volumetric fraction of liquid water (-)
 real(dp),intent(in)            :: mLayerdTheta_dPsi(:)      ! derivative in the soil water characteristic (m-1)
 real(dp),intent(in)            :: specificStorage           ! specific storage coefficient (m-1)
 real(dp),intent(in)            :: theta_sat                 ! soil porosity (-)
 ! output:
 real(dp),intent(out)           :: compress(:)               ! soil compressibility (-)
 real(dp),intent(out)           :: dCompress_dPsi(:)         ! derivative in soil compressibility w.r.t. matric head (m-1)
 integer(i4b),intent(out)       :: err                       ! error code
 character(*),intent(out)       :: message                   ! error message
 ! local variables
 character(LEN=256)             :: cmessage                  ! error message of downwind routine
 real(dp)                       :: fPart1,fPart2             ! different parts of the function
 real(dp)                       :: dPart1,dPart2             ! derivatives for different parts of the function
 integer(i4b)                   :: iLayer                    ! index of soil layer
 ! --------------------------------------------------------------
 ! initialize error control
 err=0; message='soilCmpres/'
 ! (only compute for the mixed form of Richards' equation)
 if(ixRichards==mixdform)then
  do iLayer=1,nSoil
   ! compute the compressibility term (-)
   compress(iLayer) = (specificStorage*mLayerVolFracLiqTrial(iLayer)/theta_sat) * (mLayerMatricHeadTrial(iLayer) - mLayerMatricHead(iLayer))
   ! compute the derivative for the compressibility term (m-1)
   fPart1 = specificStorage*(mLayerVolFracLiqTrial(iLayer)/theta_sat)  ! function for the 1st part (m-1)
   fPart2 = mLayerMatricHeadTrial(iLayer) - mLayerMatricHead(iLayer)   ! function for the 2nd part (m)
   dPart1 = mLayerdTheta_dPsi(iLayer)*specificStorage/theta_sat        ! derivative for the 1st part (m-2)
   dPart2 = 1._dp                                                      ! derivative for the 2nd part (-)
   dCompress_dPsi(iLayer) = fPart1*dPart2 + dPart1*fPart2              ! m-1
  end do
 else
  compress(:)       = 0._dp
  dCompress_dPsi(:) = 0._dp
 endif
 end subroutine soilCmpres





end module systemSolv_module

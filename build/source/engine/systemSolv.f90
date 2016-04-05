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
                       indx_data,      & ! intent(in):    index data
                       prog_data,      & ! intent(inout): model prognostic variables for a local HRU
                       diag_data,      & ! intent(inout): model diagnostic variables for a local HRU
                       flux_data,      & ! intent(inout): model fluxes for a local HRU
                       bvar_data,      & ! intent(in):    model variables for the local basin
                       model_decisions,& ! intent(in):    model decisions
                       ! output: model control
                       niter,          & ! number of iterations taken
                       err,message)      ! error code and error message
 ! ---------------------------------------------------------------------------------------
 ! provide access to the derived types to define the data structures
 USE data_types,only:&
                     var_i,            & ! data vector (i4b)
                     var_d,            & ! data vector (dp)
                     var_ilength,      & ! data vector with variable length dimension (i4b)
                     var_dlength,      & ! data vector with variable length dimension (dp)
                     model_options       ! defines the model decisions
 ! provide access to indices that define elements of the data structures
 USE var_lookup,only:iLookATTR           ! named variables for structure elements
 USE var_lookup,only:iLookTYPE           ! named variables for structure elements
 USE var_lookup,only:iLookPROG           ! named variables for structure elements
 USE var_lookup,only:iLookDIAG           ! named variables for structure elements
 USE var_lookup,only:iLookFLUX           ! named variables for structure elements
 USE var_lookup,only:iLookFORCE          ! named variables for structure elements
 USE var_lookup,only:iLookPARAM          ! named variables for structure elements
 USE var_lookup,only:iLookINDEX          ! named variables for structure elements
 USE var_lookup,only:iLookDECISIONS                                                 ! named variables for elements of the decision structure
 ! provide access to the numerical recipes utility modules
 USE nr_utility_module,only:arth                      ! creates a sequence of numbers (start, incr, n)
 ! structure allocations
 USE globalData,only:deriv_meta                       ! metadata on the model derivatives
 USE allocspace_module,only:allocLocal                ! allocate local data structures
 ! simulation of fluxes and residuals given a trial state vector
 USE eval8summa_module,only:eval8summa                ! simulation of fluxes and residuals given a trial state vector
 USE stateTrial_module,only:stateTrial                ! calculate the iteration increment, evaluate the new state, and refine if necessary
 ! provide access to utility modules
 USE soil_utils_module,only:volFracLiq                ! compute volumetric fraction of liquid water
 USE snow_utils_module,only:fracliquid                ! compute the fraction of liquid water at a given temperature (snow)
 USE snow_utils_module,only:templiquid                ! compute the temperature at a given fraction of liquid water (snow)
 USE snow_utils_module,only:dFracLiq_dTk              ! differentiate the freezing curve w.r.t. temperature (snow)
 USE soil_utils_module,only:dTheta_dPsi               ! derivative in the soil water characteristic (soil)
 USE soil_utils_module,only:dPsi_dTheta               ! derivative in the soil water characteristic (soil)
 USE soil_utils_module,only:dTheta_dTk                ! differentiate the freezing curve w.r.t. temperature (soil)
 USE soil_utils_module,only:matricHead                ! compute the matric head based on volumetric water content
 USE soil_utils_module,only:iceImpede                 ! compute the ice impedance factor (soil)
 USE soil_utils_module,only:dIceImpede_dTemp          ! differentiate the ice impedance factor w.r.t. temperature (soil)
 ! provide access to the flux modules
 USE getVectorz_module,only:popStateVec               ! populate the state vector
 USE getVectorz_module,only:varExtract                ! extract variables from the state vector
 USE updatState_module,only:updateSnow                ! update snow states
 USE updatState_module,only:updateSoil                ! update soil states
 USE vegnrgflux_module,only:vegnrgflux                ! compute energy fluxes over vegetation
 USE ssdnrgflux_module,only:ssdnrgflux                ! compute energy fluxes throughout the snow and soil subdomains
 USE vegliqflux_module,only:vegliqflux                ! compute liquid water fluxes through vegetation
 USE snowliqflx_module,only:snowliqflx                ! compute liquid water fluxes through snow
 USE soilliqflx_module,only:soilliqflx                ! compute liquid water fluxes through soil
 USE groundwatr_module,only:groundwatr                ! compute the baseflow flux
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
 type(var_ilength),intent(in)    :: indx_data                     ! indices for a local HRU
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
 real(dp)                        :: canopyDepth                  ! depth of the vegetation canopy (m)
 integer(i4b)                    :: iter                         ! iteration index
 integer(i4b)                    :: iLayer                       ! index of model layer
 integer(i4b)                    :: jLayer                       ! index of model layer within the full state vector (hydrology)
 integer(i4b)                    :: kLayer                       ! index of model layer within the snow-soil domain
 integer(i4b)                    :: mLayer                       ! index of model layer within the full state vector (thermodynamics)
 integer(i4b)                    :: local_ixGroundwater          ! local index for groundwater representation
 logical(lgt)                    :: printFlag                    ! flag to control printing (set to false for numerical jacobian)
 logical(lgt)                    :: printFlagInit                ! initialize flag to control printing
 logical(lgt)                    :: pauseProgress                ! flag to start looking at things more carefully
 logical(lgt)                    :: crosTempVeg                  ! flag to denoote where temperature crosses the freezing point
 real(dp)                        :: upperBoundTemp               ! temperature of the upper boundary of the snow and soil domains (K)
 real(dp),parameter              :: tempAccelerate=0.00_dp       ! factor to force initial canopy temperatures to be close to air temperature
 real(dp),parameter              :: xMinCanopyWater=0.0001_dp    ! minimum value to initialize canopy water (kg m-2)
 ! ------------------------------------------------------------------------------------------------------
 ! * model indices
 ! ------------------------------------------------------------------------------------------------------
 integer(i4b),parameter          :: nVarSnowSoil=2               ! number of state variables in the snow and soil domain (energy and total water/matric head)
 integer(i4b),parameter          :: nRHS=1                       ! number of unknown variables on the RHS of the linear system A.X=B
 integer(i4b),parameter          :: ku=3                         ! number of super-diagonal bands
 integer(i4b),parameter          :: kl=4                         ! number of sub-diagonal bands
 integer(i4b),parameter          :: ixSup3=kl+1                  ! index for the 3rd super-diagonal band
 integer(i4b),parameter          :: ixSup2=kl+2                  ! index for the 2nd super-diagonal band
 integer(i4b),parameter          :: ixSup1=kl+3                  ! index for the 1st super-diagonal band
 integer(i4b),parameter          :: ixDiag=kl+4                  ! index for the diagonal band
 integer(i4b),parameter          :: ixSub1=kl+5                  ! index for the 1st sub-diagonal band
 integer(i4b),parameter          :: ixSub2=kl+6                  ! index for the 2nd sub-diagonal band
 integer(i4b),parameter          :: ixSub3=kl+7                  ! index for the 3rd sub-diagonal band
 integer(i4b),parameter          :: ixSub4=kl+8                  ! index for the 3rd sub-diagonal band
 integer(i4b),parameter          :: nBands=2*kl+ku+1             ! length of the leading dimension of the band diagonal matrix
 integer(i4b),parameter          :: ixFullMatrix=1001            ! named variable for the full Jacobian matrix
 integer(i4b),parameter          :: ixBandMatrix=1002            ! named variable for the band diagonal matrix
 integer(i4b)                    :: ixMatrix                     ! the type of matrix used to solve the linear system A.X=B
 integer(i4b)                    :: ixSolve                      ! the type of matrix used to solve the linear system A.X=B
 integer(i4b),parameter          :: iJac1=1                      ! first layer of the Jacobian to print
 integer(i4b),parameter          :: iJac2=10                     ! last layer of the Jacobian to print
 ! ------------------------------------------------------------------------------------------------------
 ! * fluxes and derivatives
 ! ------------------------------------------------------------------------------------------------------
 ! ice content (need to keep track of this, but not part of the state vector)
 real(dp)                        :: theta                        ! liquid water equivalent of total water (liquid plus ice)
 real(dp),dimension(nLayers)     :: mLayerdTheta_dTk             ! derivative of volumetric liquid water content w.r.t. temperature (K-1)
 real(dp),dimension(nSoil)       :: dPsiLiq_dTemp                ! derivative in the liquid water matric potential w.r.t. temperature (m K-1)
 real(dp),dimension(nSnow)       :: fracLiqSnow                  ! fraction of liquid water in each snow layer (-)
 real(dp)                        :: fracLiqVeg                   ! fraction of liquid water on vegetation (-)
 real(dp)                        :: totalWaterVeg                ! total water on vegetation (kg m-2)
 real(dp)                        :: dTheta_dTkCanopy             ! derivative of volumetric liquid water content w.r.t. temperature (K-1)
 real(dp)                        :: dCanLiq_dTcanopy             ! derivative of canopy liquid storage w.r.t. temperature (kg m-2 K-1)
 ! state variables
 real(dp)                        :: scalarCanairTempTrial        ! trial value for temperature of the canopy air space (K)
 real(dp)                        :: scalarCanopyTempTrial        ! trial value for temperature of the vegetation canopy (K)
 real(dp)                        :: scalarCanopyWatTrial         ! trial value for liquid water storage in the canopy (kg m-2)
 real(dp),dimension(nLayers)     :: mLayerTempTrial              ! trial value for temperature of layers in the snow and soil domains (K)
 real(dp),dimension(nLayers)     :: mLayerVolFracWatTrial        ! trial value for volumetric fraction of total water (-)
 real(dp),dimension(nSoil)       :: mLayerMatricHeadTrial        ! trial value for matric head (m)
 ! diagnostic variables
 real(dp)                        :: scalarCanopyLiqTrial         ! trial value for mass of liquid water on the vegetation canopy (kg m-2)
 real(dp)                        :: scalarCanopyIceTrial         ! trial value for mass of ice on the vegetation canopy (kg m-2)
 real(dp),dimension(nLayers)     :: mLayerVolFracLiqTrial        ! trial value for volumetric fraction of liquid water (-)
 real(dp),dimension(nLayers)     :: mLayerVolFracIceTrial        ! trial value for volumetric fraction of ice (-)
 ! energy fluxes and derivatives for the vegetation domain

 type(var_dlength)               :: deriv_data                   ! derivatives in model fluxes w.r.t. relevant state variables 

 logical(lgt)                    :: feasible                     ! flag to define the feasibility of the solution

 real(dp)                        :: fluxVecNew(nState)           ! new flux vector
 real(qp)                        :: resVecNew(nState) ! NOTE: qp ! new residual vector


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
 real(dp)                        :: dGroundEvaporation_dCanLiq   ! derivative in ground evaporation w.r.t. canopy liquid water content (s-1)
 ! energy fluxes and derivatives for the snow and soil domains
 real(dp),dimension(nLayers)     :: ssdNetNrgFlux                ! net energy flux for each layer (J m-3 s-1)
 real(dp),dimension(0:nLayers)   :: dNrgFlux_dTempAbove          ! derivatives in the flux w.r.t. temperature in the layer above (J m-2 s-1 K-1)
 real(dp),dimension(0:nLayers)   :: dNrgFlux_dTempBelow          ! derivatives in the flux w.r.t. temperature in the layer below (J m-2 s-1 K-1)
 ! liquid water fluxes and derivatives for the vegetation domain
 real(dp)                        :: canopyNetLiqFlux             ! net liquid water flux for the vegetation canopy (kg m-2 s-1)
 real(dp)                        :: scalarCanopyLiqDeriv         ! derivative in (throughfall + canopy drainage) w.r.t. canopy liquid water (s-1)
 real(dp)                        :: scalarThroughfallRainDeriv   ! derivative in throughfall w.r.t. canopy liquid water (s-1)
 real(dp)                        :: scalarCanopyLiqDrainageDeriv ! derivative in canopy drainage w.r.t. canopy liquid water (s-1)
 real(dp)                        :: dCanopyEvaporation_dTCanair  ! derivative in canopy evaporation w.r.t. canopy air temperature (kg m-2 s-1 K-1)
 real(dp)                        :: dCanopyEvaporation_dTCanopy  ! derivative in canopy evaporation w.r.t. canopy temperature (kg m-2 s-1 K-1)
 real(dp)                        :: dCanopyEvaporation_dTGround  ! derivative in canopy evaporation w.r.t. ground temperature (kg m-2 s-1 K-1)
 real(dp)                        :: dGroundEvaporation_dTCanair  ! derivative in canopy evaporation w.r.t. canopy air temperature (kg m-2 s-1 K-1)
 real(dp)                        :: dGroundEvaporation_dTCanopy  ! derivative in canopy evaporation w.r.t. canopy temperature (kg m-2 s-1 K-1)
 real(dp)                        :: dGroundEvaporation_dTGround  ! derivative in canopy evaporation w.r.t. ground temperature (kg m-2 s-1 K-1)
 ! liquid water fluxes and derivatives for the snow domain
 real(dp),dimension(0:nSnow)     :: iLayerLiqFluxSnowDeriv       ! derivative in vertical liquid water flux at layer interfaces (m s-1)
 ! liquid water fluxes and derivatives for the soil domain
 real(dp)                        :: xMaxInfilRate                ! maximum infiltration rate (m s-1)
 real(dp)                        :: scalarInfilArea              ! fraction of unfrozen area where water can infiltrate (-)
 real(dp)                        :: scalarFrozenArea             ! fraction of area that is considered impermeable due to soil ice (-)
 real(dp)                        :: soilControl                  ! soil control on infiltration (-)
 real(dp)                        :: scalarSurfaceInfiltration    ! surface infiltration rate (m s-1) -- only computed for iter==1
 real(dp),dimension(nSoil)       :: mLayerTranspire              ! transpiration loss from each soil layer (m s-1)
 real(dp),dimension(nSoil)       :: dVolTot_dPsi0                ! derivative in total water content w.r.t. total water matric potential (m-1)
 real(dp),dimension(0:nSoil)     :: dq_dHydStateAbove            ! change in the flux in layer interfaces w.r.t. state variables in the layer above
 real(dp),dimension(0:nSoil)     :: dq_dHydStateBelow            ! change in the flux in layer interfaces w.r.t. state variables in the layer below
 real(dp),dimension(0:nSoil)     :: dq_dNrgStateAbove            ! change in the flux in layer interfaces w.r.t. state variables in the layer above
 real(dp),dimension(0:nSoil)     :: dq_dNrgStateBelow            ! change in the flux in layer interfaces w.r.t. state variables in the layer below
 real(dp),dimension(nSoil)       :: mLayerHydCond                ! hydraulic conductivity in each soil layer (m s-1)
 real(dp),dimension(nSoil)       :: dHydCond_dMatric             ! derivative in hydraulic conductivity w.r.t matric head (s-1)
 real(dp),dimension(nSoil)       :: mLayerdTheta_dPsi            ! derivative in the soil water characteristic w.r.t. psi (m-1)
 real(dp),dimension(nSoil)       :: mLayerdPsi_dTheta            ! derivative in the soil water characteristic w.r.t. theta (m)
 real(dp),dimension(nSoil)       :: dCompress_dPsi               ! derivative in compressibility w.r.t matric head (m-1)
 real(dp),dimension(nSnow)       :: snowNetLiqFlux               ! net liquid water flux for each snow layer (s-1)
 real(dp),dimension(nSoil)       :: soilNetLiqFlux               ! net liquid water flux for each soil layer (s-1)
 real(dp),allocatable            :: dBaseflow_dMatric(:,:)       ! derivative in baseflow w.r.t. matric head (s-1)  ! NOTE: allocatable, since not always needed
 integer(i4b)                    :: ixSaturation                 ! index of lowest saturated layer (NOTE: only computed on the first iteration)
 ! ------------------------------------------------------------------------------------------------------
 ! * model solver
 ! ------------------------------------------------------------------------------------------------------
 logical(lgt),parameter          :: numericalJacobian=.false.    ! flag to compute the Jacobian matrix
 logical(lgt),parameter          :: testBandDiagonal=.false.     ! flag to test the band-diagonal matrix
 logical(lgt),parameter          :: forceFullMatrix=.false.      ! flag to force the use of the full Jacobian matrix
 logical(lgt)                    :: firstFluxCall                ! flag to define the first flux call
 real(dp),dimension(nState)      :: stateVecInit                 ! initial state vector (mixed units)
 real(dp),dimension(nState)      :: stateVecTrial                ! trial state vector (mixed units)
 real(dp),dimension(nState)      :: stateVecNew                  ! new state vector (mixed units)
 real(dp),dimension(nState)      :: fluxVec0                     ! flux vector (mixed units)
 real(dp),dimension(nState)      :: fScale                       ! characteristic scale of the function evaluations (mixed units)
 real(dp),dimension(nState)      :: xScale                       ! characteristic scale of the state vector (mixed units)
 real(dp),dimension(nState)      :: dMat                         ! diagonal matrix (excludes flux derivatives)
 real(qp),dimension(nState)      :: sMul    ! NOTE: qp           ! multiplier for state vector for the residual calculations
 real(qp),dimension(nState)      :: rAdd    ! NOTE: qp           ! additional terms in the residual vector
 real(qp),dimension(nState)      :: rVec    ! NOTE: qp           ! residual vector
 real(dp),dimension(nState)      :: xInc                         ! iteration increment
 real(dp),dimension(nState)      :: grad                         ! gradient of the function vector = matmul(rVec,aJac)
 real(dp),dimension(nState,nRHS) :: rhs                          ! the nState-by-nRHS matrix of matrix B, for the linear system A.X=B
 integer(i4b),dimension(nState)  :: iPiv                         ! defines if row i of the matrix was interchanged with row iPiv(i)
 real(dp),allocatable            :: fluxVec1(:)                  ! flux vector used in the numerical Jacobian calculations (mixed units)
 real(dp),allocatable            :: aJac_test(:,:)               ! used to test the band-diagonal matrix structure
 real(dp),allocatable            :: aJac(:,:)                    ! analytical Jacobian matrix
 real(qp),allocatable            :: nJac(:,:)  ! NOTE: qp        ! numerical Jacobian matrix
 real(dp)                        :: fOld,fNew                    ! function values (-); NOTE: dimensionless because scaled
 real(dp)                        :: canopy_max                   ! absolute value of the residual in canopy water (kg m-2)
 real(dp),dimension(1)           :: energy_max                   ! maximum absolute value of the energy residual (J m-3)
 real(dp),dimension(1)           :: liquid_max                   ! maximum absolute value of the volumetric liquid water content residual (-)
 real(dp),dimension(1)           :: matric_max                   ! maximum absolute value of the matric head iteration increment (m)
 integer(i4b),dimension(1)       :: energy_loc                   ! location of maximum absolute value of the energy residual (-)
 integer(i4b),dimension(1)       :: liquid_loc                   ! location of maximum absolute value of the volumetric liquid water content residual (-)
 integer(i4b),dimension(1)       :: matric_loc                   ! location of maximum absolute value of the matric head increment (-)
 real(dp),parameter              :: absConvTol_energy=1.e-0_dp   ! convergence tolerance for energy (J m-3)
 real(dp),parameter              :: absConvTol_liquid=1.e-8_dp   ! convergence tolerance for volumetric liquid water content (-)
 real(dp),parameter              :: absConvTol_matric=1.e-3_dp   ! convergence tolerance for matric head increment in soil layers (m)
 real(dp),parameter              :: absConvTol_watbal=1.e-8_dp   ! convergence tolerance for soil water balance (m)
 real(dp),parameter              :: stepMax=1._dp                ! maximum step size (K, m, -)
 real(dp)                        :: stpmax                       ! scaled maximum step size
 real(dp),parameter              :: fScaleLiq=0.01_dp            ! func eval: characteristic scale for volumetric liquid water content (-)
 real(dp),parameter              :: fScaleMat=10._dp             ! func eval: characteristic scale for matric head (m)
 real(dp),parameter              :: fScaleNrg=1000000._dp          ! func eval: characteristic scale for energy (J m-3)
 real(dp),parameter              :: xScaleLiq=0.1_dp             ! state var: characteristic scale for volumetric liquid water content (-)
 real(dp),parameter              :: xScaleMat=10._dp             ! state var: characteristic scale for matric head (m)
 real(dp),parameter              :: xScaleTemp=1._dp             ! state var: characteristic scale for temperature (K)
 logical(lgt)                    :: converged                    ! convergence flag
 ! ------------------------------------------------------------------------------------------------------
 ! * solution constraints
 ! ------------------------------------------------------------------------------------------------------
 real(dp),dimension(nSnow)       :: mLayerTempCheck              ! updated temperatures (K) -- used to check iteration increment for snow
 real(dp),dimension(nSnow)       :: mLayerVolFracLiqCheck        ! updated volumetric liquid water content (-) -- used to check iteration increment for snow
 real(dp)                        :: cInc                         ! constrained temperature increment (K) -- simplified bi-section
 real(dp)                        :: xIncFactor                   ! scaling factor for the iteration increment (-)
 integer(i4b)                    :: iMax(1)                      ! index of maximum temperature
 logical(lgt),dimension(nSnow)   :: drainFlag                    ! flag to denote when drainage exceeds available capacity
 logical(lgt),dimension(nSoil)   :: crosFlag                     ! flag to denote temperature crossing from unfrozen to frozen (or vice-versa)
 integer(i4b)                    :: ixNrg,ixLiq                  ! index of energy and mass state variables in full state vector
 real(dp)                        :: xPsi00                       ! matric head after applying the iteration increment (m)
 real(dp)                        :: TcSoil                       ! critical point when soil begins to freeze (K)
 real(dp)                        :: critDiff                     ! temperature difference from critical (K)
 real(dp),parameter              :: epsT=1.e-7_dp               ! small interval above/below critical (K)
 ! ------------------------------------------------------------------------------------------------------
 ! * mass balance checks
 ! ------------------------------------------------------------------------------------------------------
 logical(lgt),parameter          :: checkMassBalance=.true.      ! flag to check the mass balance
 real(dp)                        :: soilWaterBalanceError        ! water balance error for soil
 real(dp)                        :: balance0,balance1            ! storage at start and end of time step
 real(dp)                        :: vertFlux                     ! change in storage due to vertical fluxes
 real(dp)                        :: tranSink,baseSink,compSink   ! change in storage sue to sink terms
 real(dp)                        :: liqError                     ! water balance error
 ! ---------------------------------------------------------------------------------------
 ! point to variables in the data structures
 ! ---------------------------------------------------------------------------------------
 associate(&

 ! model decisions
 ixRichards              => model_decisions(iLookDECISIONS%f_Richards)%iDecision   ,&  ! intent(in): [i4b] index of the form of Richards' equation
 ixGroundwater           => model_decisions(iLookDECISIONS%groundwatr)%iDecision   ,&  ! intent(in): [i4b] groundwater parameterization
 ixSpatialGroundwater    => model_decisions(iLookDECISIONS%spatial_gw)%iDecision   ,&  ! intent(in): [i4b] spatial representation of groundwater (local-column or single-basin)

 ! domain boundary conditions
 airtemp                 => forc_data%var(iLookFORCE%airtemp)                      ,&  ! intent(in): [dp] temperature of the upper boundary of the snow and soil domains (K)
 scalarRainfall          => flux_data%var(iLookFLUX%scalarRainfall)%dat(1)         ,&  ! intent(in): [dp] rainfall rate (kg m-2 s-1)
 scalarSfcMeltPond       => prog_data%var(iLookPROG%scalarSfcMeltPond)%dat(1)      ,&  ! intent(in): [dp] ponded water caused by melt of the "snow without a layer" (kg m-2)

 ! indices of model state variables
 ixCasNrg                => indx_data%var(iLookINDEX%ixCasNrg)%dat(1)             , & ! index of canopy air space energy state variable
 ixVegNrg                => indx_data%var(iLookINDEX%ixVegNrg)%dat(1)             , & ! index of canopy energy state variable
 ixVegWat                => indx_data%var(iLookINDEX%ixVegWat)%dat(1)             , & ! index of canopy hydrology state variable (mass)
 ixTopNrg                => indx_data%var(iLookINDEX%ixTopNrg)%dat(1)             , & ! index of upper-most energy state in the snow-soil subdomain
 ixTopWat                => indx_data%var(iLookINDEX%ixTopWat)%dat(1)             , & ! index of upper-most total water state in the snow-soil subdomain
 ixTopMat                => indx_data%var(iLookINDEX%ixTopMat)%dat(1)             , & ! index of upper-most matric head state in the soil subdomain
 ixSnowSoilNrg           => indx_data%var(iLookINDEX%ixSnowSoilNrg)%dat           , & ! indices for energy states in the snow-soil subdomain
 ixSnowSoilWat           => indx_data%var(iLookINDEX%ixSnowSoilWat)%dat           , & ! indices for total water states in the snow-soil subdomain
 ixSnowOnlyNrg           => indx_data%var(iLookINDEX%ixSnowOnlyNrg)%dat           , & ! indices for energy states in the snow subdomain
 ixSnowOnlyWat           => indx_data%var(iLookINDEX%ixSnowOnlyWat)%dat           , & ! indices for total water states in the snow subdomain
 ixSoilOnlyNrg           => indx_data%var(iLookINDEX%ixSoilOnlyNrg)%dat           , & ! indices for energy states in the soil subdomain
 ixSoilOnlyHyd           => indx_data%var(iLookINDEX%ixSoilOnlyHyd)%dat           , & ! indices for hydrology states in the soil subdomain
 layerType               => indx_data%var(iLookINDEX%layerType)%dat               , & ! index defining type of layer (soil or snow)

 ! type of model state variables
 ixSoilState             => indx_data%var(iLookINDEX%ixSoilState)%dat             , & ! list of indices for all soil layers
 ixLayerState            => indx_data%var(iLookINDEX%ixLayerState)%dat            , & ! list of indices for all model layers
 ixStateType             => indx_data%var(iLookINDEX%ixStateType)%dat             , & ! indices defining the type of the state (ixNrgState...)
 ixAllState              => indx_data%var(iLookINDEX%ixAllState)%dat              , & ! list of indices for all model state variables
 ixNrgOnly               => indx_data%var(iLookINDEX%ixNrgOnly)%dat               , & ! list of indices for all energy states
 ixWatOnly               => indx_data%var(iLookINDEX%ixWatOnly)%dat               , & ! list of indices for all "total water" states
 ixMatOnly               => indx_data%var(iLookINDEX%ixMatOnly)%dat               , & ! list of indices for matric head state variables
 ixMassOnly              => indx_data%var(iLookINDEX%ixMassOnly)%dat              , & ! list of indices for hydrology states (mass of water)

 ! diagnostic variables
 mLayerDepth             => prog_data%var(iLookPROG%mLayerDepth)%dat               ,&  ! intent(in): [dp(:)] depth of each layer in the snow-soil sub-domain (m)
 scalarBulkVolHeatCapVeg => diag_data%var(iLookDIAG%scalarBulkVolHeatCapVeg)%dat(1),&  ! intent(in): [dp   ] bulk volumetric heat capacity of vegetation (J m-3 K-1)
 mLayerVolHtCapBulk      => diag_data%var(iLookDIAG%mLayerVolHtCapBulk)%dat        ,&  ! intent(in): [dp(:)] bulk volumetric heat capacity in each snow and soil layer (J m-3 K-1)
 mLayerMeltFreeze        => diag_data%var(iLookDIAG%mLayerMeltFreeze)%dat          ,&  ! intent(out): [dp(:)] melt-freeze in each snow and soil layer (kg m-3)
 mLayerThetaResid        => diag_data%var(iLookDIAG%mLayerThetaResid)%dat          ,&  ! intent(out): [dp(:)] residual volumetric liquid water content in each snow layer (-)

 ! model fluxes
 mLayerCompress          => diag_data%var(iLookDIAG%mLayerCompress)%dat            ,&  ! intent(out): [dp(:)] change in storage associated with compression of the soil matrix (-)
 iLayerLiqFluxSnow       => flux_data%var(iLookFLUX%iLayerLiqFluxSnow)%dat         ,&  ! intent(out): [dp(0:)] vertical liquid water flux at snow layer interfaces (-)
 iLayerLiqFluxSoil       => flux_data%var(iLookFLUX%iLayerLiqFluxSoil)%dat         ,&  ! intent(out): [dp(0:)] vertical liquid water flux at soil layer interfaces (-)
 mLayerBaseflow          => flux_data%var(iLookFLUX%mLayerBaseflow)%dat            ,&  ! intent(out): [dp(:)] baseflow from each soil layer (m s-1)
 scalarSoilBaseflow      => flux_data%var(iLookFLUX%scalarSoilBaseflow)%dat(1)     ,&  ! intent(out): [dp] total baseflow from the soil profile (m s-1)
 scalarExfiltration      => flux_data%var(iLookFLUX%scalarExfiltration)%dat(1)     ,& ! intent(out):[dp]    exfiltration from the soil profile (m s-1)

 ! sublimation (needed to check mass balance constraints)
 scalarCanopySublimation => flux_data%var(iLookFLUX%scalarCanopySublimation)%dat(1),&  ! intent(out): [dp] sublimation of ice from the vegetation canopy (kg m-2 s-1)
 scalarSnowSublimation   => flux_data%var(iLookFLUX%scalarSnowSublimation)%dat(1)  ,&  ! intent(out): [dp] sublimation of ice from the snow surface (kg m-2 s-1)

 ! vegetation parameters
 heightCanopyTop         => mpar_data%var(iLookPARAM%heightCanopyTop)              ,&  ! intent(in): [dp] height of the top of the vegetation canopy (m)
 heightCanopyBottom      => mpar_data%var(iLookPARAM%heightCanopyBottom)           ,&  ! intent(in): [dp] height of the bottom of the vegetation canopy (m)

 ! snow parameters
 snowfrz_scale           => mpar_data%var(iLookPARAM%snowfrz_scale)                ,&  ! intent(in): [dp] scaling parameter for the snow freezing curve (K-1)

 ! soil parameters
 vGn_m                   => diag_data%var(iLookDIAG%scalarVGn_m)%dat(1)            ,&  ! intent(in): [dp] van Genutchen "m" parameter (-)
 vGn_n                   => mpar_data%var(iLookPARAM%vGn_n)                        ,&  ! intent(in): [dp] van Genutchen "n" parameter (-)
 vGn_alpha               => mpar_data%var(iLookPARAM%vGn_alpha)                    ,&  ! intent(in): [dp] van Genutchen "alpha" parameter (m-1)
 theta_sat               => mpar_data%var(iLookPARAM%theta_sat)                    ,&  ! intent(in): [dp] soil porosity (-)
 theta_res               => mpar_data%var(iLookPARAM%theta_res)                    ,&  ! intent(in): [dp] soil residual volumetric water content (-)
 specificStorage         => mpar_data%var(iLookPARAM%specificStorage)              ,&  ! intent(in): [dp] specific storage coefficient (m-1)
 fImpede                 => mpar_data%var(iLookPARAM%f_impede)                     ,&  ! intent(in): [dp] ice impedance parameter (-)

 ! model state variables (vegetation canopy)
 scalarCanairTemp        => prog_data%var(iLookPROG%scalarCanairTemp)%dat(1)       ,&  ! intent(inout): [dp] temperature of the canopy air space (K)
 scalarCanopyTemp        => prog_data%var(iLookPROG%scalarCanopyTemp)%dat(1)       ,&  ! intent(inout): [dp] temperature of the vegetation canopy (K)
 scalarCanopyIce         => prog_data%var(iLookPROG%scalarCanopyIce)%dat(1)        ,&  ! intent(inout): [dp] mass of ice on the vegetation canopy (kg m-2)
 scalarCanopyLiq         => prog_data%var(iLookPROG%scalarCanopyLiq)%dat(1)        ,&  ! intent(inout): [dp] mass of liquid water on the vegetation canopy (kg m-2)
 scalarCanopyWat         => prog_data%var(iLookPROG%scalarCanopyWat)%dat(1)        ,&  ! intent(inout): [dp] mass of total water on the vegetation canopy (kg m-2)

 ! model state variables (snow and soil domains)
 mLayerTemp              => prog_data%var(iLookPROG%mLayerTemp)%dat                ,&  ! intent(inout): [dp(:)] temperature of each snow/soil layer (K)
 mLayerVolFracIce        => prog_data%var(iLookPROG%mLayerVolFracIce)%dat          ,&  ! intent(inout): [dp(:)] volumetric fraction of ice (-)
 mLayerVolFracLiq        => prog_data%var(iLookPROG%mLayerVolFracLiq)%dat          ,&  ! intent(inout): [dp(:)] volumetric fraction of liquid water (-)
 mLayerVolFracWat        => prog_data%var(iLookPROG%mLayerVolFracWat)%dat          ,&  ! intent(inout): [dp(:)] volumetric fraction of total water (-)
 mLayerMatricHead        => prog_data%var(iLookPROG%mLayerMatricHead)%dat          ,&  ! intent(inout): [dp(:)] matric head (m)
 scalarAquiferStorage    => prog_data%var(iLookPROG%scalarAquiferStorage)%dat(1)    &  ! intent(inout): [dp   ] aquifer storage (m)

 )
 ! ---------------------------------------------------------------------------------------
 ! initialize error control
 err=0; message="systemSolv/"

 ! *****
 ! (0) PRELIMINARIES...
 ! ********************

 ! initialize the first flux call
 firstFluxCall=.true.

 ! set the flag to control printing
 printFlagInit=.false.
 printFlag=printFlagInit

 ! set the flag for pausing
 pauseProgress=.false.

 ! check dimensions
 if(printFlag)then
  print*, 'kl, ku, nBands = ', kl, ku, nBands
  print*, 'ixSup3, ixSup2, ixSup1, ixDiag = ', ixSup3, ixSup2, ixSup1, ixDiag
  print*, 'ixSub1, ixSub2, ixSub3         = ', ixSub1, ixSub2, ixSub3
 endif

 ! allocate space for the derivatives
 call allocLocal(deriv_meta(:),deriv_data,nSnow,nSoil,err,cmessage)
 if(err/=0)then; err=20; message=trim(message)//trim(cmessage); return; endif


 ! -----
 ! * initialize decisions/states and run initial checks...
 ! -------------------------------------------------------

 ! identify the matrix solution method
 ! (the type of matrix used to solve the linear system A.X=B)
 if(ixGroundwater==qbaseTopmodel .or. testBandDiagonal .or. forceFullMatrix)then
  ixMatrix=ixFullMatrix   ! full Jacobian matrix
 else
  ixMatrix=ixBandMatrix   ! band-diagonal matrix
 endif
 if(printFlag) print*, '(ixMatrix==ixFullMatrix) = ', (ixMatrix==ixFullMatrix)

 ! TEMPORARY!!
 ixSolve=ixMatrix


 ! check that dx is less that epsT
 if(dx>epsT)then
  err=20; message=trim(message)//'dx>epsT; will cause problems testing numerical derivatives'; return
 endif

 ! modify the groundwater representation for this single-column implementation
 select case(ixSpatialGroundwater)
  case(singleBasin); local_ixGroundwater = noExplicit    ! force no explicit representation of groundwater at the local scale
  case(localColumn); local_ixGroundwater = ixGroundwater ! go with the specified decision
  case default; err=20; message=trim(message)//'unable to identify spatial representation of groundwater'; return
 end select ! (modify the groundwater representation for this single-column implementation)

 ! define the scaled maximum step size (used in the line search)
 stpmax = stepMax*real(nState,dp)

 ! define canopy depth (m)
 canopyDepth = heightCanopyTop - heightCanopyBottom

 ! define the temperature of the upper boundary (K)
 upperBoundTemp = airtemp

 ! define canopy depth (m)
 canopyDepth = heightCanopyTop - heightCanopyBottom

 ! compute the total water in the vegetation canopy
 if(computeVegFlux)then
  scalarCanopyWat = scalarCanopyLiq + scalarCanopyIce   ! kg m-2
 endif

 ! compute the total water in snow
 if(nSnow>0)&
  mLayerVolFracWat(1:nSnow) = mLayerVolFracLiq(1:nSnow) + mLayerVolFracIce(1:nSnow)*(iden_ice/iden_water)

 ! compute the total water in soil
 mLayerVolFracWat(nSnow+1:nLayers) = mLayerVolFracLiq(nSnow+1:nLayers) + mLayerVolFracIce(nSnow+1:nLayers)

 ! get an initial canopy temperature if veg just starts protruding through snow on the ground
 if(computeVegFlux)then
  ! (NOTE: if canopy temperature is below absolute zero then canopy was previously buried by snow)
  if(scalarCanopyTemp < 0._dp .or. scalarCanairTemp < 0._dp)then
   message=trim(message)//'canopy temperatures below zero [previously this signified burial by snow; no longer done]'
  endif  ! (if canopy temperature undefined -- means canopy previously buried with snow)
 endif  ! (if computing vegetation fluxes -- canopy exposed)

 ! check indices
 if(iJac1 > nState .or. iJac2 > nState)then
  err=20; message=trim(message)//'index iJac1 or iJac2 is out of range'; return
 endif

 ! -----
 ! * allocate space for the model state vectors and derivative matrices...
 ! -----------------------------------------------------------------------

 ! allocate space for the baseflow derivatives
 if(ixGroundwater==qbaseTopmodel)then
  allocate(dBaseflow_dMatric(nSoil,nSoil),stat=err)
  if(err/=0)then; err=20; message=trim(message)//'unable to allocate space for the baseflow derivatives'; return; endif
 endif

 ! allocate space for the Jacobian matrix
 select case(ixSolve)
  case(ixFullMatrix); allocate(aJac(nState,nState),stat=err)
  case(ixBandMatrix); allocate(aJac(nBands,nState),stat=err)
  case default; err=20; message=trim(message)//'unable to identify option for the type of matrix'
 end select
 if(err/=0)then; err=20; message=trim(message)//'unable to allocate space for the Jacobian matrix'; return; endif

 ! allocate space for the numerical Jacobian matrix
 if(numericalJacobian)then
  allocate(fluxVec1(nState),nJac(nState,nState),stat=err)
  if(err/=0)then; err=20; message=trim(message)//'unable to allocate space for the flux vector and numerical Jacobian matrix'; return; endif
 endif  ! if calculating the numerical approximation of the Jacobian matrix

 ! allocate space for the band-diagonal matrix that is constructed from the full Jacobian matrix
 if(testBandDiagonal)then
  allocate(aJac_test(nBands,nState),stat=err)
  if(err/=0)then; err=20; message=trim(message)//'unable to allocate space for the band diagonal matrix'; return; endif
 endif

 ! -----
 ! * define components of derivative matrices that are constant over a time step (substep)...
 ! ------------------------------------------------------------------------------------------

 ! -----
 ! * initialize state vectors...
 ! -----------------------------

 ! initialize
 xInc(:)   = 0._dp  ! iteration increment

 ! initialize state vectors
 call popStateVec(&
                  ! input
                  computeVegFlux,          & ! intent(in):    flag to denote if computing energy flux over vegetation
                  canopyDepth,             & ! intent(in):    canopy depth (m)
                  prog_data,               & ! intent(in):    model prognostic variables for a local HRU
                  diag_data,               & ! intent(in):    model diagnostic variables for a local HRU
                  indx_data,               & ! intent(in):    indices defining model states and layers
                  ! output
                  stateVecInit,            & ! intent(out):   initial model state vector (mixed units)
                  fScale,                  & ! intent(out):   function scaling vector (mixed units)
                  xScale,                  & ! intent(out):   variable scaling vector (mixed units)
                  sMul,                    & ! intent(out):   multiplier for state vector (used in the residual calculations)
                  dMat,                    & ! intent(out):   diagonal of the Jacobian matrix (excludes fluxes) 
                  err,cmessage)              ! intent(out):   error control
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

 ! extract variables from the model state vector
 call varExtract(&
                 ! input
                 stateVecTrial,                             & ! intent(in):    model state vector (mixed units)
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
                 err,cmessage)                                ! intent(out):   error control
 if(err/=0)then; message=trim(message)//trim(cmessage); return; endif  ! (check for errors)

 ! -----
 ! * compute the initial function evaluation...
 ! --------------------------------------------

 ! compute flux vector and residual
 call xFluxResid(&
                 ! input
                 stateVecTrial,         & ! intent(in): full state vector (mixed units)
                 scalarCanopyLiqTrial,  & ! intent(in): liquid water storage in the canopy (kg m-2)
                 scalarCanopyIceTrial,  & ! intent(in): ice storage in the canopy (kg m-2)
                 mLayerVolFracLiqTrial, & ! intent(in): volumetric fraction of liquid water (-)
                 mLayerVolFracIceTrial, & ! intent(in): volumetric fraction of ice (-)
                 ! output
                 fluxVec0,              & ! intent(out): flux vector (mixed units)
                 rVec,                  & ! intent(out): residual vector (mixed units)
                 err,cmessage)            ! intent(out): error code and error message
 if(err/=0)then; message=trim(message)//trim(cmessage); return; endif  ! (check for errors)

 ! compute the function evaluation
 fOld=0.5_dp*norm2(real(rVec, dp)/fScale)  ! NOTE: norm2 = sqrt(sum((rVec/fScale)**2._dp))

 write(*,'(a,1x,10(e17.6,1x))') 'fluxVec0(iJac1:iJac2) = ', fluxVec0(iJac1:iJac2)
 write(*,'(a,1x,10(e17.6,1x))') 'rVec(iJac1:iJac2)     = ', rVec(iJac1:iJac2)

 ! WARNING!!!!! TESTING!!!
 firstFluxCall=.true.

 ! compute the flux and the residual vector for a given state vector
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
                 ! output
                 feasible,                & ! intent(out):   flag to denote the feasibility of the solution
                 fluxVec0,                & ! intent(out):   flux vector
                 rVec,                    & ! intent(out):   residual vector
                 fOld,                    & ! intent(out):   function evaluation
                 err,cmessage)              ! intent(out):   error control
 if(err/=0)then; message=trim(message)//trim(cmessage); return; endif  ! (check for errors)

 ! check feasibility (state vector SHOULD be feasible at this point)
 if(.not.feasible)then
  message=trim(message)//'unfeasible state vector'
  err=20; return
 endif

 write(*,'(a,1x,10(e17.6,1x))') 'fluxVec0(iJac1:iJac2) = ', fluxVec0(iJac1:iJac2)
 write(*,'(a,1x,10(e17.6,1x))') 'rVec(iJac1:iJac2)     = ', rVec(iJac1:iJac2)
 pause


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

  ! test
  if(printFlag)then
   print*, '***'
   write(*,'(a,1x,f10.2,1x,2(i4,1x),l1)') '*** new iteration: dt, iter, nstate, computeVegFlux = ', dt, iter, nstate, computeVegFlux
  endif

  ! -----
  ! * compute Jacobian...
  ! ---------------------

  ! compute terms in the Jacobian for vegetation (excluding fluxes)
  ! NOTE: energy for vegetation is computed *within* the iteration loop as it includes phase change
  if(computeVegFlux)then
   dMat(ixVegNrg) = scalarBulkVolHeatCapVeg + LH_fus*iden_water*dTheta_dTkCanopy       ! volumetric heat capacity of the vegetation (J m-3 K-1)
  endif

  ! compute additional terms for the Jacobian for the snow-soil domain (excluding fluxes)
  ! NOTE: energy for vegetation is computed *within* the iteration loop as it includes phase change
  dMat(ixSnowSoilNrg) = mLayerVolHtCapBulk(1:nLayers) + LH_fus*iden_water*mLayerdTheta_dTk(1:nLayers)

  ! compute additional terms for the Jacobian for the soil domain (excluding fluxes)
  if(ixRichards==moisture)then; err=20; message=trim(message)//'have not implemented the moisture-based form of RE yet'; return; endif
  dMat(ixSoilOnlyHyd) = dVolTot_dPsi0(1:nSoil) + dCompress_dPsi(1:nSoil)

  ! compute the analytical Jacobian matrix
  select case(ixSolve)
   case(ixFullMatrix); call analJacob(err,cmessage)
   case(ixBandMatrix); call cpactBand(err,cmessage)
   case default; err=20; message=trim(message)//'unable to identify option for the type of matrix'
  end select
  if(err/=0)then; message=trim(message)//trim(cmessage); return; endif  ! (check for errors)

  ! *** testing: compute the numerical approximation of the Jacobian matrix
  if(numericalJacobian)then
   printFlag=.false.
   call numlJacob(stateVecTrial,fluxVec0,rVec,err,cmessage)
   if(err/=0)then; message=trim(message)//trim(cmessage); return; endif  ! (check for errors)
   printFlag=printFlagInit
  endif  ! if computing the numerical Jacobian matrix

  ! -----
  ! * solve linear system...
  ! ------------------------

  ! use the lapack routines to solve the linear system A.X=B
  call lapackSolv(aJac,rVec,grad,xInc,err,cmessage)
  if(err/=0)then; message=trim(message)//trim(cmessage); return; endif  ! (check for errors)

  if(printFlag)then
   write(*,'(a,1x,10(e17.10,1x))') 'xInc = ', xInc(iJac1:iJac2)
  endif

  ! -----
  ! * impose constraints... 
  ! -----------------------

  ! constrain iteration increment
  call imposeConstraints()

  if(printFlag)then
   write(*,'(a,1x,10(e17.10,1x))') 'xInc = ', xInc(iJac1:iJac2)
  endif

  ! -----
  ! * compute model fluxes and residual
  !    NOTE: refine residual with line search...
  ! --------------------------------------------
  call lineSearch(&
                  ! input
                  .true.,                  & ! intent(in): flag to denote the need to perform line search
                  stateVecTrial,           & ! intent(in): initial state vector
                  fOld,                    & ! intent(in): function value for trial state vector (mixed units)
                  grad,                    & ! intent(in): gradient of the function vector (mixed units)
                  xInc,                    & ! intent(in): iteration increment (mixed units)
                  ! output
                  stateVecNew,             & ! intent(out): new state vector (m)
                  fluxVec0,                & ! intent(out): new flux vector (mixed units)
                  rVec,                    & ! intent(out): new residual vector (mixed units)
                  fNew,                    & ! intent(out): new function value (mixed units)
                  converged,               & ! intent(out): convergence flag
                  err,cmessage)              ! intent(out): error control
  if(err>0)then; message=trim(message)//trim(cmessage); return; endif  ! (check for errors)

  ! use full iteration increment if converged all the way to the original value
  if(err<0)then
   call lineSearch(&
                   ! input
                   .false.,                 & ! intent(in): flag to denote the need to perform line search
                   stateVecTrial,           & ! intent(in): initial state vector
                   fOld,                    & ! intent(in): function value for trial state vector (mixed units)
                   grad,                    & ! intent(in): gradient of the function vector (mixed units)
                   xInc,                    & ! intent(in): iteration increment (mixed units)
                   ! output
                   stateVecNew,             & ! intent(out): new state vector (m)
                   fluxVec0,                & ! intent(out): new flux vector (mixed units)
                   rVec,                    & ! intent(out): new residual vector (mixed units)
                   fNew,                    & ! intent(out): new function value (mixed units)
                   converged,               & ! intent(out): convergence flag
                   err,cmessage)              ! intent(out): error control
   if(err/=0)then; message=trim(message)//trim(cmessage); return; endif  ! (check for errors)
  endif

  ! update function evaluation and states
  fOld          = fNew
  stateVecTrial = stateVecNew

  ! calculate the iteration increment, evaluate the new state, and refine if necessary
  call stateTrial(&
                  ! input: model control
                  dt,                      & ! intent(in):    length of the time step (seconds)
                  nSnow,                   & ! intent(in):    number of snow layers
                  nSoil,                   & ! intent(in):    number of soil layers
                  nLayers,                 & ! intent(in):    total number of layers
                  nState,                  & ! intent(in):    total number of state variables
                  ixMatrix,                & ! intent(in):    type of matrix (full or band diagonal)
                  firstSubStep,            & ! intent(in):    flag to indicate if we are processing the first sub-step
                  firstFluxCall,           & ! intent(inout): flag to indicate if we are processing the first flux call
                  computeVegFlux,          & ! intent(in):    flag to indicate if we need to compute fluxes over vegetation
                  ! input: state vectors
                  stateVecTrial,           & ! intent(in):    trial state vector
                  rVec,                    & ! intent(in):    residual vector
                  aJac,                    & ! intent(in):    Jacobian matrix
                  fScale,                  & ! intent(in):    function scaling vector
                  xScale,                  & ! intent(in):    "variable" scaling vector, i.e., for state variables
                  sMul,                    & ! intent(in):    state vector multiplier (used in the residual calculations)
                  dMat,                    & ! intent(in):    diagonal matrix (excludes flux derivatives)
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
                  ! output
                  stateVecNew,             & ! intent(out):   new state vector
                  fluxVecNew,              & ! intent(out):   new flux vector
                  resVecNew,               & ! intent(out):   new residual vector
                  fNew,                    & ! intent(out):   new function evaluation
                  converged,               & ! intent(out):   convergence flag
                  err,cmessage)              ! intent(out):   error control
  if(err/=0)then; message=trim(message)//trim(cmessage); return; endif  ! (check for errors)




  ! exit iteration loop if converged
  if(converged) exit

  ! check convergence
  if(niter==maxiter)then; err=-20; message=trim(message)//'failed to converge'; return; endif
  !print*, 'PAUSE: iterating'; read(*,*)

 end do  ! iterating
 !print*, 'PAUSE: after iterations'; read(*,*)


 ! -----
 ! * update states and compute total volumetric melt...
 ! ----------------------------------------------------

 ! update temperatures (ensure new temperature is consistent with the fluxes)
 stateVecTrial(ixSnowSoilNrg) = stateVecInit(ixSnowSoilNrg) + (fluxVec0(ixSnowSoilNrg)*dt + real(rAdd(ixSnowSoilNrg), dp))/real(sMul(ixSnowSoilNrg), dp)

 ! update volumetric water content in the snow (ensure change in state is consistent with the fluxes)
 ! NOTE: for soil water balance is constrained within the iteration loop
 if(nSnow>0)&
 stateVecTrial(ixSnowOnlyWat) = stateVecInit(ixSnowOnlyWat) + (fluxVec0(ixSnowOnlyWat)*dt + real(rAdd(ixSnowOnlyWat), dp))

 ! compute total baseflow from the soil zone (needed for mass balance checks)
 scalarSoilBaseflow = sum(mLayerBaseflow)

 ! update states: compute liquid water and ice content from total water content
 call updatState(&
                 stateVecTrial,         & ! intent(in): full state vector (mixed units)
                 mLayerVolFracLiqTrial, & ! intent(out): volumetric fraction of liquid water (-)
                 mLayerVolFracIceTrial, & ! intent(out): volumetric fraction of ice (-)
                 scalarCanopyLiqTrial,  & ! intent(out): mass of canopy liquid (kg m-2)
                 scalarCanopyIceTrial,  & ! intent(out): mass of canopy ice (kg m-2)
                 err,cmessage)            ! intent(out): error code and error message
 if(err/=0)then; message=trim(message)//trim(cmessage); return; endif  ! (check for errors)

 ! check the mass balance for the soil domain
 ! NOTE: this should never fail since did not converge if water balance was not within tolerance=absConvTol_watbal
 if(checkMassBalance)then
  balance0 = sum( (mLayerVolFracLiq(nSnow+1:nLayers)      + mLayerVolFracIce(nSnow+1:nLayers)      )*mLayerDepth(nSnow+1:nLayers) )
  balance1 = sum( (mLayerVolFracLiqTrial(nSnow+1:nLayers) + mLayerVolFracIceTrial(nSnow+1:nLayers) )*mLayerDepth(nSnow+1:nLayers) )
  vertFlux = -(iLayerLiqFluxSoil(nSoil) - iLayerLiqFluxSoil(0))*dt  ! m s-1 --> m
  tranSink = sum(mLayerTranspire)*dt                                ! m s-1 --> m
  baseSink = sum(mLayerBaseflow)*dt                                 ! m s-1 --> m
  compSink = sum(mLayerCompress(1:nSoil) * mLayerDepth(nSnow+1:nLayers) ) ! dimensionless --> m
  liqError = balance1 - (balance0 + vertFlux + tranSink - baseSink - compSink)
  if(abs(liqError) > absConvTol_watbal*10._dp)then  ! *10 to avoid precision issues
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

 ! deallocate space for the baseflow derivatives
 if(ixGroundwater==qbaseTopmodel)then
  deallocate(dBaseflow_dMatric,stat=err)
  if(err/=0)then; err=20; message=trim(message)//'unable to deallocate space for the baseflow derivatives'; return; endif
 endif

 ! deallocate space for the Jacobian matrix
 deallocate(aJac,stat=err)
 if(err/=0)then; err=20; message=trim(message)//'unable to deallocate space for the Jacobian matrix'; return; endif

 ! deallocate space for the variables used to create the numerical Jacobian matrix
 if(numericalJacobian)then
  deallocate(fluxVec1,nJac,stat=err)
  if(err/=0)then; err=20; message=trim(message)//'unable to deallocate space for the flux vector and numerical Jacobian matrix'; return; endif
 endif

 ! de allocate space for the band-diagonal matrix that is constructed from the full Jacobian matrix
 if(testBandDiagonal)then
  deallocate(aJac_test,stat=err)
  if(err/=0)then; err=20; message=trim(message)//'unable to deallocate space for the band diagonal matrix'; return; endif
 endif

 ! end associate statement
 end associate

 contains


  ! *********************************************************************************************************
  ! internal subroutine imposeConstraints: impose solution constraints
  ! *********************************************************************************************************
  subroutine imposeConstraints()
  real(dp),parameter  :: zMaxTempIncrement=1._dp

  ! associate variables with indices of model state variables
  associate(&
  ixCasNrg                => indx_data%var(iLookINDEX%ixCasNrg)%dat(1)              ,&  ! intent(in): [i4b] index of canopy air space energy state variable
  ixVegNrg                => indx_data%var(iLookINDEX%ixVegNrg)%dat(1)              ,&  ! intent(in): [i4b] index of canopy energy state variable
  ixVegWat                => indx_data%var(iLookINDEX%ixVegWat)%dat(1)              ,&  ! intent(in): [i4b] index of canopy hydrology state variable (mass)
  ixTopNrg                => indx_data%var(iLookINDEX%ixTopNrg)%dat(1)              ,&  ! intent(in): [i4b] index of upper-most energy state in the snow-soil subdomain
  ixTopWat                => indx_data%var(iLookINDEX%ixTopWat)%dat(1)              ,&  ! intent(in): [i4b] index of upper-most total water state in the snow-soil subdomain
  ixTopMat                => indx_data%var(iLookINDEX%ixTopMat)%dat(1)              ,&  ! intent(in): [i4b] index of upper-most matric head state in the soil subdomain
  ixSnowSoilNrg           => indx_data%var(iLookINDEX%ixSnowSoilNrg)%dat            ,&  ! intent(in): [i4b(:)] indices for energy states in the snow-soil subdomain
  ixSnowSoilWat           => indx_data%var(iLookINDEX%ixSnowSoilWat)%dat            ,&  ! intent(in): [i4b(:)] indices for total water states in the snow-soil subdomain
  ixSnowOnlyNrg           => indx_data%var(iLookINDEX%ixSnowOnlyNrg)%dat            ,&  ! intent(in): [i4b(:)] indices for energy states in the snow subdomain
  ixSnowOnlyWat           => indx_data%var(iLookINDEX%ixSnowOnlyWat)%dat            ,&  ! intent(in): [i4b(:)] indices for total water states in the snow subdomain
  ixSoilOnlyNrg           => indx_data%var(iLookINDEX%ixSoilOnlyNrg)%dat            ,&  ! intent(in): [i4b(:)] indices for energy states in the soil subdomain
  ixNrgOnly               => indx_data%var(iLookINDEX%ixNrgOnly)%dat                ,&  ! intent(in): [i4b(:)] list of indices for all energy states
  layerType               => indx_data%var(iLookINDEX%layerType)%dat                 &  ! intent(in): [i4b(:)] type of each layer in the snow+soil domain (snow or soil)
  ) ! associating variables with indices of model state variables

  ! ** limit temperature increment to zMaxTempIncrement
  !if(any(abs(xInc(ixNrgOnly)) > zMaxTempIncrement))then
  ! iMax       = maxloc( abs(xInc(ixNrgOnly)) )                     ! index of maximum temperature increment
  ! xIncFactor = abs( zMaxTempIncrement/xInc(ixNrgOnly(iMax(1))) )  ! scaling factor for the iteration increment (-)
  ! xInc       = xIncFactor*xInc
  !endif

  ! vegetation
  if(computeVegFlux)then
   if(abs(xInc(ixVegNrg)) > 1._dp)then
    xIncFactor = abs(1._dp/xInc(ixVegNrg))  ! scaling factor for the iteration increment (-)
    xInc       = xIncFactor*xInc            ! scale iteration increments
   endif
  endif

  ! snow and soil
  if(any(abs(xInc(ixSnowSoilNrg)) > 1._dp))then
   iMax       = maxloc( abs(xInc(ixSnowSoilNrg)) )                   ! index of maximum temperature increment
   xIncFactor = abs( 1._dp/xInc(ixSnowSoilNrg(iMax(1))) )            ! scaling factor for the iteration increment (-)
   xInc       = xIncFactor*xInc
  endif

  ! ** impose solution constraints for vegetation
  ! (stop just above or just below the freezing point if crossing)
  if(computeVegFlux)then

   ! --------------------------------------------------------------------------------------------------------------------
   ! canopy temperatures

   ! initialize
   critDiff    = Tfreeze - stateVecTrial(ixVegNrg)
   crosTempVeg = .false.

   ! initially frozen (T < Tfreeze)
   if(critDiff > 0._dp)then
    if(xInc(ixVegNrg) > critDiff)then
     crosTempVeg = .true.
     cInc        = critDiff + epsT  ! constrained temperature increment (K)
    endif

   ! initially unfrozen (T > Tfreeze)
   else
    if(xInc(ixVegNrg) < critDiff)then
     crosTempVeg = .true.
     cInc        = critDiff - epsT  ! constrained temperature increment (K)
    endif

   endif  ! switch between frozen and unfrozen

   ! scale iterations
   if(crosTempVeg)then
    xIncFactor  = cInc/xInc(ixVegNrg)  ! scaling factor for the iteration increment (-)
    xInc        = xIncFactor*xInc      ! scale iteration increments
   endif

   ! --------------------------------------------------------------------------------------------------------------------
   ! canopy liquid water

   ! check if new value of storage will be negative
   if(stateVecTrial(ixVegWat)+xInc(ixVegWat) < 0._dp)then
    ! scale iteration increment
    cInc       = -0.5_dp*stateVecTrial(ixVegWat)                                  ! constrained iteration increment (K) -- simplified bi-section
    xIncFactor = cInc/xInc(ixVegWat)                                              ! scaling factor for the iteration increment (-)
    xInc       = xIncFactor*xInc                                                  ! new iteration increment
   endif

  endif  ! if computing fluxes through vegetation

  ! ** impose solution constraints for snow
  if(nSnow > 0)then

   ! --------------------------------------------------------------------------------------------------------------------
   ! loop through snow layers
   !checksnow: do iState=1,nSnow  ! necessary to ensure that NO layers rise above Tfreeze
    ! - get updated snow temperatures
    mLayerTempCheck = stateVecTrial(ixSnowOnlyNrg) + xInc(ixSnowOnlyNrg)
    ! - check sub-freezing temperatures for snow
    if(any(mLayerTempCheck > Tfreeze))then
     ! scale iteration increment
     iMax       = maxloc(mLayerTempCheck)                                          ! index of maximum temperature
     cInc       = 0.5_dp*(Tfreeze - stateVecTrial(ixSnowOnlyNrg(iMax(1))) )        ! constrained temperature increment (K) -- simplified bi-section
     xIncFactor = cInc/xInc(ixSnowOnlyNrg(iMax(1)))                                ! scaling factor for the iteration increment (-)
     xInc       = xIncFactor*xInc
   ! else    ! if snow temperature > freezing
   !  exit checkSnow
    endif   ! if snow temperature > freezing
   !end do checkSnow

   ! --------------------------------------------------------------------------------------------------------------------
   ! - check if drain more than what is available
   ! NOTE: change in total water is only due to liquid flux

   ! get new volumetric fraction of liquid water
   mLayerVolFracLiqCheck = mLayerVolFracLiqTrial(1:nSnow)+xInc(ixSnowOnlyWat)
   drainFlag(:) = .false.

   do iLayer=1,nSnow
    if(mLayerVolFracLiqCheck(iLayer) < 0._dp)then
     drainFlag(iLayer) = .true.
     xInc(ixSnowOnlyWat(iLayer)) = -0.5_dp*mLayerVolFracLiqTrial(iLayer)
    endif
   end do

  endif   ! if snow layers exist

  ! --------------------------------------------------------------------------------------------------------------------
  ! ** impose solution constraints for soil
  do iLayer=1,nSoil

   ! initialize crossing flag
   crosFlag(iLayer) = .false.

   ! identify indices for energy and mass state variables
   ixNrg = ixSnowSoilNrg(nSnow+iLayer)
   ixLiq = ixSnowSoilWat(nSnow+iLayer)

   ! identify the critical point when soil begins to freeze (TcSoil)
   xPsi00 = stateVecTrial(ixLiq) + xInc(ixLiq)
   TcSoil = Tfreeze + min(xPsi00,0._dp)*gravity*Tfreeze/LH_fus  ! (NOTE: J = kg m2 s-2, so LH_fus is in units of m2 s-2)

   ! get the difference from the current state and the crossing point (K)
   critDiff = TcSoil - stateVecTrial(ixNrg)

   ! * initially frozen (T < TcSoil)
   if(critDiff > 0._dp)then

    ! (check crossing above zero)
    if(xInc(ixNrg) > critDiff)then
     crosFlag(iLayer) = .true.
     xInc(ixNrg) = critDiff + epsT  ! set iteration increment to slightly above critical temperature
    endif

   ! * initially unfrozen (T > TcSoil)
   else

    ! (check crossing below zero)
    if(xInc(ixNrg) < critDiff)then
     crosFlag(iLayer) = .true.
     xInc(ixNrg) = critDiff - epsT  ! set iteration increment to slightly below critical temperature
    endif

   endif  ! (switch between initially frozen and initially unfrozen)

   ! place constraint for matric head
   if(xInc(ixLiq) > 1._dp .and. stateVecTrial(ixLiq) > 0._dp)then
    xInc(ixLiq) = 1._dp
    pauseProgress=.true.
   endif  ! if constraining matric head

  end do  ! (loop through soil layers)

  ! end association with variables with indices of model state variables
  end associate

  end subroutine imposeConstraints


  ! *********************************************************************************************************
  ! internal subroutine updatState: update model states
  ! *********************************************************************************************************
  subroutine updatState(&
                        stateVecTrial,         & ! intent(in):  full state vector (mixed units)
                        mLayerVolFracLiqTrial, & ! intent(out): volumetric fraction of liquid water (-)
                        mLayerVolFracIceTrial, & ! intent(out): volumetric fraction of ice (-)
                        scalarCanopyLiqTrial,  & ! intent(out): liquid water storage in the canopy (kg m-2)
                        scalarCanopyIceTrial,  & ! intent(out): ice storage in the canopy (kg m-2)
                        err,message)             ! intent(out): error code and error message
  ! --------------------------------------------------------------
  implicit none
  ! input
  real(dp),intent(in)            :: stateVecTrial(:)          ! full model state vector (mixed units)
  real(dp),intent(out)           :: mLayerVolFracLiqTrial(:)  ! volumetric fraction of liquid water (-)
  real(dp),intent(out)           :: mLayerVolFracIceTrial(:)  ! volumetric fraction of ice (-)
  real(dp),intent(out)           :: scalarCanopyLiqTrial      ! liquid water storage in the canopy (kg m-2)
  real(dp),intent(out)           :: scalarCanopyIceTrial      ! ice storage in the canopy (kg m-2)
  ! output
  integer(i4b),intent(out)       :: err                       ! error code
  character(*),intent(out)       :: message                   ! error message
  ! local
  character(LEN=256)             :: cmessage                  ! error message of downwind routine
  ! initialize error control
  err=0; message='updatState/'

  ! get the necessary variables from the data structures
  associate(&

  ! indices of model state variables
  ixCasNrg                => indx_data%var(iLookINDEX%ixCasNrg)%dat(1)              ,&  ! intent(in): [i4b] index of canopy air space energy state variable
  ixVegNrg                => indx_data%var(iLookINDEX%ixVegNrg)%dat(1)              ,&  ! intent(in): [i4b] index of canopy energy state variable
  ixVegWat                => indx_data%var(iLookINDEX%ixVegWat)%dat(1)              ,&  ! intent(in): [i4b] index of canopy hydrology state variable (mass)
  ixTopNrg                => indx_data%var(iLookINDEX%ixTopNrg)%dat(1)              ,&  ! intent(in): [i4b] index of upper-most energy state in the snow-soil subdomain
  ixTopWat                => indx_data%var(iLookINDEX%ixTopWat)%dat(1)              ,&  ! intent(in): [i4b] index of upper-most total water state in the snow-soil subdomain
  ixTopMat                => indx_data%var(iLookINDEX%ixTopMat)%dat(1)              ,&  ! intent(in): [i4b] index of upper-most matric head state in the soil subdomain
  ixSnowSoilNrg           => indx_data%var(iLookINDEX%ixSnowSoilNrg)%dat            ,&  ! intent(in): [i4b(:)] indices for energy states in the snow-soil subdomain
  ixSnowSoilWat           => indx_data%var(iLookINDEX%ixSnowSoilWat)%dat            ,&  ! intent(in): [i4b(:)] indices for total water states in the snow-soil subdomain
  ixSnowOnlyNrg           => indx_data%var(iLookINDEX%ixSnowOnlyNrg)%dat            ,&  ! intent(in): [i4b(:)] indices for energy states in the snow subdomain
  ixSnowOnlyWat           => indx_data%var(iLookINDEX%ixSnowOnlyWat)%dat            ,&  ! intent(in): [i4b(:)] indices for total water states in the snow subdomain
  ixSoilOnlyNrg           => indx_data%var(iLookINDEX%ixSoilOnlyNrg)%dat            ,&  ! intent(in): [i4b(:)] indices for energy states in the soil subdomain
  ixSoilOnlyHyd           => indx_data%var(iLookINDEX%ixSoilOnlyHyd)%dat            ,&  ! intent(in): [i4b(:)] indices for hydrology states in the soil subdomain
  layerType               => indx_data%var(iLookINDEX%layerType)%dat                ,&  ! intent(in): [i4b(:)] type of each layer in the snow+soil domain (snow or soil)

  ! model state variables (vegetation canopy)
  scalarCanopyIce         => prog_data%var(iLookPROG%scalarCanopyIce)%dat(1)        ,&  ! intent(inout): [dp] mass of ice on the vegetation canopy (kg m-2)
  scalarCanopyLiq         => prog_data%var(iLookPROG%scalarCanopyLiq)%dat(1)        ,&  ! intent(inout): [dp] mass of liquid water on the vegetation canopy (kg m-2)

  ! layer depth
  mLayerDepth             => prog_data%var(iLookPROG%mLayerDepth)%dat               ,&  ! intent(in): [dp] depth of each layer (m)

  ! snow parameters
  snowfrz_scale           => mpar_data%var(iLookPARAM%snowfrz_scale)                ,&  ! intent(in): [dp] scaling parameter for the snow freezing curve (K-1)

  ! soil parameters
  vGn_m                   => diag_data%var(iLookDIAG%scalarVGn_m)%dat(1)            ,&  ! intent(in): [dp] van Genutchen "m" parameter (-)
  vGn_n                   => mpar_data%var(iLookPARAM%vGn_n)                        ,&  ! intent(in): [dp] van Genutchen "n" parameter (-)
  vGn_alpha               => mpar_data%var(iLookPARAM%vGn_alpha)                    ,&  ! intent(in): [dp] van Genutchen "alpha" parameter (m-1)
  theta_sat               => mpar_data%var(iLookPARAM%theta_sat)                    ,&  ! intent(in): [dp] soil porosity (-)
  theta_res               => mpar_data%var(iLookPARAM%theta_res)                     &  ! intent(in): [dp] soil residual volumetric water content (-)
  )

  ! update states for the vegetation canopy
  if(computeVegFlux)then
   fracLiqVeg    = fracliquid(stateVecTrial(ixVegNrg),snowfrz_scale)  ! fraction of liquid water (-)
   totalWaterVeg = stateVecTrial(ixVegWat)                            ! total water (kg m-2)
   scalarCanopyLiqTrial = fracLiqVeg*totalWaterVeg                    ! mass of liquid water on the canopy (kg m-2)
   scalarCanopyIceTrial = (1._dp - fracLiqVeg)*totalWaterVeg          ! mass of ice on the canopy (kg m-2)

  ! ensure that states for the veg canopy are defined
  else
   scalarCanopyLiqTrial = scalarCanopyLiq
   scalarCanopyIceTrial = scalarCanopyIce
  endif

  ! loop through layers in the snow+soil sub-domain
  do iLayer=1,nLayers
   select case(layerType(iLayer))

    !** snow
    case(ix_snow)
     ! update states
     call updateSnow(&
                     ! input
                     stateVecTrial(ixSnowSoilNrg(iLayer)),      & ! intent(in): layer temperature (K)
                     stateVecTrial(ixSnowSoilWat(iLayer)),      & ! intent(in): volumetric fraction of total water (-)
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
                     stateVecTrial(ixSnowSoilNrg(iLayer)),      & ! intent(in): layer temperature (K)
                     stateVecTrial(ixSoilOnlyHyd(iLayer-nSnow)),& ! intent(in): matric head (m)
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

   ! sanity check
   if(mLayerVolFracIceTrial(iLayer) < -tiny(theta_sat))then
    write(message,'(a,i0,a,e20.10,a)')trim(message)//"volumetric ice content < 0; iLayer=",iLayer,"; mLayerVolFracIce =",mLayerVolFracIceTrial(iLayer),"]"
    err=10; return
   endif

  end do   ! looping through layers

  ! end association to the variables in the data structures
  end associate

  end subroutine updatState


  ! *********************************************************************************************************
  ! internal subroutine xFluxResid: compute fluxes and the residual
  ! *********************************************************************************************************
  subroutine xFluxResid(&
                        ! input
                        stateVec,                      & ! intent(in): full state vector (mixed units)
                        scalarCanopyLiqLocal,          & ! intent(in): trial value for the liquid water on the vegetation canopy (kg m-2)
                        scalarCanopyIceLocal,          & ! intent(in): trial value for the ice on the vegetation canopy (kg m-2)
                        mLayerVolFracLiqLocal,         & ! intent(in): trial value for the volumetric liquid water content in each snow and soil layer (-)
                        mLayerVolFracIceLocal,         & ! intent(in): trial value for the volumetric ice in each snow and soil layer (-)
                        ! output
                        fVec,                          & ! intent(out): flux vector (mixed units)
                        rVec,                          & ! intent(out): residual vector (mixed units)
                        err,message)                     ! intent(out): error code and error message
  ! --------------------------------------------------------------
  implicit none
  ! input variables
  real(dp),intent(in)            :: stateVec(:)               ! model state vector (mixed units)
  real(dp),intent(in)            :: scalarCanopyLiqLocal      ! trial value for mass of liquid water on the vegetation canopy (kg m-2)
  real(dp),intent(in)            :: scalarCanopyIceLocal      ! trial value for mass of ice on the vegetation canopy (kg m-2)
  real(dp),intent(in)            :: mLayerVolFracLiqLocal(:)  ! trial value for volumetric fraction of liquid water (-)
  real(dp),intent(in)            :: mLayerVolFracIceLocal(:)  ! trial value for volumetric fraction of ice (-)
  ! output variabes
  real(dp),intent(out)           :: fVec(:)                   ! flux vector (mixed units)
  real(qp),intent(out)           :: rVec(:)  ! NOTE: qp       ! residual vector (mixed units)
  integer(i4b),intent(out)       :: err                       ! error code
  character(*),intent(out)       :: message                   ! error message
  ! --------------------------------------------------------------
  ! general local variables
  character(LEN=256)             :: cmessage                  ! error message of downwind routine
  ! trial state variables (vegetation canopy)
  real(dp)                       :: scalarCanairTempTrial     ! trial value for temperature of the canopy air space (K)
  real(dp)                       :: scalarCanopyTempTrial     ! trial value for temperature of the vegetation canopy (K)
  real(dp)                       :: scalarCanopyWatTrial      ! trial value for mass of total water on the vegetation canopy (kg m-2)
  ! trial state variables (snow and soil domains)
  real(dp),dimension(nLayers)    :: mLayerTempTrial           ! trial value for temperature of each snow/soil layer (K)
  real(dp),dimension(nLayers)    :: mLayerVolFracWatTrial     ! trial value for volumetric fraction of total water (-)
  real(dp),dimension(nSoil)      :: mLayerMatricHeadTrial     ! trial value for matric head (m)
  ! temporary vectors for the soil sub-domain
  real(dp),dimension(nSoil)      :: vThetaInit                ! liquid equivalent of total water at the start of the step
  real(dp),dimension(nSoil)      :: vThetaTrial               ! liquid equivalent of total water at the current iteration
  ! initialize error control
  err=0; message='xFluxResid/'

  ! -----
  ! * associate desired variables from data structures...
  ! -----------------------------------------------------
  associate(&
  ! model decisions
  ixRichards              => model_decisions(iLookDECISIONS%f_Richards)%iDecision   ,&  ! intent(in): [i4b] index of the form of Richards' equation

  ! soil parameters
  vGn_m                   => diag_data%var(iLookDIAG%scalarVGn_m)%dat(1)            ,&  ! intent(in): [dp] van Genutchen "m" parameter (-)
  vGn_n                   => mpar_data%var(iLookPARAM%vGn_n)                        ,&  ! intent(in): [dp] van Genutchen "n" parameter (-)
  vGn_alpha               => mpar_data%var(iLookPARAM%vGn_alpha)                    ,&  ! intent(in): [dp] van Genutchen "alpha" parameter (m-1)
  theta_sat               => mpar_data%var(iLookPARAM%theta_sat)                    ,&  ! intent(in): [dp] soil porosity (-)
  theta_res               => mpar_data%var(iLookPARAM%theta_res)                    ,&  ! intent(in): [dp] soil residual volumetric water content (-)
  specificStorage         => mpar_data%var(iLookPARAM%specificStorage)              ,&  ! intent(in): [dp] specific storage coefficient (m-1)

  ! indices of model state variables
  ixCasNrg                => indx_data%var(iLookINDEX%ixCasNrg)%dat(1)              ,&  ! intent(in): [i4b] index of canopy air space energy state variable
  ixVegNrg                => indx_data%var(iLookINDEX%ixVegNrg)%dat(1)              ,&  ! intent(in): [i4b] index of canopy energy state variable
  ixVegWat                => indx_data%var(iLookINDEX%ixVegWat)%dat(1)              ,&  ! intent(in): [i4b] index of canopy hydrology state variable (mass)
  ixTopNrg                => indx_data%var(iLookINDEX%ixTopNrg)%dat(1)              ,&  ! intent(in): [i4b] index of upper-most energy state in the snow-soil subdomain
  ixTopWat                => indx_data%var(iLookINDEX%ixTopWat)%dat(1)              ,&  ! intent(in): [i4b] index of upper-most total water state in the snow-soil subdomain
  ixTopMat                => indx_data%var(iLookINDEX%ixTopMat)%dat(1)              ,&  ! intent(in): [i4b] index of upper-most matric head state in the soil subdomain
  ixSnowSoilNrg           => indx_data%var(iLookINDEX%ixSnowSoilNrg)%dat            ,&  ! intent(in): [i4b(:)] indices for energy states in the snow-soil subdomain
  ixSnowSoilWat           => indx_data%var(iLookINDEX%ixSnowSoilWat)%dat            ,&  ! intent(in): [i4b(:)] indices for total water states in the snow-soil subdomain
  ixSnowOnlyNrg           => indx_data%var(iLookINDEX%ixSnowOnlyNrg)%dat            ,&  ! intent(in): [i4b(:)] indices for energy states in the snow subdomain
  ixSnowOnlyWat           => indx_data%var(iLookINDEX%ixSnowOnlyWat)%dat            ,&  ! intent(in): [i4b(:)] indices for total water states in the snow subdomain
  ixSoilOnlyNrg           => indx_data%var(iLookINDEX%ixSoilOnlyNrg)%dat            ,&  ! intent(in): [i4b(:)] indices for energy states in the soil subdomain
  ixSoilOnlyHyd           => indx_data%var(iLookINDEX%ixSoilOnlyHyd)%dat            ,&  ! intent(in): [i4b(:)] indices for hydrology states in the soil subdomain
  layerType               => indx_data%var(iLookINDEX%layerType)%dat                ,&  ! intent(in): [i4b(:)] type of each layer in the snow+soil domain (snow or soil)

  ! model fluxes
  mLayerBaseflow          => flux_data%var(iLookFLUX%mLayerBaseflow)%dat,            &  ! intent(out): [dp(:)] baseflow from each soil layer (m s-1)
  iLayerLiqFluxSoil       => flux_data%var(iLookFLUX%iLayerLiqFluxSoil)%dat,         &  ! intent(out): [dp] liquid soil fluxes (m s-1)
  scalarSoilCompress      => diag_data%var(iLookDIAG%scalarSoilCompress)%dat(1),     &  ! intent(out): [dp] total change in storage associated with compression of the soil matrix (kg m-2)
  mLayerCompress          => diag_data%var(iLookDIAG%mLayerCompress)%dat,            &  ! intent(out): [dp(:)] change in storage associated with compression of the soil matrix (-)

  ! model state variables (vegetation canopy)
  scalarCanairTemp        => prog_data%var(iLookPROG%scalarCanairTemp)%dat(1)       ,&  ! intent(inout): [dp] temperature of the canopy air space (K)
  scalarCanopyTemp        => prog_data%var(iLookPROG%scalarCanopyTemp)%dat(1)       ,&  ! intent(inout): [dp] temperature of the vegetation canopy (K)
  scalarCanopyIce         => prog_data%var(iLookPROG%scalarCanopyIce)%dat(1)        ,&  ! intent(inout): [dp] mass of ice on the vegetation canopy (kg m-2)
  scalarCanopyLiq         => prog_data%var(iLookPROG%scalarCanopyLiq)%dat(1)        ,&  ! intent(inout): [dp] mass of liquid water on the vegetation canopy (kg m-2)
  scalarCanopyWat         => prog_data%var(iLookPROG%scalarCanopyWat)%dat(1)        ,&  ! intent(inout): [dp] mass of total water on the vegetation canopy (kg m-2)

  ! layer depth
  mLayerDepth             => prog_data%var(iLookPROG%mLayerDepth)%dat               ,&  ! intent(in): [dp(:)] depth of each layer in the snow-soil sub-domain (m)

  ! model state variables (snow and soil domains)
  mLayerTemp              => prog_data%var(iLookPROG%mLayerTemp)%dat                ,&  ! intent(inout): [dp(:)] temperature of each snow/soil layer (K)
  mLayerVolFracIce        => prog_data%var(iLookPROG%mLayerVolFracIce)%dat          ,&  ! intent(inout): [dp(:)] volumetric fraction of ice (-)
  mLayerVolFracLiq        => prog_data%var(iLookPROG%mLayerVolFracLiq)%dat          ,&  ! intent(inout): [dp(:)] volumetric fraction of liquid water (-)
  mLayerVolFracWat        => prog_data%var(iLookPROG%mLayerVolFracWat)%dat          ,&  ! intent(inout): [dp(:)] volumetric fraction of total water (-)
  mLayerMatricHead        => prog_data%var(iLookPROG%mLayerMatricHead)%dat          ,&  ! intent(inout): [dp(:)] matric head (m)
  scalarAquiferStorage    => prog_data%var(iLookPROG%scalarAquiferStorage)%dat(1)    &  ! intent(inout): [dp   ] aquifer storage (m)
  )

  ! -----
  ! * compute model fluxes...
  ! -------------------------

  ! NOTE: Need to increase cleverness and avoid copying vectors
  !  --> can we do this as an associate statement?

  ! extract the vegetation states from the state vector
  if(computeVegFlux)then
   scalarCanairTempTrial = stateVec(ixCasNrg)
   scalarCanopyTempTrial = stateVec(ixVegNrg)
   scalarCanopyWatTrial  = stateVec(ixVegWat)

  ! ensure that input values to flux routines are defined
  else
   scalarCanairTempTrial = scalarCanairTemp
   scalarCanopyTempTrial = scalarCanopyTemp
   scalarCanopyWatTrial  = scalarCanopyIce + scalarCanopyLiq
  endif

  if(printFlag)then
   write(*,'(a,1x,f20.15)') 'scalarCanairTempTrial = ', scalarCanairTempTrial
   write(*,'(a,1x,f20.15)') 'scalarCanopyTempTrial = ', scalarCanopyTempTrial
   write(*,'(a,1x,f20.15)') 'scalarCanopyWatTrial  = ', scalarCanopyWatTrial
  endif

  ! extract state variables for the snow and soil domain
  mLayerTempTrial(1:nLayers)     = stateVec(ixSnowSoilNrg)
  mLayerMatricHeadTrial(1:nSoil) = stateVec(ixSoilOnlyHyd)
  if(nSnow>0)&
   mLayerVolFracWatTrial(1:nSnow) = stateVec(ixSnowOnlyWat)

  ! compute model flux for a given state vector
  call computFlux(&
                  ! input: state variables
                  scalarCanairTempTrial,              & ! intent(in): trial value for the temperature of the canopy air space (K)
                  scalarCanopyTempTrial,              & ! intent(in): trial value for the temperature of the vegetation canopy (K)
                  mLayerTempTrial,                    & ! intent(in): trial value for the temperature of each snow and soil layer (K)
                  mLayerMatricHeadTrial,              & ! intent(in): trial value for the matric head in each soil layer (m)
                  ! input: diagnostic variables defining the liquid water and ice content
                  scalarCanopyLiqLocal,               & ! intent(in): trial value for the liquid water on the vegetation canopy (kg m-2)
                  scalarCanopyIceLocal,               & ! intent(in): trial value for the ice on the vegetation canopy (kg m-2)
                  mLayerVolFracLiqLocal,              & ! intent(in): trial value for the volumetric liquid water content in each snow and soil layer (-)
                  mLayerVolFracIceLocal,              & ! intent(in): trial value for the volumetric ice in each snow and soil layer (-)
                  ! output: flux vector
                  fVec,                               & ! intent(out): flux vector (mixed units)
                  ! output: error control
                  err,cmessage)                         ! intent(out): error code and error message
  if(err/=0)then; message=trim(message)//trim(cmessage); return; endif  ! (check for errors)

  ! compute soil compressibility (-) and its derivative w.r.t. matric head (m)
  ! NOTE: we already extracted trial matrix head and volumetric liquid water as part of the flux calculations
  call soilCmpres(&
                  ! input:
                  ixRichards,                             & ! intent(in): choice of option for Richards' equation
                  mLayerMatricHead(1:nSoil),              & ! intent(in): matric head at the start of the time step (m)
                  mLayerMatricHeadTrial(1:nSoil),         & ! intent(in): trial value of matric head (m)
                  mLayerVolFracLiqLocal(nSnow+1:nLayers), & ! intent(in): trial value for the volumetric liquid water content in each soil layer (-)
                  mLayerVolFracIceLocal(nSnow+1:nLayers), & ! intent(in): trial value for the volumetric ice content in each soil layer (-)
                  dVolTot_dPsi0,                          & ! intent(in): derivative in the soil water characteristic (m-1)
                  specificStorage,                        & ! intent(in): specific storage coefficient (m-1)
                  theta_sat,                              & ! intent(in): soil porosity (-)
                  ! output:
                  mLayerCompress,                         & ! intent(out): compressibility of the soil matrix (-)
                  dCompress_dPsi,                         & ! intent(out): derivative in compressibility w.r.t. matric head (m-1)
                  err,cmessage)                             ! intent(out): error code and error message
  if(err/=0)then; message=trim(message)//trim(cmessage); return; endif  ! (check for errors)

  ! compute the total change in storage associated with compression of the soil matrix (kg m-2)
  scalarSoilCompress = sum(mLayerCompress(1:nSoil)*mLayerDepth(nSnow+1:nLayers))*iden_water

  ! -----
  ! * compute residual vector...
  ! ----------------------------

  ! intialize additional terms on the RHS as zero
  rAdd(:) = 0._dp

  ! compute energy associated with melt freeze for the vegetation canopy
  if(computeVegFlux)then
   rAdd(ixVegNrg) = rAdd(ixVegNrg) + LH_fus*(scalarCanopyIceLocal - scalarCanopyIce)/canopyDepth   ! energy associated with melt/freeze (J m-3)
  endif

  ! compute energy associated with melt/freeze for snow
  if(nSnow>0)&
  rAdd(ixSnowOnlyNrg) = rAdd(ixSnowOnlyNrg) + LH_fus*iden_ice*(mLayerVolFracIceLocal(1:nSnow) - mLayerVolFracIce(1:nSnow))       ! energy associated with melt/freeze (J m-3)

  ! compute energy associated with melt/freeze for soil
  rAdd(ixSoilOnlyNrg) = rAdd(ixSoilOnlyNrg) + LH_fus*iden_water*(mLayerVolFracIceLocal(nSnow+1:nLayers) - mLayerVolFracIce(nSnow+1:nLayers))     ! energy associated with melt/freeze (J m-3)

  ! sink terms for water (-)
  ! NOTE: state variable is volumetric water content, so melt-freeze is not included
  ! NOTE: ground evaporation was already included in the flux at the upper boundary
  ! NOTE: rAdd(ixSnowOnlyWat)=0, and is defined in the initialization above
  rAdd(ixSoilOnlyHyd)    = rAdd(ixSoilOnlyHyd) + dt*(mLayerTranspire(1:nSoil) - mLayerBaseflow(1:nSoil) )/mLayerDepth(nSnow+1:nLayers) - mLayerCompress(1:nSoil)

  ! compute the residual vector for the vegetation canopy
  ! NOTE: sMul(ixVegWat) = 1, but include as it converts all variables to quadruple precision
  if(computeVegFlux)then
   ! --> energy balance
   rVec(ixCasNrg) = sMul(ixCasNrg)*scalarCanairTempTrial - ( (sMul(ixCasNrg)*scalarCanairTemp + fVec(ixCasNrg)*dt) + rAdd(ixCasNrg) )
   rVec(ixVegNrg) = sMul(ixVegNrg)*scalarCanopyTempTrial - ( (sMul(ixVegNrg)*scalarCanopyTemp + fVec(ixVegNrg)*dt) + rAdd(ixVegNrg) )
   ! --> mass balance
   rVec(ixVegWat) = sMul(ixVegWat)*scalarCanopyWatTrial  - ( (sMul(ixVegWat)*scalarCanopyWat  + fVec(ixVegWat)*dt) + rAdd(ixVegWat) )
   if(printFlag)then
    write(*,'(a,1x,f20.10)') 'dt                    = ', dt
    write(*,'(a,1x,f20.10)') 'rVec(ixCasNrg)        = ', rVec(ixCasNrg)
    write(*,'(a,1x,f20.10)') 'fVec(ixCasNrg)        = ', fVec(ixCasNrg)
    write(*,'(a,1x,f20.10)') 'rAdd(ixCasNrg)        = ', rAdd(ixCasNrg)
    write(*,'(a,1x,f20.10)') 'sMul(ixCasNrg)        = ', sMul(ixCasNrg)
    write(*,'(a,1x,f20.10)') 'scalarCanairTemp      = ', scalarCanairTemp
    write(*,'(a,1x,f20.10)') 'scalarCanairTempTrial = ', scalarCanairTempTrial
   endif
  endif

  ! compute the residual vector for the snow and soil sub-domains for energy
  rVec(ixSnowSoilNrg) = sMul(ixSnowSoilNrg)*mLayerTempTrial(1:nLayers) - ( (sMul(ixSnowSoilNrg)*mLayerTemp(1:nLayers)  + fVec(ixSnowSoilNrg)*dt) + rAdd(ixSnowSoilNrg) )

  ! compute the residual vector for the **snow** sub-domain for liquid water
  if(nSnow>0)&
  rVec(ixSnowOnlyWat) = mLayerVolFracWatTrial(1:nSnow) - ( (mLayerVolFracWat(1:nSnow)  + fVec(ixSnowOnlyWat)*dt) + rAdd(ixSnowOnlyWat) )

  ! compute the residual vector for the **soil** sub-domain for liquid water
  vThetaInit(1:nSoil)  = mLayerVolFracLiq(nSnow+1:nLayers)      + mLayerVolFracIce(nSnow+1:nLayers)      ! liquid equivalent of total water at the start of the step
  vThetaTrial(1:nSoil) = mLayerVolFracLiqLocal(nSnow+1:nLayers) + mLayerVolFracIceLocal(nSnow+1:nLayers) ! liquid equivalent of total water at the current iteration
  rVec(ixSoilOnlyHyd)  = vThetaTrial(1:nSoil) - ( (vThetaInit(1:nSoil) + fVec(ixSoilOnlyHyd)*dt) + rAdd(ixSoilOnlyHyd) )

  if(printFlag)then
   write(*,'(a,1x,3(e20.10,1x))')  'rVec(ixSoilOnlyHyd(1)), fVec(ixSoilOnlyHyd(1)), rAdd(ixSoilOnlyHyd(1)) = ', &
                                    rVec(ixSoilOnlyHyd(1)), fVec(ixSoilOnlyHyd(1)), rAdd(ixSoilOnlyHyd(1))
  endif

  ! compute the soil water balance error (m)
  ! NOTE: declared in the main routine so accessible in all internal routines
  soilWaterBalanceError = abs( sum(real(rVec(ixSoilOnlyHyd), dp)*mLayerDepth(nSnow+1:nSoil)) )

  ! end association to variables in the data structures
  end associate
  
  if(printFlag) write(*,'(a,1x,100(e12.5,1x))') 'rVec = ', rVec 

  end subroutine xFluxResid


  ! *********************************************************************************************************
  ! internal subroutine computFlux: compute model fluxes
  ! *********************************************************************************************************
  subroutine computFlux(&
                        ! input: state variables
                        scalarCanairTempTrial,              & ! intent(in): trial value for the temperature of the canopy air space (K)
                        scalarCanopyTempTrial,              & ! intent(in): trial value for the temperature of the vegetation canopy (K)
                        mLayerTempTrial,                    & ! intent(in): trial value for the temperature of each snow and soil layer (K)
                        mLayerMatricHeadTrial,              & ! intent(in): trial value for the matric head in each soil layer (m)
                        ! input: diagnostic variables defining the liquid water and ice content
                        scalarCanopyLiqTrial,               & ! intent(in): trial value for the liquid water on the vegetation canopy (kg m-2)
                        scalarCanopyIceTrial,               & ! intent(in): trial value for the ice on the vegetation canopy (kg m-2)
                        mLayerVolFracLiqTrial,              & ! intent(in): trial value for the volumetric liquid water content in each snow and soil layer (-)
                        mLayerVolFracIceTrial,              & ! intent(in): trial value for the volumetric ice in each snow and soil layer (-)
                        ! output: flux vector
                        fluxVec,                            & ! intent(out): flux vector (mixed units)
                        ! output: error control
                        err,message)                          ! intent(out): error code and error message
  implicit none
  ! ---------------------------------------------------------------------------------------
  ! * dummy variables
  ! ---------------------------------------------------------------------------------------
  ! input: state variables
  real(dp),intent(in)            :: scalarCanairTempTrial     ! trial value for temperature of the canopy air space (K)
  real(dp),intent(in)            :: scalarCanopyTempTrial     ! trial value for temperature of the vegetation canopy (K)
  real(dp),intent(in)            :: mLayerTempTrial(:)        ! trial value for temperature of each snow/soil layer (K)
  real(dp),intent(in)            :: mLayerMatricHeadTrial(:)  ! trial value for matric head (m)
  ! input: diagnostic variables
  real(dp),intent(in)            :: scalarCanopyLiqTrial      ! trial value for mass of liquid water on the vegetation canopy (kg m-2)
  real(dp),intent(in)            :: scalarCanopyIceTrial      ! trial value for mass of ice on the vegetation canopy (kg m-2)
  real(dp),intent(in)            :: mLayerVolFracLiqTrial(:)  ! trial value for volumetric fraction of liquid water (-)
  real(dp),intent(in)            :: mLayerVolFracIceTrial(:)  ! trial value for volumetric fraction of ice (-)
  ! output: flux vector
  real(dp),intent(out)           :: fluxVec(:)                ! model flux vector (mixed units)
  ! output: error control
  integer(i4b),intent(out)       :: err                       ! error code
  character(*),intent(out)       :: message                   ! error message
  ! ---------------------------------------------------------------------------------------
  ! * local variables
  ! ---------------------------------------------------------------------------------------
  integer(i4b)                   :: iSoil                     ! index of soil layer
  character(LEN=256)             :: cmessage                  ! error message of downwind routine
  real(dp),parameter             :: canopyTempMax=500._dp     ! expected maximum value for the canopy temperature (K)
  real(dp)                       :: xNum                      ! temporary variable: numerator
  real(dp)                       :: xDen                      ! temporary variable: denominator
  real(dp)                       :: effSat                    ! effective saturation of the soil matrix (-)
  real(dp),dimension(nSoil)      :: mLayerMatricHeadLiq       ! matric head associated with liquid water (m), f(psi0, T)
  real(dp)                       :: dPsiLiq_dEffSat           ! derivative in liquid water matric potential w.r.t. effective saturation (m)
  real(dp)                       :: dEffSat_dVolTot           ! derivative in effective saturation w.r.t. total water content (-)
  real(dp)                       :: dEffSat_dTemp             ! derivative in effective saturation w.r.t. temperature (K-1)
  real(dp),dimension(nSoil)      :: dPsiLiq_dPsi0             ! derivative in the liquid water matric potential w.r.t. the total water matric potential (-)
  ! --------------------------------------------------------------
  ! initialize error control
  err=0; message='computFlux/'

  ! *****
  ! (0) PRELIMINARIES...
  ! ********************

  ! get the necessary variables for the flux computations
  associate(&

  ! model decisions
  ixGroundwater           => model_decisions(iLookDECISIONS%groundwatr)%iDecision   ,&  ! intent(in): [i4b] groundwater parameterization

  ! domain boundary conditions
  upperBoundTemp          => forc_data%var(iLookFORCE%airtemp)                      ,&  ! intent(in): [dp] temperature of the upper boundary of the snow and soil domains (K)
  scalarRainfall          => flux_data%var(iLookFLUX%scalarRainfall)%dat(1)         ,&  ! intent(in): [dp] rainfall rate (kg m-2 s-1)
  scalarSfcMeltPond       => prog_data%var(iLookPROG%scalarSfcMeltPond)%dat(1)      ,&  ! intent(in): [dp] ponded water caused by melt of the "snow without a layer" (kg m-2)

  ! layer depth
  mLayerDepth             => prog_data%var(iLookPROG%mLayerDepth)%dat               ,&  ! intent(in): [dp(:)] depth of each layer in the snow-soil sub-domain (m)

  ! indices of model state variables
  ixCasNrg                => indx_data%var(iLookINDEX%ixCasNrg)%dat(1)              , & ! intent(in): [i4b] index of canopy air space energy state variable
  ixVegNrg                => indx_data%var(iLookINDEX%ixVegNrg)%dat(1)              , & ! intent(in): [i4b] index of canopy energy state variable
  ixVegWat                => indx_data%var(iLookINDEX%ixVegWat)%dat(1)              , & ! intent(in): [i4b] index of canopy hydrology state variable (mass)
  ixTopNrg                => indx_data%var(iLookINDEX%ixTopNrg)%dat(1)              , & ! intent(in): [i4b] index of upper-most energy state in the snow-soil subdomain
  ixTopWat                => indx_data%var(iLookINDEX%ixTopWat)%dat(1)              , & ! intent(in): [i4b] index of upper-most total water state in the snow-soil subdomain
  ixTopMat                => indx_data%var(iLookINDEX%ixTopMat)%dat(1)              , & ! intent(in): [i4b] index of upper-most matric head state in the soil subdomain
  ixSnowSoilNrg           => indx_data%var(iLookINDEX%ixSnowSoilNrg)%dat            , & ! intent(in): [i4b(:)] indices for energy states in the snow-soil subdomain
  ixSnowSoilWat           => indx_data%var(iLookINDEX%ixSnowSoilWat)%dat            , & ! intent(in): [i4b(:)] indices for total water states in the snow-soil subdomain
  ixSnowOnlyNrg           => indx_data%var(iLookINDEX%ixSnowOnlyNrg)%dat            , & ! intent(in): [i4b(:)] indices for energy states in the snow subdomain
  ixSnowOnlyWat           => indx_data%var(iLookINDEX%ixSnowOnlyWat)%dat            , & ! intent(in): [i4b(:)] indices for total water states in the snow subdomain
  ixSoilOnlyNrg           => indx_data%var(iLookINDEX%ixSoilOnlyNrg)%dat            , & ! intent(in): [i4b(:)] indices for energy states in the soil subdomain
  ixSoilOnlyHyd           => indx_data%var(iLookINDEX%ixSoilOnlyHyd)%dat            ,&  ! intent(in): [i4b(:)] indices for hydrology states in the soil subdomain
  layerType               => indx_data%var(iLookINDEX%layerType)%dat                , & ! intent(in): [i4b(:)] type of each layer in the snow+soil domain (snow or soil)

  ! snow parameters
  snowfrz_scale           => mpar_data%var(iLookPARAM%snowfrz_scale)                ,&  ! intent(in): [dp] scaling parameter for the snow freezing curve (K-1)

  ! soil parameters
  vGn_m                   => diag_data%var(iLookDIAG%scalarVGn_m)%dat(1)            ,&  ! intent(in): [dp] van Genutchen "m" parameter (-)
  vGn_n                   => mpar_data%var(iLookPARAM%vGn_n)                        ,&  ! intent(in): [dp] van Genutchen "n" parameter (-)
  vGn_alpha               => mpar_data%var(iLookPARAM%vGn_alpha)                    ,&  ! intent(in): [dp] van Genutchen "alpha" parameter (m-1)
  theta_sat               => mpar_data%var(iLookPARAM%theta_sat)                    ,&  ! intent(in): [dp] soil porosity (-)
  theta_res               => mpar_data%var(iLookPARAM%theta_res)                    ,&  ! intent(in): [dp] soil residual volumetric water content (-)

  ! model diagnostic variables
  scalarThroughfallRain   => flux_data%var(iLookFLUX%scalarThroughfallRain)%dat(1)  ,&  ! intent(out): [dp] rain that reaches the ground without ever touching the canopy (kg m-2 s-1)
  scalarCanopyLiqDrainage => flux_data%var(iLookFLUX%scalarCanopyLiqDrainage)%dat(1),&  ! intent(out): [dp] drainage of liquid water from the vegetation canopy (kg m-2 s-1)
  scalarSurfaceRunoff     => flux_data%var(iLookFLUX%scalarSurfaceRunoff)%dat(1)    ,&  ! intent(out): [dp] surface runoff (m s-1)
  scalarRainPlusMelt      => flux_data%var(iLookFLUX%scalarRainPlusMelt)%dat(1)     ,&  ! intent(out): [dp] rain plus melt (m s-1)
  scalarExfiltration      => flux_data%var(iLookFLUX%scalarExfiltration)%dat(1)     ,&  ! intent(out): [dp] exfiltration from the soil profile (m s-1)
  mLayerColumnOutflow     => flux_data%var(iLookFLUX%mLayerColumnOutflow)%dat       ,&  ! intent(out): [dp(:)] column outflow from each soil layer (m3 s-1)

  ! soil fluxes
  iLayerNrgFlux           => flux_data%var(iLookFLUX%iLayerNrgFlux)%dat             ,&  ! intent(out): [dp(0:)] vertical energy flux at the interface of snow and soil layers
  iLayerLiqFluxSnow       => flux_data%var(iLookFLUX%iLayerLiqFluxSnow)%dat         ,&  ! intent(out): [dp(0:)] vertical liquid water flux at snow layer interfaces (-)
  iLayerLiqFluxSoil       => flux_data%var(iLookFLUX%iLayerLiqFluxSoil)%dat         ,&  ! intent(out): [dp(0:)] vertical liquid water flux at soil layer interfaces (-)
  mLayerBaseflow          => flux_data%var(iLookFLUX%mLayerBaseflow)%dat            ,&  ! intent(out): [dp(:)] baseflow from each soil layer (m s-1)

  ! aquifer fluxes
  scalarAquiferTranspire  => flux_data%var(iLookFLUX%scalarAquiferTranspire)%dat(1) ,&  ! intent(out): [dp] transpiration loss from the aquifer (m s-1
  scalarAquiferRecharge   => flux_data%var(iLookFLUX%scalarAquiferRecharge)%dat(1)  ,&  ! intent(out): [dp] recharge to the aquifer (m s-1)
  scalarAquiferBaseflow   => flux_data%var(iLookFLUX%scalarAquiferBaseflow)%dat(1)   &  ! intent(out): [dp] total baseflow from the aquifer (m s-1)

  )  ! association to data in structures

  ! check that canopy temperature is reasonable
  if(scalarCanopyTempTrial > canopyTempMax)then
   print*, 'scalarCanopyTempTrial = ', scalarCanopyTempTrial
   message=trim(message)//'canopy temperature is > expected maximum'
   err=20; return
  endif

  ! * vegetation domain: compute derivative of volumetric liquid water content w.r.t. temperature (K-1)
  if(computeVegFlux)then
   if(scalarCanopyIceTrial > verySmall)then
    theta = (scalarCanopyIceTrial + scalarCanopyLiqTrial)/(canopyDepth*iden_water)
    dTheta_dTkCanopy = dFracLiq_dTk(scalarCanopyTempTrial,snowfrz_scale)*theta   ! K-1
    dCanLiq_dTcanopy = dTheta_dTkCanopy*iden_water*canopyDepth                   ! kg m-2 K-1
   else
    dTheta_dTkCanopy = 0._dp
    dCanLiq_dTcanopy = 0._dp
   endif
  endif

  ! * snow+soil domain: compute derivative of volumetric liquid water content w.r.t. temperature (K-1)
  do iLayer=1,nLayers  ! loop through all snow and soil layers
   select case(layerType(iLayer))
    case(ix_snow) ! (snow layers)
     theta = mLayerVolFracIceTrial(iLayer)*(iden_ice/iden_water) + mLayerVolFracLiqTrial(iLayer)
     mLayerdTheta_dTk(iLayer) = dFracLiq_dTk(mLayerTempTrial(iLayer),snowfrz_scale)*theta
    case(ix_soil) ! (soil layers)
     if(mLayerVolFracIceTrial(iLayer)>verySmall)then
      mLayerdTheta_dTk(iLayer)        = dTheta_dTk(mLayerTempTrial(iLayer),theta_res,theta_sat,vGn_alpha,vGn_n,vGn_m)  ! assume no volume expansion
     else
      mLayerdTheta_dTk(iLayer)        = 0._dp
     endif
    case default; err=40; message=trim(message)//"cannot identify the layer as snow or soil"; return
   endselect
  end do  ! (looping through snow+soil layers)

  ! * compute the matric head associated with liquid water
  do iSoil=1,nSoil  ! loop through soil layers

   ! - compute derivative in total water content w.r.t. total water matric potential (m-1)
   dVolTot_dPsi0(iSoil) = dTheta_dPsi(mLayerMatricHeadTrial(iSoil),vGn_alpha,theta_res,theta_sat,vGn_n,vGn_m)  ! valid for both frozen and unfrozen conditions

   ! ** partially frozen soil
   if(mLayerVolFracIceTrial(nSnow+iSoil) > verySmall .and. mLayerMatricHeadTrial(iSoil) < 0._dp)then  ! check that ice exists and that the soil is unsaturated
    ! - compute effective saturation
    ! NOTE: include ice content as part of the solid porosity - major effect of ice is to reduce the pore size; ensure that effSat=1 at saturation
    ! (from Zhao et al., J. Hydrol., 1997: Numerical analysis of simultaneous heat and mass transfer...)
    xNum   = mLayerVolFracLiqTrial(nSnow+iSoil) - theta_res
    xDen   = theta_sat - mLayerVolFracIceTrial(nSnow+iSoil) - theta_res
    effSat = xNum/xDen          ! effective saturation
    ! - matric head associated with liquid water
    mLayerMatricHeadLiq(iSoil) = matricHead(effSat,vGn_alpha,0._dp,1._dp,vGn_n,vGn_m)  ! argument is effective saturation, so theta_res=0 and theta_sat=1
    ! - compute derivative in the liquid water matric potential w.r.t. the total water matric potential
    dPsiLiq_dEffSat      = dPsi_dTheta(effSat,vGn_alpha,0._dp,1._dp,vGn_n,vGn_m) ! derivative in liquid water matric potential w.r.t. effective saturation (m)
    dEffSat_dVolTot      = xNum/(xDen**2._dp) ! derivative in effective saturation w.r.t. total water content (-)
    dPsiLiq_dPsi0(iSoil) = dVolTot_dPsi0(iSoil)*dPsiLiq_dEffSat*dEffSat_dVolTot
    ! compute the derivative in the liquid water matric potential w.r.t. temperature (m K-1)
    dEffSat_dTemp        = -mLayerdTheta_dTk(nSnow+iSoil)*xNum/(xDen**2._dp) + mLayerdTheta_dTk(nSnow+iSoil)/xDen
    dPsiLiq_dTemp(iSoil) = dPsiLiq_dEffSat*dEffSat_dTemp
   ! ** unfrozen soil
   else   ! (no ice)
    dPsiLiq_dPsi0(iSoil)       = 1._dp  ! derivative=1 because values are identical
    dPsiLiq_dTemp(iSoil)       = 0._dp  ! derivative=0 because no impact of temperature for unfrozen conditions
    mLayerMatricHeadLiq(iSoil) = mLayerMatricHeadTrial(iSoil) ! liquid water matric potential is equal to the total water matic potential when there is no ice
   endif  ! (if ice exists)

  end do  ! (looping through soil layers)

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
                  firstSubStep,                           & ! intent(in): flag to indicate if we are processing the first sub-step
                  firstFluxCall,                          & ! intent(in): flag to indicate if we are processing the first flux call
                  computeVegFlux,                         & ! intent(in): flag to indicate if we need to compute fluxes over vegetation
                  ! input: model state variables
                  upperBoundTemp,                         & ! intent(in): temperature of the upper boundary (K) --> NOTE: use air temperature
                  scalarCanairTempTrial,                  & ! intent(in): trial value of the canopy air space temperature (K)
                  scalarCanopyTempTrial,                  & ! intent(in): trial value of canopy temperature (K)
                  mLayerTempTrial(1),                     & ! intent(in): trial value of ground temperature (K)
                  scalarCanopyIceTrial,                   & ! intent(in): trial value of mass of ice on the vegetation canopy (kg m-2)
                  scalarCanopyLiqTrial,                   & ! intent(in): trial value of mass of liquid water on the vegetation canopy (kg m-2)
                  ! input: model derivatives
                  dCanLiq_dTcanopy,                       & ! intent(in): derivative in canopy liquid storage w.r.t. canopy temperature (kg m-2 K-1)
                  ! input/output: data structures
                  type_data,                              & ! intent(in):    type of vegetation and soil
                  attr_data,                              & ! intent(in):    spatial attributes
                  forc_data,                              & ! intent(in):    model forcing data
                  mpar_data,                              & ! intent(in):    model parameters
                  indx_data,                              & ! intent(in):    index data
                  prog_data,                              & ! intent(in):    model prognostic variables for a local HRU
                  diag_data,                              & ! intent(inout): model diagnostic variables for a local HRU
                  flux_data,                              & ! intent(inout): model fluxes for a local HRU
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
                  ! output: liquid water flux derivarives (canopy evap)
                  dCanopyEvaporation_dCanLiq,             & ! intent(out): derivative in canopy evaporation w.r.t. canopy liquid water content (s-1)
                  dCanopyEvaporation_dTCanair,            & ! intent(out): derivative in canopy evaporation w.r.t. canopy air temperature (kg m-2 s-1 K-1)
                  dCanopyEvaporation_dTCanopy,            & ! intent(out): derivative in canopy evaporation w.r.t. canopy temperature (kg m-2 s-1 K-1)
                  dCanopyEvaporation_dTGround,            & ! intent(out): derivative in canopy evaporation w.r.t. ground temperature (kg m-2 s-1 K-1)
                  ! output: liquid water flux derivarives (ground evap)
                  dGroundEvaporation_dCanLiq,             & ! intent(out): derivative in ground evaporation w.r.t. canopy liquid water content (s-1)
                  dGroundEvaporation_dTCanair,            & ! intent(out): derivative in ground evaporation w.r.t. canopy air temperature (kg m-2 s-1 K-1)
                  dGroundEvaporation_dTCanopy,            & ! intent(out): derivative in ground evaporation w.r.t. canopy temperature (kg m-2 s-1 K-1)
                  dGroundEvaporation_dTGround,            & ! intent(out): derivative in ground evaporation w.r.t. ground temperature (kg m-2 s-1 K-1)
                  ! output: cross derivative terms
                  dCanopyNetFlux_dCanLiq,                 & ! intent(out): derivative in net canopy fluxes w.r.t. canopy liquid water content (J kg-1 s-1)
                  dGroundNetFlux_dCanLiq,                 & ! intent(out): derivative in net ground fluxes w.r.t. canopy liquid water content (J kg-1 s-1)
                  ! output: error control
                  err,cmessage)                             ! intent(out): error control
  if(err/=0)then; message=trim(message)//trim(cmessage); return; endif  ! (check for errors)

  ! check fluxes
  if(printFlag)then
   write(*,'(a,1x,f30.20)') 'canairNetNrgFlux = ', canairNetNrgFlux
   write(*,'(a,1x,f30.20)') 'canopyNetNrgFlux = ', canopyNetNrgFlux
   write(*,'(a,1x,f30.20)') 'groundNetNrgFlux = ', groundNetNrgFlux
   write(*,'(a,1x,f30.20)') 'dGroundNetFlux_dGroundTemp = ', dGroundNetFlux_dGroundTemp
  endif

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
                  ! input-output: data structures
                  mpar_data,                              & ! intent(in):    model parameters
                  indx_data,                              & ! intent(in):    model indices
                  prog_data,                              & ! intent(in):    model prognostic variables for a local HRU
                  diag_data,                              & ! intent(in):    model diagnostic variables for a local HRU
                  flux_data,                              & ! intent(inout): model fluxes for a local HRU
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
   if(printFlag)then
    if(iLayer < 3) write(*,'(a,1x,i4,1x,10(f25.15,1x))') 'iLayer, iLayerNrgFlux(iLayer-1:iLayer), ssdNetNrgFlux(iLayer)   = ', iLayer, iLayerNrgFlux(iLayer-1:iLayer), ssdNetNrgFlux(iLayer)
   endif
  end do

  ! *****
  ! (3) CALCULATE THE LIQUID FLUX THROUGH VEGETATION...
  ! ***************************************************
  call vegLiqFlux(&
                  ! input
                  computeVegFlux,                         & ! intent(in): flag to denote if computing energy flux over vegetation
                  scalarCanopyLiqTrial,                   & ! intent(in): trial mass of liquid water on the vegetation canopy at the current iteration (kg m-2)
                  scalarRainfall,                         & ! intent(in): rainfall rate (kg m-2 s-1)
                  ! input-output: data structures
                  mpar_data,                              & ! intent(in): model parameters
                  diag_data,                              & ! intent(in): local HRU diagnostic model variables
                  ! output
                  scalarThroughfallRain,                  & ! intent(out): rain that reaches the ground without ever touching the canopy (kg m-2 s-1)
                  scalarCanopyLiqDrainage,                & ! intent(out): drainage of liquid water from the vegetation canopy (kg m-2 s-1)
                  scalarThroughfallRainDeriv,             & ! intent(out): derivative in throughfall w.r.t. canopy liquid water (s-1)
                  scalarCanopyLiqDrainageDeriv,           & ! intent(out): derivative in canopy drainage w.r.t. canopy liquid water (s-1)
                  err,cmessage)                             ! intent(out): error control
  if(err/=0)then; message=trim(message)//trim(cmessage); return; endif
  ! calculate the net liquid water flux for the vegetation canopy
  canopyNetLiqFlux = scalarRainfall + scalarCanopyEvaporation - scalarThroughfallRain - scalarCanopyLiqDrainage
  ! calculate the total derivative in the downward liquid flux
  scalarCanopyLiqDeriv = scalarThroughfallRainDeriv + scalarCanopyLiqDrainageDeriv
  ! test
  if(printFlag)then
   print*, 'scalarRainfall = ', scalarRainfall
   print*, 'scalarThroughfallRain   = ', scalarThroughfallRain
   print*, 'scalarCanopyEvaporation = ', scalarCanopyEvaporation
   print*, 'scalarCanopyLiqDrainage = ', scalarCanopyLiqDrainage
  endif

  ! *****
  ! (4) CALCULATE THE LIQUID FLUX THROUGH SNOW...
  ! *********************************************

  if(nSnow > 0)then
   ! compute liquid fluxes
   call snowLiqFlx(&
                   ! input: model control
                   nSnow,                                 & ! intent(in): number of snow layers
                   firstFluxCall,                         & ! intent(in): the first flux call (compute variables that are constant over the iterations)
                   ! input: forcing for the snow domain
                   scalarThroughfallRain,                 & ! intent(in): rain that reaches the snow surface without ever touching vegetation (kg m-2 s-1)
                   scalarCanopyLiqDrainage,               & ! intent(in): liquid drainage from the vegetation canopy (kg m-2 s-1)
                   ! input: model state vector
                   mLayerVolFracLiqTrial(1:nSnow),        & ! intent(in): trial value of volumetric fraction of liquid water at the current iteration (-)
                   ! input-output: data structures
                   mpar_data,                             & ! intent(in):    model parameters
                   prog_data,                             & ! intent(in):    model prognostic variables for a local HRU
                   diag_data,                             & ! intent(inout): model diagnostic variables for a local HRU
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
   scalarRainPlusMelt = (scalarThroughfallRain + scalarCanopyLiqDrainage)/iden_water &  ! liquid flux from the canopy (m s-1)
                         + (scalarSfcMeltPond/dt)/iden_water  ! melt of the snow without a layer (m s-1)
  endif

  ! *****
  ! (5) CALCULATE THE LIQUID FLUX THROUGH SOIL...
  ! *********************************************
  call soilLiqFlx(&
                  ! input: model control
                  nSoil,                                  & ! intent(in): number of soil layers
                  firstFluxCall,                          & ! intent(in): flag indicating first call
                  .true.,                                 & ! intent(in): flag indicating if derivatives are desired
                  ! input: trial state variables
                  mLayerTempTrial(nSnow+1:nLayers),       & ! intent(in): trial temperature at the current iteration (K)
                  mLayerMatricHeadLiq(1:nSoil),           & ! intent(in): liquid water matric potential (m)
                  mLayerVolFracLiqTrial(nSnow+1:nLayers), & ! intent(in): volumetric fraction of liquid water (-)
                  mLayerVolFracIceTrial(nSnow+1:nLayers), & ! intent(in): volumetric fraction of ice (-)
                  ! input: pre-computed deriavatives
                  mLayerdTheta_dTk(nSnow+1:nLayers),      & ! intent(in): derivative in volumetric liquid water content w.r.t. temperature (K-1)
                  dPsiLiq_dTemp(1:nSoil),                 & ! intent(in): derivative in liquid water matric potential w.r.t. temperature (m K-1)
                  ! input: fluxes
                  scalarCanopyTranspiration,              & ! intent(in): canopy transpiration (kg m-2 s-1)
                  scalarGroundEvaporation,                & ! intent(in): ground evaporation (kg m-2 s-1)
                  scalarRainPlusMelt,                     & ! intent(in): rain plus melt (m s-1)
                  ! input-output: data structures
                  mpar_data,                              & ! intent(in):    model parameters
                  indx_data,                              & ! intent(in):    model indices
                  prog_data,                              & ! intent(in):    model prognostic variables for a local HRU
                  diag_data,                              & ! intent(in):    model diagnostic variables for a local HRU
                  flux_data,                              & ! intent(in):    model fluxes for a local HRU
                  ! output: diagnostic variables for surface runoff
                  xMaxInfilRate,                          & ! intent(inout): maximum infiltration rate (m s-1)
                  scalarInfilArea,                        & ! intent(inout): fraction of unfrozen area where water can infiltrate (-)
                  scalarFrozenArea,                       & ! intent(inout): fraction of area that is considered impermeable due to soil ice (-)
                  scalarSurfaceRunoff,                    & ! intent(out): surface runoff (m s-1)
                  ! output: diagnostic variables for model layers
                  mLayerdTheta_dPsi,                      & ! intent(out): derivative in the soil water characteristic w.r.t. psi (m-1)
                  mLayerdPsi_dTheta,                      & ! intent(out): derivative in the soil water characteristic w.r.t. theta (m)
                  dHydCond_dMatric,                       & ! intent(out): derivative in hydraulic conductivity w.r.t matric head (s-1)
                  ! output: fluxes
                  scalarSurfaceInfiltration,              & ! intent(out): surface infiltration rate (m s-1) -- only computed for iter==1
                  iLayerLiqFluxSoil,                      & ! intent(out): liquid fluxes at layer interfaces (m s-1)
                  mLayerTranspire,                        & ! intent(out): transpiration loss from each soil layer (m s-1)
                  mLayerHydCond,                          & ! intent(out): hydraulic conductivity in each layer (m s-1)
                  ! output: derivatives in fluxes w.r.t. state variables -- matric head or volumetric lquid water -- in the layer above and layer below (m s-1 or s-1)
                  dq_dHydStateAbove,                      & ! intent(out): derivatives in the flux w.r.t. matric head in the layer above (s-1)
                  dq_dHydStateBelow,                      & ! intent(out): derivatives in the flux w.r.t. matric head in the layer below (s-1)
                  ! output: derivatives in fluxes w.r.t. energy state variables -- now just temperature -- in the layer above and layer below (m s-1 K-1)
                  dq_dNrgStateAbove,                      & ! intent(out): derivatives in the flux w.r.t. temperature in the layer above (m s-1 K-1)
                  dq_dNrgStateBelow,                      & ! intent(out): derivatives in the flux w.r.t. temperature in the layer below (m s-1 K-1)
                  ! output: error control
                  err,cmessage)                             ! intent(out): error control
  if(err/=0)then; message=trim(message)//trim(cmessage); return; endif
  ! calculate net liquid water fluxes for each soil layer (s-1)
  do iLayer=1,nSoil
   soilNetLiqFlux(iLayer) = -(iLayerLiqFluxSoil(iLayer) - iLayerLiqFluxSoil(iLayer-1))/mLayerDepth(iLayer+nSnow)
  end do
  ! calculate the soil control on infiltration
  if(nSnow==0) then
   ! * case of infiltration into soil
   if(xMaxInfilRate > scalarRainPlusMelt)then  ! infiltration is not rate-limited
    soilControl = (1._dp - scalarFrozenArea)*scalarInfilArea
   else
    soilControl = 0._dp  ! (scalarRainPlusMelt exceeds maximum infiltration rate
   endif
  else
   ! * case of infiltration into snow
   soilControl = 1._dp
  endif

  ! expand derivatives to the total water matric potential
  dq_dHydStateAbove(1:nSoil)   = dq_dHydStateAbove(1:nSoil)  *dPsiLiq_dPsi0(1:nSoil)
  dq_dHydStateBelow(0:nSoil-1) = dq_dHydStateBelow(0:nSoil-1)*dPsiLiq_dPsi0(1:nSoil)

  ! *****
  ! (6) CALCULATE THE GROUNDWATER FLOW...
  ! *************************************

  ! set baseflow fluxes to zero if the baseflow routine is not used
  if(local_ixGroundwater/=qbaseTopmodel)then
   ! (diagnostic variables in the data structures)
   scalarExfiltration     = 0._dp  ! exfiltration from the soil profile (m s-1)
   mLayerColumnOutflow(:) = 0._dp  ! column outflow from each soil layer (m3 s-1)
   ! (variables needed for the numerical solution)
   mLayerBaseflow(:)      = 0._dp  ! baseflow from each soil layer (m s-1)

  ! compute the basdeflow flux
  else ! local_ixGroundwater==qbaseTopmodel
   call groundwatr(&
                   ! input: model control
                   nSnow,                                   & ! intent(in):    number of snow layers
                   nSoil,                                   & ! intent(in):    number of soil layers
                   nLayers,                                 & ! intent(in):    total number of layers
                   firstFluxCall,                           & ! intent(in):    logical flag to compute index of the lowest saturated layer
                   ! input: state and diagnostic variables
                   mLayerdTheta_dPsi,                       & ! intent(in):    derivative in the soil water characteristic w.r.t. matric head in each layer (m-1)
                   mLayerMatricHeadLiq,                     & ! intent(in):    liquid water matric potential (m)
                   mLayerVolFracLiqTrial(nSnow+1:nLayers),  & ! intent(in):    volumetric fraction of liquid water (-)
                   mLayerVolFracIceTrial(nSnow+1:nLayers),  & ! intent(in):    volumetric fraction of ice (-)
                   ! input: data structures
                   attr_data,                               & ! intent(in):    model attributes
                   mpar_data,                               & ! intent(in):    model parameters
                   prog_data,                               & ! intent(in):    model prognostic variables for a local HRU
                   diag_data,                               & ! intent(in):    model diagnostic variables for a local HRU
                   flux_data,                               & ! intent(inout): model fluxes for a local HRU
                   ! output
                   ixSaturation,                            & ! intent(inout) index of lowest saturated layer (NOTE: only computed on the first iteration)
                   mLayerBaseflow,                          & ! intent(out): baseflow from each soil layer (m s-1)
                   dBaseflow_dMatric,                       & ! intent(out): derivative in baseflow w.r.t. matric head (s-1)
                   err,cmessage)                              ! intent(out): error control
   if(err/=0)then; message=trim(message)//trim(cmessage); return; endif
  endif  ! computing baseflow flux

  ! *****
  ! (7) CALCUALTE FLUXES FOR THE DEEP AQUIFER...
  ! ********************************************

  ! identify modeling decision
  if(ixGroundwater==bigBucket)then
   ! deep aquifer is not yet transfered from old code structure
   message=trim(message)//'bigBucket groundwater parameterization is not yet transfered from old code structure'
   err=20; return
  else
   ! if no quifer, then fluxes are zero
   scalarAquiferTranspire = 0._dp  ! transpiration loss from the aquifer (m s-1
   scalarAquiferRecharge  = 0._dp  ! recharge to the aquifer (m s-1)
   scalarAquiferBaseflow  = 0._dp  ! total baseflow from the aquifer (m s-1)
  endif

  ! *****
  ! (X) WRAP UP...
  ! **************

  ! define model flux vector for the vegetation sub-domain
  if(computeVegFlux)then
   fluxVec(ixCasNrg) = canairNetNrgFlux/canopyDepth
   fluxVec(ixVegNrg) = canopyNetNrgFlux/canopyDepth
   fluxVec(ixVegWat) = canopyNetLiqFlux   ! NOTE: solid fluxes are handled separately
  endif

  ! define the model flux vector for the snow and soil sub-domains
  fluxVec(ixSnowSoilNrg) = ssdNetNrgFlux(1:nLayers)
  fluxVec(ixSoilOnlyHyd) = soilNetLiqFlux(1:nSoil)
  if(nSnow>0)&
  fluxVec(ixSnowOnlyWat) = snowNetLiqFlux(1:nSnow)

  ! end association to variables in the data structures
  end associate

  ! set the first flux call to false
  firstFluxCall=.false.

  end subroutine computFlux


  ! *********************************************************************************************************
  ! internal subroutine cpactBand: compute the compact band-diagonal matric
  ! *********************************************************************************************************
  subroutine cpactBand(err,message)
  ! dummy variables
  integer(i4b),intent(out)       :: err                     ! error code
  character(*),intent(out)       :: message                 ! error message
  ! --------------------------------------------------------------
  ! initialize error control
  err=0; message='cpactBand/'

  ! associate variables from data structures
  associate(&

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

            ! layer depth
            mLayerDepth   => prog_data%var(iLookPROG%mLayerDepth)%dat) ! intent(in): [dp(:)] depth of each layer in the snow-soil sub-domain (m)

  ! initialize the Jacobian
  ! NOTE: this needs to be done every time, since Jacobian matrix is modified in the solver
  aJac(:,:) = 0._dp  ! analytical Jacobian matrix

  ! -----
  ! * energy and liquid fluxes over vegetation...
  ! ---------------------------------------------
  if(computeVegFlux)then  ! (derivatives only defined when vegetation protrudes over the surface)

   ! liquid water fluxes for vegetation canopy (-)
   aJac(ixDiag,ixVegWat) = -fracLiqVeg*(dCanopyEvaporation_dCanLiq - scalarCanopyLiqDeriv)*dt + 1._dp     ! ixVegWat: CORRECT

   ! cross-derivative terms w.r.t. system temperatures (kg m-2 K-1)
   aJac(ixSub2,ixCasNrg) = -dCanopyEvaporation_dTCanair*dt                                                        ! ixCasNrg: CORRECT
   aJac(ixSub1,ixVegNrg) = -dCanopyEvaporation_dTCanopy*dt + dt*scalarCanopyLiqDeriv*dCanLiq_dTcanopy     ! ixVegNrg: CORRECT
   aJac(ixSup1,ixTopNrg) = -dCanopyEvaporation_dTGround*dt                                                        ! ixTopNrg: CORRECT

   ! cross-derivative terms w.r.t. canopy water (kg-1 m2)
   aJac(ixSub2,ixVegWat) = (dt/mLayerDepth(1))*(-soilControl*fracLiqVeg*scalarCanopyLiqDeriv)/iden_water  ! ixVegWat: CORRECT

   ! cross-derivative terms w.r.t. canopy temperature (K-1)
   aJac(ixSub3,ixVegNrg) = (dt/mLayerDepth(1))*(-soilControl*scalarCanopyLiqDeriv*dCanLiq_dTcanopy)/iden_water    ! ixVegNrg: CORRECT
   !print*, 'soilControl, scalarCanopyLiqDeriv, dCanLiq_dTcanopy = ', soilControl, scalarCanopyLiqDeriv, dCanLiq_dTcanopy

   ! cross-derivative terms w.r.t. canopy liquid water (J m-1 kg-1)
   ! NOTE: dIce/dLiq = (1 - fracLiqVeg); dIce*LH_fus/canopyDepth = J m-3; dLiq = kg m-2
   aJac(ixSup1,ixVegWat) = (dt/canopyDepth)   *(-dCanopyNetFlux_dCanLiq) - (1._dp - fracLiqVeg)*LH_fus/canopyDepth   ! dF/dLiq    ! ixVegWat: CORRECT
   aJac(ixSub1,ixVegWat) = (dt/mLayerDepth(1))*(-dGroundNetFlux_dCanLiq)                                          ! ixVegWat: CORRECT

   ! energy fluxes with the canopy air space (J m-3 K-1)
   aJac(ixDiag,ixCasNrg) = (dt/canopyDepth)*(-dCanairNetFlux_dCanairTemp) + dMat(ixCasNrg)                        ! ixCasNrg: CORRECT
   aJac(ixSup1,ixVegNrg) = (dt/canopyDepth)*(-dCanairNetFlux_dCanopyTemp)                                         ! ixVegNrg: CORRECT
   aJac(ixSup3,ixTopNrg) = (dt/canopyDepth)*(-dCanairNetFlux_dGroundTemp)                                         ! ixTopNrg: CORRECT

   ! energy fluxes with the vegetation canopy (J m-3 K-1)
   aJac(ixSub1,ixCasNrg) = (dt/canopyDepth)*(-dCanopyNetFlux_dCanairTemp)                                         ! ixCasNrg: CORRECT
   aJac(ixDiag,ixVegNrg) = (dt/canopyDepth)*(-dCanopyNetFlux_dCanopyTemp) + dMat(ixVegNrg)                        ! ixVegNrg: CORRECT
   aJac(ixSup2,ixTopNrg) = (dt/canopyDepth)*(-dCanopyNetFlux_dGroundTemp)                                         ! ixTopNrg: CORRECT

   ! energy fluxes with the surface (J m-3 K-1)
   aJac(ixSub3,ixCasNrg) = (dt/mLayerDepth(1))*(-dGroundNetFlux_dCanairTemp)                                      ! ixCasNrg: CORRECT
   aJac(ixSub2,ixVegNrg) = (dt/mLayerDepth(1))*(-dGroundNetFlux_dCanopyTemp)                                      ! ixVegNrg: CORRECT

  endif  ! if there is a need to compute energy fluxes within vegetation

  ! -----
  ! * energy fluxes for the snow-soil domain...
  ! -------------------------------------------
  do iLayer=1,nLayers  ! loop through layers in the snow-soil domain
   ! (define layer indices)
   jLayer = ixSnowSoilNrg(iLayer)   ! layer index within the full state vector
   ! (define the compact band-diagonal matrix)
   if(iLayer > 1)       aJac(ixSup2,jLayer) = (dt/mLayerDepth(iLayer-1))*( dNrgFlux_dTempBelow(iLayer-1) )
                        aJac(ixDiag,jLayer) = (dt/mLayerDepth(iLayer))  *(-dNrgFlux_dTempBelow(iLayer-1) + dNrgFlux_dTempAbove(iLayer)) + dMat(jLayer)
   if(iLayer < nLayers) aJac(ixSub2,jLayer) = (dt/mLayerDepth(iLayer+1))*(-dNrgFlux_dTempAbove(iLayer  ) )
  end do  ! (looping through layers in the snow-soil system)

  ! -----
  ! * liquid water fluxes for the snow domain...
  ! --------------------------------------------
  do iLayer=1,nSnow
   ! - define layer indices
   jLayer = ixSnowOnlyWat(iLayer)   ! layer index within the full state vector
   mLayer = ixSnowSoilNrg(iLayer)   ! energy layer index within the full state vector
   ! - compute the diagonal
   aJac(ixDiag,jLayer) = (dt/mLayerDepth(iLayer))*iLayerLiqFluxSnowDeriv(iLayer)*fracLiqSnow(iLayer) + dMat(jLayer)
   ! - compute cross-derivative terms for the current layer
   ! NOTE: increase in volumetric liquid water content balanced by a decrease in volumetric ice content
   aJac(ixSub1,mLayer) = (dt/mLayerDepth(iLayer))*iLayerLiqFluxSnowDeriv(iLayer)*mLayerdTheta_dTk(iLayer)  ! (dVol/dT)
   aJac(ixSup1,jLayer) = -(1._dp - fracLiqSnow(iLayer))*LH_fus*iden_water     ! (dF/dLiq)
   ! - compute cross-derivative terms for the layer below (w.r.t. state in the current layer)
   if(iLayer < nSnow)then
    aJac(ixSub3,mLayer) = -(dt/mLayerDepth(iLayer+1))*iLayerLiqFluxSnowDeriv(iLayer)*mLayerdTheta_dTk(iLayer)        ! dVol(below)/dT(above) -- K-1
    aJac(ixSub2,jLayer) = (dt/mLayerDepth(iLayer+1))*iLayerLiqFluxSnowDeriv(iLayer)*fracLiqSnow(iLayer)              ! dVol(below)/dLiq(above) -- (-)
   endif
  end do  ! (looping through snow layers)

  ! -----
  ! * liquid water fluxes for the soil domain...
  ! --------------------------------------------
  do iLayer=1,nSoil    ! loop through layers in the soil domain
   ! - define layer indices
   jLayer = ixSoilOnlyHyd(iLayer)  ! layer index within the full state vector
   kLayer = iLayer+nSnow           ! layer index within the full snow-soil vector
   ! - compute the Jacobian
   if(kLayer > nSnow+1) aJac(ixSup2,jLayer) = (dt/mLayerDepth(kLayer-1))*( dq_dHydStateBelow(iLayer-1))
                        aJac(ixDiag,jLayer) = (dt/mLayerDepth(kLayer))  *(-dq_dHydStateBelow(iLayer-1) + dq_dHydStateAbove(iLayer)) + dMat(jLayer)
   if(kLayer < nLayers) aJac(ixSub2,jLayer) = (dt/mLayerDepth(kLayer+1))*(-dq_dHydStateAbove(iLayer))
  end do  ! (looping through soil layers)

  ! -----
  ! * derivative in liquid water fluxes w.r.t. temperature for the soil domain...
  ! -----------------------------------------------------------------------------
  do iLayer=1,nSoil    ! loop through layers in the soil domain

   ! - define layer indices
   kLayer = iLayer+nSnow                ! layer index within the full snow-soil vector
   jLayer = ixSoilOnlyHyd(iLayer)       ! hydrology layer index within the full state vector
   mLayer = ixSnowSoilNrg(kLayer)       ! thermodynamics layer index within the full state vector

   ! - compute the Jacobian for the layer itself
   aJac(ixSub1,mLayer) = (dt/mLayerDepth(kLayer))*(-dq_dNrgStateBelow(iLayer-1) + dq_dNrgStateAbove(iLayer))   ! dVol/dT (K-1) -- flux depends on ice impedance

   ! - include derivatives w.r.t. ground evaporation
   if(nSnow==0 .and. iLayer==1)then  ! upper-most soil layer
    if(computeVegFlux)then
     aJac(ixSub4,ixCasNrg) = (dt/mLayerDepth(kLayer))*(-dGroundEvaporation_dTCanair/iden_water) ! dVol/dT (K-1)
     aJac(ixSub3,ixVegNrg) = (dt/mLayerDepth(kLayer))*(-dGroundEvaporation_dTCanopy/iden_water) ! dVol/dT (K-1)
     aJac(ixSub2,ixVegWat) = (dt/mLayerDepth(kLayer))*(-dGroundEvaporation_dCanLiq/iden_water)  ! dVol/dLiq (kg m-2)-1
    endif
    aJac(ixSub1,ixTopNrg)   = (dt/mLayerDepth(kLayer))*(-dGroundEvaporation_dTGround/iden_water) + aJac(ixSub1,ixTopNrg) ! dVol/dT (K-1)
   endif

   ! melt-freeze: compute derivative in energy with respect to mass
   if(mLayerVolFracIceTrial(kLayer) > tiny(dt))then
    aJac(ixSup1,jLayer) = -dVolTot_dPsi0(iLayer)*LH_fus*iden_water    ! dNrg/dMat (J m-3 m-1) -- dMat changes volumetric water, and hence ice content
   else
    aJac(ixSup1,jLayer) = 0._dp
   endif

   ! - compute the Jacobian for neighboring layers (dVol/dT)
   if(kLayer > nSnow+1) aJac(ixSup1,mLayer) = (dt/mLayerDepth(kLayer-1))*( dq_dNrgStateBelow(iLayer-1))   ! K-1
   if(kLayer < nLayers) aJac(ixSub3,mLayer) = (dt/mLayerDepth(kLayer+1))*(-dq_dNrgStateAbove(iLayer))     ! K-1

  end do  ! (looping through soil layers)

  ! end association to variables in the data structures
  end associate

  if(printFlag)then
   print*, '** in cpact: banded analytical Jacobian:'
   write(*,'(a4,1x,100(i17,1x))') 'xCol', (iLayer, iLayer=iJac1,iJac2)
   do iLayer=kl+1,nBands
    write(*,'(i4,1x,100(e17.10,1x))') iLayer, (aJac(iLayer,jLayer),jLayer=iJac1,iJac2)
   end do
  endif


  end subroutine cpactBand

  ! *********************************************************************************************************
  ! internal subroutine analJacob: compute the Jacobian matrix (analytical)
  ! *********************************************************************************************************
  subroutine analJacob(err,message)
  implicit none
  ! dummy variables
  integer(i4b),intent(out)       :: err                     ! error code
  character(*),intent(out)       :: message                 ! error message
  ! local variables
  integer(i4b)                   :: pLayer,qLayer           ! indices of model layers
  ! --------------------------------------------------------------
  ! initialize error control
  err=0; message='analJacob/'

  ! associate variables from data structures
  associate(&

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

            ! groundwater variables
            ixGroundwater => model_decisions(iLookDECISIONS%groundwatr)%iDecision,&  ! intent(in): [i4b] groundwater parameterization
            mLayerDepth   => prog_data%var(iLookPROG%mLayerDepth)%dat             )  ! intent(in): [dp(:)] depth of each layer in the snow-soil sub-domain (m)

  ! initialize the Jacobian
  ! NOTE: this needs to be done every time, since Jacobian matrix is modified in the solver
  aJac(:,:) = 0._dp  ! analytical Jacobian matrix

  ! -----
  ! * energy and liquid fluxes over vegetation...
  ! ---------------------------------------------
  if(computeVegFlux)then  ! (derivatives only defined when vegetation protrudes over the surface)

   ! liquid water fluxes for vegetation canopy (-)
   aJac(ixVegWat,ixVegWat) = -fracLiqVeg*(dCanopyEvaporation_dCanLiq - scalarCanopyLiqDeriv)*dt + 1._dp

   ! cross-derivative terms w.r.t. system temperatures (kg m-2 K-1)
   aJac(ixVegWat,ixCasNrg) = -dCanopyEvaporation_dTCanair*dt
   aJac(ixVegWat,ixVegNrg) = -dCanopyEvaporation_dTCanopy*dt + dt*scalarCanopyLiqDeriv*dCanLiq_dTcanopy
   aJac(ixVegWat,ixTopNrg) = -dCanopyEvaporation_dTGround*dt

   ! cross-derivative terms w.r.t. canopy water (kg-1 m2)
   aJac(ixTopWat,ixVegWat) = (dt/mLayerDepth(1))*(-soilControl*fracLiqVeg*scalarCanopyLiqDeriv)/iden_water

   ! cross-derivative terms w.r.t. canopy temperature (K-1)
   aJac(ixTopWat,ixVegNrg) = (dt/mLayerDepth(1))*(-soilControl*scalarCanopyLiqDeriv*dCanLiq_dTcanopy)/iden_water

   ! cross-derivative terms w.r.t. canopy liquid water (J m-1 kg-1)
   ! NOTE: dIce/dLiq = (1 - fracLiqVeg); dIce*LH_fus/canopyDepth = J m-3; dLiq = kg m-2
   aJac(ixVegNrg,ixVegWat) = (dt/canopyDepth)   *(-dCanopyNetFlux_dCanLiq) - (1._dp - fracLiqVeg)*LH_fus/canopyDepth   ! dF/dLiq
   aJac(ixTopNrg,ixVegWat) = (dt/mLayerDepth(1))*(-dGroundNetFlux_dCanLiq)

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

  ! -----
  ! * energy fluxes for the snow-soil domain...
  ! -------------------------------------------
  do iLayer=1,nLayers  ! loop through layers in the snow-soil domain
   ! - define layer indices
   jLayer = ixSnowSoilNrg(iLayer)
   ! - compute the Jacobian
   aJac(jLayer,jLayer)   = (dt/mLayerDepth(iLayer))*(-dNrgFlux_dTempBelow(iLayer-1) + dNrgFlux_dTempAbove(iLayer)) + dMat(jLayer)
   if(iLayer > 1)       aJac(jLayer-nVarSnowSoil,jLayer) = (dt/mLayerDepth(iLayer-1))*( dNrgFlux_dTempBelow(iLayer-1) )
   if(iLayer < nLayers) aJac(jLayer+nVarSnowSoil,jLayer) = (dt/mLayerDepth(iLayer+1))*(-dNrgFlux_dTempAbove(iLayer  ) )
  end do  ! (looping through layers in the snow-soil system)

  ! -----
  ! * liquid water fluxes for the snow domain...
  ! --------------------------------------------
  do iLayer=1,nSnow
   ! - define layer indices
   jLayer = ixSnowOnlyWat(iLayer)   ! hydrology layer index within the full state vector
   mLayer = ixSnowSoilNrg(iLayer)   ! energy layer index within the full state vector
   ! - compute the Jacobian
   aJac(jLayer,jLayer) = (dt/mLayerDepth(iLayer))*iLayerLiqFluxSnowDeriv(iLayer)*fracLiqSnow(iLayer) + dMat(jLayer)
   if(iLayer > 1)     aJac(jLayer-nVarSnowSoil,jLayer) = 0._dp  ! sub-diagonal: no dependence on other layers
   ! - compute cross-derivative terms for the current layer
   ! NOTE: increase in volumetric liquid water content balanced by a decrease in volumetric ice content
   aJac(mLayer,jLayer) = -(1._dp - fracLiqSnow(iLayer))*LH_fus*iden_water     ! (dF/dLiq)
   aJac(jLayer,mLayer) = (dt/mLayerDepth(iLayer))*iLayerLiqFluxSnowDeriv(iLayer)*mLayerdTheta_dTk(iLayer)  ! (dVol/dT)
   ! - compute cross-derivative terms for the layer below (w.r.t. state in the current layer)
   if(iLayer < nSnow)then
    aJac(jLayer+nVarSnowSoil,mLayer) = -(dt/mLayerDepth(iLayer+1))*iLayerLiqFluxSnowDeriv(iLayer)*mLayerdTheta_dTk(iLayer)        ! dVol(below)/dT(above) -- K-1
    aJac(jLayer+nVarSnowSoil,jLayer) = -(dt/mLayerDepth(iLayer+1))*iLayerLiqFluxSnowDeriv(iLayer)*fracLiqSnow(iLayer)             ! dVol(below)/dLiq(above) -- (-)
   endif
  end do

  ! -----
  ! * liquid water fluxes for the soil domain...
  ! --------------------------------------------
  do iLayer=1,nSoil    ! loop through layers in the soil domain

   ! - define layer indices
   jLayer = ixSoilOnlyHyd(iLayer)  ! layer index within the full state vector
   kLayer = iLayer+nSnow           ! layer index within the full snow-soil vector

   ! - compute the Jacobian
   ! all terms *excluding* baseflow
   aJac(jLayer,jLayer) = (dt/mLayerDepth(kLayer))*(-dq_dHydStateBelow(iLayer-1) + dq_dHydStateAbove(iLayer)) + dMat(jLayer)
   if(kLayer > nSnow+1) aJac(jLayer-nVarSnowSoil,jLayer) = (dt/mLayerDepth(kLayer-1))*( dq_dHydStateBelow(iLayer-1))
   if(kLayer < nLayers) aJac(jLayer+nVarSnowSoil,jLayer) = (dt/mLayerDepth(kLayer+1))*(-dq_dHydStateAbove(iLayer))

   ! include terms for baseflow
   if(ixGroundwater==qbaseTopmodel)then
    do pLayer=1,nSoil
     qLayer = ixSoilOnlyHyd(pLayer)  ! layer index within the full state vector
     aJac(jLayer,qLayer) = aJac(jLayer,qLayer) + (dt/mLayerDepth(kLayer))*dBaseflow_dMatric(iLayer,pLayer)
    end do
   endif

  end do  ! (looping through soil layers)

  ! -----
  ! * derivative in liquid water fluxes w.r.t. temperature for the soil domain...
  ! -----------------------------------------------------------------------------
  do iLayer=1,nSoil    ! loop through layers in the soil domain

   ! - define layer indices
   kLayer = iLayer+nSnow                ! layer index within the full snow-soil vector
   jLayer = ixSoilOnlyHyd(iLayer)       ! hydrology layer index within the full state vector
   mLayer = ixSnowSoilNrg(kLayer)       ! thermodynamics layer index within the full state vector

   ! - compute the Jacobian for the layer itself
   aJac(jLayer,mLayer) = (dt/mLayerDepth(kLayer))*(-dq_dNrgStateBelow(iLayer-1) + dq_dNrgStateAbove(iLayer))   ! dVol/dT (K-1) -- flux depends on ice impedance

   ! - include derivatives w.r.t. ground evaporation
   if(nSnow==0 .and. iLayer==1)then  ! upper-most soil layer
    if(computeVegFlux)then
     aJac(jLayer,ixVegWat) = (dt/mLayerDepth(kLayer))*(-dGroundEvaporation_dCanLiq/iden_water)  ! dVol/dLiq (kg m-2)-1
     aJac(jLayer,ixCasNrg) = (dt/mLayerDepth(kLayer))*(-dGroundEvaporation_dTCanair/iden_water) ! dVol/dT (K-1)
     aJac(jLayer,ixVegNrg) = (dt/mLayerDepth(kLayer))*(-dGroundEvaporation_dTCanopy/iden_water) ! dVol/dT (K-1)
    endif
    aJac(jLayer,ixTopNrg) = (dt/mLayerDepth(kLayer))*(-dGroundEvaporation_dTGround/iden_water) + aJac(jLayer,ixTopNrg) ! dVol/dT (K-1)
   endif

   ! melt-freeze: compute derivative in energy with respect to mass
   if(mLayerVolFracIceTrial(iLayer+nSnow) > tiny(dt))then
    aJac(mLayer,jLayer) = -dVolTot_dPsi0(iLayer)*LH_fus*iden_water    ! dNrg/dMat (J m-3 m-1) -- dMat changes volumetric water, and hence ice content
   else
    aJac(mLayer,jLayer) = 0._dp
   endif

   ! - compute the Jacobian for neighboring layers
   if(kLayer > nSnow+1) aJac(jLayer-nVarSnowSoil,mLayer) = (dt/mLayerDepth(kLayer-1))*( dq_dNrgStateBelow(iLayer-1))   ! K-1
   if(kLayer < nLayers) aJac(jLayer+nVarSnowSoil,mLayer) = (dt/mLayerDepth(kLayer+1))*(-dq_dNrgStateAbove(iLayer))     ! K-1

  end do  ! (looping through soil layers)

  ! print the Jacobian
  if(printFlag)then
   print*, '** analytical Jacobian:'
   write(*,'(a4,1x,100(i12,1x))') 'xCol', (iLayer, iLayer=iJac1,iJac2)
   do iLayer=iJac1,iJac2; write(*,'(i4,1x,100(e12.5,1x))') iLayer, aJac(iJac1:iJac2,iLayer); end do
  endif
  !pause 'testing analytical jacobian'

  ! end the association to data structures
  end associate

  end subroutine analJacob


  ! *********************************************************************************************************
  ! internal subroutine numlJacob: compute the Jacobian matrix (numerical)
  ! *********************************************************************************************************
  subroutine numlJacob(stateVec,fluxVec,resVec,err,message)
  implicit none
  ! dummy
  real(dp),intent(in)            :: stateVec(:)             ! model state vector (mixed units)
  real(dp),intent(in)            :: fluxVec(:)              ! model flux vector (mixed units)
  real(qp),intent(in)            :: resVec(:)  ! qp         ! model residual vector (mixed units)
  integer(i4b),intent(out)       :: err                     ! error code
  character(*),intent(out)       :: message                 ! error message
  ! local
  character(len=256)             :: cmessage                ! error message of downwind routine
  real(dp),dimension(nState)     :: stateVecPerturbed       ! perturbed state vector
  real(dp),dimension(nState)     :: fluxVecJac              ! flux vector
  real(qp),dimension(nState)     :: resVecJac   ! qp        ! residual vector (mixed units)
  integer(i4b)                   :: iJac                    ! index of row of the Jacobian matrix
  integer(i4b),parameter         :: iTry=-999               ! index of trial model state variable (used for testing)
  ! trial state variables (vegetation canopy)
  real(dp)                       :: scalarCanairTempLocal   ! trial value for temperature of the canopy air space (K)
  real(dp)                       :: scalarCanopyTempLocal   ! trial value for temperature of the vegetation canopy (K)
  real(dp)                       :: scalarCanopyLiqLocal    ! trial value for mass of liquid water on the vegetation canopy (kg m-2)
  real(dp)                       :: scalarCanopyIceLocal    ! trial value for mass of ice on the vegetation canopy (kg m-2)
  ! trial state variables (snow and soil domains)
  real(dp),dimension(nLayers)    :: mLayerTempLocal         ! trial value for temperature of each snow/soil layer (K)
  real(dp),dimension(nLayers)    :: mLayerVolFracLiqLocal   ! trial value for volumetric fraction of liquid water (-)
  real(dp),dimension(nLayers)    :: mLayerVolFracIceLocal   ! trial value for volumetric fraction of ice (-)
  real(dp),dimension(nSoil)      :: mLayerMatricHeadLocal   ! trial value for matric head (m)
  ! model control -- swith between flux-based form and residual-based form of the numerical Jacobian
  integer(i4b),parameter         :: ixNumFlux=1001          ! named variable for the flux-based form of the numerical Jacobian
  integer(i4b),parameter         :: ixNumRes=1002           ! named variable for the residual-based form of the numerical Jacobian
  integer(i4b)                   :: ixNumType=ixNumRes      ! method used to calculate the numerical Jacobian
  ! --------------------------------------------------------------
  ! initialize error control
  err=0; message='numlJacob/'

  ! associate local variables with the indices of model state variables
  associate(&
            ixCasNrg      => indx_data%var(iLookINDEX%ixCasNrg)%dat(1)  , & ! index of canopy air space energy state variable
            ixVegNrg      => indx_data%var(iLookINDEX%ixVegNrg)%dat(1)  , & ! index of canopy energy state variable
            ixVegWat      => indx_data%var(iLookINDEX%ixVegWat)%dat(1)  , & ! index of canopy hydrology state variable (mass)
            ixSnowSoilNrg => indx_data%var(iLookINDEX%ixSnowSoilNrg)%dat, & ! indices for energy states in the snow-soil subdomain
            ixSoilOnlyHyd => indx_data%var(iLookINDEX%ixSoilOnlyHyd)%dat  & ! indices for hydrology states in the soil subdomain
  ) ! association to local variables with the indices of model state variables

  ! get a copy of the state vector to perturb
  stateVecPerturbed(:) = stateVec(:)

  ! get a copy of the canopy ice storage
  scalarCanopyLiqLocal  = scalarCanopyLiqTrial
  scalarCanopyIceLocal  = scalarCanopyIceTrial

  ! get a copy of the volumetric liquid water and ice content
  mLayerVolFracLiqLocal = mLayerVolFracLiqTrial
  mLayerVolFracIceLocal = mLayerVolFracIceTrial

  ! loop through state variables
  do iJac=1,nState

   ! define printFlag
   if(iJac==iTry) printFlag=.true.
   if(iJac/=iTry) printFlag=.false.

   ! (perturb state vector)
   stateVecPerturbed(iJac) = stateVec(iJac) + dx

   ! (use constitutive functions to compute unknown terms removed from the state equations...)
   call updatState(&
                   stateVecPerturbed,     & ! intent(in):  full state vector (mixed units)
                   mLayerVolFracLiqLocal, & ! intent(out): volumetric fraction of liquid water (-)
                   mLayerVolFracIceLocal, & ! intent(out): volumetric fraction of ice (-)
                   scalarCanopyLiqLocal,  & ! intent(out): mass of canopy liquid (kg m-2)
                   scalarCanopyIceLocal,  & ! intent(out): mass of canopy ice (kg m-2)
                   err,cmessage)            ! intent(out): error code and error message
   if(err/=0)then; message=trim(message)//trim(cmessage); return; endif  ! (check for errors)

   ! **
   ! ** residual-based calculation of the numerical Jacobian
   if(ixNumType==ixNumRes)then ! switch between the residual-based form and flux-based form

    ! (compute residual vector)
    call xFluxResid(&
                    ! input
                    stateVecPerturbed,             & ! intent(in): full state vector (mixed units)
                    scalarCanopyLiqLocal,          & ! intent(in): trial value for the liquid water on the vegetation canopy (kg m-2)
                    scalarCanopyIceLocal,          & ! intent(in): trial value for the ice on the vegetation canopy (kg m-2)
                    mLayerVolFracLiqLocal,         & ! intent(in): trial value for the volumetric liquid water content in each snow and soil layer (-)
                    mLayerVolFracIceLocal,         & ! intent(in): trial value for the volumetric ice in each snow and soil layer (-)
                    ! output
                    fluxVecJac,                    & ! intent(out): flux vector (mixed units)
                    resVecJac,                     & ! intent(out): residual vector (mixed units)
                    err,cmessage)                    ! intent(out): error code and error message
    if(err/=0)then; message=trim(message)//trim(cmessage); return; endif  ! (check for errors)

    ! (compute the row of the Jacobian matrix)
    nJac(:,iJac) = (resVecJac - resVec)/dx

   ! **
   ! ** flux-based calculation of the numerical Jacobian
   else

    ! NOTE: Need to increase cleverness and avoid copying vectors
    !  --> can we do this as an associate statement?

    ! extract the vegetation states from the state vector
    if(computeVegFlux)then
     scalarCanairTempLocal = stateVecPerturbed(ixCasNrg)
     scalarCanopyTempLocal = stateVecPerturbed(ixVegNrg)
    endif

    ! extract state variables for the snow and soil domain
    mLayerTempLocal(1:nLayers)     = stateVecPerturbed(ixSnowSoilNrg)
    mLayerMatricHeadLocal(1:nSoil) = stateVecPerturbed(ixSoilOnlyHyd)

    ! (compute fluxes)
    call computFlux(&
                    ! input: state variables
                    scalarCanairTempLocal,              & ! intent(in): trial value for the temperature of the canopy air space (K)
                    scalarCanopyTempLocal,              & ! intent(in): trial value for the temperature of the vegetation canopy (K)
                    mLayerTempLocal,                    & ! intent(in): trial value for the temperature of each snow and soil layer (K)
                    mLayerMatricHeadLocal,              & ! intent(in): trial value for the matric head in each soil layer (m)
                    ! input: diagnostic variables defining the liquid water and ice content
                    scalarCanopyLiqLocal,               & ! intent(in): trial value for the liquid water on the vegetation canopy (kg m-2)
                    scalarCanopyIceLocal,               & ! intent(in): trial value for the ice on the vegetation canopy (kg m-2)
                    mLayerVolFracLiqLocal,              & ! intent(in): trial value for the volumetric liquid water content in each snow and soil layer (-)
                    mLayerVolFracIceLocal,              & ! intent(in): trial value for the volumetric ice in each snow and soil layer (-)
                    ! output: flux vector
                    fluxVecJac,                         & ! intent(out): flux vector (mixed units)
                    ! output: error control
                    err,cmessage)                         ! intent(out): error code and error message
    if(err/=0)then; message=trim(message)//trim(cmessage); return; endif  ! (check for errors)

    ! (compute the row of the Jacobian matrix)
    nJac(:,iJac) = -dt*(fluxVecJac(:) - fluxVec(:))/dx

    ! (add in the diagonal matrix)
    nJac(iJac,iJac) = nJac(iJac,iJac) + dMat(iJac)

   endif

   ! (set the state back to the input value)
   stateVecPerturbed(iJac) = stateVec(iJac)

   ! (set the liquid water content back to the input value)
   mLayerVolFracLiqLocal = mLayerVolFracLiqTrial

   ! (set the ice content back to the input value)
   scalarCanopyIceLocal  = scalarCanopyIceTrial
   mLayerVolFracIceLocal = mLayerVolFracIceTrial

  end do  ! (looping through state variables)

  ! print the Jacobian
  print*, '** numerical Jacobian:', ixNumType==ixNumRes
  write(*,'(a4,1x,100(i12,1x))') 'xCol', (iLayer, iLayer=iJac1,iJac2)
  do iJac=iJac1,iJac2; write(*,'(i4,1x,100(e12.5,1x))') iJac, nJac(iJac1:iJac2,iJac); end do
  !pause 'testing Jacobian'

  ! end association to local variables with the indices of model state variables
  end associate

  end subroutine numlJacob


  ! *********************************************************************************************************
  ! internal subroutine lapackSolv: use the lapack routines to solve the linear system A.X=B
  ! *********************************************************************************************************
  subroutine lapackSolv(aJac,rVec,grad,xInc,err,message)
  implicit none
  ! dummy
  real(dp),intent(inout)         :: aJac(:,:)     ! input = the Jacobian matrix A; output = decomposed matrix
  real(qp),intent(in)            :: rVec(:)  ! qp ! the residual vector B
  real(dp),intent(out)           :: grad(:)       ! gradient of the function vector
  real(dp),intent(out)           :: xInc(:)       ! the solution vector X
  integer(i4b),intent(out)       :: err           ! error code
  character(*),intent(out)       :: message       ! error message
  ! local
  integer(i4b)                   :: iJac                    ! index of row of the Jacobian matrix
  integer(i4b)                   :: iState,jState,kState    ! indices of state variables
  ! initialize error control
  select case(ixSolve)
   case(ixFullMatrix); message='lapackSolv/dgesv/'
   case(ixBandMatrix); message='lapackSolv/dgbsv/'
   case default; err=20; message=trim(message)//'unable to identify option for the type of matrix'
  end select

  ! --------------------------------------------------------------
  ! * scale variables
  ! --------------------------------------------------------------

  if(printFlag)then
   write(*,'(a,1x,10(e17.10,1x))') 'fScale = ', fScale(iJac1:iJac2)
  endif

  ! select the option used to solve the linear system A.X=B
  select case(ixSolve)

   ! * full Jacobian matrix
   case(ixFullMatrix)
    do iJac=1,nState
     aJac(iJac,1:nState) = aJac(iJac,1:nState)/fscale(iJac)
    end do

    ! ** test band diagonal matrix
    if(testBandDiagonal)then

     aJac_test(:,:)=0._dp
     ! form band-diagonal matrix
     do iState=1,nState
      do jState=max(1,iState-ku),min(nState,iState+kl)
       aJac_test(kl + ku + 1 + jState - iState, iState) = aJac(jState,iState)
      end do
     end do
     print*, '** test banded analytical Jacobian:'
     write(*,'(a4,1x,100(i11,1x))') 'xCol', (iLayer, iLayer=iJac1,iJac2)
     do iLayer=kl+1,nBands
      write(*,'(i4,1x,100(e17.10,1x))') iLayer, (aJac_test(iLayer,iJac),iJac=iJac1,iJac2)
     end do

    endif  ! (if desire to test band-diagonal matric

   ! * band-diagonal matrix
   case(ixBandMatrix)
    do iJac=1,nState   ! (loop through state variables)
     do iState=kl+1,nBands  ! (loop through elements of the band-diagonal matrix)
      kState = iState + iJac - kl - ku - 1
      if(kState<1 .or. kState>nState)cycle
      aJac(iState,iJac) = aJac(iState,iJac)/fscale(kState)
     end do  ! looping through elements of the band-diagonal matric
    end do  ! looping through state variables

    ! check
    if(printFlag)then
     print*, '** banded analytical Jacobian:'
     write(*,'(a4,1x,100(i11,1x))') 'xCol', (iLayer, iLayer=iJac1,iJac2)
     do iLayer=kl+1,nBands
      write(*,'(i4,1x,100(e17.10,1x))') iLayer, (aJac(iLayer,iJac),iJac=iJac1,iJac2)
     end do
    endif

  end select  ! (option to solve the linear system A.X=B)


  ! form the rhs matrix
  ! NOTE: scale the residual vector
  rhs(1:nState,1) = -real(rVec(1:nState), dp)/fScale(1:nState)

  ! --------------------------------------------------------------
  ! * compute the gradient of the function vector
  ! --------------------------------------------------------------

  ! compute the gradient of the function vector
  select case(ixSolve)

   ! full Jacobian matrix
   case(ixFullMatrix)
   grad = matmul(-rhs(1:nState,1),aJac)

   ! band-diagonal matrix
   case(ixBandMatrix)
   do iJac=1,nState  ! (loop through state variables)

    grad(iJac) = 0._dp
    do iState=kl+1,nBands  ! (loop through elements of the band-diagonal matrix)

     ! identify indices in the band-diagonal matrix
     kState = iJac + iState-2*kl
     if(kState < 1 .or. kState > nState)cycle

     ! calculate gradient (long-hand matrix multiplication)
     grad(iJac) = grad(iJac) - aJac(iState,iJac)*rhs(kState,1)

    end do  ! looping through elements of the band-diagonal matric
   end do  ! looping through state variables

  end select  ! (option to solve the linear system A.X=B)

  ! --------------------------------------------------------------
  ! * solve the linear system A.X=B
  ! --------------------------------------------------------------

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

  ! --------------------------------------------------------------
  ! * wrap-up
  ! --------------------------------------------------------------

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
  xInc(1:nState) = rhs(1:nState,1)

  end subroutine lapackSolv


  ! *********************************************************************************************************
  ! internal subroutine lineSearch: perform the line search
  ! *********************************************************************************************************
  ! Routine modified extentively from Numerical Recipes in Fortran (Press et al. 1998) to
  !  1) Make use of local variables for the flux and residual calculations;
  !  2) Scale function evaluations and state vectors;
  !  3) Return error code and message;
  !  4) Additonal comments.
  ! ************************************************************************************************
  subroutine lineSearch(&
                        ! input
                        doLineSearch,            & ! intent(in): flag to denote the need to perform line search
                        xOld,                    & ! intent(in): initial state vector
                        fOld,                    & ! intent(in): function value for trial state vector (mixed units)
                        g,                       & ! intent(in): gradient of the function vector (mixed units)
                        p,                       & ! intent(in): iteration increment (mixed units)
                        ! output
                        x,                       & ! intent(out): new state vector (m)
                        fVec,                    & ! intent(out): new flux vector (mixed units)
                        rVec,                    & ! intent(out): new residual vector (mixed units)
                        f,                       & ! intent(out): new function value (mixed units)
                        converged,               & ! intent(out): convergence flag
                        err,message)               ! intent(out): error control
  IMPLICIT NONE
  ! input variables
  logical(lgt),intent(in) :: doLineSearch
  REAL(DP), DIMENSION(:), INTENT(IN) :: xOld,g
  REAL(DP), DIMENSION(:), INTENT(INOUT) :: p
  REAL(DP), INTENT(IN) :: fOld
  ! output variables
  REAL(DP), DIMENSION(:), INTENT(OUT) :: x,fVec
  REAL(QP), DIMENSION(:), INTENT(OUT) :: rVec
  REAL(DP), INTENT(OUT) :: f
  logical(lgt)              :: converged   ! convergence flag
  integer(i4b),intent(out)  :: err         ! error code
  character(*),intent(out)  :: message     ! error message
  ! local variables
  character(LEN=256)            :: cmessage                 ! error message of downwind routine
  ! variables for the line search
  REAL(DP), PARAMETER :: ALF=1.0e-4_dp,TOLX=epsilon(x),xTolInc=1.0e-4_dp
  INTEGER(I4B) :: ndum,iterLS
  integer(i4b),parameter :: maxiterLS=5
  REAL(DP) :: a,alam,alam2,alamin,b,disc,f2,fold2,pabs,rhs1,rhs2,slope,&
      tmplam
  ! initialize error control
  err=0; message="lineSearch/"

  ! check arguments
  if ( all((/size(g),size(p),size(x)/) == size(xold)) ) then
   ndum=size(xold)
  else
   err=20; message=trim(message)//"sizeMismatch"; return
  endif

  ! define step size and initialize tolerances
  if(doLineSearch)then
   pabs=norm2(p/xscale)   ! NOTE: norm2 is the Euclidean norm
   if (pabs > stpmax) p(:)=p(:)*stpmax/pabs  ! reduce step if it is too big
   slope=dot_product(g,p)
   alamin=TOLX/maxval(abs(p(:))/max(abs(xold(:)),xscale))  ! minimum lambda
  endif
  alam=1.0_dp

  ! backtrack
  do iterLS=1,maxIterLS

   ! update the state vector
   x(:)=xold(:)+alam*p(:)

   ! use constitutive functions to compute unknown terms removed from the state equations...
   call updatState(&
                   x,                     & ! intent(in): full state vector (mixed units)
                   mLayerVolFracLiqTrial, & ! intent(out): volumetric fraction of liquid water (-)
                   mLayerVolFracIceTrial, & ! intent(out): volumetric fraction of ice (-)
                   scalarCanopyLiqTrial,  & ! intent(out): mass of canopy liquid (kg m-2)
                   scalarCanopyIceTrial,  & ! intent(out): mass of canopy ice (kg m-2)
                   err,cmessage)            ! intent(out): error code and error message
   if(err/=0)then; message=trim(message)//trim(cmessage); return; endif  ! (check for errors)

   ! compute flux vector and residual
   call xFluxResid(&
                   ! input
                   x,                     & ! intent(in): full state vector (mixed units)
                   scalarCanopyLiqTrial,  & ! intent(in): trial value for the liquid water on the vegetation canopy (kg m-2)
                   scalarCanopyIceTrial,  & ! intent(in): trial value for the ice on the vegetation canopy (kg m-2)
                   mLayerVolFracLiqTrial, & ! intent(in): trial value for the volumetric liquid water content in each snow and soil layer (-)
                   mLayerVolFracIceTrial, & ! intent(in): trial value for the volumetric ice in each snow and soil layer (-)
                   ! output
                   fVec,                  & ! intent(out): flux vector (mixed units)
                   rVec,                  & ! intent(out): residual vector (mixed units)
                   err,cmessage)            ! intent(out): error code and error message
   if(err/=0)then; message=trim(message)//trim(cmessage); return; endif  ! (check for errors)

   ! compute the function evaluation
   f=0.5_dp*norm2(real(rVec, dp)/fScale)  ! NOTE: norm2 = sqrt(sum((rVec/fScale)**2._dp))

   ! check
   if(printFlag)then
    print*, '***'
    write(*,'(a,1x,100(e14.5,1x))')  trim(message)//': alam, fOld, f       = ', alam, fOld, f
    write(*,'(a,1x,100(f20.8,1x))')  trim(message)//': x(iJac1:iJac2)      = ', x(iJac1:iJac2)
    write(*,'(a,1x,100(f20.12,1x))') trim(message)//': p(iJac1:iJac2)      = ', p(iJac1:iJac2)
    write(*,'(a,1x,100(e20.5,1x))')  trim(message)//': rVec(iJac1:iJac2)   = ', rVec(iJac1:iJac2)
   endif

   ! return if not doing the line search
   if(.not.doLineSearch)then
    converged = checkConv(rVec,p,x,soilWaterBalanceError)
    return
   endif

   ! check convergence
   ! NOTE: this must be after the first flux call
   converged = checkConv(rVec,p,x,soilWaterBalanceError)
   if(converged) return

   ! check if backtracked all the way to the original value
   if (iterLS==maxIterLS) then   !if (alam < alamin) then
    x(:)=xold(:)
    err=-10; message=trim(message)//'warning: check convergence'
    RETURN

   ! check if improved the solution sufficiently
   else if (f <= fold+ALF*alam*slope) then
    RETURN

   ! build another trial vector
   else
    if (alam == 1.0_dp) then
     tmplam=-slope/(2.0_dp*(f-fold-slope))
     if (tmplam > 0.5_dp*alam) tmplam=0.5_dp*alam
    else
     rhs1=f-fold-alam*slope
     rhs2=f2-fold2-alam2*slope
     a=(rhs1/alam**2-rhs2/alam2**2)/(alam-alam2)
     b=(-alam2*rhs1/alam**2+alam*rhs2/alam2**2)/&
         (alam-alam2)
     if (a == 0.0_dp) then
      tmplam=-slope/(2.0_dp*b)
     else
      disc=b*b-3.0_dp*a*slope
      if (disc < 0.0_dp)then; err=-10; message=trim(message)//'warning: roundoff problem in lnsrch'; return; endif
      tmplam=(-b+sqrt(disc))/(3.0_dp*a)
     end if
     if (tmplam > 0.5_dp*alam) tmplam=0.5_dp*alam
    end if
   end if
   alam2=alam
   f2=f
   fold2=fold
   alam=max(tmplam,0.1_dp*alam)

  end do
  END SUBROUTINE lineSearch


  ! *********************************************************************************************************
  ! internal function checkConv: check convergence based on the residual vector
  ! *********************************************************************************************************
  function checkConv(rVec,xInc,xVec,soilWatbalErr)
  implicit none
  ! dummies
  real(qp),intent(in)       :: rVec(:) ! qp    ! residual vector (mixed units)
  real(dp),intent(in)       :: xInc(:)         ! iteration increment (mixed units)
  real(dp),intent(in)       :: xVec(:)         ! state vector (mixed units)
  real(dp),intent(in)       :: soilWatbalErr   ! soil water balance error (m)
  logical(lgt)              :: checkConv       ! flag to denote convergence
  ! locals
  real(dp),dimension(nSoil) :: psiScale        ! scaling factor for matric head
  real(dp),parameter        :: xSmall=1.e-0_dp ! a small offset
  logical(lgt)              :: canopyConv      ! flag for canopy water balance convergence
  logical(lgt)              :: watbalConv      ! flag for soil water balance convergence
  logical(lgt)              :: liquidConv      ! flag for residual convergence
  logical(lgt)              :: matricConv      ! flag for matric head convergence
  logical(lgt)              :: energyConv      ! flag for energy convergence
  ! make association with variables in the data structures
  associate(&
   ixCasNrg      => indx_data%var(iLookINDEX%ixCasNrg)%dat(1)  , & ! index of canopy air space energy state variable
   ixVegNrg      => indx_data%var(iLookINDEX%ixVegNrg)%dat(1)  , & ! index of canopy energy state variable
   ixVegWat      => indx_data%var(iLookINDEX%ixVegWat)%dat(1)  , & ! index of canopy hydrology state variable (mass)
   ixSnowSoilNrg => indx_data%var(iLookINDEX%ixSnowSoilNrg)%dat, & ! indices for energy states in the snow-soil subdomain
   ixSnowSoilWat => indx_data%var(iLookINDEX%ixSnowSoilWat)%dat, & ! indices for total water states in the snow-soil subdomain
   ixSoilOnlyHyd => indx_data%var(iLookINDEX%ixSoilOnlyHyd)%dat  & ! indices for hydrology states in the soil subdomain
  ) ! association with variables in the data structures

  ! check convergence based on the residuals for energy (J m-3)
  if(computeVegFlux)then
   canopy_max = real(abs(rVec(ixVegWat)), dp)*iden_water
   energy_max = real(maxval(abs( (/rVec(ixCasNrg), rVec(ixVegNrg), rVec(ixSnowSoilNrg)/) ) ), dp)
   energy_loc =      maxloc(abs( (/rVec(ixCasNrg), rVec(ixVegNrg), rVec(ixSnowSoilNrg)/) ) )
   canopyConv = (canopy_max    < absConvTol_watbal)  ! absolute error in canopy water balance (m)
  else
   energy_max = real(maxval(abs( rVec(ixSnowSoilNrg) ) ), dp)
   energy_loc =      maxloc(abs( rVec(ixSnowSoilNrg) ) )
   canopyConv = .true. ! don't check canopy convergence if not computing the vegetation flux (canopy is buried by snow or canopy non-existent)
  endif

  ! check convergence based on the residuals for volumetric liquid water content (-)
  liquid_max = real(maxval(abs( rVec(ixSnowSoilWat) ) ), dp)
  liquid_loc =      maxloc(abs( rVec(ixSnowSoilWat) ) )

  ! check convergence based on the iteration increment for matric head
  ! NOTE: scale by matric head to avoid unnecessairly tight convergence when there is no water
  psiScale   = abs(xVec(ixSoilOnlyHyd)) + xSmall ! avoid divide by zero
  matric_max = maxval(abs( xInc(ixSoilOnlyHyd)/psiScale ) )
  matric_loc = maxloc(abs( xInc(ixSoilOnlyHyd)/psiScale ) )

  ! convergence check
  watbalConv = (soilWatbalErr < absConvTol_watbal)  ! absolute error in total soil water balance (m)
  matricConv = (matric_max(1) < absConvTol_matric)  ! NOTE: based on iteration increment
  liquidConv = (liquid_max(1) < absConvTol_liquid)  ! (based on the residual)
  energyConv = (energy_max(1) < absConvTol_energy)  ! (based on the residual)

  ! print progress towards solution
  if(printFlag)then
   print*, 'iter, dt = ', iter, dt
   write(*,'(a,1x,4(e15.5,1x),3(i4,1x),5(L1,1x))') 'fNew, matric_max(1), liquid_max(1), energy_max(1), matric_loc(1), liquid_loc(1), energy_loc(1), matricConv, liquidConv, energyConv, watbalConv, canopyConv = ', &
                                                    fNew, matric_max(1), liquid_max(1), energy_max(1), matric_loc(1), liquid_loc(1), energy_loc(1), matricConv, liquidConv, energyConv, watbalConv, canopyConv
  endif

  ! final convergence check
  checkConv = (canopyConv .and. watbalConv .and. matricConv .and. liquidConv .and. energyConv)

  ! end association with variables in the data structures
  end associate

  end function checkConv

 end subroutine systemSolv


 ! **********************************************************************************************************
 ! private subroutine soilCmpres: compute soil compressibility (-) and its derivative w.r.t matric head (m-1)
 ! **********************************************************************************************************
 subroutine soilCmpres(&
                       ! input:
                       ixRichards,                         & ! intent(in): choice of option for Richards' equation
                       mLayerMatricHead,                   & ! intent(in): matric head at the start of the time step (m)
                       mLayerMatricHeadTrial,              & ! intent(in): trial value of matric head (m)
                       mLayerVolFracLiqTrial,              & ! intent(in): trial value for the volumetric liquid water content in each soil layer (-)
                       mLayerVolFracIceTrial,              & ! intent(in): trial value for the volumetric ice content in each soil layer (-)
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
 real(dp),intent(in)            :: mLayerVolFracIceTrial(:)  ! trial value for volumetric fraction of ice (-)
 real(dp),intent(in)            :: mLayerdTheta_dPsi(:)      ! derivative in the soil water characteristic (m-1)
 real(dp),intent(in)            :: specificStorage           ! specific storage coefficient (m-1)
 real(dp),intent(in)            :: theta_sat                 ! soil porosity (-)
 ! output:
 real(dp),intent(out)           :: compress(:)               ! soil compressibility (-)
 real(dp),intent(out)           :: dCompress_dPsi(:)         ! derivative in soil compressibility w.r.t. matric head (m-1)
 integer(i4b),intent(out)       :: err                       ! error code
 character(*),intent(out)       :: message                   ! error message
 ! local variables
 real(dp)                       :: volFracWat                ! total volumetric fraction of water (-)
 real(dp)                       :: fPart1,fPart2             ! different parts of the function
 real(dp)                       :: dPart1,dPart2             ! derivatives for different parts of the function
 integer(i4b)                   :: iLayer                    ! index of soil layer
 ! --------------------------------------------------------------
 ! initialize error control
 err=0; message='soilCmpres/'
 ! (only compute for the mixed form of Richards' equation)
 if(ixRichards==mixdform)then
  do iLayer=1,size(mLayerMatricHead)
   ! compute the total volumetric fraction of water (-)
   volFracWat = mLayerVolFracLiqTrial(iLayer) + mLayerVolFracIceTrial(iLayer)
   ! compute the compressibility term (-)
   compress(iLayer) = (specificStorage*volFracWat/theta_sat) * (mLayerMatricHeadTrial(iLayer) - mLayerMatricHead(iLayer))
   ! compute the derivative for the compressibility term (m-1)
   fPart1 = specificStorage*(volFracWat/theta_sat)  ! function for the 1st part (m-1)
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

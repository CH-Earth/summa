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
 USE summaSolve_module,only:summaSolve                ! calculate the iteration increment, evaluate the new state, and refine if necessary
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


 USE matrixOper_module, only: lapackSolv
 USE matrixOper_module, only: scaleMatrices



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
 integer(i4b),parameter          :: iJac1=1,iJac2=10             ! desired indices for printing
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
 real(qp),dimension(nState)      :: sMul          ! NOTE: qp     ! multiplier for state vector for the residual calculations
 real(qp),dimension(nState)      :: rAdd          ! NOTE: qp     ! additional terms in the residual vector
 real(qp),dimension(nState)      :: rVec          ! NOTE: qp     ! residual vector
 real(dp),dimension(nState)      :: rVecScaled                   ! residual vector (scaled)
 real(dp),dimension(nState)      :: newtStepScaled               ! full newton step (scaled)
 real(dp),dimension(nState)      :: xInc                         ! iteration increment
 real(dp),dimension(nState)      :: grad                         ! gradient of the function vector = matmul(rVec,aJac)
 real(dp),dimension(nState,nRHS) :: rhs                          ! the nState-by-nRHS matrix of matrix B, for the linear system A.X=B
 integer(i4b),dimension(nState)  :: iPiv                         ! defines if row i of the matrix was interchanged with row iPiv(i)
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

 ! -----
 ! * initialize...
 ! ---------------

 ! initialize the first flux call
 firstFluxCall=.true.

 ! define canopy depth
 if(computeVegFlux)then
  canopyDepth = heightCanopyTop - heightCanopyBottom
 else
  canopyDepth = realMissing
 endif

 ! compute the total water content in the vegetation canopy
 scalarCanopyWat = scalarCanopyLiq + scalarCanopyIce  ! kg m-2

 ! compute the total water content in snow and soil
 ! NOTE: no ice expansion allowed for soil
 if(nSnow>0)& 
 mLayerVolFracWat(1:      nSnow)   = mLayerVolFracLiq(1:      nSnow)   + mLayerVolFracIce(1:      nSnow)*(iden_ice/iden_water)
 mLayerVolFracWat(nSnow+1:nLayers) = mLayerVolFracLiq(nSnow+1:nLayers) + mLayerVolFracIce(nSnow+1:nLayers)

 ! identify the matrix solution method
 ! (the type of matrix used to solve the linear system A.X=B)
 if(ixGroundwater==qbaseTopmodel .or. testBandDiagonal .or. forceFullMatrix)then
  ixMatrix=ixFullMatrix   ! full Jacobian matrix
 else
  ixMatrix=ixBandMatrix   ! band-diagonal matrix
 endif

 ! allocate space for the derivative structure
 call allocLocal(deriv_meta(:),deriv_data,nSnow,nSoil,err,cmessage)
 if(err/=0)then; err=20; message=trim(message)//trim(cmessage); return; endif

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
                 canopyDepth,             & ! intent(in):    canopy depth (m)
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
                  nSnow,                   & ! intent(in):    number of snow layers
                  nSoil,                   & ! intent(in):    number of soil layers
                  nLayers,                 & ! intent(in):    total number of layers
                  nBands,                  & ! intent(in):    total number of bands in the band-diagonal matrix (nBands=nState for the full matrix)
                  nState,                  & ! intent(in):    total number of state variables
                  ixMatrix,                & ! intent(in):    type of matrix (full or band diagonal)
                  firstSubStep,            & ! intent(in):    flag to indicate if we are processing the first sub-step
                  firstFluxCall,           & ! intent(inout): flag to indicate if we are processing the first flux call
                  computeVegFlux,          & ! intent(in):    flag to indicate if we need to compute fluxes over vegetation
                  canopyDepth,             & ! intent(in):    depth of the vegetation canopy (m)
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
                  ! output
                  stateVecNew,             & ! intent(out):   new state vector
                  fluxVecNew,              & ! intent(out):   new flux vector
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

 ! end associate statement
 end associate

 contains


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

 end subroutine systemSolv

end module systemSolv_module

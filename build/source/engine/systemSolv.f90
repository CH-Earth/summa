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
USE data_struc,only:ix_soil,ix_snow ! named variables for snow and soil

! access the global print flag
USE data_struc,only:globalPrintFlag

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
real(dp),parameter  :: verySmall=tiny(1.0_dp)     ! a very small number
real(dp),parameter  :: veryBig=1.e+20_dp          ! a very big number
real(dp),parameter  :: dx = 1.e-8_dp             ! finite difference increment
contains


 ! **********************************************************************************************************
 ! public subroutine systemSolv: run the coupled energy-mass model for one timestep
 ! **********************************************************************************************************
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
                       indx_data,      & ! intent(in):    index data
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
                     var_ilength,      & ! data vector with variable length dimension (i4b)
                     var_dlength,      & ! data vector with variable length dimension (dp)
                     model_options       ! defines the model decisions
 ! provide access to indices that define elements of the data structures
 USE var_lookup,only:iLookATTR,iLookTYPE,iLookPARAM,iLookFORCE,iLookMVAR,iLookINDEX ! named variables for structure elements
 USE var_lookup,only:iLookDECISIONS                                                 ! named variables for elements of the decision structure
 ! provide access to the numerical recipes modules
 USE nr_utility_module,only:arth                      ! creates a sequence of numbers (start, incr, n)
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
 real(dp),intent(in)             :: dt                            ! time step (seconds)
 integer(i4b),intent(in)         :: maxiter                       ! maximum number of iterations
 logical(lgt),intent(in)         :: firstSubStep                  ! flag to indicate if we are processing the first sub-step
 logical(lgt),intent(in)         :: computeVegFlux                ! flag to indicate if we are computing fluxes over vegetation (.false. means veg is buried with snow)
 ! input/output: data structures
 type(var_i),intent(in)          :: type_data                     ! type of vegetation and soil
 type(var_d),intent(in)          :: attr_data                     ! spatial attributes
 type(var_d),intent(in)          :: forc_data                     ! model forcing data
 type(var_d),intent(in)          :: mpar_data                     ! model parameters
 type(var_ilength),intent(in)    :: indx_data
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
 integer(i4b),dimension(nLayers) :: layerType                    ! intent(in): type of layer in the snow+soil domain (snow or soil)
 ! domain boundary conditions
 real(dp)                        :: upperBoundTemp               ! intent(in): temperature of the upper boundary of the snow and soil domains (K)
 real(dp)                        :: scalarRainfall               ! intent(in): rainfall (kg m-2 s-1)
 real(dp)                        :: scalarSfcMeltPond            ! intent(in): ponded water caused by melt of the "snow without a layer" (kg m-2)
 ! diagnostic variables
 real(dp),dimension(nLayers)     :: mLayerDepth                  ! intent(in): depth of each layer in the snow-soil sub-domain (m)
 real(dp)                        :: scalarBulkVolHeatCapVeg      ! intent(in): bulk volumetric heat capacity of vegetation (J m-3 K-1)
 real(dp),dimension(nLayers)     :: mLayerVolHtCapBulk           ! intent(in): bulk volumetric heat capacity in each snow and soil layer (J m-3 K-1)
 real(dp),dimension(nLayers)     :: mLayerMeltFreeze             ! intent(out): melt-freeze in each snow and soil layer (kg m-3)
 real(dp),dimension(nSnow)       :: mLayerThetaResid             ! intent(out): residual volumetric liquid water content in each snow layer (-)
 ! model fluxes
 real(dp)                        :: scalarSurfaceRunoff          ! intent(out): surface runoff (m s-1)
 real(dp),dimension(0:nSnow)     :: iLayerLiqFluxSnow            ! intent(out): vertical liquid water flux at layer interfaces (m s-1)
 real(dp),dimension(0:nSoil)     :: iLayerLiqFluxSoil            ! intent(out): liquid flux at soil layer interfaces (m s-1)
 real(dp),dimension(nSoil)       :: mLayerColumnOutflow          ! intent(out): column outflow from each soil layer (m3 s-1)
 real(dp),dimension(nSoil)       :: mLayerBaseflow               ! intent(out): baseflow from each soil layer -- only compute at the start of the step (m s-1)
 real(dp),dimension(nSoil)       :: mLayerCompress               ! intent(out): change in storage associated with compression of the soil matrix (-)
 real(dp)                        :: scalarCanopySublimation      ! intent(out): sublimation of ice from the vegetation canopy (kg m-2 s-1)
 real(dp)                        :: scalarSnowSublimation        ! intent(out): sublimation of ice from the snow surface (kg m-2 s-1)
 real(dp)                        :: scalarExfiltration           ! intent(out): exfiltration from the soil profile (m s-1)
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
 real(dp)                        :: fImpede                      ! intent(in): ice impedance parameter (-)
 ! snow parameters
 real(dp)                        :: snowfrz_scale                ! intent(in): scaling parameter for the snow freezing curve (K-1)
 ! model state variables (vegetation canopy)
 real(dp)                        :: scalarCanairTemp             ! intent(inout): temperature of the canopy air space (K)
 real(dp)                        :: scalarCanopyTemp             ! intent(inout): temperature of the vegetation canopy (K)
 real(dp)                        :: scalarCanopyIce              ! intent(inout): mass of ice on the vegetation canopy (kg m-2)
 real(dp)                        :: scalarCanopyLiq              ! intent(inout): mass of liquid water on the vegetation canopy (kg m-2)
 real(dp)                        :: scalarCanopyWat              ! intent(inout): mass of total water on the vegetation canopy (kg m-2)
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
 real(dp),parameter              :: xMinCanopyWater=0.0001_dp    ! minimum value to initialize canopy water (kg m-2)
 ! ------------------------------------------------------------------------------------------------------
 ! * model indices
 ! ------------------------------------------------------------------------------------------------------
 integer(i4b)                    :: iPos                         ! position in vector desire to print
 integer(i4b),parameter          :: nVegNrg=2                    ! number of energy state variables for vegetation
 integer(i4b),parameter          :: nVegLiq=1                    ! number of hydrology state variables for vegetation
 integer(i4b)                    :: nVegState                    ! number of vegetation state variables (defines position of snow-soil states in the state vector)
 integer(i4b)                    :: nState                       ! total number of model state variables
 integer(i4b),parameter          :: ixCasNrg=1                   ! index of the canopy air space state variable
 integer(i4b),parameter          :: ixVegNrg=2                   ! index of the canopy energy state variable
 integer(i4b),parameter          :: ixVegWat=3                   ! index of the canopy total water state variable
 integer(i4b)                    :: ixTopNrg                     ! index of the upper-most energy state variable in the snow-soil subdomain
 integer(i4b)                    :: ixTopLiq                     ! index of the upper-most liquid water state variable in the snow subdomain
 integer(i4b)                    :: ixTopMat                     ! index of the upper-most matric head state variable in the soil subdomain
 integer(i4b),dimension(nLayers) :: ixSnowSoilNrg                ! indices for energy state variables in the snow-soil subdomain
 integer(i4b),dimension(nLayers) :: ixSnowSoilWat                ! indices for total water state variables in the snow-soil subdomain
 integer(i4b),dimension(nSnow)   :: ixSnowOnlyNrg                ! indices for energy state variables in the snow subdomain
 integer(i4b),dimension(nSnow)   :: ixSnowOnlyWat                ! indices for total water state variables in the snow subdomain
 integer(i4b),dimension(nSoil)   :: ixSoilOnlyNrg                ! indices for energy state variables in the soil subdomain
 integer(i4b),dimension(nSoil)   :: ixSoilOnlyMat                ! indices for matric head state variables in the soil subdomain
 integer(i4b),parameter          :: nVarSnowSoil=2               ! number of state variables in the snow and soil domain (energy and liquid water/matric head)
 integer(i4b),parameter          :: nRHS=1                       ! number of unknown variables on the RHS of the linear system A.X=B
 integer(i4b),parameter          :: ku=3                         ! number of super-diagonal bands
 integer(i4b),parameter          :: kl=3                         ! number of sub-diagonal bands
 integer(i4b),parameter          :: ixSup3=kl+1                  ! index for the 3rd super-diagonal band
 integer(i4b),parameter          :: ixSup2=kl+2                  ! index for the 2nd super-diagonal band
 integer(i4b),parameter          :: ixSup1=kl+3                  ! index for the 1st super-diagonal band
 integer(i4b),parameter          :: ixDiag=kl+4                  ! index for the diagonal band
 integer(i4b),parameter          :: ixSub1=kl+5                  ! index for the 1st sub-diagonal band
 integer(i4b),parameter          :: ixSub2=kl+6                  ! index for the 2nd sub-diagonal band
 integer(i4b),parameter          :: ixSub3=kl+7                  ! index for the 3rd sub-diagonal band
 integer(i4b),parameter          :: nBands=2*kl+ku+1             ! length of tyhe leading dimension of the band diagonal matrix
 integer(i4b),parameter          :: ixFullMatrix=1001            ! named variable for the full Jacobian matrix
 integer(i4b),parameter          :: ixBandMatrix=1002            ! named variable for the band diagonal matrix
 integer(i4b)                    :: ixSolve                      ! the type of matrix used to solve the linear system A.X=B
 integer(i4b),parameter          :: iJac1=1                      ! first layer of the Jacobian to print
 integer(i4b),parameter          :: iJac2=10                     ! last layer of the Jacobian to print
 !integer(i4b),parameter          :: iJac1=457                    ! first layer of the Jacobian to print
 !integer(i4b),parameter          :: iJac2=466                    ! last layer of the Jacobian to print
 ! ------------------------------------------------------------------------------------------------------
 ! * fluxes and derivatives
 ! ------------------------------------------------------------------------------------------------------
 ! ice content (need to keep track of this, but not part of the state vector)
 real(dp)                        :: theta                        ! liquid water equivalent of total water (liquid plus ice)
 real(dp)                        :: scalarCanopyLiqTrial         ! trial value for mass of liquid water on the vegetation canopy (kg m-2)
 real(dp)                        :: scalarCanopyIceTrial         ! trial value for mass of ice on the vegetation canopy (kg m-2)
 real(dp),dimension(nLayers)     :: mLayerdTheta_dTk             ! derivative of volumetric liquid water content w.r.t. temperature (K-1)
 real(dp),dimension(nSoil)       :: dPsiLiq_dTemp                ! derivative in the liquid water matric potential w.r.t. temperature (m K-1)
 real(dp),dimension(nSnow)       :: fracLiqSnow                  ! fraction of liquid water in each snow layer (-)
 real(dp)                        :: fracLiqVeg                   ! fraction of liquid water on vegetation (-)
 real(dp)                        :: totalWaterVeg                ! total water on vegetation (kg m-2)
 real(dp)                        :: dTheta_dTkCanopy             ! derivative of volumetric liquid water content w.r.t. temperature (K-1)
 real(dp)                        :: dCanLiq_dTcanopy             ! derivative of canopy liquid storage w.r.t. temperature (kg m-2 K-1)
 ! volumetric liquid water content (need to keep track of this, but not part of the state vector for snow and soil)
 real(dp),dimension(nSnow)       :: mLayerVolFracWat             ! initial value of mass fraction of total water (-)
 real(dp),dimension(nLayers)     :: mLayerVolFracIceTrial        ! trial value for volumetric fraction of ice (-)
 real(dp),dimension(nLayers)     :: mLayerVolFracLiqTrial        ! trial value for volumetric fraction of liquid water (-)
 real(dp)                        :: fLiq0,fLiq1                  ! fraction of liquid water -- used to compute numerical derivatives (-)
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
 real(dp),dimension(0:nSnow)     :: iLayerLiqFluxSnowDeriv       ! derivative in vertical liquid water flux at layer interfaces (m s-1)
 real(dp)                        :: scalarRainPlusMelt           ! surface water input to the soil zone (m s-1)
 ! liquid water fluxes and derivatives for the soil domain
 real(dp)                        :: xMaxInfilRate                ! maximum infiltration rate (m s-1)
 real(dp)                        :: scalarInfilArea              ! fraction of unfrozen area where water can infiltrate (-)
 real(dp)                        :: scalarFrozenArea             ! fraction of area that is considered impermeable due to soil ice (-)
 real(dp)                        :: scalarSoilBaseflow           ! total baseflow from the soil profile (m s-1)
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
 real(dp),dimension(nSoil)       :: mLayerdIceImpede_dT          ! derivative in the ice impedance factor w.r.t. temperature (K-1)
 real(dp),dimension(nSoil)       :: dCompress_dPsi               ! derivative in compressibility w.r.t matric head (m-1)
 real(dp),dimension(nSnow)       :: snowNetLiqFlux               ! net liquid water flux for each snow layer (s-1)
 real(dp),dimension(nSoil)       :: soilNetLiqFlux               ! net liquid water flux for each soil layer (s-1)
 real(dp),allocatable            :: dBaseflow_dMatric(:,:)       ! derivative in baseflow w.r.t. matric head (s-1)
 real(dp),dimension(nSoil)       :: mLayerMatricHeadDiff         ! iteration increment for the matric head (m)
 integer(i4b)                    :: ixSaturation                 ! index of lowest saturated layer (NOTE: only computed on the first iteration)
 ! liquid water fluxes and derivatives for the aquifer
 real(dp)                        :: scalarAquiferTranspire       ! transpiration loss from the aquifer at the start-of-step (m s-1)
 real(dp)                        :: scalarAquiferRecharge        ! recharge to the aquifer (m s-1)
 real(dp)                        :: scalarAquiferBaseflow        ! total baseflow from the aquifer (m s-1)
 ! ------------------------------------------------------------------------------------------------------
 ! * model solver
 ! ------------------------------------------------------------------------------------------------------
 logical(lgt),parameter          :: numericalJacobian=.false.     ! flag to compute the Jacobian matrix
 logical(lgt),parameter          :: testBandDiagonal=.false.     ! flag to test the band-diagonal matrix
 logical(lgt)                    :: firstFluxCall                ! flag to define the first flux call
 real(dp),allocatable            :: stateVecInit(:)              ! initial state vector (mixed units)
 real(dp),allocatable            :: stateVecTrial(:)             ! trial state vector (mixed units)
 real(dp),allocatable            :: stateVecNew(:)               ! new state vector (mixed units)
 real(dp),allocatable            :: fluxVec0(:)                  ! flux vector (mixed units)
 real(dp),allocatable            :: fluxVec1(:)                  ! flux vector used in the numerical Jacobian calculations (mixed units)
 real(dp),allocatable            :: fScale(:)                    ! characteristic scale of the function evaluations (mixed units)
 real(dp),allocatable            :: xScale(:)                    ! characteristic scale of the state vector (mixed units)
 real(dp),allocatable            :: aJac_test(:,:)               ! used to test the band-diagonal matrix structure
 real(dp),allocatable            :: aJac(:,:)                    ! analytical Jacobian matrix
 real(qp),allocatable            :: nJac(:,:)  ! NOTE: qp        ! numerical Jacobian matrix
 real(dp),allocatable            :: dMat(:)                      ! diagonal matrix (excludes flux derivatives)
 real(qp),allocatable            :: sMul(:)    ! NOTE: qp        ! multiplier for state vector for the residual calculations
 real(qp),allocatable            :: rAdd(:)    ! NOTE: qp        ! additional terms in the residual vector
 real(qp),allocatable            :: rVec(:)    ! NOTE: qp        ! residual vector
 real(dp),allocatable            :: xInc(:)                      ! iteration increment
 real(dp),allocatable            :: grad(:)                      ! gradient of the function vector = matmul(rVec,aJac)
 real(dp),allocatable            :: rhs(:,:)                     ! the nState-by-nRHS matrix of matrix B, for the linear system A.X=B
 integer(i4b),allocatable        :: iPiv(:)                      ! defines if row i of the matrix was interchanged with row iPiv(i)
 real(dp)                        :: fOld,fNew                    ! function values (-); NOTE: dimensionless because scaled
 real(dp)                        :: canopy_max                   ! absolute value of the residual in canopy water (kg m-2)
 real(dp),dimension(1)           :: energy_max                   ! maximum absolute value of the energy residual (J m-3)
 real(dp),dimension(1)           :: liquid_max                   ! maximum absolute value of the volumetric liquid water content residual (-)
 real(dp),dimension(1)           :: matric_max                   ! maximum absolute value of the matric head iteration increment (m)
 integer(i4b),dimension(1)       :: energy_loc                   ! location of maximum absolute value of the energy residual (-)
 integer(i4b),dimension(1)       :: liquid_loc                   ! location of maximum absolute value of the volumetric liquid water content residual (-)
 integer(i4b),dimension(1)       :: matric_loc                   ! location of maximum absolute value of the matric head increment (-)
 real(dp),parameter              :: absConvTol_energy=1.e-0_dp   ! convergence tolerance for energy (J m-3)
 real(dp),parameter              :: absConvTol_liquid=1.e-6_dp   ! convergence tolerance for volumetric liquid water content (-)
 real(dp),parameter              :: absConvTol_matric=1.e-3_dp   ! convergence tolerance for matric head increment in soil layers (m)
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
 real(dp)                        :: xIncScale                    ! scaling factor for the iteration increment (-)
 integer(i4b)                    :: iMin(1)                      ! index of most excessive drainage
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
 logical(lgt),parameter          :: checkMassBalance=.false.     ! flag to check the mass balance
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
 upperBoundTemp          => forc_data%var(iLookFORCE%airtemp)                      ,&  ! intent(in): [dp] temperature of the upper boundary of the snow and soil domains (K)
 scalarRainfall          => mvar_data%var(iLookMVAR%scalarRainfall)%dat(1)         ,&  ! intent(in): [dp] rainfall rate (kg m-2 s-1)
 scalarSfcMeltPond       => mvar_data%var(iLookMVAR%scalarSfcMeltPond)%dat(1)      ,&  ! intent(in): [dp] ponded water caused by melt of the "snow without a layer" (kg m-2)

 ! diagnostic variables
 mLayerDepth             => mvar_data%var(iLookMVAR%mLayerDepth)%dat               ,&  ! intent(in): [dp(:)] depth of each layer in the snow-soil sub-domain (m)
 scalarBulkVolHeatCapVeg => mvar_data%var(iLookMVAR%scalarBulkVolHeatCapVeg)%dat(1),&  ! intent(in): [dp   ] bulk volumetric heat capacity of vegetation (J m-3 K-1)
 mLayerVolHtCapBulk      => mvar_data%var(iLookMVAR%mLayerVolHtCapBulk)%dat        ,&  ! intent(in): [dp(:)] bulk volumetric heat capacity in each snow and soil layer (J m-3 K-1)
 mLayerMeltFreeze        => mvar_data%var(iLookMVAR%mLayerMeltFreeze)%dat          ,&  ! intent(out): [dp(:)] melt-freeze in each snow and soil layer (kg m-3)
 mLayerThetaResid        => mvar_data%var(iLookMVAR%mLayerThetaResid)%dat          ,&  ! intent(out): [dp(:)] residual volumetric liquid water content in each snow layer (-)

 ! model fluxes
 iLayerLiqFluxSnow       => mvar_data%var(iLookMVAR%iLayerLiqFluxSnow)%dat         ,&  ! intent(out): [dp(0:)] vertical liquid water flux at snow layer interfaces (-)
 iLayerLiqFluxSoil       => mvar_data%var(iLookMVAR%iLayerLiqFluxSoil)%dat         ,&  ! intent(out): [dp(0:)] vertical liquid water flux at soil layer interfaces (-)
 mLayerBaseflow          => mvar_data%var(iLookMVAR%mLayerBaseflow)%dat            ,&  ! intent(out): [dp(:)] baseflow from each soil layer (m s-1)
 mLayerCompress          => mvar_data%var(iLookMVAR%mLayerCompress)%dat            ,&  ! intent(out): [dp(:)] change in storage associated with compression of the soil matrix (-)
 scalarSoilBaseflow      => mvar_data%var(iLookMVAR%scalarSoilBaseflow)%dat(1)     ,&  ! intent(out): [dp] total baseflow from the soil profile (m s-1)
 scalarExfiltration      => mvar_data%var(iLookMVAR%scalarExfiltration)%dat(1)     ,& ! intent(out):[dp]    exfiltration from the soil profile (m s-1)

 ! sublimation (needed to check mass balance constraints)
 scalarCanopySublimation => mvar_data%var(iLookMVAR%scalarCanopySublimation)%dat(1),&  ! intent(out): [dp] sublimation of ice from the vegetation canopy (kg m-2 s-1)
 scalarSnowSublimation   => mvar_data%var(iLookMVAR%scalarSnowSublimation)%dat(1)  ,&  ! intent(out): [dp] sublimation of ice from the snow surface (kg m-2 s-1)

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
 fImpede                 => mpar_data%var(iLookPARAM%f_impede)                     ,&  ! intent(in): [dp] ice impedance parameter (-)

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

 !print*, 'kl, ku, nBands = ', kl, ku, nBands
 !print*, 'ixSup3, ixSup2, ixSup1 = ', ixSup3, ixSup2, ixSup1
 !print*, 'ixDiag, ixSub1, ixSub2 = ', ixDiag, ixSub1, ixSub2
 !pause

 ! initialize the first flux call
 firstFluxCall=.true.

 ! set the flag to control printing
 printFlagInit=.false.
 printFlag=printFlagInit

 ! set the flag for pausing
 pauseProgress=.false.

 ! identify the matrix solution method
 ! (the type of matrix used to solve the linear system A.X=B)
 if(ixGroundwater==qbaseTopmodel)then
  ixSolve=ixFullMatrix   ! full Jacobian matrix
 else
  ixSolve=ixBandMatrix   ! band-diagonal matrix
 endif
 if(globalPrintFlag) print*, '(ixSolve==ixFullMatrix) = ', (ixSolve==ixFullMatrix)

 ! print states
 !do iLayer=1,nLayers
 ! write(*,'(a10,1x,2(f12.7,1x),f10.3,1x,f17.6,1x,f16.6,1x,f16.6)') 'soil', mvar_data%var(iLookMVAR%iLayerHeight)%dat(iLayer-1), mLayerDepth(iLayer), &
 !  mLayerTemp(iLayer), mLayerVolFracIce(iLayer), mLayerVolFracLiq(iLayer), mLayerMatricHead(iLayer)
 !end do

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

 ! check indices
 if(iJac1 > nState .or. iJac2 > nState)then
  err=20; message=trim(message)//'index iJac1 or iJac2 is out of range'; return
 endif

 ! allocate space for the state vectors
 allocate(stateVecInit(nState),stateVecTrial(nState),stateVecNew(nState),stat=err)
 if(err/=0)then; err=20; message=trim(message)//'unable to allocate space for the state vector'; return; endif

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

 ! allocate space for the band-diagonal matrix that is constructed from the full Jacobian matrix
 if(testBandDiagonal)then
  allocate(aJac_test(nBands,nState),stat=err)
  if(err/=0)then; err=20; message=trim(message)//'unable to allocate space for the band diagonal matrix'; return; endif
 endif

 ! allocate space for the flux vectors and Jacobian matrix
 allocate(dMat(nState),sMul(nState),rAdd(nState),fScale(nState),xScale(nState),fluxVec0(nState),grad(nState),rVec(nState),rhs(nState,nRHS),iPiv(nState),xInc(nState),stat=err)
 if(err/=0)then; err=20; message=trim(message)//'unable to allocate space for the solution vectors'; return; endif

 ! define variables to calculate the numerical Jacobian matrix
 if(numericalJacobian)then
  ! (allocate space for the flux vector and Jacobian matrix
  allocate(fluxVec1(nState),nJac(nState,nState),stat=err)
  if(err/=0)then; err=20; message=trim(message)//'unable to allocate space for the flux vector and numerical Jacobian matrix'; return; endif
 endif  ! if calculating the numerical approximation of the Jacobian matrix

 ! define the index of the top layer
 ixTopNrg = nVegState + 1                       ! energy
 ixTopLiq = nVegState + 2                       ! total water (only snow)
 ixTopMat = nVegState + nSnow*nVarSnowSoil + 2  ! matric head (only soil)

 ! define the indices within the snow-soil domain
 ixSnowSoilNrg = arth(ixTopNrg,nVarSnowSoil,nLayers)  ! energy
 ixSnowSoilWat = arth(ixTopLiq,nVarSnowSoil,nLayers)  ! total water

 ! define indices just for the snow and soil domains
 ixSoilOnlyNrg = arth(ixTopNrg + nSnow*nVarSnowSoil,nVarSnowSoil,nSoil)    ! matric head
 ixSoilOnlyMat = arth(ixTopMat,nVarSnowSoil,nSoil)    ! matric head

 if(nSnow>0)then  ! (liquid water in snow only defined if snow layers exist)
  ixSnowOnlyNrg = arth(ixTopNrg,nVarSnowSoil,nSnow)    ! energy
  ixSnowOnlyWat = arth(ixTopLiq,nVarSnowSoil,nSnow)    ! total water
 endif
 !print*, 'nLayers       = ', nLayers
 !print*, 'nVegState     = ', nVegState
 !print*, 'nSnow, nSoil  = ', nSnow, nSoil
 !print*, 'ixSnowSoilNrg = ', ixSnowSoilNrg
 !print*, 'ixSoilOnlyNrg = ', ixSoilOnlyNrg
 !print*, 'ixSoilOnlyMat = ', ixSoilOnlyMat
 !print*, 'ixSnowOnlyNrg = ', ixSnowOnlyNrg
 !print*, 'ixSnowOnlyWat = ', ixSnowOnlyWat

 ! define the scaled maximum step size (used in the line search)
 stpmax = stepMax*real(nState,dp)

 ! define additional vectors used in the residual calculations
 sMul(:) = 1._dp  ! multiplier for the state vector
 rAdd(:) = 0._dp  ! additional terms in the residual calculations (phase change, compressibility, etc.)

 ! define the multiplier for the state vector for residual calculations (vegetation canopy)
 if(computeVegFlux)then
  sMul(ixCasNrg) = Cp_air*iden_air          ! volumetric heat capacity of air (J m-3 K-1)
  sMul(ixVegNrg) = scalarBulkVolHeatCapVeg  ! volumetric heat capacity of the vegetation (J m-3 K-1)
  sMul(ixVegWat) = 1._dp                    ! nothing else on the left hand side
 endif

 ! define the multiplier for the state vector for residual calculations (snow-soil domain)
 sMul(ixSnowSoilNrg) = mLayerVolHtCapBulk(1:nLayers)
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

 ! initialize
 xInc(:)   = 0._dp  ! iteration increment

 ! compute the total water in the vegetation canopy
 if(computeVegFlux)then
  scalarCanopyWat = scalarCanopyLiq + scalarCanopyIce   ! kg m-2
 endif

 ! compute the total water in snow
 if(nSnow>0)&
  mLayerVolFracWat(1:nSnow) = mLayerVolFracLiq(1:nSnow) + mLayerVolFracIce(1:nSnow)*(iden_ice/iden_water)

 ! define the scaling for the function evaluation -- vegetation
 if(computeVegFlux)then
  fScale(ixCasNrg) = fScaleNrg   ! (J m-3)
  fScale(ixVegNrg) = fScaleNrg   ! (J m-3)
  fScale(ixVegWat) = fScaleLiq*canopyDepth*iden_water  ! (kg m-2)
 endif

 ! define the scaling for the function evaluation -- snow and soil
 fScale(ixSnowSoilNrg) = fScaleNrg  ! (J m-3)
 fScale(ixSnowSoilWat) = fScaleLiq  ! (-)

 ! define scaling for the state vector -- vegetation
 if(computeVegFlux)then
  xScale(ixCasNrg) = xScaleTemp   ! (K)
  xScale(ixVegNrg) = xScaleTemp   ! (K)
  xScale(ixVegWat) = xScaleLiq*canopyDepth*iden_water  ! (kg m-2)
 endif

 ! define the scaling for the function evaluation -- snow and soil
 xScale(ixSnowSoilNrg) = xScaleTemp  ! (K)
 xScale(ixSnowOnlyWat) = xScaleLiq   ! (-)
 xScale(ixSoilOnlyMat) = xScaleMat   ! (m)

 ! build the state vector for the vegetation canopy
 if(computeVegFlux)then
  stateVecInit(ixCasNrg) = scalarCanairTemp
  stateVecInit(ixVegNrg) = scalarCanopyTemp
  stateVecInit(ixVegWat) = scalarCanopyWat  ! kg m-2
 endif

 ! build the state vector for the snow and soil domain
 stateVecInit(ixSnowSoilNrg) = mLayerTemp(1:nLayers)
 stateVecInit(ixSoilOnlyMat) = mLayerMatricHead(1:nSoil)
 if(nSnow>0)&
 stateVecInit(ixSnowOnlyWat) = mLayerVolFracWat(1:nSnow)

 ! initialize the trial state vectors
 stateVecTrial = stateVecInit

 ! need to intialize canopy water at a positive value
 if(computeVegFlux)then
  if(scalarCanopyWat < xMinCanopyWater) stateVecTrial(ixVegWat) = scalarCanopyWat + xMinCanopyWater
 endif

 ! initialize the volumetric fraction of liquid water and ice in the vegetation canopy
 !print*, 'scalarCanopyIce = ', scalarCanopyIce
 !scalarCanopyLiqTrial = scalarCanopyLiq
 !scalarCanopyIceTrial = scalarCanopyIce

 ! initialize the volumetric fraction of liquid water and ice in snow and soil layers
 !mLayerVolFracLiqTrial(1:nLayers) = mLayerVolFracLiq(1:nLayers)  ! additional state variable for all layers
 !mLayerVolFracIceTrial(1:nLayers) = mLayerVolFracIce(1:nLayers)  ! additional state variable for all layers

 ! initialize the function variable
 fOld=veryBig  ! initialize to a very big number

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
  !print*, '***'
  !print*, '***'
  !print*, '***'
  !print*, '***'
  !print*, '***'
  !print*, '***'
  !print*, '***'
  !write(*,'(a,1x,f10.2,1x,2(i4,1x),l1)') '*** new iteration: dt, iter, nstate, computeVegFlux = ', dt, iter, nstate, computeVegFlux
  !write(*,'(a,1x,10(e15.5,1x))') 'stateVecInit(1:10)  = ', stateVecInit(1:10)
  !write(*,'(a,1x,10(e15.5,1x))') 'stateVecTrial(1:10) = ', stateVecTrial(1:10)
  !write(*,'(a,1x,10(e15.5,1x))') 'xInc(1:10)          = ', xInc(1:10)

  ! -----
  ! * compute model fluxes and residual
  !    NOTE: refine residual with line search...
  ! --------------------------------------------
  call lineSearch(&
                  ! input
                  (iter>1),                & ! intent(in): flag to denote the need to perform line search
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

  ! exit iteration loop if converged
  if(converged) exit

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
  !do iLayer=1,nLayers
  ! write(*,'(a,1x,i4,1x,100(e15.5,1x))') 'iLayer, LH_fus*iden_water*mLayerdTheta_dTk(iLayer) = ', iLayer, LH_fus*iden_water*mLayerdTheta_dTk(iLayer)
  !end do
  !write(*,'(a,1x,100(e15.5,1x))')'dMat(ixSoilOnlyNrg) = ', dMat(ixSoilOnlyNrg)

  ! compute additional terms for the Jacobian for the soil domain (excluding fluxes)
  if(ixRichards==moisture)then; err=20; message=trim(message)//'have not implemented the moisture-based form of RE yet'; return; endif
  dMat(ixSoilOnlyMat) = dVolTot_dPsi0(1:nSoil) + dCompress_dPsi(1:nSoil)
  !print*, 'dVolTot_dPsi0(1:nSoil) = ', dVolTot_dPsi0(1:nSoil)
  !print*, 'dCompress_dPsi(1:nSoil)    = ', dCompress_dPsi(1:nSoil)

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

  !if(computeVegFlux .and. printFlag)then
  ! write(*,'(a,1x,10(e15.5,1x))') 'xInc(ixCasNrg)      = ', xInc(ixCasNrg)
  ! write(*,'(a,1x,10(e15.5,1x))') 'xInc(ixVegNrg)      = ', xInc(ixVegNrg)
  ! write(*,'(a,1x,10(e15.5,1x))') 'xInc(ixVegWat)      = ', xInc(ixVegWat)
  !endif

  !write(*,'(a,1x,10(e15.5,1x))') 'rVec(ixSoilOnlyMat) = ', rVec(ixSoilOnlyMat)
  !write(*,'(a,1x,10(e15.5,1x))') 'grad(ixSoilOnlyMat) = ', grad(ixSoilOnlyMat)

  !if(printFlag)then
  ! write(*,'(a,1x,10(e15.5,1x))') 'xInc(ixSoilOnlyMat) = ', xInc(ixSoilOnlyMat)
  ! write(*,'(a,1x,10(e15.5,1x))') 'xInc(ixSnowOnlyWat) = ', xInc(ixSnowOnlyWat)
  ! write(*,'(a,1x,10(e15.5,1x))') 'xInc(ixSnowSoilNrg) = ', xInc(ixSnowSoilNrg)
  ! pause
  ! if(pauseProgress) pause
  !endif

  ! print iteration increment
  !write(*,'(a,1x,10(f20.12,1x))') 'xInc(iJac1:iJac2)              = ', xInc(iJac1:iJac2)

  ! -----
  ! * impose solution constraints...
  ! --------------------------------

  ! ** limit temperature increment to 1K

  ! vegetation
  if(computeVegFlux)then
   if(abs(xInc(ixVegNrg)) > 1._dp)then
    !write(*,'(a,1x,10(f20.12,1x))') 'before scale: xInc(iJac1:iJac2) = ', xInc(iJac1:iJac2)
    xIncScale = abs(1._dp/xInc(ixVegNrg))  ! scaling factor for the iteration increment (-)
    xInc      = xIncScale*xInc             ! scale iteration increments
    !write(*,'(a,1x,10(f20.12,1x))') 'after scale: xInc(iJac1:iJac2) = ', xInc(iJac1:iJac2)
   endif
  endif

  ! snow and soil
  if(any(abs(xInc(ixSnowSoilNrg)) > 1._dp))then
   !write(*,'(a,1x,10(f20.12,1x))') 'before scale: xInc(iJac1:iJac2) = ', xInc(iJac1:iJac2)
   iMax      = maxloc( abs(xInc(ixSnowSoilNrg)) )                   ! index of maximum temperature increment
   xIncScale = abs( 1._dp/xInc(ixSnowSoilNrg(iMax(1))) )            ! scaling factor for the iteration increment (-)
   xInc      = xIncScale*xInc
   !write(*,'(a,1x,10(f20.12,1x))') 'after scale: xInc(iJac1:iJac2) = ', xInc(iJac1:iJac2)
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
    xIncScale   = cInc/xInc(ixVegNrg)  ! scaling factor for the iteration increment (-)
    xInc        = xIncScale*xInc       ! scale iteration increments
   endif

   !print*, 'crosTempVeg = ', crosTempVeg

   ! --------------------------------------------------------------------------------------------------------------------
   ! canopy liquid water

   ! check if new value of storage will be negative
   if(stateVecTrial(ixVegWat)+xInc(ixVegWat) < 0._dp)then
    ! scale iteration increment
    cInc      = -0.5_dp*stateVecTrial(ixVegWat)                                  ! constrained iteration increment (K) -- simplified bi-section
    xIncScale = cInc/xInc(ixVegWat)                                              ! scaling factor for the iteration increment (-)
    xInc      = xIncScale*xInc                                                   ! new iteration increment
    !print*, 'canopy liquid water constraint'
   endif

  endif  ! if computing fluxes through vegetation

  ! ** impose solution constraints for snow
  if(nSnow > 0)then

   ! --------------------------------------------------------------------------------------------------------------------
   ! get new temperatures
   mLayerTempCheck = stateVecTrial(ixSnowOnlyNrg) + xInc(ixSnowOnlyNrg)

   ! - check sub-freezing temperatures for snow
   if(any(mLayerTempCheck > Tfreeze))then
    ! scale iteration increment
    iMax      = maxloc(mLayerTempCheck)                                          ! index of maximum temperature
    cInc      = 0.5_dp*(Tfreeze - stateVecTrial(ixSnowOnlyNrg(iMax(1))) )        ! constrained temperature increment (K) -- simplified bi-section
    xIncScale = cInc/xInc(ixSnowOnlyNrg(iMax(1)))                                ! scaling factor for the iteration increment (-)
    xInc      = xIncScale*xInc
    !print*, 'stateVecTrial(ixSnowOnlyNrg(iMax(1))), mLayerTempCheck(iMax(1)), cInc, xIncScale = ', &
    !         stateVecTrial(ixSnowOnlyNrg(iMax(1))), mLayerTempCheck(iMax(1)), cInc, xIncScale
   endif   ! if snow temperature > freezing

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
    !write(*,'(a,1x,i4,1x,l1,1x,10(f15.8,1x))') 'iLayer, drainFlag(iLayer), xInc(ixSnowOnlyWat(iLayer)), mLayerVolFracLiqTrial(iLayer), mLayerThetaResid(iLayer) = ',&
    !                                            iLayer, drainFlag(iLayer), xInc(ixSnowOnlyWat(iLayer)), mLayerVolFracLiqTrial(iLayer), mLayerThetaResid(iLayer)
   end do

   ! check if the iteration increment removes all the water
   !if(any(mLayerVolFracLiqCheck < 0._dp))then
   ! ! print original iteration increment
   ! do iLayer=1,nSnow
   !  write(*,'(a,1x,i4,1x,10(f15.8,1x))') 'iLayer, xInc(ixSnowOnlyWat(iLayer)) = ', iLayer, xInc(ixSnowOnlyWat(iLayer))
   ! end do
   ! ! scale iteration increment
   ! iMin      = minloc(mLayerVolFracLiqCheck)                 ! index of the most excessive drainage
   ! cInc      = -0.5_dp*mLayerVolFracLiqTrial(iMin(1))        ! constrained drainage increment (-) -- simplified bi-secion
   ! xIncScale = cInc/xInc(ixSnowOnlyWat(iMin(1)))        ! scaling factor for the iteration increment (-)
   ! xInc      = xIncScale*xInc
   ! drainFlag(iMin(1)) = .true.
   ! ! print results
   ! do iLayer=1,nSnow
   !  write(*,'(a,1x,i4,1x,l1,1x,10(f15.8,1x))') 'iLayer, drainFlag(iLayer), xInc(ixSnowOnlyWat(iLayer)), mLayerVolFracLiqTrial(iLayer), mLayerThetaResid(iLayer) = ',&
   !                                              iLayer, drainFlag(iLayer), xInc(ixSnowOnlyWat(iLayer)), mLayerVolFracLiqTrial(iLayer), mLayerThetaResid(iLayer)
   ! end do
   ! !pause
   !endif   ! if iteration increment removes all the water

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

   !write(*,'(i4,3x,a20,1x,f20.10)') iLayer, ' - ', xInc(ixNrg)

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
   !write(*,'(i4,1x,l1,1x,2(f20.10,1x))') iLayer, crosFlag(iLayer), TcSoil, xInc(ixNrg)

   ! place constraint for matric head
   if(xInc(ixLiq) > 1._dp .and. stateVecTrial(ixLiq) > 0._dp)then
    xInc(ixLiq) = 1._dp
    pauseProgress=.true.
   endif  ! if constraining matric head

  end do  ! (loop through soil layers

  !print*, ' SWE = ', sum( (mLayerVolFracLiqTrial(1:nSnow)*iden_water + mLayerVolFracIceTrial(1:nSnow)*iden_ice) * mLayerDepth(1:nSnow) )

  ! check convergence
  if(niter==maxiter)then; err=-20; message=trim(message)//'failed to converge'; return; endif
  !pause 'iterating'


 end do  ! iterating
 !pause 'after iterations'

 ! check that we got baseflow
 !print*, 'mLayerBaseflow(:) = ', mLayerBaseflow(:)

 ! -----
 ! * update states and compute total volumetric melt...
 ! ----------------------------------------------------

 ! update temperatures (ensure new temperature is consistent with the fluxes)
 stateVecTrial(ixSnowSoilNrg) = stateVecInit(ixSnowSoilNrg) + (fluxVec0(ixSnowSoilNrg)*dt + rAdd(ixSnowSoilNrg))/sMul(ixSnowSoilNrg)

 ! update volumetric liquid water content in the snow (ensure change in state is consistent with the fluxes)
 ! NOTE: should really update fluxes based on the iteration increment (as for soil below)
 if(nSnow>0)&
 stateVecTrial(ixSnowOnlyWat) = stateVecInit(ixSnowOnlyWat) + (fluxVec0(ixSnowOnlyWat)*dt + rAdd(ixSnowOnlyWat))

 ! update fluxes in the soil layers
 ! (extract iteration increment for the matric head)
 mLayerMatricHeadDiff = xInc(ixSoilOnlyMat)
 ! (loop through soil layers)
 do iLayer=1,nSoil
  ! (update vertical fluxes)
  if(iLayer==0)then;         iLayerLiqFluxSoil(iLayer) = iLayerLiqFluxSoil(iLayer) + dq_dHydStateBelow(iLayer)*mLayerMatricHeadDiff(iLayer+1)
  elseif(iLayer==nSoil)then; iLayerLiqFluxSoil(iLayer) = iLayerLiqFluxSoil(iLayer) + dq_dHydStateAbove(iLayer)*mLayerMatricHeadDiff(iLayer)
  else;                      iLayerLiqFluxSoil(iLayer) = iLayerLiqFluxSoil(iLayer) + dq_dHydStateAbove(iLayer)*mLayerMatricHeadDiff(iLayer) &
                                                                                   + dq_dHydStateBelow(iLayer)*mLayerMatricHeadDiff(iLayer+1)
  endif
  !print*, '**'
  ! (update the baseflow fluxes)
  !do jLayer=1,nSoil
  ! !write(*,'(a,1x,2(i4,1x),10(e20.10,1x))') 'iLayer, jLayer, dBaseflow_dMatric(jLayer,iLayer), dBaseflow_dMatric(iLayer,jLayer), dVolTot_dPsi0(iLayer), mLayerMatricHeadDiff(iLayer) = ', &
  ! !                                          iLayer, jLayer, dBaseflow_dMatric(jLayer,iLayer), dBaseflow_dMatric(iLayer,jLayer), dVolTot_dPsi0(iLayer), mLayerMatricHeadDiff(iLayer)
  ! mLayerBaseflow(iLayer) = mLayerBaseflow(iLayer) + dBaseflow_dMatric(jLayer,iLayer)*mLayerMatricHeadDiff(iLayer)
  !end do  ! looping through soil layers
 end do  ! looping through soil layers

 ! compute total baseflow from the soil zone (needed for mass balance checks)
 scalarSoilBaseflow = sum(mLayerBaseflow)
 !write(*,'(a,1x,e20.10)') 'scalarSoilBaseflow = ', scalarSoilBaseflow

 ! update states
 call updatState(&
                 stateVecTrial,         & ! intent(in): full state vector (mixed units)
                 mLayerVolFracLiqTrial, & ! intent(out): volumetric fraction of liquid water (-)
                 mLayerVolFracIceTrial, & ! intent(out): volumetric fraction of ice (-)
                 scalarCanopyLiqTrial,  & ! intent(out): mass of canopy liquid (kg m-2)
                 scalarCanopyIceTrial,  & ! intent(out): mass of canopy ice (kg m-2)
                 err,cmessage)            ! intent(out): error code and error message
 if(err/=0)then; message=trim(message)//trim(cmessage); return; endif  ! (check for errors)

 ! check the mass balance for each soil layer
 if(checkMassBalance)then
  do iLayer=1,nSoil
   balance0 = (mLayerVolFracLiq(nSnow+iLayer)      + mLayerVolFracIce(nSnow+iLayer)      )*mLayerDepth(nSnow+iLayer)   ! dimensionless --> m
   balance1 = (mLayerVolFracLiqTrial(nSnow+iLayer) + mLayerVolFracIceTrial(nSnow+iLayer) )*mLayerDepth(nSnow+iLayer)   ! dimensionless --> m
   vertFlux = dt * -(iLayerLiqFluxSoil(iLayer) - iLayerLiqFluxSoil(iLayer-1))  ! m s-1 --> m
   tranSink = dt*mLayerTranspire(iLayer)                                       ! m s-1 --> m
   baseSink = dt*mLayerBaseflow(iLayer)                                        ! m s-1 --> m
   compSink = mLayerCompress(iLayer)*mLayerDepth(nSnow+iLayer)                 ! dimensionless --> m
   liqError = balance1 - (balance0 + vertFlux + tranSink - baseSink - compSink)
   write(*,'(a,1x,e20.10,1x,10(f20.10,1x))') 'liqError, balance0, balance1, vertFlux, tranSink, baseSink, compSink = ', &
                                              liqError, balance0, balance1, vertFlux, tranSink, baseSink, compSink
  end do  ! looping through soil layers
 endif  ! checking mass balance

 ! compute the melt in each snow and soil layer
 if(nSnow>0) mLayerMeltFreeze(      1:nSnow  ) = -(mLayerVolFracIceTrial(      1:nSnow  ) - mLayerVolFracIce(      1:nSnow  ))*iden_ice
             mLayerMeltFreeze(nSnow+1:nLayers) = -(mLayerVolFracIceTrial(nSnow+1:nLayers) - mLayerVolFracIce(nSnow+1:nLayers))*iden_water
 !write(*,'(a,1x,10(f20.10,1x))') 'mLayerVolFracIce(1),      mLayerVolFracLiq(1)      = ', mLayerVolFracIce(1)*iden_ice, mLayerVolFracLiq(1)*iden_water
 !write(*,'(a,1x,10(f20.10,1x))') 'mLayerVolFracIceTrial(1), mLayerVolFracLiqTrial(1) = ', mLayerVolFracIceTrial(1)*iden_ice, mLayerVolFracLiqTrial(1)*iden_water
 !write(*,'(a,1x,10(f20.10,1x))') 'mLayerMeltFreeze(      1:nSnow  ) = ', mLayerMeltFreeze(      1:nSnow  )

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
 mLayerMatricHead(1:nSoil) = stateVecTrial(ixSoilOnlyMat)

 ! save the volumetric liquid water and ice content
 mLayerVolFracLiq = mLayerVolFracLiqTrial  ! computed in updatState
 mLayerVolFracIce = mLayerVolFracIceTrial  ! computed in updatState
 !print*, 'mLayerVolFracLiq = ', mLayerVolFracLiq

 ! ==========================================================================================================================================
 ! ==========================================================================================================================================
 ! ==========================================================================================================================================
 ! ==========================================================================================================================================

 ! deallocate space for the state vectors etc.
 deallocate(stateVecInit,stateVecTrial,stateVecNew,dMat,sMul,rAdd,fScale,xScale,fluxVec0,aJac,grad,rVec,rhs,iPiv,xInc,stat=err)
 if(err/=0)then; err=20; message=trim(message)//'unable to deallocate space for the state/flux vectors and analytical Jacobian matrix'; return; endif

 ! deallocate space for the baseflow derivatives
 if(ixGroundwater==qbaseTopmodel)then
  deallocate(dBaseflow_dMatric,stat=err)
  if(err/=0)then; err=20; message=trim(message)//'unable to deallocate space for the baseflow derivatives'; return; endif
 endif

 ! deallocate space for the variables used to create the numerical Jacobian matrix
 if(numericalJacobian)then
  deallocate(fluxVec1,nJac,stat=err)
  if(err/=0)then; err=20; message=trim(message)//'unable to allocate space for the flux vector and numerical Jacobian matrix'; return; endif
 endif

 if(testBandDiagonal)then
  deallocate(aJac_test,stat=err)
  if(err/=0)then; err=20; message=trim(message)//'unable to allocate space for the band diagonal matrix'; return; endif
 endif

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
  real(dp),dimension(nSoil)      :: mLayerPsiLiq              ! liquid water matric potential (m)
  ! initialize error control
  err=0; message='updatState/'

  ! get the necessary variables from the data structures
  associate(&

  ! layer type (snow or soil)
  layerType               => indx_data%var(iLookINDEX%layerType)%dat                ,&  ! intent(in): [i4b(:)] type of each layer in the snow+soil domain (snow or soil)

  ! layer depth
  mLayerDepth             => mvar_data%var(iLookMVAR%mLayerDepth)%dat               ,&  ! intent(in): [dp] depth of each layer (m)

  ! snow parameters
  snowfrz_scale           => mpar_data%var(iLookPARAM%snowfrz_scale)                ,&  ! intent(in): [dp] scaling parameter for the snow freezing curve (K-1)

  ! model state variables (vegetation canopy)
  scalarCanopyIce         => mvar_data%var(iLookMVAR%scalarCanopyIce)%dat(1)        ,&  ! intent(inout): [dp] mass of ice on the vegetation canopy (kg m-2)
  scalarCanopyLiq         => mvar_data%var(iLookMVAR%scalarCanopyLiq)%dat(1)        ,&  ! intent(inout): [dp] mass of liquid water on the vegetation canopy (kg m-2)

  ! soil parameters
  vGn_alpha               => mpar_data%var(iLookPARAM%vGn_alpha)                    ,&  ! intent(in): [dp] van Genutchen "alpha" parameter (m-1)
  vGn_n                   => mpar_data%var(iLookPARAM%vGn_n)                        ,&  ! intent(in): [dp] van Genutchen "n" parameter (-)
  vGn_m                   => mvar_data%var(iLookMVAR%scalarVGn_m)%dat(1)            ,&  ! intent(in): [dp] van Genutchen "m" parameter (-)
  theta_sat               => mpar_data%var(iLookPARAM%theta_sat)                    ,&  ! intent(in): [dp] soil porosity (-)
  theta_res               => mpar_data%var(iLookPARAM%theta_res)                     &  ! intent(in): [dp] soil residual volumetric water content (-)
  )

  ! update states for the vegetation canopy
  if(computeVegFlux)then
   fracLiqVeg    = fracliquid(stateVecTrial(ixVegNrg),snowfrz_scale)  ! fraction of liquid water (-)
   totalWaterVeg = stateVecTrial(ixVegWat)                            ! total water (kg m-2)
   scalarCanopyLiqTrial = fracLiqVeg*totalWaterVeg                    ! mass of liquid water on the canopy (kg m-2)
   scalarCanopyIceTrial = (1._dp - fracLiqVeg)*totalWaterVeg          ! mass of ice on the canopy (kg m-2)
   !write(*,'(a,1x,10(f20.15,1x))') 'fracLiqVeg, totalWaterVeg, stateVecTrial(ixVegWat), scalarCanopyLiqTrial, scalarCanopyIceTrial = ', &
   !                                 fracLiqVeg, totalWaterVeg, stateVecTrial(ixVegWat), scalarCanopyLiqTrial, scalarCanopyIceTrial

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

     !if(printflag .and. iLayer < 5)then
     ! if(iLayer==1) write(*,'(a,1x)')     'iLayer, ixSnowOnlyWat(iLayer), mLayerDepth(iLayer), mLayerVolFracLiqTrial(iLayer), mLayerVolFracIceTrial(iLayer), stateVecTrial(ixSnowSoilNrg(iLayer)), stateVecTrial(ixSnowSoilWat(iLayer)) = '
     ! write(*,'(2(i4,1x),10(f20.10,1x))')  iLayer, ixSnowOnlyWat(iLayer), mLayerDepth(iLayer), mLayerVolFracLiqTrial(iLayer), mLayerVolFracIceTrial(iLayer), stateVecTrial(ixSnowSoilNrg(iLayer)), stateVecTrial(ixSnowSoilWat(iLayer))
     !endif

    !** soil
    case(ix_soil)

     !if(printflag)&
     !write(*,'(a,1x,2(i4,1x),2(f20.10,1x))') 'iLayer, ixSnowSoilWat(iLayer), mLayerVolFracLiqTrial(iLayer), mLayerVolFracIceTrial(iLayer) = ', &
     !                                         iLayer, ixSnowSoilWat(iLayer), mLayerVolFracLiqTrial(iLayer), mLayerVolFracIceTrial(iLayer)

     call updateSoil(&
                     ! input
                     stateVecTrial(ixSnowSoilNrg(iLayer)),      & ! intent(in): layer temperature (K)
                     stateVecTrial(ixSoilOnlyMat(iLayer-nSnow)),& ! intent(in): matric head (m)
                     vGn_alpha,vGn_n,theta_sat,theta_res,vGn_m, & ! intent(in): van Genutchen soil parameters
                     ! output
                     mLayerPsiLiq(iLayer-nSnow),                & ! intent(out): liquid water matric potential
                     mLayerVolFracLiqTrial(iLayer),             & ! intent(out): volumetric fraction of liquid water (-)
                     mLayerVolFracIceTrial(iLayer),             & ! intent(out): volumetric fraction of ice (-)
                     err,cmessage)                                ! intent(out): error control
     if(err/=0)then; message=trim(message)//trim(cmessage); return; endif  ! (check for errors)

     !if(printFlag)&
     !if(iLayer==1)     write(*,'(a)')        'stateVecTrial(ixSnowSoilNrg(iLayer)), stateVecTrial(ixSoilOnlyMat(iLayer-nSnow)), mLayerPsiLiq(iLayer-nSnow), mLayerVolFracLiqTrial(iLayer), mLayerVolFracIceTrial(iLayer) = '
     !write(*,'(i4,1x,10(f20.10,1x))') iLayer, stateVecTrial(ixSnowSoilNrg(iLayer)), stateVecTrial(ixSoilOnlyMat(iLayer-nSnow)), mLayerPsiLiq(iLayer-nSnow), mLayerVolFracLiqTrial(iLayer), mLayerVolFracIceTrial(iLayer)

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
  real(qp),intent(out)           :: rVec(:)                   ! residual vector (mixed units)
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
  real(dp),dimension(nSnow)      :: mLayerVolFracWatTrial     ! trial value for volumetric fraction of total water (-)
  real(dp),dimension(nSoil)      :: mLayerMatricHeadTrial     ! trial value for matric head (m)
  ! temporary vectors for the soil sub-domain
  real(dp),dimension(nSoil)      :: vThetaInit                ! liquid equivalent of total water at the start of the step
  real(dp),dimension(nSoil)      :: vThetaTrial               ! liquid equivalent of total water at the current iteration
  ! variables for testing
  real(dp)                       :: xCompress                 ! compression in a given layer (m)
  real(dp)                       :: xFlux0,xFlux1             ! fluxes at the layer boundaries (m)
  real(dp)                       :: xBalance                  ! water balance (m)
  ! initialize error control
  err=0; message='xFluxResid/'

  ! -----
  ! * associate desired variables from data structures...
  ! -----------------------------------------------------
  associate(&
  ! model decisions
  ixRichards              => model_decisions(iLookDECISIONS%f_Richards)%iDecision   ,&  ! intent(in): [i4b] index of the form of Richards' equation

  ! soil parameters
  vGn_alpha               => mpar_data%var(iLookPARAM%vGn_alpha)                    ,&  ! intent(in): [dp] van Genutchen "alpha" parameter (m-1)
  vGn_n                   => mpar_data%var(iLookPARAM%vGn_n)                        ,&  ! intent(in): [dp] van Genutchen "n" parameter (-)
  vGn_m                   => mvar_data%var(iLookMVAR%scalarVGn_m)%dat(1)            ,&  ! intent(in): [dp] van Genutchen "m" parameter (-)
  theta_sat               => mpar_data%var(iLookPARAM%theta_sat)                    ,&  ! intent(in): [dp] soil porosity (-)
  theta_res               => mpar_data%var(iLookPARAM%theta_res)                    ,&  ! intent(in): [dp] soil residual volumetric water content (-)
  specificStorage         => mpar_data%var(iLookPARAM%specificStorage)              ,&  ! intent(in): [dp] specific storage coefficient (m-1)

  ! layer depth
  mLayerDepth             => mvar_data%var(iLookMVAR%mLayerDepth)%dat               ,&  ! intent(in): [dp(:)] depth of each layer in the snow-soil sub-domain (m)

  ! model fluxes
  mLayerBaseflow          => mvar_data%var(iLookMVAR%mLayerBaseflow)%dat,            &  ! intent(out): [dp(:)] baseflow from each soil layer (m s-1)
  mLayerCompress          => mvar_data%var(iLookMVAR%mLayerCompress)%dat,            &  ! intent(out): [dp(:)] change in storage associated with compression of the soil matrix (-)
  scalarSoilCompress      => mvar_data%var(iLookMVAR%scalarSoilCompress)%dat(1),     &  ! intent(out): [dp] total change in storage associated with compression of the soil matrix (kg m-2)
  iLayerLiqFluxSoil       => mvar_data%var(iLookMVAR%iLayerLiqFluxSoil)%dat,         &  ! intent(out): [dp] liquid soil fluxes (m s-1)

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

  ! extract state variables for the snow and soil domain
  mLayerTempTrial(1:nLayers)     = stateVec(ixSnowSoilNrg)
  mLayerMatricHeadTrial(1:nSoil) = stateVec(ixSoilOnlyMat)
  if(nSnow>0)&
   mLayerVolFracWatTrial(1:nSnow) = stateVec(ixSnowOnlyWat)
  ! (test)
  !if(printFlag)then
  ! write(*,'(a,1x,f20.10)') 'iden_water*mLayerVolFracWatTrial(1:nSnow)*mLayerDepth(1:nSnow) = ', &
  !                           iden_water*mLayerVolFracWatTrial(1:nSnow)*mLayerDepth(1:nSnow)
  !endif

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

  !if(printFlag)then
   !write(*,'(a,1x,100(e25.15,1x))') 'mLayerVolFracLiqLocal(1), mLayerVolFracLiqLocal(1) - mLayerVolFracLiq(1) = ', mLayerVolFracLiqLocal(1), mLayerVolFracLiqLocal(1) - mLayerVolFracLiq(1)
   !write(*,'(a,1x,100(e25.15,1x))') 'mLayerVolFracIceLocal(1), mLayerVolFracIceLocal(1) - mLayerVolFracIce(1) = ', mLayerVolFracIceLocal(1), (mLayerVolFracIceLocal(1) - mLayerVolFracIce(1))*(iden_ice/iden_water)
   !write(*,'(a,1x,100(e25.15,1x))') 'mLayerMatricHeadTrial = ', mLayerMatricHeadTrial
   !write(*,'(a,1x,10(e15.5,1x))') 'mLayerMatricHeadTrial(1:13) = ', mLayerMatricHeadTrial(1:13)
   !write(*,'(a,1x,10(e15.5,1x))') 'fVec(ixSnowSoilNrg) = ', fVec(ixSnowSoilNrg)
   !write(*,'(a,1x,10(e15.5,1x))') 'fVec(ixSnowSoilWat) = ', fVec(ixSnowSoilWat)
  !endif

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
  !print*, 'scalarSoilCompress = ', scalarSoilCompress

  ! -----
  ! * compute residual vector...
  ! ----------------------------

  !write(*,'(a,1x,10(f25.15,1x))') 'scalarCanopyLiqLocal - scalarCanopyLiq = ', scalarCanopyLiqLocal - scalarCanopyLiq
  !write(*,'(a,1x,10(f25.15,1x))') 'scalarCanopyIceLocal - scalarCanopyIce = ', scalarCanopyIceLocal - scalarCanopyIce
  !write(*,'(a,1x,10(f25.15,1x))') 'change in energy (J m-3)                  = ', LH_fus*(scalarCanopyIceLocal - scalarCanopyIce)/canopyDepth  ! J m-3

  ! intialize additional terms on the RHS as zero
  rAdd(:) = 0._dp

  ! compute energy associated with melt freeze for the vegetation canopy
  if(computeVegFlux)then
   rAdd(ixVegNrg) = rAdd(ixVegNrg) + LH_fus*(scalarCanopyIceLocal - scalarCanopyIce)/canopyDepth   ! energy associated with melt/freeze (J m-3)
   !if(printFlag)then
   ! print*, 'rAdd(ixVegNrg), scalarCanopyIceLocal, scalarCanopyIce = ', rAdd(ixVegNrg), scalarCanopyIceLocal, scalarCanopyIce
   ! pause
   !endif
  endif

  ! compute energy associated with melt/freeze for snow
  if(nSnow>0)&
  rAdd(ixSnowOnlyNrg) = rAdd(ixSnowOnlyNrg) + LH_fus*iden_ice*(mLayerVolFracIceLocal(1:nSnow) - mLayerVolFracIce(1:nSnow))       ! energy associated with melt/freeze (J m-3)
  if(printFlag)then
   !write(*,'(a,1x,10(e20.10,1x))') 'rAdd(ixSnowOnlyNrg) = ', rAdd(ixSnowOnlyNrg)
   !write(*,'(a,1x,10(e20.10,1x))') 'mLayerVolFracIce(1:5) = ', mLayerVolFracIce(1:5)
   !write(*,'(a,1x,10(e20.10,1x))') 'mLayerVolFracIceLocal(1:5) = ', mLayerVolFracIceLocal(1:5)
   !write(*,'(a,1x,10(e20.10,1x))') 'delIce = ', mLayerVolFracIceLocal(1:5) - mLayerVolFracIce(1:5)
  endif

  ! compute energy associated with melt/freeze for soil
  rAdd(ixSoilOnlyNrg) = rAdd(ixSoilOnlyNrg) + LH_fus*iden_water*(mLayerVolFracIceLocal(nSnow+1:nLayers) - mLayerVolFracIce(nSnow+1:nLayers))     ! energy associated with melt/freeze (J m-3)
  !if(printFlag)then
  ! write(*,'(a,1x,10(e30.15,1x))') 'vIce01 = ', mLayerVolFracIce(nSnow+1:nLayers)
  ! write(*,'(a,1x,10(e30.15,1x))') 'vIce02 = ', mLayerVolFracIceLocal(nSnow+1:nLayers)
  ! write(*,'(a,1x,10(e15.5,1x))') 'delIce = ', mLayerVolFracIceLocal(nSnow+1:nLayers) - mLayerVolFracIce(nSnow+1:nLayers)
  ! write(*,'(a,1x,10(e15.5,1x))') 'delNrg = ', LH_fus*iden_water*(mLayerVolFracIceLocal(nSnow+1:nLayers) - mLayerVolFracIce(nSnow+1:nLayers))
  ! print*, 'precision(mLayerVolFracIce) = ', precision(mLayerVolFracIce)
  !endif

  ! sink terms (-)
  ! NOTE: state variable is volumetric water content, so melt-freeze is not included
  ! NOTE: ground evaporation was already included in the flux at the upper boundary
  !print*, 'mLayerTranspire(1:nSoil) = ', mLayerTranspire(1:nSoil)
  !print*, 'mLayerBaseflow(1:nSoil)  = ', mLayerBaseflow(1:nSoil)
  !print*, 'mLayerCompress(1:nSoil)  = ', mLayerCompress(1:nSoil)

  !mLayerCompress(:)  = 0._dp
  !mLayerTranspire(:) = 0._dp
  !mLayerBaseflow(:) = 0._dp
  rAdd(ixSoilOnlyMat)    = rAdd(ixSoilOnlyMat) + dt*(mLayerTranspire(1:nSoil) - mLayerBaseflow(1:nSoil) )/mLayerDepth(nSnow+1:nLayers) - mLayerCompress(1:nSoil)
  !print*, 'rAdd(ixSoilOnlyMat)      = ', rAdd(ixSoilOnlyMat)

  ! liquid water equivalent of melt/freeze for snow layers (-)
  ! NOTE: state equation for soil is based on the total equivalent liquid water content (liquid plus ice)
  !if(nSnow>0)&
  !rAdd(ixSnowOnlyWat) = rAdd(ixSnowOnlyWat) - (iden_ice/iden_water)*(mLayerVolFracIceLocal(1:nSnow) - mLayerVolFracIce(1:nSnow)) ! liquid water equivalent of melt/freeze (-)
  !if(printFlag)then
  ! write(*,'(a,1x,10(e20.10,1x))') 'rAdd(ixSnowOnlyWat) = ', rAdd(ixSnowOnlyWat)
  !endif

  ! compute the residual vector for the vegetation canopy
  ! NOTE: sMul(ixVegWat) = 1, but include as it converts all variables to quadruple precision
  if(computeVegFlux)then
   ! --> energy balance
   rVec(ixCasNrg) = sMul(ixCasNrg)*scalarCanairTempTrial - ( (sMul(ixCasNrg)*scalarCanairTemp + fVec(ixCasNrg)*dt) + rAdd(ixCasNrg) )
   rVec(ixVegNrg) = sMul(ixVegNrg)*scalarCanopyTempTrial - ( (sMul(ixVegNrg)*scalarCanopyTemp + fVec(ixVegNrg)*dt) + rAdd(ixVegNrg) )
   ! --> mass balance
   rVec(ixVegWat) = sMul(ixVegWat)*scalarCanopyWatTrial  - ( (sMul(ixVegWat)*scalarCanopyWat  + fVec(ixVegWat)*dt) + rAdd(ixVegWat) )
  endif
  !write(*,'(a,1x,2(e20.10,1x))') 'rVec(ixVegWat), fVec(ixVegWat) = ', rVec(ixVegWat), fVec(ixVegWat)

  ! compute the residual vector for the snow and soil sub-domains for energy
  !rAdd(ixSnowSoilNrg) = 0._dp
  !fVec(ixSnowSoilNrg) = 0._dp
  !print*, 'fVec(ixSnowSoilNrg) = ', fVec(ixSnowSoilNrg)
  rVec(ixSnowSoilNrg) = sMul(ixSnowSoilNrg)*mLayerTempTrial(1:nLayers) - ( (sMul(ixSnowSoilNrg)*mLayerTemp(1:nLayers)  + fVec(ixSnowSoilNrg)*dt) + rAdd(ixSnowSoilNrg) )
  !if(printFlag)then
  ! write(*,'(a,1x,10(e25.15,1x))') 'fVec(1:2)*dt = ', fVec(1:2)*dt
  ! write(*,'(a,1x,10(e20.10,1x))') 'rAdd(1:8)    = ', rAdd(1:8)
  !endif

  !if(printFlag)then
  ! do iLayer=nSnow+1,nLayers
  !  jLayer = ixSnowSoilNrg(iLayer)
  !  write(*,'(a,1x,2(i4,1x),10(e20.10,1x))') 'iLayer, jLayer, fVec(jLayer), sMul(jLayer), rAdd(jLayer), rVec(jLayer), mLayerVolFracIceLocal(iLayer) = ', &
  !                                            iLayer, jLayer, fVec(jLayer), sMul(jLayer), rAdd(jLayer), rVec(jLayer), mLayerVolFracIceLocal(iLayer)
  ! end do
  !endif

  ! compute the residual vector for the **snow** sub-domain for liquid water
  if(nSnow>0)&
  rVec(ixSnowOnlyWat) = mLayerVolFracWatTrial(1:nSnow) - ( (mLayerVolFracWat(1:nSnow)  + fVec(ixSnowOnlyWat)*dt) + rAdd(ixSnowOnlyWat) )
  !if(printFlag)then
  ! do iLayer=1,min(nSnow,5)
  !  jLayer = ixSnowOnlyWat(iLayer)
  !  write(*,'(a,1x,2(i4,1x),10(e20.10,1x))') 'iLayer, jLayer, fVec(jLayer), sMul(jLayer), rAdd(jLayer), rVec(jLayer), mLayerVolFracIceLocal(iLayer) = ', &
  !                                            iLayer, jLayer, fVec(jLayer), sMul(jLayer), rAdd(jLayer), rVec(jLayer), mLayerVolFracIceLocal(iLayer)
  ! end do
  !endif

  ! compute the residual vector for the **soil** sub-domain for liquid water
  !fVec(ixSoilOnlyMat) = 0._dp
  vThetaInit(1:nSoil)  = mLayerVolFracLiq(nSnow+1:nLayers)      + mLayerVolFracIce(nSnow+1:nLayers)      ! liquid equivalent of total water at the start of the step
  vThetaTrial(1:nSoil) = mLayerVolFracLiqLocal(nSnow+1:nLayers) + mLayerVolFracIceLocal(nSnow+1:nLayers) ! liquid equivalent of total water at the current iteration
  rVec(ixSoilOnlyMat)  = vThetaTrial(1:nSoil) - ( (vThetaInit(1:nSoil) + fVec(ixSoilOnlyMat)*dt) + rAdd(ixSoilOnlyMat) )

  !do iLayer=1,nSoil
  ! xCompress = mLayerCompress(iLayer)*mLayerDepth(iLayer)  ! m
  ! xFlux0    = iLayerLiqFluxSoil(iLayer-1)*dt              ! m
  ! xFlux1    = iLayerLiqFluxSoil(iLayer)*dt                ! m
  ! if(iLayer==1) write(*,'(a)')                     'iLayer, vThetaTrial(iLayer), vThetaInit(iLayer), (xFlux1 - xFlux0), xCompress, fVec(ixSoilOnlyMat(iLayer))*dt, rAdd(ixSoilOnlyMat(iLayer)), rVec(ixSoilOnlyMat(iLayer)) = '
  !               write(*,'(1x,i4,1x,10(e20.10,1x))') iLayer, vThetaTrial(iLayer), vThetaInit(iLayer), (xFlux1 - xFlux0), xCompress, fVec(ixSoilOnlyMat(iLayer))*dt, rAdd(ixSoilOnlyMat(iLayer)), rVec(ixSoilOnlyMat(iLayer))
  !end do
  !xBalance = sum((vThetaTrial(1:nSoil) - vThetaInit(1:nSoil))*mLayerDepth(nSnow+1:nSoil))
  !write(*,'(a,e20.10)') 'xBalance                                                       = ', xBalance
  !write(*,'(a,e20.10)') 'sum(rAdd(ixSoilOnlyMat)*mLayerDepth(nSnow+1:nSoil))            = ', sum(rAdd(ixSoilOnlyMat)*mLayerDepth(nSnow+1:nSoil))
  !write(*,'(a,e20.10)') 'sum(rAdd(ixSoilOnlyMat)*mLayerDepth(nSnow+1:nSoil)) - xBalance = ', sum(rAdd(ixSoilOnlyMat)*mLayerDepth(nSnow+1:nSoil)) - xBalance
  !write(*,'(a,e20.10)') 'sum(fVec(ixSoilOnlyMat)*mLayerDepth(nSnow+1:nSoil))*dt         = ', sum(fVec(ixSoilOnlyMat)*mLayerDepth(nSnow+1:nSoil))*dt
  !write(*,'(a,e20.10)') 'iLayerLiqFluxSoil(0)*dt                                        = ', iLayerLiqFluxSoil(0)*dt


  !if(printFlag)then
  ! write(*,'(a,1x,10(e20.10,1x))') 'vThetaInit(1:nSoil)  = ', vThetaInit(1:nSoil)
  ! write(*,'(a,1x,10(e20.10,1x))') 'vThetaTrial(1:nSoil) = ', vThetaTrial(1:nSoil)
  ! write(*,'(a,1x,10(e20.10,1x))') 'fVec(ixSnowSoilWat)  = ', fVec(ixSnowSoilWat)
  ! write(*,'(a,1x,10(e20.10,1x))') 'rAdd(ixSnowSoilWat)  = ', rAdd(ixSnowSoilWat)
  ! write(*,'(a,1x,10(e20.10,1x))') 'rVec(ixSnowSoilWat)  = ', rVec(ixSnowSoilWat)
  !endif

  ! test
  !write(*,'(a,1x,10(e15.5,1x))') 'rVec(1:10) = ', rVec(1:10)
  !write(*,'(a,1x,10(e15.5,1x))') 'rAdd(1:10) = ', rAdd(1:10)
  !write(*,'(a,1x,10(e15.5,1x))') 'fVec(1:10) = ', fVec(1:10)

  !print*, '***'
  !write(*,'(a,1x,10(e15.5,1x))') 'mLayerVolFracLiqTrial(1:10) = ', mLayerVolFracLiqTrial(1:10)

  ! end association to variables in the data structures
  end associate

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
  !real(dp)                       :: vTheta1,volIce1,effSat1,psiLiq1           ! test derivatives (-)
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

  ! layer type (snow or soil)
  layerType               => indx_data%var(iLookINDEX%layerType)%dat                ,&  ! intent(in): [i4b(:)] type of each layer in the snow+soil domain (snow or soil)

  ! layer depth
  mLayerDepth             => mvar_data%var(iLookMVAR%mLayerDepth)%dat               ,&  ! intent(in): [dp(:)] depth of each layer in the snow-soil sub-domain (m)

  ! snow parameters
  snowfrz_scale           => mpar_data%var(iLookPARAM%snowfrz_scale)                ,&  ! intent(in): [dp] scaling parameter for the snow freezing curve (K-1)

  ! soil parameters
  vGn_alpha               => mpar_data%var(iLookPARAM%vGn_alpha)                    ,&  ! intent(in): [dp] van Genutchen "alpha" parameter (m-1)
  vGn_n                   => mpar_data%var(iLookPARAM%vGn_n)                        ,&  ! intent(in): [dp] van Genutchen "n" parameter (-)
  vGn_m                   => mvar_data%var(iLookMVAR%scalarVGn_m)%dat(1)            ,&  ! intent(in): [dp] van Genutchen "m" parameter (-)
  theta_sat               => mpar_data%var(iLookPARAM%theta_sat)                    ,&  ! intent(in): [dp] soil porosity (-)
  theta_res               => mpar_data%var(iLookPARAM%theta_res)                    ,&  ! intent(in): [dp] soil residual volumetric water content (-)

  ! model diagnostic variables
  scalarThroughfallRain   => mvar_data%var(iLookMVAR%scalarThroughfallRain)%dat(1)  ,&  ! intent(out): [dp] rain that reaches the ground without ever touching the canopy (kg m-2 s-1)
  scalarCanopyLiqDrainage => mvar_data%var(iLookMVAR%scalarCanopyLiqDrainage)%dat(1),&  ! intent(out): [dp] drainage of liquid water from the vegetation canopy (kg m-2 s-1)
  scalarSurfaceRunoff     => mvar_data%var(iLookMVAR%scalarSurfaceRunoff)%dat(1)    ,&  ! intent(out): [dp] surface runoff (m s-1)
  scalarRainPlusMelt      => mvar_data%var(iLookMVAR%scalarRainPlusMelt)%dat(1)     ,&  ! intent(out): [dp] rain plus melt (m s-1)
  scalarExfiltration      => mvar_data%var(iLookMVAR%scalarExfiltration)%dat(1)     ,&  ! intent(out): [dp] exfiltration from the soil profile (m s-1)
  mLayerColumnOutflow     => mvar_data%var(iLookMVAR%mLayerColumnOutflow)%dat       ,&  ! intent(out): [dp(:)] column outflow from each soil layer (m3 s-1)

  ! soil fluxes
  iLayerLiqFluxSnow       => mvar_data%var(iLookMVAR%iLayerLiqFluxSnow)%dat         ,&  ! intent(out): [dp(0:)] vertical liquid water flux at snow layer interfaces (-)
  iLayerLiqFluxSoil       => mvar_data%var(iLookMVAR%iLayerLiqFluxSoil)%dat         ,&  ! intent(out): [dp(0:)] vertical liquid water flux at soil layer interfaces (-)
  mLayerBaseflow          => mvar_data%var(iLookMVAR%mLayerBaseflow)%dat            ,&  ! intent(out): [dp(:)] baseflow from each soil layer (m s-1)

  ! aquifer fluxes
  scalarAquiferTranspire  => mvar_data%var(iLookMVAR%scalarAquiferTranspire)%dat(1) ,&  ! intent(out): [dp] transpiration loss from the aquifer (m s-1
  scalarAquiferRecharge   => mvar_data%var(iLookMVAR%scalarAquiferRecharge)%dat(1)  ,&  ! intent(out): [dp] recharge to the aquifer (m s-1)
  scalarAquiferBaseflow   => mvar_data%var(iLookMVAR%scalarAquiferBaseflow)%dat(1)   &  ! intent(out): [dp] total baseflow from the aquifer (m s-1)

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
     ! compute numerical derivative (testing)
     !if(printFlag .and. iLayer==1)then
     ! fLiq0 = fracliquid(mLayerTempTrial(iLayer),snowfrz_scale)
     ! fLiq1 = fracliquid(mLayerTempTrial(iLayer)+dx,snowfrz_scale)
     ! print*, 'testderivative', (fLiq1*theta - fLiq0*theta)/dx, mLayerdTheta_dTk(iLayer)
     !endif
    case(ix_soil) ! (soil layers)
     if(mLayerVolFracIceTrial(iLayer)>verySmall)then
      mLayerdTheta_dTk(iLayer)        = dTheta_dTk(mLayerTempTrial(iLayer),theta_res,theta_sat,vGn_alpha,vGn_n,vGn_m)  ! assume no volume expansion
     else
      mLayerdTheta_dTk(iLayer)        = 0._dp
     endif
    case default; err=40; message=trim(message)//"cannot identify the layer as snow or soil"; return
   endselect
   ! (check)
   !if(printFlag .and. iLayer==58)&
   !write(*,'(a,1x,i4,1x,2(f20.10,1x))') 'mLayerdTheta_dTk(iLayer) = ', iLayer, mLayerVolFracIceTrial(iLayer), mLayerdTheta_dTk(iLayer)
  end do  ! (looping through snow+soil layers)

  ! * compute the matric head associated with liquid water
  do iSoil=1,nSoil  ! loop through soil layers

   ! - compute derivative in total water content w.r.t. total water matric potential (m-1)
   dVolTot_dPsi0(iSoil) = dTheta_dPsi(mLayerMatricHeadTrial(iSoil),vGn_alpha,theta_res,theta_sat,vGn_n,vGn_m)  ! valid for both frozen and unfrozen conditions
   !if(printFlag.and.iSoil<5) print*, 'iSoil, dVolTot_dPsi0 = ', iSoil, dVolTot_dPsi0(iSoil)
   !if(printFlag.and.iSoil<5) print*, 'mLayerVolFracIceTrial(nSnow+iSoil), mLayerVolFracLiqTrial(nSnow+iSoil) = ', mLayerVolFracIceTrial(nSnow+iSoil), mLayerVolFracLiqTrial(nSnow+iSoil)

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
    !if(printFlag) print*, 'mLayerMatricHeadLiq(iSoil) = ', mLayerMatricHeadLiq(iSoil)
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

  ! print volumetric ice content
  !if(printFlag)then
  ! write(*,'(a,1x,100(f20.12,1x))') 'mLayerMatricHeadTrial(1:nSoil)         = ', mLayerMatricHeadTrial(1:nSoil)
  ! write(*,'(a,1x,100(f20.12,1x))') 'mLayerVolFracLiqTrial(nSnow+1:nLayers) = ', mLayerVolFracLiqTrial(nSnow+1:nLayers)
  ! write(*,'(a,1x,100(f20.12,1x))') 'mLayerVolFracIceTrial(nSnow+1:nLayers) = ', mLayerVolFracIceTrial(nSnow+1:nLayers)
  !endif

  ! *****
  ! (1) CALCULATE ENERGY FLUXES OVER VEGETATION...
  ! **********************************************

  call vegNrgFlux(&
                  ! input: model control
                  iter,                                   & ! intent(in): iteration index
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
  if(printFlag)then
   write(*,'(a,1x,f30.20)') 'canairNetNrgFlux = ', canairNetNrgFlux
   write(*,'(a,1x,f30.20)') 'canopyNetNrgFlux = ', canopyNetNrgFlux
   write(*,'(a,1x,f30.20)') 'groundNetNrgFlux = ', groundNetNrgFlux
   write(*,'(a,1x,f30.20)') 'dGroundNetFlux_dGroundTemp = ', dGroundNetFlux_dGroundTemp
  endif

  !if(printFlag)then
   !print*, 'in systemSolv: scalarGroundEvaporation = ', scalarGroundEvaporation
   !print*, 'in systemSolv: scalarCanopyEvaporation = ', scalarCanopyEvaporation
   !print*, 'in systemSolv: dCanopyEvaporation_dCanLiq = ', dCanopyEvaporation_dCanLiq
  !endif

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

  !iLayerNrgFlux(0) = 0._dp

  do iLayer=1,nLayers
   ssdNetNrgFlux(iLayer) = -(iLayerNrgFlux(iLayer) - iLayerNrgFlux(iLayer-1))/mLayerDepth(iLayer)
   if(printFlag)then
    if(iLayer < 3) write(*,'(a,1x,i4,1x,10(f25.15,1x))') 'iLayer, iLayerNrgFlux(iLayer-1:iLayer), ssdNetNrgFlux(iLayer)   = ', iLayer, iLayerNrgFlux(iLayer-1:iLayer), ssdNetNrgFlux(iLayer)
   endif
  end do
  !print*, 'iLayerNrgFlux = ', iLayerNrgFlux
  !print*, 'ssdNetNrgFlux = ', ssdNetNrgFlux


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
  !if(printFlag)then
  ! write(*,'(a,1x,10(e20.10,1x))') 'fracLiqSnow(1), mLayerVolFracLiqTrial(1), iLayerLiqFluxSnow(1), iLayerLiqFluxSnowDeriv(1) = ', &
  !                                  fracLiqSnow(1), mLayerVolFracLiqTrial(1), iLayerLiqFluxSnow(1), iLayerLiqFluxSnowDeriv(1)
  !endif
  !if(printFlag)then
  ! print*, 'scalarRainPlusMelt = ', scalarRainPlusMelt
  !endif

  ! *****
  ! (5) CALCULATE THE LIQUID FLUX THROUGH SOIL...
  ! *********************************************
  call soilLiqFlx(&
                  ! input: model control
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
   !if(printFlag .and. iLayer==1)&
   !write(*,'(a,1x,i4,1x,10(f25.15,1x))') 'iLayer, soilNetLiqFlux(iLayer), scalarSurfaceInfiltration, iLayerLiqFluxSoil(iLayer-1), iLayerLiqFluxSoil(iLayer) = ', &
   !                                       iLayer, soilNetLiqFlux(iLayer), scalarSurfaceInfiltration, iLayerLiqFluxSoil(iLayer-1), iLayerLiqFluxSoil(iLayer)
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

  !do iLayer=1,nSoil
   !vTheta1 = volFracLiq(mLayerMatricHeadTrial(iLayer)+dx,vGn_alpha,theta_res,theta_sat,vGn_n,vGn_m)
   !volIce1 = vTheta1 - mLayerVolFracLiqTrial(nSnow+iLayer)
   !effSat1 = (mLayerVolFracLiqTrial(nSnow+iLayer) - theta_res) / (theta_sat - volIce1 - theta_res)
   !psiLiq1 = matricHead(effSat1,vGn_alpha,0._dp,1._dp,vGn_n,vGn_m)  ! use effective saturation, so theta_res=0 and theta_sat=1
   !print*, 'numerical derivative = ', (psiLiq1 -  mLayerMatricHeadLiq(iLayer))/dx
   !print*, 'analytical derivative = ', dPsiLiq_dPsi0(iLayer)
  !end do
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
                   dt,                                      & ! intent(in):    length of the model time step (s)
                   firstFluxCall,                           & ! intent(in):    logical flag to compute index of the lowest saturated layer
                   ! input: state and diagnostic variables
                   mLayerdTheta_dPsi,                       & ! intent(in):    derivative in the soil water characteristic w.r.t. matric head in each layer (m-1)
                   mLayerMatricHeadLiq,                     & ! intent(in):    liquid water matric potential (m)
                   mLayerVolFracLiqTrial(nSnow+1:nLayers),  & ! intent(in):    volumetric fraction of liquid water (-)
                   mLayerVolFracIceTrial(nSnow+1:nLayers),  & ! intent(in):    volumetric fraction of ice (-)
                   ! input: data structures
                   attr_data,                               & ! intent(in):    model attributes
                   mpar_data,                               & ! intent(in):    model parameters
                   mvar_data,                               & ! intent(inout): model variables for a local HRU
                   ! output
                   ixSaturation,                            & ! intent(inout) index of lowest saturated layer (NOTE: only computed on the first iteration)
                   mLayerBaseflow,                          & ! intent(out): baseflow from each soil layer (m s-1)
                   dBaseflow_dMatric,                       & ! intent(out): derivative in baseflow w.r.t. matric head (s-1)
                   err,cmessage)                              ! intent(out): error control
   if(err/=0)then; message=trim(message)//trim(cmessage); return; endif
   !write(*,'(a,1x,10(e20.10,1x))') 'iter, mLayerBaseflow(:) = ', mLayerBaseflow(:)
   !pause 'computing baseflow fluxes'

   ! check
   if(printFlag)then
    ! check baseflow
    write(*,'(a,1x,10(e30.20,1x))') 'baseflow: ', mLayerBaseflow(:)
    ! check baseflow derivatives
    !do iLayer=1,nSoil
    ! write(*,'(a,1x,i4,1x,100(e20.10,1x))') 'dBaseflow: ', iLayer, dBaseflow_dMatric(:,iLayer)
    !end do
   endif

  endif

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
  fluxVec(ixSoilOnlyMat) = soilNetLiqFlux(1:nSoil)
  if(nSnow>0)&
  fluxVec(ixSnowOnlyWat) = snowNetLiqFlux(1:nSnow)

  ! print progress
  !print*, '**'
  !if(printFlag)then
  ! write(*,'(a,1x,100(f15.9,1x))') 'stateVec(:) = ', stateVec(:)
  ! write(*,'(a,1x,100(e15.9,1x))') 'fluxVec(:)  = ', fluxVec(:)
  !endif

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
  associate(mLayerDepth => mvar_data%var(iLookMVAR%mLayerDepth)%dat) ! intent(in): [dp(:)] depth of each layer in the snow-soil sub-domain (m)

  ! initialize the Jacobian
  ! NOTE: this needs to be done every time, since Jacobian matrix is modified in the solver
  aJac(:,:) = 0._dp  ! analytical Jacobian matrix

  ! -----
  ! * energy and liquid fluxes over vegetation...
  ! ---------------------------------------------
  if(computeVegFlux)then  ! (derivatives only defined when vegetation protrudes over the surface)

   ! liquid water fluxes for vegetation canopy (-)
   aJac(ixDiag,ixVegWat) = -fracLiqVeg*(dCanopyEvaporation_dCanLiq - scalarCanopyLiqDrainageDeriv)*dt + 1._dp     ! ixVegWat: CORRECT

   ! cross-derivative terms w.r.t. system temperatures (kg m-2 K-1)
   aJac(ixSub2,ixCasNrg) = -dCanopyEvaporation_dTCanair*dt                                                        ! ixCasNrg: CORRECT
   aJac(ixSub1,ixVegNrg) = -dCanopyEvaporation_dTCanopy*dt + dt*scalarCanopyLiqDrainageDeriv*dCanLiq_dTcanopy     ! ixVegNrg: CORRECT
   aJac(ixSup1,ixTopNrg) = -dCanopyEvaporation_dTGround*dt                                                        ! ixTopNrg: CORRECT

   ! cross-derivative terms w.r.t. canopy water (kg-1 m2)
   aJac(ixSub2,ixVegWat) = (dt/mLayerDepth(1))*(-soilControl*fracLiqVeg*scalarCanopyLiqDrainageDeriv)/iden_water  ! ixVegWat: CORRECT

   ! cross-derivative terms w.r.t. canopy temperature (K-1)
   aJac(ixSub3,ixVegNrg) = (dt/mLayerDepth(1))*(-soilControl*scalarCanopyLiqDrainageDeriv*dCanLiq_dTcanopy)/iden_water    ! ixVegNrg: CORRECT
   !print*, 'soilControl, scalarCanopyLiqDrainageDeriv, dCanLiq_dTcanopy = ', soilControl, scalarCanopyLiqDrainageDeriv, dCanLiq_dTcanopy

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

   ! test
   !print*, 'aJac(ixSub2,ixVegWat) = ', aJac(ixSub2,ixVegWat)
   !print*, 'aJac(ixSub3,ixVegNrg) = ', aJac(ixSub3,ixVegNrg)

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
   aJac(ixSub1,mLayer) = (dt/mLayerDepth(iLayer))*iLayerLiqFluxSnowDeriv(iLayer)*mLayerdTheta_dTk(iLayer)  ! (dVol/dT) -- increase in volumetric liquid water content balanced by a decrease in volumetric ice content
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
   jLayer = ixSoilOnlyMat(iLayer)  ! layer index within the full state vector
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
   jLayer = ixSoilOnlyMat(iLayer)       ! hydrology layer index within the full state vector
   mLayer = ixSnowSoilNrg(kLayer)       ! thermodynamics layer index within the full state vector
   !write(*,'(a,1x,10(i4,1x))') 'iLayer, jLayer, jLayer-nVarSnowSoil, jLayer+nVarSnowSoil, kLayer, mLayer = ', &
   !                             iLayer, jLayer, jLayer-nVarSnowSoil, jLayer+nVarSnowSoil, kLayer, mLayer
   ! - compute the Jacobian for the layer itself
   aJac(ixSub1,mLayer) = (dt/mLayerDepth(kLayer))*(-dq_dNrgStateBelow(iLayer-1) + dq_dNrgStateAbove(iLayer))   ! dVol/dT (K-1) -- flux depends on ice impedance
   if(mLayerVolFracIceTrial(kLayer) > tiny(dt))then
    aJac(ixSup1,jLayer) = -dVolTot_dPsi0(iLayer)*LH_fus*iden_water    ! dNrg/dMat (J m-3 m-1) -- dMat changes volumetric water, and hence ice content
   else
    aJac(ixSup1,jLayer) = 0._dp
   endif

   ! - compute the Jacobian for neighboring layers (dVol/dT)
   if(kLayer > nSnow+1) aJac(ixSup1,mLayer) = (dt/mLayerDepth(kLayer-1))*( dq_dNrgStateBelow(iLayer-1))   ! K-1
   if(kLayer < nLayers) aJac(ixSub3,mLayer) = (dt/mLayerDepth(kLayer+1))*(-dq_dNrgStateAbove(iLayer))     ! K-1

  end do  ! (looping through soil layers)

  ! -----
  ! * testing.....
  ! --------------
  !print*, '** analytical Jacobian:'
  !write(*,'(a4,1x,100(i11,1x))') 'xCol', (iLayer, iLayer=iJac1,iJac2)
  !do iLayer=kl+1,nBands; write(*,'(i4,1x,100(e11.5,1x))') iLayer, aJac(iLayer,iJac1:iJac2); end do
  !pause


  ! end association to variables in the data structures
  end associate

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
  associate(mLayerDepth => mvar_data%var(iLookMVAR%mLayerDepth)%dat) ! intent(in): [dp(:)] depth of each layer in the snow-soil sub-domain (m)

  ! initialize the Jacobian
  ! NOTE: this needs to be done every time, since Jacobian matrix is modified in the solver
  aJac(:,:) = 0._dp  ! analytical Jacobian matrix

  ! -----
  ! * energy and liquid fluxes over vegetation...
  ! ---------------------------------------------
  if(computeVegFlux)then  ! (derivatives only defined when vegetation protrudes over the surface)

   ! liquid water fluxes for vegetation canopy (-)
   aJac(ixVegWat,ixVegWat) = -fracLiqVeg*(dCanopyEvaporation_dCanLiq - scalarCanopyLiqDrainageDeriv)*dt + 1._dp

   ! cross-derivative terms w.r.t. system temperatures (kg m-2 K-1)
   aJac(ixVegWat,ixCasNrg) = -dCanopyEvaporation_dTCanair*dt
   aJac(ixVegWat,ixVegNrg) = -dCanopyEvaporation_dTCanopy*dt + dt*scalarCanopyLiqDrainageDeriv*dCanLiq_dTcanopy
   aJac(ixVegWat,ixTopNrg) = -dCanopyEvaporation_dTGround*dt

   ! cross-derivative terms w.r.t. canopy water (kg-1 m2)
   aJac(ixTopLiq,ixVegWat) = (dt/mLayerDepth(1))*(-soilControl*fracLiqVeg*scalarCanopyLiqDrainageDeriv)/iden_water

   ! cross-derivative terms w.r.t. canopy temperature (K-1)
   aJac(ixTopLiq,ixVegNrg) = (dt/mLayerDepth(1))*(-soilControl*scalarCanopyLiqDrainageDeriv*dCanLiq_dTcanopy)/iden_water
   !print*, 'soilControl, scalarCanopyLiqDrainageDeriv, dCanLiq_dTcanopy = ', soilControl, scalarCanopyLiqDrainageDeriv, dCanLiq_dTcanopy

   ! cross-derivative terms w.r.t. canopy liquid water (J m-1 kg-1)
   ! NOTE: dIce/dLiq = (1 - fracLiqVeg); dIce*LH_fus/canopyDepth = J m-3; dLiq = kg m-2
   aJac(ixVegNrg,ixVegWat) = (dt/canopyDepth)   *(-dCanopyNetFlux_dCanLiq) - (1._dp - fracLiqVeg)*LH_fus/canopyDepth   ! dF/dLiq
   aJac(ixTopNrg,ixVegWat) = (dt/mLayerDepth(1))*(-dGroundNetFlux_dCanLiq)
   !print*, '(dt/canopyDepth)   *(-dCanopyNetFlux_dCanLiq) = ', (dt/canopyDepth)   *(-dCanopyNetFlux_dCanLiq)
   !print*, '(1._dp - fracLiqVeg)*LH_fus/canopyDepth = ', (1._dp - fracLiqVeg)*LH_fus/canopyDepth

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

   !print*, 'aJac(ixVegWat,ixVegNrg) = ', aJac(ixVegWat,ixVegNrg)
   !print*, 'aJac(ixVegNrg,ixVegWat) = ', aJac(ixVegNrg,ixVegWat)
   !print*, 'aJac(ixTopNrg,ixVegWat) = ', aJac(ixTopNrg,ixVegWat)
   !print*, 'aJac(ixTopNrg,ixVegNrg) = ', aJac(ixTopNrg,ixVegNrg)
   !print*, 'aJac(ixTopLiq,ixVegNrg) = ', aJac(ixTopLiq,ixVegNrg)
   !print*, 'aJac(ixTopLiq,ixVegWat) = ', aJac(ixTopLiq,ixVegWat)

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
   aJac(mLayer,jLayer) = -(1._dp - fracLiqSnow(iLayer))*LH_fus*iden_water     ! (dF/dLiq)
   aJac(jLayer,mLayer) = (dt/mLayerDepth(iLayer))*iLayerLiqFluxSnowDeriv(iLayer)*mLayerdTheta_dTk(iLayer)  ! (dVol/dT) -- increase in volumetric liquid water content balanced by a decrease in volumetric ice content
   ! - compute cross-derivative terms for the layer below (w.r.t. state in the current layer)
   if(iLayer < nSnow)then
    aJac(jLayer+nVarSnowSoil,mLayer) = -(dt/mLayerDepth(iLayer+1))*iLayerLiqFluxSnowDeriv(iLayer)*mLayerdTheta_dTk(iLayer)        ! dVol(below)/dT(above) -- K-1
    aJac(jLayer+nVarSnowSoil,jLayer) = -(dt/mLayerDepth(iLayer+1))*iLayerLiqFluxSnowDeriv(iLayer)*fracLiqSnow(iLayer)             ! dVol(below)/dLiq(above) -- (-)
    !print*, 'aJac(jLayer+nVarSnowSoil,jLayer) = ', aJac(jLayer+nVarSnowSoil,jLayer)
   endif
  end do

  ! -----
  ! * liquid water fluxes for the soil domain...
  ! --------------------------------------------
  do iLayer=1,nSoil    ! loop through layers in the soil domain

   ! - define layer indices
   jLayer = ixSoilOnlyMat(iLayer)  ! layer index within the full state vector
   kLayer = iLayer+nSnow           ! layer index within the full snow-soil vector

   ! - compute the Jacobian
   ! all terms *excluding* baseflow
   aJac(jLayer,jLayer) = (dt/mLayerDepth(kLayer))*(-dq_dHydStateBelow(iLayer-1) + dq_dHydStateAbove(iLayer)) + dMat(jLayer)
   if(kLayer > nSnow+1) aJac(jLayer-nVarSnowSoil,jLayer) = (dt/mLayerDepth(kLayer-1))*( dq_dHydStateBelow(iLayer-1))
   if(kLayer < nLayers) aJac(jLayer+nVarSnowSoil,jLayer) = (dt/mLayerDepth(kLayer+1))*(-dq_dHydStateAbove(iLayer))

   ! include terms for baseflow
   do pLayer=1,nSoil
    qLayer = ixSoilOnlyMat(pLayer)  ! layer index within the full state vector
    aJac(jLayer,qLayer) = aJac(jLayer,qLayer) + (dt/mLayerDepth(kLayer))*dBaseflow_dMatric(iLayer,pLayer)
   end do

  end do  ! (looping through soil layers)

  ! -----
  ! * derivative in liquid water fluxes w.r.t. temperature for the soil domain...
  ! -----------------------------------------------------------------------------
  do iLayer=1,nSoil    ! loop through layers in the soil domain
   ! - define layer indices
   kLayer = iLayer+nSnow                ! layer index within the full snow-soil vector
   jLayer = ixSoilOnlyMat(iLayer)       ! hydrology layer index within the full state vector
   mLayer = ixSnowSoilNrg(kLayer)       ! thermodynamics layer index within the full state vector
   ! - compute the Jacobian for the layer itself
   aJac(jLayer,mLayer) = (dt/mLayerDepth(kLayer))*(-dq_dNrgStateBelow(iLayer-1) + dq_dNrgStateAbove(iLayer))   ! dVol/dT (K-1) -- flux depends on ice impedance
   if(mLayerVolFracIceTrial(iLayer+nSnow) > tiny(dt))then
    aJac(mLayer,jLayer) = -dVolTot_dPsi0(iLayer)*LH_fus*iden_water    ! dNrg/dMat (J m-3 m-1) -- dMat changes volumetric water, and hence ice content
   else
    aJac(mLayer,jLayer) = 0._dp
   endif
   !if(iLayer==1) write(*,'(a)')        'iLayer, jLayer, jLayer-nVarSnowSoil, jLayer+nVarSnowSoil, kLayer, mLayer, aJac(mLayer,jLayer), aJac(jLayer,mLayer), dq_dNrgStateBelow(iLayer-1), dq_dNrgStateAbove(iLayer) = '
   !write(*,'(6(i4,1x),10(e20.10,1x))')  iLayer, jLayer, jLayer-nVarSnowSoil, jLayer+nVarSnowSoil, kLayer, mLayer, aJac(mLayer,jLayer), aJac(jLayer,mLayer), dq_dNrgStateBelow(iLayer-1), dq_dNrgStateAbove(iLayer)

   ! - compute the Jacobian for neighboring layers
   if(kLayer > nSnow+1) aJac(jLayer-nVarSnowSoil,mLayer) = (dt/mLayerDepth(kLayer-1))*( dq_dNrgStateBelow(iLayer-1))   ! K-1
   if(kLayer < nLayers) aJac(jLayer+nVarSnowSoil,mLayer) = (dt/mLayerDepth(kLayer+1))*(-dq_dNrgStateAbove(iLayer))     ! K-1

  end do  ! (looping through soil layers)

  ! print the Jacobian
  if(globalPrintFlag)then
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
  real(qp),intent(in)            :: resVec(:)               ! model residual vector (mixed units)
  integer(i4b),intent(out)       :: err                     ! error code
  character(*),intent(out)       :: message                 ! error message
  ! local
  character(len=256)             :: cmessage                ! error message of downwind routine
  real(dp),dimension(nState)     :: stateVecPerturbed       ! perturbed state vector
  real(dp),dimension(nState)     :: fluxVecJac              ! flux vector
  real(qp),dimension(nState)     :: resVecJac               ! residual vector (mixed units)
  integer(i4b)                   :: iJac                    ! index of row of the Jacobian matrix
  integer(i4b),parameter         :: iTry=-999               ! index of trial model state variable (used for testing)
  integer(i4b)                   :: ixDesire                ! index of a desired layer (used for testing)
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
  !integer(i4b)                   :: ixNumType=ixNumFlux      ! method used to calculate the numerical Jacobian
  ! --------------------------------------------------------------
  ! initialize error control
  err=0; message='numlJacob/'

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

   ! print progress
   !print*, '*** iJac = ', iJac

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


    if(iJac==iTry)then

     !write(*,'(a)')             'fluxVecJac(iTry), fluxVec(iTry), (fluxVecJac(iTry) - fluxVec(iTry))/dx, (LH_fus*iden_water*(mLayerVolFracIceLocal(15)-mLayerVolFracIceTrial(15)))/dx = '
     !write(*,'(100(e25.15,1x))') fluxVecJac(iTry), fluxVec(iTry), (fluxVecJac(iTry) - fluxVec(iTry))/dx, (LH_fus*iden_water*(mLayerVolFracIceLocal(15)-mLayerVolFracIceTrial(15)))/dx
     !write(*,'(a,1x,100(e25.15,1x))') 'fluxVecJac(4), fluxVec(4), (fluxVecJac(4) - fluxVec(4))/dx  = ', &
     !                                  fluxVecJac(4), fluxVec(4), (fluxVecJac(4) - fluxVec(4))/dx
     write(*,'(a,1x,i4,1x,100(e25.15,1x))') 'test: iJac; resVec(iJac1:iJac2):    ', iJac, resVec(iJac1:iJac2)
     write(*,'(a,1x,i4,1x,100(e25.15,1x))') 'test: iJac; resVecJac(iJac1:iJac2): ', iJac, resVecJac(iJac1:iJac2)
     write(*,'(a,1x,i4,1x,100(e25.15,1x))') 'test: iJac; resVecJac(iJac1:iJac2) - resVec(iJac1:iJac2): ', iJac, resVecJac(iJac1:iJac2) - resVec(iJac1:iJac2)
    endif


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
    mLayerMatricHeadLocal(1:nSoil) = stateVecPerturbed(ixSoilOnlyMat)

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

    !if(iJac==iTry)then
    ! write(*,'(a,1x,100(e25.15,1x))') 'mLayerMatricHeadLocal = ', mLayerMatricHeadLocal
    ! write(*,'(a,1x,i4,1x,100(e25.15,1x))') 'test: iJac; fluxVec:    ', iJac, fluxVec(iJac1:iJac2)
    ! write(*,'(a,1x,i4,1x,100(e25.15,1x))') 'test: iJac; fluxVecJac: ', iJac, fluxVecJac(iJac1:iJac2)
    !endif

    ! (compute the row of the Jacobian matrix)
    nJac(:,iJac) = -dt*(fluxVecJac(:) - fluxVec(:))/dx
    !if(iJac==iTry)then
    ! write(*,'(a,1x,i4,1x,100(e15.5,1x))') 'test: iJac; nJac: ', iJac, nJac(iJac1:iJac2,iJac)
    !endif

    ! (add in the diagonal matrix)
    nJac(iJac,iJac) = nJac(iJac,iJac) + dMat(iJac)

   endif

   ! (print progress)
   if(iJac==iTry)then
    write(*,'(a,1x,3(f20.10,1x))') 'stateVec(iJac), stateVecPerturbed(iJac), dx = ', stateVec(iJac), stateVecPerturbed(iJac), dx
    write(*,'(a,1x,2(e20.10,1x))') 'fluxVec(iJac), fluxVecJac(iJac) = ', fluxVec(iJac), fluxVecJac(iJac)
    write(*,'(a,1x,2(e20.10,1x))') 'fluxVec(iJac+1), fluxVecJac(iJac+1) = ', fluxVec(iJac+1), fluxVecJac(iJac+1)
    write(*,'(a,1x,i4,1x,100(e15.5,1x))') 'iJac; nJac: ', iJac, nJac(iJac1:iJac2,iJac)
   endif

   ! (test)
   !if(iJac<10) write(*,'(a,1x,10(e15.5,1x))') 'fluxVecJac(1:10) = ', fluxVecJac(1:10)
   !if(iJac==iTry) write(*,'(a,1x,i4,1x,10(f20.14,1x))') 'iTry, stateVec(iTry), stateVecPerturbed(iTry) = ', iTry, stateVec(iTry), stateVecPerturbed(iTry)
   !if(iJac==iTry) write(*,'(a,1x,i4,1x,10(f20.8,1x))'), 'iTry, -dt*(fluxVecJac(iTry) - fluxVec(iTry))/dx = ', iTry, -dt*(fluxVecJac(iTry) - fluxVec(iTry))/dx
   !if(iJac==iTry) pause ' in numerical Jacobian calculations'

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

  end subroutine numlJacob


  ! *********************************************************************************************************
  ! internal subroutine lapackSolv: use the lapack routines to solve the linear system A.X=B
  ! *********************************************************************************************************
  subroutine lapackSolv(aJac,rVec,grad,xInc,err,message)
  implicit none
  ! dummy
  real(dp),intent(inout)         :: aJac(:,:)     ! input = the Jacobian matrix A; output = decomposed matrix
  real(qp),intent(in)            :: rVec(:)       ! the residual vector B
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

  ! select the option used to solve the linear system A.X=B
  select case(ixSolve)

   ! * full Jacobian matrix
   case(ixFullMatrix)
    do iJac=1,nState
     aJac(iJac,1:nState) = aJac(iJac,1:nState)/fscale(iJac)
    end do
    !print*, '** analytical Jacobian:'
    !write(*,'(a4,1x,100(i12,1x))') 'xCol', (iLayer, iLayer=iJac1,iJac2)
    !do iLayer=iJac1,iJac2; write(*,'(i4,1x,100(e12.5,1x))') iLayer, aJac(iJac1:iJac2,iLayer); end do

    ! ** test band diagonal matrix
    if(testBandDiagonal)then

     aJac_test(:,:)=0._dp
     ! form band-diagonal matrix
     do iState=1,nState
      do jState=max(1,iState-ku),min(nState,iState+kl)
       aJac_test(kl + ku + 1 + jState - iState, iState) = aJac(jState,iState)
       if(iState<6 .or. jState<6) write(*,'(2(i4,1x),e11.5)') jState,iState,aJac(jState,iState)
      end do
     end do
     print*, '** test banded analytical Jacobian:'
     write(*,'(a4,1x,100(i11,1x))') 'xCol', (iLayer, iLayer=iJac1,iJac2)
     do iLayer=kl+1,nBands; write(*,'(i4,1x,100(e11.5,1x))') iLayer, aJac_test(iLayer,iJac1:iJac2); end do
     pause

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
    !print*, '** analytical Jacobian:'
    !write(*,'(a4,1x,100(i11,1x))') 'xCol', (iLayer, iLayer=iJac1,iJac2)
    !do iLayer=kl+1,nBands; write(*,'(i4,1x,100(e11.5,1x))') iLayer, aJac(iLayer,iJac1:iJac2); end do
    !pause

  end select  ! (option to solve the linear system A.X=B)


  ! form the rhs matrix
  ! NOTE: scale the residual vector
  rhs(1:nState,1) = -rVec(1:nState)/fScale(1:nState)

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
     !if(iJac>99)&
     !write(*,'(a,1x,3(i4,1x),e15.5)') 'iJac,iState,kState = ', iJac,iState,kState,aJac(iState,iJac)

     ! calculate gradient (long-hand matrix multiplication)
     grad(iJac) = grad(iJac) - aJac(iState,iJac)*rhs(kState,1)

    end do  ! looping through elements of the band-diagonal matric
   end do  ! looping through state variables

  end select  ! (option to solve the linear system A.X=B)

  !print*, '(ixSolve == ixFullMatrix) =  ', (ixSolve == ixFullMatrix)
  !do iLayer=100,nState
  ! write(*,'(i4,1x,f20.10)') iLayer, grad(iLayer)
  !end do
  !pause

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

  !print*, 'xInc: (ixSolve == ixFullMatrix) =  ', (ixSolve == ixFullMatrix)
  !do iLayer=1,nState
  ! write(*,'(i4,1x,f20.10)') iLayer, xInc(iLayer)
  !end do
  !pause

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
  INTEGER(I4B) :: ndum,iterLS,iMax(1)
  integer(i4b),parameter :: maxiterLS=5
  REAL(DP) :: a,alam,alam2,alamin,b,disc,f2,fold2,pabs,rhs1,rhs2,slope,&
      tmplam
  ! NOTE: these variables are only used for testing
  !real(dp),dimension(size(xOld)) :: rVecOld
  !integer(i4b) :: iCheck
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
   f=0.5_dp*norm2(rVec/fScale)  ! NOTE: norm2 = sqrt(sum((rVec/fScale)**2._dp))

   ! check
   if(globalPrintFlag)then
    print*, '***'
    write(*,'(a,1x,100(e14.5,1x))')  trim(message)//': alam, fOld, f       = ', alam, fOld, f
    write(*,'(a,1x,100(f20.8,1x))')  trim(message)//': x(iJac1:iJac2)      = ', x(iJac1:iJac2)
    write(*,'(a,1x,100(f20.12,1x))') trim(message)//': p(iJac1:iJac2)      = ', p(iJac1:iJac2)
    write(*,'(a,1x,100(e20.5,1x))')  trim(message)//': rVec(iJac1:iJac2)   = ', rVec(iJac1:iJac2)
   endif

   ! check
   !if(iterLS>1)then   !.and. printFlag)then
   ! do iCheck=1,size(xOld)
   !  write(*,'(i4,1x,10(e20.10,1x))') iCheck, rVec(iCheck), rVecOld(iCheck), fScale(iCheck)
   ! end do
   !endif
   !rVecOld = rVec

   ! return if not doing the line search
   if(.not.doLineSearch)then
    converged = checkConv(rVec,p,x)
    return
   endif

   ! check convergence
   ! NOTE: this must be after the first flux call
   converged = checkConv(rVec,p,x)
   if(converged) return

   ! check if backtracked all the way to the original value
   if (iterLS==maxIterLS) then   !if (alam < alamin) then
    x(:)=xold(:)
    !print*, '*****************************************************************************'
    !print*, '*****************************************************************************'
    !print*, '*****************************************************************************'
    !print*, '*****************************************************************************'
    !print*, '** backtrack'
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
  function checkConv(rVec,xInc,xVec)
  implicit none
  ! dummies
  real(qp),intent(in)       :: rVec(:)         ! residual vector (mixed units)
  real(dp),intent(in)       :: xInc(:)         ! iteration increment (mixed units)
  real(dp),intent(in)       :: xVec(:)         ! state vector (mixed units)
  logical(lgt)              :: checkConv       ! flag to denote convergence
  ! locals
  real(dp),dimension(nSoil) :: psiScale        ! scaling factor for matric head
  real(dp),parameter        :: xSmall=1.e-0_dp ! a small offset
  logical(lgt)              :: liquidConv      ! flag for residual convergence
  logical(lgt)              :: matricConv      ! flag for matric head convergence
  logical(lgt)              :: energyConv      ! flag for energy convergence

  ! check convergence based on the residuals for energy (J m-3)
  if(computeVegFlux)then
   !canopy_max = abs(rVec(ixVegWat))
   energy_max = maxval(abs( (/rVec(ixCasNrg), rVec(ixVegNrg), rVec(ixSnowSoilNrg)/) ) )
   energy_loc = maxloc(abs( (/rVec(ixCasNrg), rVec(ixVegNrg), rVec(ixSnowSoilNrg)/) ) )
  else
   energy_max = maxval(abs( rVec(ixSnowSoilNrg) ) )
   energy_loc = maxloc(abs( rVec(ixSnowSoilNrg) ) )
  endif

  ! check convergence based on the residuals for volumetric liquid water content (-)
  liquid_max = maxval(abs( rVec(ixSnowSoilWat) ) )
  liquid_loc = maxloc(abs( rVec(ixSnowSoilWat) ) )

  ! check convergence based on the iteration increment for matric head
  ! NOTE: scale by matric head to avoid unnecessairly tight convergence when there is no water
  psiScale   = abs(xVec(ixSoilOnlyMat)) + xSmall ! avoid divide by zero
  matric_max = maxval(abs( xInc(ixSoilOnlyMat)/psiScale ) )
  matric_loc = maxloc(abs( xInc(ixSoilOnlyMat)/psiScale ) )

  ! convergence check
  matricConv = (matric_max(1) < absConvTol_matric)  ! NOTE: based on iteration increment
  liquidConv = (liquid_max(1) < absConvTol_liquid)  ! (based on the residual)
  energyConv = (energy_max(1) < absConvTol_energy)  ! (based on the residual)

  ! print progress towards solution
  if(globalPrintFlag)then
   print*, 'iter, dt = ', iter, dt
   write(*,'(a,1x,4(e15.5,1x),3(i4,1x),3(L1,1x))') 'fNew, matric_max(1), liquid_max(1), energy_max(1), matric_loc(1), liquid_loc(1), energy_loc(1), matricConv, liquidConv, energyConv = ', &
                                                    fNew, matric_max(1), liquid_max(1), energy_max(1), matric_loc(1), liquid_loc(1), energy_loc(1), matricConv, liquidConv, energyConv
  endif

  ! final convergence check
  checkConv = (matricConv .and. liquidConv .and. energyConv)

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
 character(LEN=256)             :: cmessage                  ! error message of downwind routine
 real(dp)                       :: volFracWat                ! total volumetric fraction of water (-)
 real(dp)                       :: fPart1,fPart2             ! different parts of the function
 real(dp)                       :: dPart1,dPart2             ! derivatives for different parts of the function
 integer(i4b)                   :: iLayer                    ! index of soil layer
 ! --------------------------------------------------------------
 ! initialize error control
 err=0; message='soilCmpres/'
 ! (only compute for the mixed form of Richards' equation)
 if(ixRichards==mixdform)then
  do iLayer=1,nSoil
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

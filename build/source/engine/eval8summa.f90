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

module eval8summa_module
! data types
USE nrtype
! access the global print flag
USE data_struc,only:globalPrintFlag
! constants
USE multiconst,only:&
                    Tfreeze,      & ! temperature at freezing              (K)
                    LH_fus,       & ! latent heat of fusion                (J kg-1)
                    LH_vap,       & ! latent heat of vaporization          (J kg-1)
                    LH_sub,       & ! latent heat of sublimation           (J kg-1)
                    Cp_air,       & ! specific heat of air                 (J kg-1 K-1)
                    iden_air,     & ! intrinsic density of air             (kg m-3)
                    iden_ice,     & ! intrinsic density of ice             (kg m-3)
                    iden_water      ! intrinsic density of liquid water    (kg m-3)
! layer types
USE data_struc,only:ix_soil,ix_snow ! named variables for snow and soil
! provide access to the derived types to define the data structures
USE data_struc,only:&
                    var_i,        & ! data vector (i4b)
                    var_d,        & ! data vector (dp)
                    var_ilength,  & ! data vector with variable length dimension (i4b)
                    var_dlength,  & ! data vector with variable length dimension (dp)
                    model_options   ! defines the model decisions
! provide access to indices that define elements of the data structures
USE var_lookup,only:iLookATTR,iLookTYPE,iLookPARAM,iLookFORCE,iLookMVAR,iLookINDEX ! named variables for structure elements
USE var_lookup,only:iLookDECISIONS                                                 ! named variables for elements of the decision structure
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
implicit none
private
! control parameters
real(dp),parameter  :: valueMissing=-9999._dp     ! missing value
real(dp),parameter  :: verySmall=tiny(1.0_dp)     ! a very small number
real(dp),parameter  :: veryBig=1.e+20_dp          ! a very big number
real(dp),parameter  :: dx = 1.e-8_dp              ! finite difference increment
contains

 ! **********************************************************************************************************
 ! public subroutine eval8summa: compute the residual vector and the Jacobian matrix
 ! **********************************************************************************************************
 subroutine eval8summa(&
                       ! input
                       dt,                      & ! intent(in):    length of the time step (seconds)
                       firstSubStep,            & ! intent(in):    flag to indicate if we are processing the first sub-step
                       firstFluxCall,           & ! intent(in):    flag to indicate if we are processing the first flux call
                       nSnow,                   & ! intent(in):    number of snow layers
                       nSoil,                   & ! intent(in):    number of soil layers
                       nLayers,                 & ! intent(in):    total number of layers
                       nState,                  & ! intent(in):    total number of state variables
                       stateVec,                & ! intent(in):    model state vector
                       sMul,                    & ! intent(inout): multiplier for state vector (used in the residual calculations)
                       dMat,                    & ! intent(inout): diagonal of the Jacobian matrix (excludes fluxes)
                       ! input/output: data structures
                       type_data,               & ! intent(in):    type of vegetation and soil
                       attr_data,               & ! intent(in):    spatial attributes
                       forc_data,               & ! intent(in):    model forcing data
                       mpar_data,               & ! intent(in):    model parameters
                       indx_data,               & ! intent(in):    index data
                       mvar_data,               & ! intent(inout): model variables for a local HRU
                       bvar_data,               & ! intent(in):    model variables for the local basin
                       model_decisions,         & ! intent(in):    model decisions
                       ! output
                       fluxVec,                 & ! intent(out):   flux vector
                       resVec,                  & ! intent(out):   residual vector
                       aJac,                    & ! intent(out):   analytical Jacobian matrix (either band or full)
                       err,message)               ! intent(out):   error control
 ! --------------------------------------------------------------------------------------------------------------------------------
 USE getVectorz_module,only:varExtract            ! extract variables from the state vector
 implicit none
 ! --------------------------------------------------------------------------------------------------------------------------------
 ! --------------------------------------------------------------------------------------------------------------------------------
 ! input
 real(dp),intent(in)             :: dt            ! length of the time step (seconds)
 logical(lgt),intent(in)         :: firstSubStep  ! flag to indicate if we are processing the first sub-step
 logical(lgt),intent(inout)      :: firstFluxCall ! flag to indicate if we are processing the first flux call
 integer(i4b),intent(in)         :: nSnow         ! number of snow layers
 integer(i4b),intent(in)         :: nSoil         ! number of soil layers
 integer(i4b),intent(in)         :: nLayers       ! total number of layers
 integer(i4b),intent(in)         :: nState        ! total number of state variables
 real(dp),intent(in)             :: stateVec(:)   ! model state vector (mixed units)
 real(dp),intent(inout)          :: sMul(:)       ! multiplier for state vector (used in the residual calculations)
 real(dp),intent(inout)          :: dMat(:)       ! diagonal of the Jacobian matrix (excludes fluxes)
 ! input/output: data structures
 type(var_i),intent(in)          :: type_data     ! type of vegetation and soil
 type(var_d),intent(in)          :: attr_data     ! spatial attributes
 type(var_d),intent(in)          :: forc_data     ! model forcing data
 type(var_d),intent(in)          :: mpar_data     ! model parameters
 type(var_ilength),intent(in)    :: indx_data     ! indices defining model states and layers
 type(var_dlength),intent(inout) :: mvar_data     ! model variables for a local HRU
 type(var_dlength),intent(in)    :: bvar_data     ! model variables for the local basin
 type(model_options),intent(in)  :: model_decisions(:)  ! model decisions
 ! output
 real(dp),intent(out)            :: fluxVec(:)    ! flux vector
 real(dp),intent(out)            :: resVec(:)     ! residual vector
 real(dp),intent(out)            :: aJac(:,:)     ! analytical Jacobian matrix (either band or full)
 integer(i4b),intent(out)        :: err           ! error code
 character(*),intent(out)        :: message       ! error message
 ! --------------------------------------------------------------------------------------------------------------------------------
 ! local variables
 ! --------------------------------------------------------------------------------------------------------------------------------
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
 ! cross-derivative terms associated with melt-freeze 
 real(dp)                        :: theta                        ! liquid water equivalent of total water (liquid plus ice)
 real(dp),dimension(nLayers)     :: mLayerdTheta_dTk             ! derivative of volumetric liquid water content w.r.t. temperature (K-1)
 real(dp),dimension(nSoil)       :: dPsiLiq_dTemp                ! derivative in the liquid water matric potential w.r.t. temperature (m K-1)
 real(dp),dimension(nSnow)       :: fracLiqSnow                  ! fraction of liquid water in each snow layer (-)
 real(dp)                        :: fracLiqVeg                   ! fraction of liquid water on vegetation (-)
 real(dp)                        :: totalWaterVeg                ! total water on vegetation (kg m-2)
 real(dp)                        :: dTheta_dTkCanopy             ! derivative of volumetric liquid water content w.r.t. temperature (K-1)
 real(dp)                        :: dCanLiq_dTcanopy             ! derivative of canopy liquid storage w.r.t. temperature (kg m-2 K-1)
 ! Jacobian matrices
 logical(lgt),parameter          :: numericalJacobian=.false.    ! flag to compute the Jacobian matrix
 logical(lgt),parameter          :: testBandDiagonal=.false.     ! flag to test the band-diagonal matrix
 logical(lgt),parameter          :: forceFullMatrix=.false.      ! flag to force the use of the full Jacobian matrix
 integer(i4b),parameter          :: ixFullMatrix=1001            ! named variable for the full Jacobian matrix
 integer(i4b),parameter          :: ixBandMatrix=1002            ! named variable for the band diagonal matrix
 integer(i4b)                    :: ixSolve                      ! the type of matrix used to solve the linear system A.X=B
 integer(i4b),parameter          :: ku=3                         ! number of super-diagonal bands
 integer(i4b),parameter          :: kl=4                         ! number of sub-diagonal bands
 integer(i4b),parameter          :: nBands=2*kl+ku+1             ! length of the leading dimension of the band diagonal matrix
 ! general local variables
 logical(lgt)                    :: computeVegFlux               ! flag to compute the vegetation flux
 real(dp)                        :: canopyDepth                  ! canopy depth (m)
 real(dp)                        :: soilWaterBalanceError        ! water balance error for soil
 integer(i4b)                    :: local_ixGroundwater          ! local index for groundwater representation
 logical(lgt)                    :: printFlag                    ! flag to control printing (set to false for numerical jacobian)
 logical(lgt)                    :: printFlagInit                ! initialize flag to control printing
 integer(i4b)                    :: iState                       ! index of matrix row
 integer(i4b)                    :: jState                       ! index of matrix column
 integer(i4b)                    :: iLayer                       ! index of model layer
 integer(i4b)                    :: jLayer                       ! index of model layer within the full state vector (hydrology)
 integer(i4b)                    :: kLayer                       ! index of model layer within the snow-soil domain
 integer(i4b)                    :: mLayer                       ! index of model layer within the full state vector (thermodynamics)
 integer(i4b),parameter          :: nVarSnowSoil=2               ! number of state variables in the snow and soil domain (energy and total water/matric head)
 integer(i4b),parameter          :: iJac1=1                      ! first layer of the Jacobian to print
 integer(i4b),parameter          :: iJac2=10                     ! last layer of the Jacobian to print
 character(len=256)              :: cMessage                     ! error message of downwind routine
 ! --------------------------------------------------------------------------------------------------------------------------------
 ! association to variables in the data structures
 ! --------------------------------------------------------------------------------------------------------------------------------
 associate(&
 ! model decisions
 ixRichards              => model_decisions(iLookDECISIONS%f_Richards)%iDecision   ,&  ! intent(in): [i4b] index of the form of Richards' equation
 ixGroundwater           => model_decisions(iLookDECISIONS%groundwatr)%iDecision   ,&  ! intent(in): [i4b] groundwater parameterization
 ixSpatialGroundwater    => model_decisions(iLookDECISIONS%spatial_gw)%iDecision   ,&  ! intent(in): [i4b] spatial representation of groundwater (local-column or single-basin)

 ! vegetation parameters
 nVegState               => indx_data%var(iLookINDEX%nVegState)%dat(1)             ,&  ! intent(in): [i4b]    number of vegetation state variables
 heightCanopyTop         => mpar_data%var(iLookPARAM%heightCanopyTop)              ,&  ! intent(in): [dp] height of the top of the vegetation canopy (m)
 heightCanopyBottom      => mpar_data%var(iLookPARAM%heightCanopyBottom)           ,&  ! intent(in): [dp] height of the bottom of the vegetation canopy (m)

 ! diagnostic variables
 scalarBulkVolHeatCapVeg => mvar_data%var(iLookMVAR%scalarBulkVolHeatCapVeg)%dat(1),&  ! intent(in): [dp   ] bulk volumetric heat capacity of vegetation (J m-3 K-1)
 mLayerVolHtCapBulk      => mvar_data%var(iLookMVAR%mLayerVolHtCapBulk)%dat        ,&  ! intent(in): [dp(:)] bulk volumetric heat capacity in each snow and soil layer (J m-3 K-1)

 ! indices for specific state variables
 ixVegNrg                => indx_data%var(iLookINDEX%ixVegNrg)%dat(1)              ,&  ! intent(in): [i4b] index of canopy energy state variable
 ixSnowSoilNrg           => indx_data%var(iLookINDEX%ixSnowSoilNrg)%dat            ,&  ! intent(in): [i4b(:)] indices for energy states in the snow-soil subdomain
 ixSoilOnlyHyd           => indx_data%var(iLookINDEX%ixSoilOnlyHyd)%dat             &  ! intent(in): [i4b(:)] indices for hydrology states in the soil subdomain

 ) ! association to variables in the data structures

 ! -----
 ! * preliminaries...
 ! ------------------

 ! NOTE: these preliminaries could be done outside this routine, which is called EVERY iteration!!

 ! set the flag to control printing
 printFlagInit=globalPrintFlag
 printFlag=printFlagInit

 ! identify the matrix solution method
 ! (the type of matrix used to solve the linear system A.X=B)
 if(ixGroundwater==qbaseTopmodel .or. testBandDiagonal .or. forceFullMatrix .or. numericalJacobian)then
  ixSolve=ixFullMatrix   ! full Jacobian matrix
 else
  ixSolve=ixBandMatrix   ! band-diagonal matrix
 endif
 if(printFlag) print*, '(ixSolve==ixFullMatrix) = ', (ixSolve==ixFullMatrix)

 ! define the need to compute the fluxes over vegetation
 computeVegFlux = (nVegState==3)
 if(.not.computeVegFlux .and. nVegState/=0)then
  message=trim(message)//'expect the number of vegetation states to be zero when not computing the vegetation flux'
  err=20; return
 endif

 ! compute the canopy depth (m)
 canopyDepth = heightCanopyTop - heightCanopyBottom

 ! modify the groundwater representation for this single-column implementation
 select case(ixSpatialGroundwater)
  case(singleBasin); local_ixGroundwater = noExplicit    ! force no explicit representation of groundwater at the local scale
  case(localColumn); local_ixGroundwater = ixGroundwater ! go with the specified decision
  case default; err=20; message=trim(message)//'unable to identify spatial representation of groundwater'; return
 end select ! (modify the groundwater representation for this single-column implementation)

 ! -----
 ! * compute the residual vector...
 ! --------------------------------

 ! compute residual vector
 call xFluxResid(&
                 ! input
                 stateVec,      & ! intent(in): full state vector (mixed units)
                 ! output
                 fluxVec,       & ! intent(out): flux vector (mixed units)
                 resVec,        & ! intent(out): residual vector (mixed units)
                 err,cmessage)    ! intent(out): error code and error message
 if(err/=0)then; message=trim(message)//trim(cmessage); return; endif  ! (check for errors)

 ! -----
 ! * compute the Jacobian matrix...
 ! --------------------------------

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
  call numlJacob(stateVec,fluxVec,resVec,err,cmessage)
  if(err/=0)then; message=trim(message)//trim(cmessage); return; endif  ! (check for errors)
  printFlag=printFlagInit
 endif  ! if computing the numerical Jacobian matrix

 ! end associtation to variables in the data structures
 end associate


 ! start of internal subroutines 
 ! **********************************************************************************************************
 ! **********************************************************************************************************
 ! **********************************************************************************************************
 ! **********************************************************************************************************
 ! **********************************************************************************************************
 ! **********************************************************************************************************
 ! **********************************************************************************************************
 ! **********************************************************************************************************
 ! **********************************************************************************************************
 ! **********************************************************************************************************
 ! **********************************************************************************************************
 contains

  ! **********************************************************************************************************
  ! internal subroutine xFluxResid: compute the model fluxes and the residual vector
  ! **********************************************************************************************************
  subroutine xFluxResid(&
                        ! input
                        stateVec,                       & ! intent(in): full state vector (mixed units)
                        ! output
                        fVec,                           & ! intent(out): flux vector (mixed units)
                        rVec,                           & ! intent(out): residual vector (mixed units)
                        err,message)                      ! intent(out): error code and error message
  ! --------------------------------------------------------------
  implicit none
  ! input variables
  real(dp),intent(in)            :: stateVec(:)           ! model state vector (mixed units)
  ! output variabes
  real(dp),intent(out)           :: fVec(:)               ! flux vector (mixed units)
  real(dp),intent(out)           :: rVec(:)               ! residual vector (mixed units)
  integer(i4b),intent(out)       :: err                   ! error code
  character(*),intent(out)       :: message               ! error message
  ! --------------------------------------------------------------
  ! general local variables
  real(dp),dimension(nState)     :: rAdd                  ! additional terms in the residual vector
  real(dp)                       :: scalarCanopyWatInit   ! initial value for liquid water storage in the canopy (kg m-2)
  real(dp),dimension(nSnow)      :: mLayerVolFracWatInit  ! initial value for volumetric fraction of total water (-)
  real(dp),dimension(nSoil)      :: vThetaInit            ! liquid equivalent of total water at the start of the step
  real(dp),dimension(nSoil)      :: vThetaTrial           ! liquid equivalent of total water at the current iteration
  character(LEN=256)             :: cmessage              ! error message of downwind routine
  ! --------------------------------------------------------------
  ! initialize error control
  err=0; message='xFluxResid/'

  ! -----
  ! * associate desired variables from data structures...
  ! -----------------------------------------------------
  associate(&
  ! model decisions
  ixRichards              => model_decisions(iLookDECISIONS%f_Richards)%iDecision   ,&  ! intent(in): [i4b] index of the form of Richards' equation

  ! soil parameters
  theta_sat               => mpar_data%var(iLookPARAM%theta_sat)                    ,&  ! intent(in): [dp] soil porosity (-)
  specificStorage         => mpar_data%var(iLookPARAM%specificStorage)              ,&  ! intent(in): [dp] specific storage coefficient (m-1)

  ! indices of model state variables
  ixCasNrg                => indx_data%var(iLookINDEX%ixCasNrg)%dat(1)              ,&  ! intent(in): [i4b] index of canopy air space energy state variable
  ixVegNrg                => indx_data%var(iLookINDEX%ixVegNrg)%dat(1)              ,&  ! intent(in): [i4b] index of canopy energy state variable
  ixVegWat                => indx_data%var(iLookINDEX%ixVegWat)%dat(1)              ,&  ! intent(in): [i4b] index of canopy hydrology state variable (mass)
  ixSnowSoilNrg           => indx_data%var(iLookINDEX%ixSnowSoilNrg)%dat            ,&  ! intent(in): [i4b(:)] indices for energy states in the snow-soil subdomain
  ixSnowSoilWat           => indx_data%var(iLookINDEX%ixSnowSoilWat)%dat            ,&  ! intent(in): [i4b(:)] indices for total water states in the snow-soil subdomain
  ixSnowOnlyNrg           => indx_data%var(iLookINDEX%ixSnowOnlyNrg)%dat            ,&  ! intent(in): [i4b(:)] indices for energy states in the snow subdomain
  ixSnowOnlyWat           => indx_data%var(iLookINDEX%ixSnowOnlyWat)%dat            ,&  ! intent(in): [i4b(:)] indices for total water states in the snow subdomain
  ixSoilOnlyNrg           => indx_data%var(iLookINDEX%ixSoilOnlyNrg)%dat            ,&  ! intent(in): [i4b(:)] indices for energy states in the soil subdomain
  ixSoilOnlyHyd           => indx_data%var(iLookINDEX%ixSoilOnlyHyd)%dat            ,&  ! intent(in): [i4b(:)] indices for hydrology states in the soil subdomain

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

  ! extract summa variables from the model state vector
  call varExtract(&
                  ! input
                  stateVec,                      & ! intent(in):    model state vector (mixed units)
                  mpar_data,                     & ! intent(in):    model parameters
                  mvar_data,                     & ! intent(in):    model variables for a local HRU
                  indx_data,                     & ! intent(in):    indices defining model states and layers
                  ! output: variables for the vegetation canopy
                  fracLiqVeg,                    & ! intent(out):   fraction of liquid water on the vegetation canopy (-)
                  scalarCanairTempTrial,         & ! intent(out):   trial value of canopy air temperature (K)
                  scalarCanopyTempTrial,         & ! intent(out):   trial value of canopy temperature (K)
                  scalarCanopyWatTrial,          & ! intent(out):   trial value of canopy total water (kg m-2)
                  scalarCanopyLiqTrial,          & ! intent(out):   trial value of canopy liquid water (kg m-2)
                  scalarCanopyIceTrial,          & ! intent(out):   trial value of canopy ice content (kg m-2)
                  ! output: variables for the snow-soil domain
                  fracLiqSnow,                   & ! intent(out):   volumetric fraction of water in each snow layer (-)
                  mLayerTempTrial,               & ! intent(out):   trial vector of layer temperature (K)
                  mLayerVolFracWatTrial,         & ! intent(out):   trial vector of volumetric total water content (-)
                  mLayerVolFracLiqTrial,         & ! intent(out):   trial vector of volumetric liquid water content (-)
                  mLayerVolFracIceTrial,         & ! intent(out):   trial vector of volumetric ice water content (-)
                  mLayerMatricHeadTrial,         & ! intent(out):   trial vector of matric head (m)
                  ! output: error control
                  err,cmessage)                    ! intent(out):   error control
  if(err/=0)then; message=trim(message)//trim(cmessage); return; endif  ! (check for errors)

  ! compute model flux for a given state vector
  call computFlux(&
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
                  mLayerVolFracLiqTrial(nSnow+1:nLayers), & ! intent(in): trial value for the volumetric liquid water content in each soil layer (-)
                  mLayerVolFracIceTrial(nSnow+1:nLayers), & ! intent(in): trial value for the volumetric ice content in each soil layer (-)
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
   rAdd(ixVegNrg) = rAdd(ixVegNrg) + LH_fus*(scalarCanopyIceTrial - scalarCanopyIce)/canopyDepth   ! energy associated with melt/freeze (J m-3)
  endif

  ! compute energy associated with melt/freeze for snow
  if(nSnow>0)&
  rAdd(ixSnowOnlyNrg) = rAdd(ixSnowOnlyNrg) + LH_fus*iden_ice*(mLayerVolFracIceTrial(1:nSnow) - mLayerVolFracIce(1:nSnow))       ! energy associated with melt/freeze (J m-3)

  ! compute energy associated with melt/freeze for soil
  rAdd(ixSoilOnlyNrg) = rAdd(ixSoilOnlyNrg) + LH_fus*iden_water*(mLayerVolFracIceTrial(nSnow+1:nLayers) - mLayerVolFracIce(nSnow+1:nLayers))     ! energy associated with melt/freeze (J m-3)

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
   scalarCanopyWatInit = scalarCanopyLiq + scalarCanopyIce   ! kg m-2
   rVec(ixVegWat) = sMul(ixVegWat)*scalarCanopyWatTrial  - ( (sMul(ixVegWat)*scalarCanopyWatInit  + fVec(ixVegWat)*dt) + rAdd(ixVegWat) )
  endif

  ! compute the residual vector for the snow and soil sub-domains for energy
  rVec(ixSnowSoilNrg) = sMul(ixSnowSoilNrg)*mLayerTempTrial(1:nLayers) - ( (sMul(ixSnowSoilNrg)*mLayerTemp(1:nLayers)  + fVec(ixSnowSoilNrg)*dt) + rAdd(ixSnowSoilNrg) )

  ! compute the residual vector for the **snow** sub-domain for liquid water
  if(nSnow>0)then
   mLayerVolFracWatInit(1:nSnow) = mLayerVolFracLiq(1:nSnow) + mLayerVolFracIce(1:nSnow)*(iden_ice/iden_water)
   rVec(ixSnowOnlyWat) = mLayerVolFracWatTrial(1:nSnow) - ( (mLayerVolFracWatInit(1:nSnow)  + fVec(ixSnowOnlyWat)*dt) + rAdd(ixSnowOnlyWat) )
  endif

  ! compute the residual vector for the **soil** sub-domain for liquid water
  vThetaInit(1:nSoil)  = mLayerVolFracLiq(nSnow+1:nLayers)      + mLayerVolFracIce(nSnow+1:nLayers)      ! liquid equivalent of total water at the start of the step
  vThetaTrial(1:nSoil) = mLayerVolFracLiqTrial(nSnow+1:nLayers) + mLayerVolFracIceTrial(nSnow+1:nLayers) ! liquid equivalent of total water at the current iteration
  rVec(ixSoilOnlyHyd)  = vThetaTrial(1:nSoil) - ( (vThetaInit(1:nSoil) + fVec(ixSoilOnlyHyd)*dt) + rAdd(ixSoilOnlyHyd) )

  ! compute the soil water balance error (m)
  ! NOTE: declared in the main routine so accessible in all internal routines
  soilWaterBalanceError = abs( sum(real(rVec(ixSoilOnlyHyd), dp)*mLayerDepth(nSnow+1:nSoil)) )

  ! end association to variables in the data structures
  end associate

  if(printFlag) write(*,'(a,1x,100(e12.5,1x))') 'rVec = ', rVec

  end subroutine xFluxResid


  ! **********************************************************************************************************
  ! **********************************************************************************************************


  ! **********************************************************************************************************
  ! internal subroutine computFlux: compute the model fluxes
  ! **********************************************************************************************************
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
  ! ---------------------------------------------------------------------------------------
  ! utility modules
  USE snow_utils_module,only:dFracLiq_dTk              ! differentiate the freezing curve w.r.t. temperature (snow)
  USE soil_utils_module,only:dTheta_dTk                ! differentiate the freezing curve w.r.t. temperature (soil)
  USE soil_utils_module,only:dTheta_dPsi               ! derivative in the soil water characteristic (soil)
  USE soil_utils_module,only:dPsi_dTheta               ! derivative in the soil water characteristic (soil)
  USE soil_utils_module,only:matricHead                ! compute the matric head based on volumetric water content
  ! flux modules
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
  integer(i4b)                   :: iLayer                    ! index of model layer
  real(dp),parameter             :: canopyTempMax=500._dp     ! expected maximum value for the canopy temperature (K)
  real(dp)                       :: xNum                      ! temporary variable: numerator
  real(dp)                       :: xDen                      ! temporary variable: denominator
  real(dp)                       :: effSat                    ! effective saturation of the soil matrix (-)
  real(dp),dimension(nSoil)      :: mLayerMatricHeadLiq       ! matric head associated with liquid water (m), f(psi0, T)
  real(dp)                       :: dPsiLiq_dEffSat           ! derivative in liquid water matric potential w.r.t. effective saturation (m)
  real(dp)                       :: dEffSat_dVolTot           ! derivative in effective saturation w.r.t. total water content (-)
  real(dp)                       :: dEffSat_dTemp             ! derivative in effective saturation w.r.t. temperature (K-1)
  real(dp),dimension(nSoil)      :: dPsiLiq_dPsi0             ! derivative in the liquid water matric potential w.r.t. the total water matric potential (-)
  character(LEN=256)             :: cmessage                  ! error message of downwind routine
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
  scalarRainfall          => mvar_data%var(iLookMVAR%scalarRainfall)%dat(1)         ,&  ! intent(in): [dp] rainfall rate (kg m-2 s-1)
  scalarSfcMeltPond       => mvar_data%var(iLookMVAR%scalarSfcMeltPond)%dat(1)      ,&  ! intent(in): [dp] ponded water caused by melt of the "snow without a layer" (kg m-2)

  ! indices of model state variables
  ixCasNrg                => indx_data%var(iLookINDEX%ixCasNrg)%dat(1)              , & ! intent(in): [i4b] index of canopy air space energy state variable
  ixVegNrg                => indx_data%var(iLookINDEX%ixVegNrg)%dat(1)              , & ! intent(in): [i4b] index of canopy energy state variable
  ixVegWat                => indx_data%var(iLookINDEX%ixVegWat)%dat(1)              , & ! intent(in): [i4b] index of canopy hydrology state variable (mass)
  ixSnowSoilNrg           => indx_data%var(iLookINDEX%ixSnowSoilNrg)%dat            , & ! intent(in): [i4b(:)] indices for energy states in the snow-soil subdomain
  ixSnowSoilWat           => indx_data%var(iLookINDEX%ixSnowSoilWat)%dat            , & ! intent(in): [i4b(:)] indices for total water states in the snow-soil subdomain
  ixSnowOnlyNrg           => indx_data%var(iLookINDEX%ixSnowOnlyNrg)%dat            , & ! intent(in): [i4b(:)] indices for energy states in the snow subdomain
  ixSnowOnlyWat           => indx_data%var(iLookINDEX%ixSnowOnlyWat)%dat            , & ! intent(in): [i4b(:)] indices for total water states in the snow subdomain
  ixSoilOnlyNrg           => indx_data%var(iLookINDEX%ixSoilOnlyNrg)%dat            , & ! intent(in): [i4b(:)] indices for energy states in the soil subdomain
  ixSoilOnlyHyd           => indx_data%var(iLookINDEX%ixSoilOnlyHyd)%dat            ,&  ! intent(in): [i4b(:)] indices for hydrology states in the soil subdomain
  layerType               => indx_data%var(iLookINDEX%layerType)%dat                , & ! intent(in): [i4b(:)] type of each layer in the snow+soil domain (snow or soil)

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
  iLayerNrgFlux           => mvar_data%var(iLookMVAR%iLayerNrgFlux)%dat             ,&  ! intent(out): [dp(0:)] vertical energy flux at the interface of snow and soil layers
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
  if(firstFluxCall)then
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
                   mvar_data,                               & ! intent(inout): model variables for a local HRU
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


  ! **********************************************************************************************************
  ! **********************************************************************************************************


  ! *********************************************************************************************************
  ! internal subroutine cpactBand: compute the compact band-diagonal matric
  ! *********************************************************************************************************
  subroutine cpactBand(err,message)
  implicit none
  ! dummy variables
  integer(i4b),intent(out)       :: err                     ! error code
  character(*),intent(out)       :: message                 ! error message
  ! named variables for elements of the compact band-diagonal matrix
  integer(i4b),parameter         :: ixSup3=kl+1             ! index for the 3rd super-diagonal band
  integer(i4b),parameter         :: ixSup2=kl+2             ! index for the 2nd super-diagonal band
  integer(i4b),parameter         :: ixSup1=kl+3             ! index for the 1st super-diagonal band
  integer(i4b),parameter         :: ixDiag=kl+4             ! index for the diagonal band
  integer(i4b),parameter         :: ixSub1=kl+5             ! index for the 1st sub-diagonal band
  integer(i4b),parameter         :: ixSub2=kl+6             ! index for the 2nd sub-diagonal band
  integer(i4b),parameter         :: ixSub3=kl+7             ! index for the 3rd sub-diagonal band
  integer(i4b),parameter         :: ixSub4=kl+8             ! index for the 3rd sub-diagonal band
  ! local variables
  ! -------------------------------------------------------------
  ! initialize error control
  err=0; message='cpactBand/'

  ! associate variables from data structures
  associate(&

            ! indices of model state variables
            ixCasNrg      => indx_data%var(iLookINDEX%ixCasNrg)%dat(1)  , & ! index of canopy air space energy state variable
            ixVegNrg      => indx_data%var(iLookINDEX%ixVegNrg)%dat(1)  , & ! index of canopy energy state variable
            ixVegWat      => indx_data%var(iLookINDEX%ixVegWat)%dat(1)  , & ! index of canopy hydrology state variable (mass)
            ixTopNrg      => indx_data%var(iLookINDEX%ixTopNrg)%dat(1)  , & ! index of upper-most energy state in the snow-soil subdomain
            ixSnowSoilNrg => indx_data%var(iLookINDEX%ixSnowSoilNrg)%dat, & ! indices for energy states in the snow-soil subdomain
            ixSnowOnlyWat => indx_data%var(iLookINDEX%ixSnowOnlyWat)%dat, & ! indices for total water states in the snow subdomain
            ixSoilOnlyHyd => indx_data%var(iLookINDEX%ixSoilOnlyHyd)%dat, & ! indices for hydrology states in the soil subdomain

            ! layer depth
            mLayerDepth => mvar_data%var(iLookMVAR%mLayerDepth)%dat) ! intent(in): [dp(:)] depth of each layer in the snow-soil sub-domain (m)

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


  ! **********************************************************************************************************
  ! **********************************************************************************************************


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
            mLayerDepth   => mvar_data%var(iLookMVAR%mLayerDepth)%dat             )  ! intent(in): [dp(:)] depth of each layer in the snow-soil sub-domain (m)

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


  ! **********************************************************************************************************
  ! **********************************************************************************************************


  ! *********************************************************************************************************
  ! internal subroutine numlJacob: compute the Jacobian matrix (numerical)
  ! *********************************************************************************************************
  subroutine numlJacob(stateVec,fluxVec,resVec,err,message)
  implicit none
  ! dummy
  real(dp),intent(in)            :: stateVec(:)             ! model state vector (mixed units)
  real(dp),intent(in)            :: fluxVec(:)              ! model flux vector (mixed units)
  real(dp),intent(in)            :: resVec(:)               ! model residual vector (mixed units)
  integer(i4b),intent(out)       :: err                     ! error code
  character(*),intent(out)       :: message                 ! error message
  ! local
  character(len=256)             :: cmessage                ! error message of downwind routine
  real(dp),dimension(nState)     :: stateVecPerturbed       ! perturbed state vector
  real(dp),dimension(nState)     :: fluxVecJac              ! flux vector
  real(dp),dimension(nState)     :: resVecJac               ! residual vector (mixed units)
  real(dp),dimension(nState,nState) :: nJac                 ! the numerical Jacobian matrix
  integer(i4b)                   :: iJac                    ! index of row of the Jacobian matrix
  integer(i4b),parameter         :: iTry=-999               ! index of trial model state variable (used for testing)
  integer(i4b),parameter         :: ixNumFlux=1001          ! named variable for the flux-based form of the numerical Jacobian
  integer(i4b),parameter         :: ixNumRes=1002           ! named variable for the residual-based form of the numerical Jacobian
  integer(i4b)                   :: ixNumType=ixNumRes      ! method used to calculate the numerical Jacobian
  ! --------------------------------------------------------------
  ! initialize error control
  err=0; message='numlJacob/'

  ! get a copy of the state vector to perturb
  stateVecPerturbed(:) = stateVec(:)

  ! loop through state variables
  do iJac=1,nState

   ! define printFlag
   if(iJac==iTry) printFlag=.true.
   if(iJac/=iTry) printFlag=.false.

   ! (perturb state vector)
   stateVecPerturbed(iJac) = stateVec(iJac) + dx

   ! **
   ! ** residual-based calculation of the numerical Jacobian
   if(ixNumType==ixNumRes)then ! switch between the residual-based form and flux-based form

    ! (compute residual vector)
    call xFluxResid(&
                    ! input
                    stateVecPerturbed,              & ! intent(in): full state vector (mixed units)
                    ! output
                    fluxVecJac,                     & ! intent(out): flux vector (mixed units)
                    resVecJac,                      & ! intent(out): residual vector (mixed units)
                    err,cmessage)                     ! intent(out): error code and error message
    if(err/=0)then; message=trim(message)//trim(cmessage); return; endif  ! (check for errors)

    ! (compute the row of the Jacobian matrix)
    nJac(:,iJac) = (resVecJac - resVec)/dx

   ! **
   ! ** flux-based calculation of the numerical Jacobian
   else

    ! extract summa variables from the model state vector
    call varExtract(&
                    ! input
                    stateVecPerturbed,             & ! intent(in):    model state vector (mixed units)
                    mpar_data,                     & ! intent(in):    model parameters
                    mvar_data,                     & ! intent(in):    model variables for a local HRU
                    indx_data,                     & ! intent(in):    indices defining model states and layers
                    ! output: variables for the vegetation canopy
                    fracLiqVeg,                    & ! intent(out):   fraction of liquid water on the vegetation canopy (-)
                    scalarCanairTempTrial,         & ! intent(out):   trial value of canopy air temperature (K)
                    scalarCanopyTempTrial,         & ! intent(out):   trial value of canopy temperature (K)
                    scalarCanopyWatTrial,          & ! intent(out):   trial value of canopy total water (kg m-2)
                    scalarCanopyLiqTrial,          & ! intent(out):   trial value of canopy liquid water (kg m-2)
                    scalarCanopyIceTrial,          & ! intent(out):   trial value of canopy ice content (kg m-2)
                    ! output: variables for the snow-soil domain
                    fracLiqSnow,                   & ! intent(out):   volumetric fraction of water in each snow layer (-)
                    mLayerTempTrial,               & ! intent(out):   trial vector of layer temperature (K)
                    mLayerVolFracWatTrial,         & ! intent(out):   trial vector of volumetric total water content (-)
                    mLayerVolFracLiqTrial,         & ! intent(out):   trial vector of volumetric liquid water content (-)
                    mLayerVolFracIceTrial,         & ! intent(out):   trial vector of volumetric ice water content (-)
                    mLayerMatricHeadTrial,         & ! intent(out):   trial vector of matric head (m)
                    ! output: error control
                    err,cmessage)                    ! intent(out):   error control
    if(err/=0)then; message=trim(message)//trim(cmessage); return; endif  ! (check for errors)

    ! compute model flux for a given state vector
    call computFlux(&
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

  end do  ! (looping through state variables)

  ! print the Jacobian
  print*, '** numerical Jacobian:', ixNumType==ixNumRes
  write(*,'(a4,1x,100(i12,1x))') 'xCol', (iLayer, iLayer=iJac1,iJac2)
  do iJac=iJac1,iJac2; write(*,'(i4,1x,100(e12.5,1x))') iJac, nJac(iJac1:iJac2,iJac); end do
  !pause 'testing Jacobian'

  end subroutine numlJacob


  ! **********************************************************************************************************
  ! **********************************************************************************************************


  subroutine testBandMatrix()
  real(dp)  :: aJac_test(nState,nState)
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
   write(*,'(i4,1x,100(e17.10,1x))') iLayer, (aJac_test(iLayer,iState),iState=iJac1,iJac2)
  end do
  end subroutine testBandMatrix

 end subroutine eval8summa


 ! **********************************************************************************************************
 ! **********************************************************************************************************


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



end module eval8summa_module

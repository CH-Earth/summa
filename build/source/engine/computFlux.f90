! SUMMA - Structure for Unifying Multiple Modeling Alternatives
! Copyright (C) 2014-2020 NCAR/RAL; University of Saskatchewan; University of Washington
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

module computFlux_module

! data types
USE nrtype

! provide access to the derived types to define the data structures
USE data_types,only:&
                    var_i,              & ! data vector (i4b)
                    var_d,              & ! data vector (rkind)
                    var_ilength,        & ! data vector with variable length dimension (i4b)
                    var_dlength,        & ! data vector with variable length dimension (rkind)
                    model_options,      & ! defines the model decisions
                    in_type_vegNrgFlux,out_type_vegNrgFlux,                   & ! arguments for vegNrgFlux call
                    in_type_ssdNrgFlux,io_type_ssdNrgFlux,out_type_ssdNrgFlux,& ! arguments for ssdNrgFlux call
                    in_type_vegLiqFlux,out_type_vegLiqFlux,                   & ! arguments for vegLiqFlux call
                    in_type_snowLiqFlx,io_type_snowLiqFlx,out_type_snowLiqFlx,& ! arguments for snowLiqFlx call                
                    in_type_soilLiqFlx,io_type_soilLiqFlx,out_type_soilLiqFlx,& ! arguments for soilLiqFlx call
                    in_type_groundwatr,io_type_groundwatr,out_type_groundwatr,& ! arguments for groundwatr call
                    in_type_bigAquifer,io_type_bigAquifer,out_type_bigAquifer   ! arguments for bigAquifer call

! indices that define elements of the data structures
USE var_lookup,only:iLookDECISIONS  ! named variables for elements of the decision structure
USE var_lookup,only:iLookPARAM      ! named variables for structure elements
USE var_lookup,only:iLookFORCE      ! named variables for structure elements
USE var_lookup,only:iLookPROG       ! named variables for structure elements
USE var_lookup,only:iLookINDEX      ! named variables for structure elements
USE var_lookup,only:iLookDIAG       ! named variables for structure elements
USE var_lookup,only:iLookFLUX       ! named variables for structure elements
USE var_lookup,only:iLookDERIV      ! named variables for structure elements
USE var_lookup,only:iLookROUTINE    ! named variables for structure elements
USE var_lookup,only:iLookOP         ! named variables for structure elements

! missing values
USE globalData,only:integerMissing  ! missing integer
USE globalData,only:realMissing     ! missing real number

! layer types
USE globalData,only:iname_snow      ! named variables for snow
USE globalData,only:iname_soil      ! named variables for soil

! access the global print flag
USE globalData,only:globalPrintFlag

! control parameters
USE globalData,only:verySmall       ! a very small number
USE globalData,only:veryBig         ! a very big number
USE globalData,only:dx              ! finite difference increment

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
USE mDecisions_module,only:       &
 localColumn,                     & ! separate groundwater representation in each local soil column
 singleBasin                        ! single groundwater store over the entire basin

! look-up values for the choice of groundwater parameterization
USE mDecisions_module,only:       &
 qbaseTopmodel,                   & ! TOPMODEL-ish baseflow parameterization
 bigBucket,                       & ! a big bucket (lumped aquifer model)
 noExplicit                         ! no explicit groundwater parameterization

! look-up values for the form of Richards' equation
USE mDecisions_module,only:       &
 moisture,                        & ! moisture-based form of Richards' equation
 mixdform                           ! mixed form of Richards' equation

! look-up values for the choice of boundary conditions for hydrology
USE mDecisions_module,only:       &
 prescribedHead,                  & ! prescribed head (volumetric liquid water content for mixed form of Richards' eqn)
 funcBottomHead,                  & ! function of matric head in the lower-most layer
 freeDrainage,                    & ! free drainage
 liquidFlux,                      & ! liquid water flux
 zeroFlux                           ! zero flux

implicit none
private
public::computFlux
public::soilCmpres
public::soilCmpresPrime

contains
! *********************************************************************************************************
! public subroutine computFlux: compute model fluxes
! *********************************************************************************************************
subroutine computFlux(&
                      ! input-output: model control
                      nSnow,                    & ! intent(in):    number of snow layers
                      nSoil,                    & ! intent(in):    number of soil layers
                      nLayers,                  & ! intent(in):    total number of layers
                      firstSubStep,             & ! intent(in):    flag to indicate if we are processing the first sub-step
                      firstFluxCall,            & ! intent(inout): flag to denote the first flux call
                      firstSplitOper,           & ! intent(in):    flag to indicate if we are processing the first flux call in a splitting operation
                      computeVegFlux,           & ! intent(in):    flag to indicate if we need to compute fluxes over vegetation
                      scalarSolution,           & ! intent(in):    flag to indicate the scalar solution
                      checkLWBalance,           & ! intent(in):    flag to check longwave balance
                      drainageMeltPond,         & ! intent(in):    drainage from the surface melt pond (kg m-2 s-1)
                      ! input: state variables
                      scalarCanairTempTrial,    & ! intent(in):    trial value for the temperature of the canopy air space (K)
                      scalarCanopyTempTrial,    & ! intent(in):    trial value for the temperature of the vegetation canopy (K)
                      mLayerTempTrial,          & ! intent(in):    trial value for the temperature of each snow and soil layer (K)
                      mLayerMatricHeadLiqTrial, & ! intent(in):    trial value for the liquid water matric potential in each soil layer (m)
                      mLayerMatricHeadTrial,    & ! intent(in):    trial vector of total water matric potential (m)
                      scalarAquiferStorageTrial,& ! intent(in):    trial value of storage of water in the aquifer (m)
                      ! input: diagnostic variables defining the liquid water and ice content
                      scalarCanopyLiqTrial,     & ! intent(in):    trial value for the liquid water on the vegetation canopy (kg m-2)
                      scalarCanopyIceTrial,     & ! intent(in):    trial value for the ice on the vegetation canopy (kg m-2)
                      mLayerVolFracLiqTrial,    & ! intent(in):    trial value for the volumetric liquid water content in each snow and soil layer (-)
                      mLayerVolFracIceTrial,    & ! intent(in):    trial value for the volumetric ice in each snow and soil layer (-)
                      ! input: data structures
                      model_decisions,          & ! intent(in):    model decisions
                      type_data,                & ! intent(in):    type of vegetation and soil
                      attr_data,                & ! intent(in):    spatial attributes
                      mpar_data,                & ! intent(in):    model parameters
                      forc_data,                & ! intent(in):    model forcing data
                      bvar_data,                & ! intent(in):    average model variables for the entire basin
                      prog_data,                & ! intent(in):    model prognostic variables for a local HRU
                      indx_data,                & ! intent(in):    index data
                      ! input-output: data structures
                      diag_data,                & ! intent(inout): model diagnostic variables for a local HRU
                      flux_data,                & ! intent(inout): model fluxes for a local HRU
                      deriv_data,               & ! intent(inout): derivatives in model fluxes w.r.t. relevant state variables
                      ! input-output: flux vector and baseflow derivatives
                      ixSaturation,             & ! intent(inout): index of the lowest saturated layer (NOTE: only computed on the first iteration)
                      dBaseflow_dMatric,        & ! intent(out):   derivative in baseflow w.r.t. matric head (s-1)
                      fluxVec,                  & ! intent(out):   flux vector (mixed units)
                      ! output: error control
                      err,message)                ! intent(out):   error code and error message
  ! provide access to flux subroutines
  USE vegNrgFlux_module,only:vegNrgFlux           ! compute energy fluxes over vegetation
  USE ssdNrgFlux_module,only:ssdNrgFlux           ! compute energy fluxes throughout the snow and soil subdomains
  USE vegLiqFlux_module,only:vegLiqFlux           ! compute liquid water fluxes through vegetation
  USE snowLiqFlx_module,only:snowLiqflx           ! compute liquid water fluxes through snow
  USE soilLiqFlx_module,only:soilLiqflx           ! compute liquid water fluxes through soil
  USE groundwatr_module,only:groundwatr           ! compute the baseflow flux
  USE bigAquifer_module,only:bigAquifer           ! compute fluxes for the big aquifer
  implicit none
  ! ---------------------------------------------------------------------------------------
  ! * dummy variables
  ! ---------------------------------------------------------------------------------------
  ! input-output: control
  integer(i4b),intent(in)            :: nSnow                       ! number of snow layers
  integer(i4b),intent(in)            :: nSoil                       ! number of soil layers
  integer(i4b),intent(in)            :: nLayers                     ! total number of layers
  logical(lgt),intent(in)            :: firstSubStep                ! flag to indicate if we are processing the first sub-step
  logical(lgt),intent(inout)         :: firstFluxCall               ! flag to indicate if we are processing the first flux call
  logical(lgt),intent(in)            :: firstSplitOper              ! flag to indicate if we are processing the first flux call in a splitting operation
  logical(lgt),intent(in)            :: computeVegFlux              ! flag to indicate if computing fluxes over vegetation
  logical(lgt),intent(in)            :: scalarSolution              ! flag to denote if implementing the scalar solution
  logical(lgt),intent(in)            :: checkLWBalance              ! flag to check longwave balance
  real(rkind),intent(in)             :: drainageMeltPond            ! drainage from the surface melt pond (kg m-2 s-1)
  ! input: state variables
  real(rkind),intent(in)             :: scalarCanairTempTrial       ! trial value for temperature of the canopy air space (K)
  real(rkind),intent(in)             :: scalarCanopyTempTrial       ! trial value for temperature of the vegetation canopy (K)
  real(rkind),intent(in)             :: mLayerTempTrial(:)          ! trial value for temperature of each snow/soil layer (K)
  real(rkind),intent(in)             :: mLayerMatricHeadLiqTrial(:) ! trial value for the liquid water matric potential (m)
  real(rkind),intent(in)             :: mLayerMatricHeadTrial(:)    ! trial value for the total water matric potential (m)
  real(rkind),intent(in)             :: scalarAquiferStorageTrial   ! trial value of aquifer storage (m)
  ! input: diagnostic variables
  real(rkind),intent(in)             :: scalarCanopyLiqTrial        ! trial value for mass of liquid water on the vegetation canopy (kg m-2)
  real(rkind),intent(in)             :: scalarCanopyIceTrial        ! trial value for mass of ice on the vegetation canopy (kg m-2)
  real(rkind),intent(in)             :: mLayerVolFracLiqTrial(:)    ! trial value for volumetric fraction of liquid water (-)
  real(rkind),intent(in)             :: mLayerVolFracIceTrial(:)    ! trial value for volumetric fraction of ice (-)
  ! input: data structures
  type(model_options),intent(in)     :: model_decisions(:)          ! model decisions
  type(var_i),        intent(in)     :: type_data                   ! type of vegetation and soil
  type(var_d),        intent(in)     :: attr_data                   ! spatial attributes
  type(var_dlength),  intent(in)     :: mpar_data                   ! model parameters
  type(var_d),        intent(in)     :: forc_data                   ! model forcing data
  type(var_dlength),  intent(in)     :: bvar_data                   ! model variables for the local basin
  type(var_dlength),  intent(in)     :: prog_data                   ! prognostic variables for a local HRU
  type(var_ilength),  intent(in)     :: indx_data                   ! indices defining model states and layers
  ! input-output: data structures
  type(var_dlength),intent(inout)    :: diag_data                   ! diagnostic variables for a local HRU
  type(var_dlength),intent(inout)    :: flux_data                   ! model fluxes for a local HRU
  type(var_dlength),intent(inout)    :: deriv_data                  ! derivatives in model fluxes w.r.t. relevant state variables
  ! input-output: flux vector and baseflow derivatives
  integer(i4b),intent(inout)         :: ixSaturation                ! index of the lowest saturated layer (NOTE: only computed on the first iteration)
  real(rkind),intent(out)            :: dBaseflow_dMatric(:,:)      ! derivative in baseflow w.r.t. matric head (s-1)
  real(rkind),intent(out)            :: fluxVec(:)                  ! model flux vector (mixed units)
  ! output: error control
  integer(i4b),intent(out)           :: err                         ! error code
  character(*),intent(out)           :: message                     ! error message
  ! ---------------------------------------------------------------------------------------
  ! * local variables
  ! ---------------------------------------------------------------------------------------
  integer(i4b)                       :: local_ixGroundwater         ! local index for groundwater representation
  integer(i4b)                       :: iLayer                      ! index of model layers
  logical(lgt)                       :: doVegNrgFlux                ! flag to compute the energy flux over vegetation
  real(rkind),dimension(nSoil)       :: dHydCond_dMatric            ! derivative in hydraulic conductivity w.r.t matric head (s-1)
  character(LEN=256)                 :: cmessage                    ! error message of downwind routine
  real(rkind)                        :: above_soilLiqFluxDeriv      ! derivative in layer above soil (canopy or snow) liquid flux w.r.t. liquid water
  real(rkind)                        :: above_soildLiq_dTk          ! derivative of layer above soil (canopy or snow) liquid flux w.r.t. temperature
  real(rkind)                        :: above_soilFracLiq           ! fraction of liquid water layer above soil (canopy or snow) (-)
  ! data structures for flux subroutine arguments (derived types defined in data_types module)
  !      ** intent(in) arguments **       ||       ** intent(inout) arguments **        ||      ** intent(out) arguments **
  type(in_type_vegNrgFlux) :: in_vegNrgFlux;                                            type(out_type_vegNrgFlux) :: out_vegNrgFlux ! vegNrgFlux arguments
  type(in_type_ssdNrgFlux) :: in_ssdNrgFlux; type(io_type_ssdNrgFlux) :: io_ssdNrgFlux; type(out_type_ssdNrgFlux) :: out_ssdNrgFlux ! ssdNrgFlux arguments
  type(in_type_vegLiqFlux) :: in_vegLiqFlux;                                            type(out_type_vegLiqFlux) :: out_vegLiqFlux ! vegLiqFlux arguments
  type(in_type_snowLiqFlx) :: in_snowLiqFlx; type(io_type_snowLiqFlx) :: io_snowLiqFlx; type(out_type_snowLiqFlx) :: out_snowLiqFlx ! snowLiqFlx arguments
  type(in_type_soilLiqFlx) :: in_soilLiqFlx; type(io_type_soilLiqFlx) :: io_soilLiqFlx; type(out_type_soilLiqFlx) :: out_soilLiqFlx ! soilLiqFlx arguments
  type(in_type_groundwatr) :: in_groundwatr; type(io_type_groundwatr) :: io_groundwatr; type(out_type_groundwatr) :: out_groundwatr ! groundwatr arguments
  type(in_type_bigAquifer) :: in_bigAquifer; type(io_type_bigAquifer) :: io_bigAquifer; type(out_type_bigAquifer) :: out_bigAquifer ! bigAquifer arguments
  ! --------------------------------------------------------------
  ! initialize error control
  err=0; message='computFlux/'

  ! *****
  ! * PRELIMINARIES...
  ! ********************
  ! get the necessary variables for the flux computations
  associate(&
    ! model decisions
    ixGroundwater                => model_decisions(iLookDECISIONS%groundwatr)%iDecision            ,& ! intent(in): [i4b]    groundwater parameterization
    ixSpatialGroundwater         => model_decisions(iLookDECISIONS%spatial_gw)%iDecision            ,& ! intent(in): [i4b]    spatial representation of groundwater (local-column or single-basin)
    ! canopy and layer depth
    canopyDepth                  => diag_data%var(iLookDIAG%scalarCanopyDepth)%dat(1)               ,& ! intent(in): [dp   ]  canopy depth (m)
    ! indices of model state variables for the vegetation subdomain
    ixCasNrg                     => indx_data%var(iLookINDEX%ixCasNrg)%dat(1)                       ,& ! intent(in): [i4b]    index of canopy air space energy state variable
    ixVegNrg                     => indx_data%var(iLookINDEX%ixVegNrg)%dat(1)                       ,& ! intent(in): [i4b]    index of canopy energy state variable
    ixVegHyd                     => indx_data%var(iLookINDEX%ixVegHyd)%dat(1)                       ,& ! intent(in): [i4b]    index of canopy hydrology state variable (mass)
    ixTopNrg                     => indx_data%var(iLookINDEX%ixTopNrg)%dat(1)                       ,& ! intent(in): [i4b]    index of upper-most energy state in the snow+soil subdomain
    ixAqWat                      => indx_data%var(iLookINDEX%ixAqWat)%dat(1)                        ,& ! intent(in): [i4b]    index of water storage in the aquifer
    ! indices of model state variables for the snow+soil domain
    ixSnowSoilNrg                => indx_data%var(iLookINDEX%ixSnowSoilNrg)%dat                     ,& ! intent(in): [i4b(:)] indices for energy states in the snow+soil subdomain
    ixSnowSoilHyd                => indx_data%var(iLookINDEX%ixSnowSoilHyd)%dat                     ,& ! intent(in): [i4b(:)] indices for hydrology states in the snow+soil subdomain
    layerType                    => indx_data%var(iLookINDEX%layerType)%dat                         ,& ! intent(in): [i4b(:)] type of layer (iname_soil or iname_snow)
    ! number of state variables of a specific type
    nSnowSoilNrg                 => indx_data%var(iLookINDEX%nSnowSoilNrg )%dat(1)                  ,& ! intent(in): [i4b]    number of energy state variables in the snow+soil domain
    nSnowSoilHyd                 => indx_data%var(iLookINDEX%nSnowSoilHyd )%dat(1)                  ,& ! intent(in): [i4b]    number of hydrology variables in the snow+soil domain
    nSnowOnlyHyd                 => indx_data%var(iLookINDEX%nSnowOnlyHyd )%dat(1)                  ,& ! intent(in): [i4b]    number of hydrology variables in the snow domain
    nSoilOnlyHyd                 => indx_data%var(iLookINDEX%nSoilOnlyHyd )%dat(1)                  ,& ! intent(in): [i4b]    number of hydrology variables in the soil domain
    ! derivatives
    mLayerdTheta_dTk             => deriv_data%var(iLookDERIV%mLayerdTheta_dTk)%dat                 ,&  ! intent(in):  [dp(:)] derivative of volumetric liquid water content w.r.t. temperature
    ! number of flux calls
    numFluxCalls                 => diag_data%var(iLookDIAG%numFluxCalls)%dat(1)                    ,&  ! intent(out): [dp] number of flux calls (-)
    ! net fluxes over the vegetation domain
    scalarCanairNetNrgFlux       => flux_data%var(iLookFLUX%scalarCanairNetNrgFlux)%dat(1)          ,&  ! intent(out): [dp] net energy flux for the canopy air space        (W m-2)
    scalarCanopyNetNrgFlux       => flux_data%var(iLookFLUX%scalarCanopyNetNrgFlux)%dat(1)          ,&  ! intent(out): [dp] net energy flux for the vegetation canopy       (W m-2)
    scalarCanopyNetLiqFlux       => flux_data%var(iLookFLUX%scalarCanopyNetLiqFlux)%dat(1)          ,&  ! intent(out): [dp] net liquid water flux for the vegetation canopy (kg m-2 s-1)
    ! net fluxes over the snow+soil domain
    mLayerNrgFlux                => flux_data%var(iLookFLUX%mLayerNrgFlux)%dat                      ,&  ! intent(out): [dp] net energy flux for each layer within the snow+soil domain (J m-3 s-1)
    mLayerLiqFluxSnow            => flux_data%var(iLookFLUX%mLayerLiqFluxSnow)%dat                  ,&  ! intent(out): [dp] net liquid water flux for each snow layer (s-1)
    mLayerLiqFluxSoil            => flux_data%var(iLookFLUX%mLayerLiqFluxSoil)%dat                  ,&  ! intent(out): [dp] net liquid water flux for each soil layer (s-1)
    ! fluxes for the snow+soil domain
    iLayerLiqFluxSnow            => flux_data%var(iLookFLUX%iLayerLiqFluxSnow)%dat                  ,&  ! intent(out): [dp(0:)] vertical liquid water flux at snow layer interfaces (-)
    iLayerLiqFluxSoil            => flux_data%var(iLookFLUX%iLayerLiqFluxSoil)%dat                  ,&  ! intent(out): [dp(0:)] vertical liquid water flux at soil layer interfaces (-)
    mLayerBaseflow               => flux_data%var(iLookFLUX%mLayerBaseflow)%dat                     ,&  ! intent(out): [dp(:)]  baseflow from each soil layer (m s-1)
    scalarSoilDrainage           => flux_data%var(iLookFLUX%scalarSoilDrainage)%dat(1)              ,&  ! intent(out): [dp]     drainage from the soil profile (m s-1)
    scalarSoilBaseflow           => flux_data%var(iLookFLUX%scalarSoilBaseflow)%dat(1)              ,&  ! intent(out): [dp]     total baseflow from the soil profile (m s-1)
    ! infiltration
    scalarFracLiqVeg             => diag_data%var(iLookDIAG%scalarFracLiqVeg)%dat(1)                ,& ! intent(inout): [dp] fraction of liquid water on vegetation (-)
    mLayerFracLiqSnow            => diag_data%var(iLookDIAG%mLayerFracLiqSnow)%dat                  ,& ! intent(inout): [dp(:)] fraction of liquid water in each snow layer (-)
    ! boundary fluxes in the soil domain
    scalarThroughfallRain        => flux_data%var(iLookFLUX%scalarThroughfallRain)%dat(1)           ,&  ! intent(out): [dp] rain that reaches the ground without ever touching the canopy (kg m-2 s-1)
    scalarCanopyLiqDrainage      => flux_data%var(iLookFLUX%scalarCanopyLiqDrainage)%dat(1)         ,&  ! intent(out): [dp] drainage of liquid water from the vegetation canopy (kg m-2 s-1)
    scalarRainPlusMelt           => flux_data%var(iLookFLUX%scalarRainPlusMelt)%dat(1)              ,&  ! intent(out): [dp] rain plus melt (m s-1)
    scalarSurfaceRunoff          => flux_data%var(iLookFLUX%scalarSurfaceRunoff)%dat(1)             ,&  ! intent(out): [dp] surface runoff (m s-1)
    scalarExfiltration           => flux_data%var(iLookFLUX%scalarExfiltration)%dat(1)              ,&  ! intent(out): [dp] exfiltration from the soil profile (m s-1)
    mLayerColumnOutflow          => flux_data%var(iLookFLUX%mLayerColumnOutflow)%dat                ,&  ! intent(out): [dp(:)] column outflow from each soil layer (m3 s-1)
    ! fluxes for the aquifer
    scalarAquiferTranspire       => flux_data%var(iLookFLUX%scalarAquiferTranspire)%dat(1)          ,&  ! intent(out): [dp] transpiration loss from the aquifer (m s-1
    scalarAquiferRecharge        => flux_data%var(iLookFLUX%scalarAquiferRecharge)%dat(1)           ,&  ! intent(out): [dp] recharge to the aquifer (m s-1)
    scalarAquiferBaseflow        => flux_data%var(iLookFLUX%scalarAquiferBaseflow)%dat(1)           ,&  ! intent(out): [dp] total baseflow from the aquifer (m s-1)
    ! total runoff
    scalarTotalRunoff            => flux_data%var(iLookFLUX%scalarTotalRunoff)%dat(1)               ,&  ! intent(out): [dp] total runoff (m s-1)
    ! derivatives in canopy water w.r.t canopy temperature
    dCanLiq_dTcanopy             => deriv_data%var(iLookDERIV%dCanLiq_dTcanopy            )%dat(1)  ,&  ! intent(out): [dp] derivative of canopy liquid storage w.r.t. temperature
    ! derivatives in canopy liquid fluxes w.r.t. canopy water
    scalarCanopyLiqDeriv         => deriv_data%var(iLookDERIV%scalarCanopyLiqDeriv        )%dat(1)  ,&  ! intent(out): [dp] derivative in (throughfall + drainage) w.r.t. canopy liquid water
    ! derivative in liquid water fluxes at the interface of snow layers w.r.t. volumetric liquid water content in the layer above
    iLayerLiqFluxSnowDeriv       => deriv_data%var(iLookDERIV%iLayerLiqFluxSnowDeriv      )%dat     ,&  ! intent(out): [dp(:)] derivative in vertical liquid water flux at layer interfaces
    ! derivative in baseflow flux w.r.t. aquifer storage
    dBaseflow_dAquifer           => deriv_data%var(iLookDERIV%dBaseflow_dAquifer          )%dat(1)   &  ! intent(out): [dp(:)] derivative in baseflow flux w.r.t. aquifer storage (s-1)
    )  ! end association to data in structures

    numFluxCalls = numFluxCalls+1 ! increment the number of flux calls

    ! modify the groundwater representation for this single-column implementation
    select case(ixSpatialGroundwater)
      case(singleBasin); local_ixGroundwater = noExplicit    ! force no explicit representation of groundwater at the local scale
      case(localColumn); local_ixGroundwater = ixGroundwater ! go with the specified decision
      case default; err=20; message=trim(message)//'unable to identify spatial representation of groundwater'; return
    end select ! end modify the groundwater representation for this single-column implementation

    ! initialize liquid water fluxes throughout the snow and soil domains
    ! NOTE: used in the energy routines, which is called before the hydrology routines
    if (firstFluxCall) then
      if (nSnow>0) iLayerLiqFluxSnow(0:nSnow) = 0._rkind
      iLayerLiqFluxSoil(0:nSoil) = 0._rkind
    end if

    ! *** CALCULATE ENERGY FLUXES OVER VEGETATION ***
    ! identify the need to calculate the energy flux over vegetation
    doVegNrgFlux = (ixCasNrg/=integerMissing .or. ixVegNrg/=integerMissing .or. ixTopNrg/=integerMissing)
    if (doVegNrgFlux) then ! if necessary, calculate the energy fluxes over vegetation
      call subTools(iLookOP%pre,iLookROUTINE%vegNrgFlux)  ! pre-processing for call to vegNrgFlux
      call vegNrgFlux(in_vegNrgFlux,type_data,forc_data,mpar_data,indx_data,prog_data,diag_data,flux_data,bvar_data,model_decisions,out_vegNrgFlux)
      call subTools(iLookOP%post,iLookROUTINE%vegNrgFlux) ! post-processing for call to vegNrgFlux
    end if ! end if calculating the energy fluxes over vegetation

    ! *** CALCULATE ENERGY FLUXES THROUGH THE SNOW-SOIL DOMAIN ***
    if (nSnowSoilNrg>0) then ! if necessary, calculate energy fluxes at layer interfaces through the snow and soil domain
      call subTools(iLookOP%pre,iLookROUTINE%ssdNrgFlux)  ! pre-processing for call to ssdNrgFlux
      call ssdNrgFlux(in_ssdNrgFlux,mpar_data,indx_data,prog_data,diag_data,flux_data,io_ssdNrgFlux,out_ssdNrgFlux)
      call subTools(iLookOP%post,iLookROUTINE%ssdNrgFlux) ! post-processing for call to ssdNrgFlux
    end if  ! end if computing energy fluxes throughout the snow+soil domain

    ! *** CALCULATE THE LIQUID FLUX THROUGH VEGETATION ***
    if (ixVegHyd/=integerMissing) then ! if necessary, calculate liquid water fluxes through vegetation
      call in_vegLiqFlux%initialize(computeVegFlux,scalarCanopyLiqTrial,flux_data)
      call vegLiqFlux(in_vegLiqFlux,mpar_data,diag_data,out_vegLiqFlux)
      call out_vegLiqFlux%finalize(globalPrintFlag,scalarCanopyLiqTrial,flux_data,deriv_data,message,err,cmessage)
    end if  ! end if computing the liquid water fluxes through vegetation

    ! *** CALCULATE THE LIQUID FLUX THROUGH SNOW ***
    if (nSnowOnlyHyd>0) then ! if necessary, compute liquid fluxes through snow
      call initialize_snowLiqFlx
      call snowLiqFlx(in_snowLiqFlx,indx_data,mpar_data,prog_data,diag_data,io_snowLiqFlx,out_snowLiqFlx)
      call finalize_snowLiqFlx
    else
      ! define forcing for the soil domain for the case of no snow layers
      ! NOTE: in case where nSnowOnlyHyd==0 AND snow layers exist, then scalarRainPlusMelt is taken from the previous flux evaluation
      if (nSnow==0) then !no snow layers
        scalarRainPlusMelt = (scalarThroughfallRain + scalarCanopyLiqDrainage)/iden_water &  ! liquid flux from the canopy (m s-1)
                            + drainageMeltPond/iden_water  ! melt of the snow without a layer (m s-1)
        if (ixVegHyd/=integerMissing) then
          ! save canopy derivatives
          above_soilLiqFluxDeriv = scalarCanopyLiqDeriv/iden_water ! derivative in (throughfall + drainage) w.r.t. canopy liquid water
          above_soildLiq_dTk     = dCanLiq_dTcanopy     ! derivative of canopy liquid storage w.r.t. temperature
          above_soilFracLiq      = scalarFracLiqVeg     ! fraction of liquid water in canopy (-)
        else
          above_soilLiqFluxDeriv = 0._rkind
          above_soildLiq_dTk     = 0._rkind
          above_soilFracLiq      = 0._rkind
        end if
      else ! snow layers, take from previous flux calculation
        above_soilLiqFluxDeriv = iLayerLiqFluxSnowDeriv(nSnow) ! derivative in vertical liquid water flux at bottom snow layer interface
        above_soildLiq_dTk     = mLayerdTheta_dTk(nSnow)  ! derivative in volumetric liquid water content in bottom snow layer w.r.t. temperature
        above_soilFracLiq      = mLayerFracLiqSnow(nSnow) ! fraction of liquid water in bottom snow layer (-)
      end if  ! snow layers or not
    end if ! if calculating the liquid flux through snow

    ! *** CALCULATE THE LIQUID FLUX THROUGH SOIL ***
    if (nSoilOnlyHyd>0) then ! if necessary, calculate the liquid flux through soil
      call subTools(iLookOP%pre,iLookROUTINE%soilLiqFlx)  ! pre-processing for call to soilLiqFlx
      call soilLiqFlx(in_soilLiqFlx,mpar_data,indx_data,prog_data,diag_data,flux_data,io_soilLiqFlx,out_soilLiqFlx)
      call subTools(iLookOP%post,iLookROUTINE%soilLiqFlx) ! post-processing for call to soilLiqFlx
    end if  ! end if calculating the liquid flux through soil

    ! *** CALCULATE THE GROUNDWATER FLOW ***
    if (nSoilOnlyHyd>0) then ! check if computing soil hydrology
      if (local_ixGroundwater/=qbaseTopmodel) then ! set baseflow fluxes to zero if the topmodel baseflow routine is not used
        ! diagnostic variables in the data structures
        scalarExfiltration     = 0._rkind  ! exfiltration from the soil profile (m s-1)
        mLayerColumnOutflow(:) = 0._rkind  ! column outflow from each soil layer (m3 s-1)
        ! variables needed for the numerical solution
        mLayerBaseflow(:)      = 0._rkind  ! baseflow from each soil layer (m s-1)
      else ! compute the baseflow flux for topmodel-ish shallow groundwater
        call subTools(iLookOP%pre,iLookROUTINE%groundwatr)  ! pre-processing for call to groundwatr
        call groundwatr(in_groundwatr,attr_data,mpar_data,prog_data,diag_data,flux_data,io_groundwatr,out_groundwatr)
        call subTools(iLookOP%post,iLookROUTINE%groundwatr) ! post-processing for call to groundwatr
      end if  ! computing baseflow flux
      scalarSoilBaseflow = sum(mLayerBaseflow) ! compute total baseflow from the soil zone (needed for mass balance checks)
      ! compute total runoff
      ! (Note: scalarSoilBaseflow is zero if topmodel is not used)
      ! (Note: scalarSoilBaseflow may need to re-envisioned in topmodel formulation if parts of it flow into neighboring soil rather than exfiltrate)
      scalarTotalRunoff  = scalarSurfaceRunoff + scalarSoilDrainage + scalarSoilBaseflow
    end if  ! end if computing soil hydrology


    ! *** CALCULATE FLUXES FOR THE DEEP AQUIFER ***
    if (ixAqWat/=integerMissing) then ! check if computing aquifer fluxes
      if (local_ixGroundwater==bigBucket) then ! compute fluxes for the big bucket
        call subTools(iLookOP%pre,iLookROUTINE%bigAquifer)  ! pre-processing for call to bigAquifer
        call bigAquifer(in_bigAquifer,mpar_data,diag_data,io_bigAquifer,out_bigAquifer)
        call subTools(iLookOP%post,iLookROUTINE%bigAquifer) ! post-processing for call to bigAquifer
      else ! if no aquifer, then fluxes are zero
        scalarAquiferTranspire = 0._rkind  ! transpiration loss from the aquifer (m s-1)
        scalarAquiferRecharge  = 0._rkind  ! recharge to the aquifer (m s-1)
        scalarAquiferBaseflow  = 0._rkind  ! total baseflow from the aquifer (m s-1)
        dBaseflow_dAquifer     = 0._rkind  ! change in baseflow flux w.r.t. aquifer storage (s-1)
      end if ! end check aquifer model decision
    end if  ! if computing aquifer fluxes

    ! *** WRAP UP ***
    ! define model flux vector for the vegetation sub-domain
    if (ixCasNrg/=integerMissing) fluxVec(ixCasNrg) = scalarCanairNetNrgFlux/canopyDepth
    if (ixVegNrg/=integerMissing) fluxVec(ixVegNrg) = scalarCanopyNetNrgFlux/canopyDepth
    if (ixVegHyd/=integerMissing) fluxVec(ixVegHyd) = scalarCanopyNetLiqFlux   ! NOTE: solid fluxes are handled separately
    if (nSnowSoilNrg>0) then ! if necessary, populate the flux vector for energy
      do concurrent (iLayer=1:nLayers,ixSnowSoilNrg(iLayer)/=integerMissing)   ! loop through non-missing energy state variables in the snow+soil domain
        fluxVec( ixSnowSoilNrg(iLayer) ) = mLayerNrgFlux(iLayer)
      end do
    end if
    ! populate the flux vector for hydrology
    ! NOTE: ixVolFracWat  and ixVolFracLiq can also include states in the soil domain, hence enable primary variable switching
    if (nSnowSoilHyd>0) then  ! check if any hydrology states exist
      do iLayer=1,nLayers     ! loop through non-missing energy state variables in the snow+soil domain
        if (ixSnowSoilHyd(iLayer)/=integerMissing) then   ! check if a given hydrology state exists
          select case(layerType(iLayer))
            case(iname_snow); fluxVec(ixSnowSoilHyd(iLayer)) = mLayerLiqFluxSnow(iLayer)
            case(iname_soil); fluxVec(ixSnowSoilHyd(iLayer)) = mLayerLiqFluxSoil(iLayer-nSnow)
            case default; err=20; message=trim(message)//'expect layerType to be either iname_snow or iname_soil'; return
          end select
        end if  ! end if a given hydrology state exists
      end do
    end if  ! end if any hydrology states exist
    ! compute the flux vector for the aquifer
    if (ixAqWat/=integerMissing) fluxVec(ixAqWat) = scalarAquiferTranspire + scalarAquiferRecharge - scalarAquiferBaseflow

    firstFluxCall=.false. ! set the first flux call to false

  end associate  ! end association to variables in the data structures

contains

 subroutine subTools(op,sub)
  implicit none
  integer(i4b),intent(in) :: op  ! index of requested operation: iLookOP%pre for pre-processing or iLookOP%post for post-processing
  integer(i4b),intent(in) :: sub ! index of subroutine in computFlux (e.g., iLookROUTINE%vegNrgFlux for vegNrgFlux routine)

  associate(&
    ! domain boundary conditions
    upperBoundTemp               => forc_data%var(iLookFORCE%airtemp)                               ,& ! intent(in): [dp]     temperature of the upper boundary of the snow and soil domains (K)
    scalarRainfall               => flux_data%var(iLookFLUX%scalarRainfall)%dat(1)                  ,& ! intent(in): [dp]     rainfall rate (kg m-2 s-1)
    ! canopy and layer depth
    canopyDepth                  => diag_data%var(iLookDIAG%scalarCanopyDepth)%dat(1)               ,& ! intent(in): [dp   ]  canopy depth (m)
    mLayerDepth                  => prog_data%var(iLookPROG%mLayerDepth)%dat                        ,& ! intent(in): [dp(:)]  depth of each layer in the snow-soil sub-domain (m)
    ! derivatives
    dPsiLiq_dPsi0                => deriv_data%var(iLookDERIV%dPsiLiq_dPsi0   )%dat                 ,&  ! intent(in):  [dp(:)] derivative in liquid water matric pot w.r.t. the total water matric pot (-)
    dPsiLiq_dTemp                => deriv_data%var(iLookDERIV%dPsiLiq_dTemp   )%dat                 ,&  ! intent(in):  [dp(:)] derivative in the liquid water matric potential w.r.t. temperature
    mLayerdTheta_dTk             => deriv_data%var(iLookDERIV%mLayerdTheta_dTk)%dat                 ,&  ! intent(in):  [dp(:)] derivative of volumetric liquid water content w.r.t. temperature
    dTheta_dTkCanopy             => deriv_data%var(iLookDERIV%dTheta_dTkCanopy)%dat(1)              ,&  ! intent(in):  [dp]    derivative of volumetric liquid water content w.r.t. temperature
    ! net fluxes over the vegetation domain
    scalarCanairNetNrgFlux       => flux_data%var(iLookFLUX%scalarCanairNetNrgFlux)%dat(1)          ,&  ! intent(out): [dp] net energy flux for the canopy air space        (W m-2)
    scalarCanopyNetNrgFlux       => flux_data%var(iLookFLUX%scalarCanopyNetNrgFlux)%dat(1)          ,&  ! intent(out): [dp] net energy flux for the vegetation canopy       (W m-2)
    scalarGroundNetNrgFlux       => flux_data%var(iLookFLUX%scalarGroundNetNrgFlux)%dat(1)          ,&  ! intent(out): [dp] net energy flux for the ground surface          (W m-2)
    scalarCanopyNetLiqFlux       => flux_data%var(iLookFLUX%scalarCanopyNetLiqFlux)%dat(1)          ,&  ! intent(out): [dp] net liquid water flux for the vegetation canopy (kg m-2 s-1)
    ! net fluxes over the snow+soil domain
    mLayerNrgFlux                => flux_data%var(iLookFLUX%mLayerNrgFlux)%dat                      ,&  ! intent(out): [dp] net energy flux for each layer within the snow+soil domain (J m-3 s-1)
    mLayerLiqFluxSnow            => flux_data%var(iLookFLUX%mLayerLiqFluxSnow)%dat                  ,&  ! intent(out): [dp] net liquid water flux for each snow layer (s-1)
    mLayerLiqFluxSoil            => flux_data%var(iLookFLUX%mLayerLiqFluxSoil)%dat                  ,&  ! intent(out): [dp] net liquid water flux for each soil layer (s-1)
    ! evaporative fluxes
    scalarCanopyTranspiration    => flux_data%var(iLookFLUX%scalarCanopyTranspiration)%dat(1)       ,&  ! intent(out): [dp]    canopy transpiration (kg m-2 s-1)
    scalarCanopyEvaporation      => flux_data%var(iLookFLUX%scalarCanopyEvaporation)%dat(1)         ,&  ! intent(out): [dp]    canopy evaporation/condensation (kg m-2 s-1)
    scalarGroundEvaporation      => flux_data%var(iLookFLUX%scalarGroundEvaporation)%dat(1)         ,&  ! intent(out): [dp]    ground evaporation/condensation -- below canopy or non-vegetated (kg m-2 s-1)
    mLayerTranspire              => flux_data%var(iLookFLUX%mLayerTranspire)%dat                    ,&  ! intent(out): [dp(:)] transpiration loss from each soil layer (m s-1)
    ! fluxes for the snow+soil domain
    iLayerNrgFlux                => flux_data%var(iLookFLUX%iLayerNrgFlux)%dat                      ,&  ! intent(out): [dp(0:)] vertical energy flux at the interface of snow and soil layers
    iLayerLiqFluxSnow            => flux_data%var(iLookFLUX%iLayerLiqFluxSnow)%dat                  ,&  ! intent(out): [dp(0:)] vertical liquid water flux at snow layer interfaces (-)
    iLayerLiqFluxSoil            => flux_data%var(iLookFLUX%iLayerLiqFluxSoil)%dat                  ,&  ! intent(out): [dp(0:)] vertical liquid water flux at soil layer interfaces (-)
    mLayerHydCond                => flux_data%var(iLookFLUX%mLayerHydCond)%dat                      ,&  ! intent(out): [dp(:)]  hydraulic conductivity in each soil layer (m s-1)
    mLayerBaseflow               => flux_data%var(iLookFLUX%mLayerBaseflow)%dat                     ,&  ! intent(out): [dp(:)]  baseflow from each soil layer (m s-1)
    scalarSnowDrainage           => flux_data%var(iLookFLUX%scalarSnowDrainage)%dat(1)              ,&  ! intent(out): [dp]     drainage from the snow profile (m s-1)
    scalarSoilDrainage           => flux_data%var(iLookFLUX%scalarSoilDrainage)%dat(1)              ,&  ! intent(out): [dp]     drainage from the soil profile (m s-1)
    ! infiltration
    scalarInfilArea              => diag_data%var(iLookDIAG%scalarInfilArea   )%dat(1)              ,&  ! intent(out): [dp] fraction of unfrozen area where water can infiltrate (-)
    scalarFrozenArea             => diag_data%var(iLookDIAG%scalarFrozenArea  )%dat(1)              ,&  ! intent(out): [dp] fraction of area that is considered impermeable due to soil ice (-)
    scalarSoilControl            => diag_data%var(iLookDIAG%scalarSoilControl )%dat(1)              ,&  ! intent(out): [dp] soil control on infiltration, zero or one
    scalarMaxInfilRate           => flux_data%var(iLookFLUX%scalarMaxInfilRate)%dat(1)              ,&  ! intent(out): [dp] maximum infiltration rate (m s-1)
    scalarInfiltration           => flux_data%var(iLookFLUX%scalarInfiltration)%dat(1)              ,&  ! intent(out): [dp] infiltration of water into the soil profile (m s-1)
    mLayerFracLiqSnow            => diag_data%var(iLookDIAG%mLayerFracLiqSnow)%dat                  ,& ! intent(inout): [dp(:)] fraction of liquid water in each snow layer (-)
    ! boundary fluxes in the soil domain
    scalarThroughfallRain        => flux_data%var(iLookFLUX%scalarThroughfallRain)%dat(1)           ,&  ! intent(out): [dp] rain that reaches the ground without ever touching the canopy (kg m-2 s-1)
    scalarCanopyLiqDrainage      => flux_data%var(iLookFLUX%scalarCanopyLiqDrainage)%dat(1)         ,&  ! intent(out): [dp] drainage of liquid water from the vegetation canopy (kg m-2 s-1)
    scalarRainPlusMelt           => flux_data%var(iLookFLUX%scalarRainPlusMelt)%dat(1)              ,&  ! intent(out): [dp] rain plus melt (m s-1)
    scalarSurfaceRunoff          => flux_data%var(iLookFLUX%scalarSurfaceRunoff)%dat(1)             ,&  ! intent(out): [dp] surface runoff (m s-1)
    ! fluxes for the aquifer
    scalarAquiferTranspire       => flux_data%var(iLookFLUX%scalarAquiferTranspire)%dat(1)          ,&  ! intent(out): [dp] transpiration loss from the aquifer (m s-1
    scalarAquiferRecharge        => flux_data%var(iLookFLUX%scalarAquiferRecharge)%dat(1)           ,&  ! intent(out): [dp] recharge to the aquifer (m s-1)
    scalarAquiferBaseflow        => flux_data%var(iLookFLUX%scalarAquiferBaseflow)%dat(1)           ,&  ! intent(out): [dp] total baseflow from the aquifer (m s-1)
    ! total runoff
    scalarTotalRunoff            => flux_data%var(iLookFLUX%scalarTotalRunoff)%dat(1)               ,&  ! intent(out): [dp] total runoff (m s-1)
    ! derivatives in net vegetation energy fluxes w.r.t. relevant state variables
    dCanairNetFlux_dCanairTemp   => deriv_data%var(iLookDERIV%dCanairNetFlux_dCanairTemp  )%dat(1)  ,&  ! intent(out): [dp] derivative in net canopy air space flux w.r.t. canopy air temperature
    dCanairNetFlux_dCanopyTemp   => deriv_data%var(iLookDERIV%dCanairNetFlux_dCanopyTemp  )%dat(1)  ,&  ! intent(out): [dp] derivative in net canopy air space flux w.r.t. canopy temperature
    dCanairNetFlux_dGroundTemp   => deriv_data%var(iLookDERIV%dCanairNetFlux_dGroundTemp  )%dat(1)  ,&  ! intent(out): [dp] derivative in net canopy air space flux w.r.t. ground temperature
    dCanopyNetFlux_dCanairTemp   => deriv_data%var(iLookDERIV%dCanopyNetFlux_dCanairTemp  )%dat(1)  ,&  ! intent(out): [dp] derivative in net canopy flux w.r.t. canopy air temperature
    dCanopyNetFlux_dCanopyTemp   => deriv_data%var(iLookDERIV%dCanopyNetFlux_dCanopyTemp  )%dat(1)  ,&  ! intent(out): [dp] derivative in net canopy flux w.r.t. canopy temperature
    dCanopyNetFlux_dGroundTemp   => deriv_data%var(iLookDERIV%dCanopyNetFlux_dGroundTemp  )%dat(1)  ,&  ! intent(out): [dp] derivative in net canopy flux w.r.t. ground temperature
    dCanopyNetFlux_dCanWat       => deriv_data%var(iLookDERIV%dCanopyNetFlux_dCanWat      )%dat(1)  ,&  ! intent(out): [dp] derivative in net canopy fluxes w.r.t. canopy total water content
    dGroundNetFlux_dCanairTemp   => deriv_data%var(iLookDERIV%dGroundNetFlux_dCanairTemp  )%dat(1)  ,&  ! intent(out): [dp] derivative in net ground flux w.r.t. canopy air temperature
    dGroundNetFlux_dCanopyTemp   => deriv_data%var(iLookDERIV%dGroundNetFlux_dCanopyTemp  )%dat(1)  ,&  ! intent(out): [dp] derivative in net ground flux w.r.t. canopy temperature
    dGroundNetFlux_dGroundTemp   => deriv_data%var(iLookDERIV%dGroundNetFlux_dGroundTemp  )%dat(1)  ,&  ! intent(out): [dp] derivative in net ground flux w.r.t. ground temperature
    dGroundNetFlux_dCanWat       => deriv_data%var(iLookDERIV%dGroundNetFlux_dCanWat      )%dat(1)  ,&  ! intent(out): [dp] derivative in net ground fluxes w.r.t. canopy total water content
    ! derivatives in evaporative fluxes w.r.t. relevant state variables
    dCanopyEvaporation_dTCanair  => deriv_data%var(iLookDERIV%dCanopyEvaporation_dTCanair )%dat(1)  ,&  ! intent(out): [dp] derivative in canopy evaporation w.r.t. canopy air temperature
    dCanopyEvaporation_dTCanopy  => deriv_data%var(iLookDERIV%dCanopyEvaporation_dTCanopy )%dat(1)  ,&  ! intent(out): [dp] derivative in canopy evaporation w.r.t. canopy temperature
    dCanopyEvaporation_dTGround  => deriv_data%var(iLookDERIV%dCanopyEvaporation_dTGround )%dat(1)  ,&  ! intent(out): [dp] derivative in canopy evaporation w.r.t. ground temperature
    dCanopyEvaporation_dCanWat   => deriv_data%var(iLookDERIV%dCanopyEvaporation_dCanWat  )%dat(1)  ,&  ! intent(out): [dp] derivative in canopy evaporation w.r.t. canopy total water content
    dGroundEvaporation_dTCanair  => deriv_data%var(iLookDERIV%dGroundEvaporation_dTCanair )%dat(1)  ,&  ! intent(out): [dp] derivative in ground evaporation w.r.t. canopy air temperature
    dGroundEvaporation_dTCanopy  => deriv_data%var(iLookDERIV%dGroundEvaporation_dTCanopy )%dat(1)  ,&  ! intent(out): [dp] derivative in ground evaporation w.r.t. canopy temperature
    dGroundEvaporation_dTGround  => deriv_data%var(iLookDERIV%dGroundEvaporation_dTGround )%dat(1)  ,&  ! intent(out): [dp] derivative in ground evaporation w.r.t. ground temperature
    dGroundEvaporation_dCanWat   => deriv_data%var(iLookDERIV%dGroundEvaporation_dCanWat  )%dat(1)  ,&  ! intent(out): [dp] derivative in ground evaporation w.r.t. canopy total water content
    ! derivatives in transpiration
    dCanopyTrans_dTCanair        => deriv_data%var(iLookDERIV%dCanopyTrans_dTCanair       )%dat(1)  ,&  ! intent(out): [dp] derivative in canopy transpiration w.r.t. canopy air temperature (kg m-2 s-1 K-1)
    dCanopyTrans_dTCanopy        => deriv_data%var(iLookDERIV%dCanopyTrans_dTCanopy       )%dat(1)  ,&  ! intent(out): [dp] derivative in canopy transpiration w.r.t. canopy temperature (kg m-2 s-1 K-1)
    dCanopyTrans_dTGround        => deriv_data%var(iLookDERIV%dCanopyTrans_dTGround       )%dat(1)  ,&  ! intent(out): [dp] derivative in canopy transpiration w.r.t. ground temperature (kg m-2 s-1 K-1)
    dCanopyTrans_dCanWat         => deriv_data%var(iLookDERIV%dCanopyTrans_dCanWat        )%dat(1)  ,&  ! intent(out): [dp] derivative in canopy transpiration w.r.t. canopy total water content (s-1)
    ! derivatives in canopy water w.r.t canopy temperature
    dCanLiq_dTcanopy             => deriv_data%var(iLookDERIV%dCanLiq_dTcanopy            )%dat(1)  ,&  ! intent(out): [dp] derivative of canopy liquid storage w.r.t. temperature
    ! derivatives in canopy liquid fluxes w.r.t. canopy water
    scalarCanopyLiqDeriv         => deriv_data%var(iLookDERIV%scalarCanopyLiqDeriv        )%dat(1)  ,&  ! intent(out): [dp] derivative in (throughfall + drainage) w.r.t. canopy liquid water
    scalarThroughfallRainDeriv   => deriv_data%var(iLookDERIV%scalarThroughfallRainDeriv  )%dat(1)  ,&  ! intent(out): [dp] derivative in throughfall w.r.t. canopy liquid water
    scalarCanopyLiqDrainageDeriv => deriv_data%var(iLookDERIV%scalarCanopyLiqDrainageDeriv)%dat(1)  ,&  ! intent(out): [dp] derivative in canopy drainage w.r.t. canopy liquid water
    ! derivatives in energy fluxes at the interface of snow+soil layers w.r.t. temperature in layers above and below
    dNrgFlux_dTempAbove          => deriv_data%var(iLookDERIV%dNrgFlux_dTempAbove         )%dat     ,&  ! intent(out): [dp(:)] derivatives in the flux w.r.t. temperature in the layer above
    dNrgFlux_dTempBelow          => deriv_data%var(iLookDERIV%dNrgFlux_dTempBelow         )%dat     ,&  ! intent(out): [dp(:)] derivatives in the flux w.r.t. temperature in the layer below
    dThermalC_dWatAbove          => deriv_data%var(iLookDERIV%dThermalC_dWatAbove         )%dat     ,&  ! intent(in):  [dp(:)] derivative in the thermal conductivity w.r.t. water state in the layer above
    dThermalC_dWatBelow          => deriv_data%var(iLookDERIV%dThermalC_dWatBelow         )%dat     ,&  ! intent(in):  [dp(:)] derivative in the thermal conductivity w.r.t. water state in the layer above
    ! derivatives in energy fluxes at the interface of snow+soil layers w.r.t. water state in layers above and below
    dNrgFlux_dWatAbove           => deriv_data%var(iLookDERIV%dNrgFlux_dWatAbove          )%dat     ,&  ! intent(out):  [dp(:)] derivatives in the flux w.r.t. water state in the layer above
    dNrgFlux_dWatBelow           => deriv_data%var(iLookDERIV%dNrgFlux_dWatBelow          )%dat     ,&  ! intent(out): [dp(:)] derivatives in the flux w.r.t. water state in the layer below
    dThermalC_dTempAbove         => deriv_data%var(iLookDERIV%dThermalC_dTempAbove        )%dat     ,&  ! intent(in):  [dp(:)] derivative in the thermal conductivity w.r.t. energy state in the layer above
    dThermalC_dTempBelow         => deriv_data%var(iLookDERIV%dThermalC_dTempBelow        )%dat     ,&  ! intent(in):  [dp(:)] derivative in the thermal conductivity w.r.t. energy state in the layer above
    ! derivative in liquid water fluxes at the interface of snow layers w.r.t. volumetric liquid water content in the layer above
    iLayerLiqFluxSnowDeriv       => deriv_data%var(iLookDERIV%iLayerLiqFluxSnowDeriv      )%dat     ,&  ! intent(out): [dp(:)] derivative in vertical liquid water flux at layer interfaces
    ! derivative in liquid water fluxes for the soil domain w.r.t hydrology state variables
    dq_dHydStateAbove            => deriv_data%var(iLookDERIV%dq_dHydStateAbove           )%dat     ,&  ! intent(out): [dp(:)] change in flux at layer interfaces w.r.t. states in the layer above
    dq_dHydStateBelow            => deriv_data%var(iLookDERIV%dq_dHydStateBelow           )%dat     ,&  ! intent(out): [dp(:)] change in flux at layer interfaces w.r.t. states in the layer below
    dq_dHydStateLayerSurfVec     => deriv_data%var(iLookDERIV%dq_dHydStateLayerSurfVec    )%dat     ,&  ! intent(out): [dp(:)] change in the flux in soil surface interface w.r.t. state variables in layers
    mLayerdTheta_dPsi            => deriv_data%var(iLookDERIV%mLayerdTheta_dPsi           )%dat     ,&  ! intent(out): [dp(:)] derivative in the soil water characteristic w.r.t. psi
    mLayerdPsi_dTheta            => deriv_data%var(iLookDERIV%mLayerdPsi_dTheta           )%dat     ,&  ! intent(out): [dp(:)] derivative in the soil water characteristic w.r.t. theta
    ! derivative in baseflow flux w.r.t. aquifer storage
    dBaseflow_dAquifer           => deriv_data%var(iLookDERIV%dBaseflow_dAquifer          )%dat(1)  ,&  ! intent(out): [dp(:)] derivative in baseflow flux w.r.t. aquifer storage (s-1)
    ! derivative in liquid water fluxes for the soil domain w.r.t energy state variables
    dq_dNrgStateAbove            => deriv_data%var(iLookDERIV%dq_dNrgStateAbove           )%dat     ,&  ! intent(out): [dp(:)] change in flux at layer interfaces w.r.t. states in the layer above
    dq_dNrgStateBelow            => deriv_data%var(iLookDERIV%dq_dNrgStateBelow           )%dat     ,&  ! intent(out): [dp(:)] change in flux at layer interfaces w.r.t. states in the layer below
    dq_dNrgStateLayerSurfVec     => deriv_data%var(iLookDERIV%dq_dNrgStateLayerSurfVec    )%dat     ,&  ! intent(out): [dp(:)] change in the flux in soil surface interface w.r.t. state variables in layers
    ! derivatives in soil transpiration w.r.t. canopy state variables
    mLayerdTrans_dTCanair        => deriv_data%var(iLookDERIV%mLayerdTrans_dTCanair       )%dat     ,&  ! intent(out): derivatives in the soil layer transpiration flux w.r.t. canopy air temperature
    mLayerdTrans_dTCanopy        => deriv_data%var(iLookDERIV%mLayerdTrans_dTCanopy       )%dat     ,&  ! intent(out): derivatives in the soil layer transpiration flux w.r.t. canopy temperature
    mLayerdTrans_dTGround        => deriv_data%var(iLookDERIV%mLayerdTrans_dTGround       )%dat     ,&  ! intent(out): derivatives in the soil layer transpiration flux w.r.t. ground temperature
    mLayerdTrans_dCanWat         => deriv_data%var(iLookDERIV%mLayerdTrans_dCanWat        )%dat     ,&  ! intent(out): derivatives in the soil layer transpiration flux w.r.t. canopy total water
    ! derivatives in aquifer transpiration w.r.t. canopy state variables
    dAquiferTrans_dTCanair       => deriv_data%var(iLookDERIV%dAquiferTrans_dTCanair      )%dat(1)  ,&  ! intent(out): derivatives in the aquifer transpiration flux w.r.t. canopy air temperature
    dAquiferTrans_dTCanopy       => deriv_data%var(iLookDERIV%dAquiferTrans_dTCanopy      )%dat(1)  ,&  ! intent(out): derivatives in the aquifer transpiration flux w.r.t. canopy temperature
    dAquiferTrans_dTGround       => deriv_data%var(iLookDERIV%dAquiferTrans_dTGround      )%dat(1)  ,&  ! intent(out): derivatives in the aquifer transpiration flux w.r.t. ground temperature
    dAquiferTrans_dCanWat        => deriv_data%var(iLookDERIV%dAquiferTrans_dCanWat       )%dat(1)   &  ! intent(out): derivatives in the aquifer transpiration flux w.r.t. canopy total water
  )  ! end associating to data in structures
  ! initialize error control
  err=0; message='subTools/'

  ! Validate operation argument
  if ((op/=iLookOP%pre).and.(op/=iLookOP%post)) then
   err=20
   cmessage="Error in subTools: invalid op argument requested."
   message=trim(message)//trim(cmessage); return
  end if

  select case(sub)
   case(iLookROUTINE%vegNrgFlux) ! vegNrgFlux
    if (op==iLookOP%pre) then ! pre-processing
     dCanLiq_dTcanopy = dTheta_dTkCanopy*iden_water*canopyDepth     ! derivative in canopy liquid storage w.r.t. canopy temperature (kg m-2 K-1)
     ! intent(in) arguments
     in_vegNrgFlux % firstSubStep=firstSubStep                      ! intent(in): flag to indicate if we are processing the first sub-step
     in_vegNrgFlux % firstFluxCall=firstFluxCall                    ! intent(in): flag to indicate if we are processing the first flux call
     in_vegNrgFlux % computeVegFlux=computeVegFlux                  ! intent(in): flag to indicate if we need to compute fluxes over vegetation
     in_vegNrgFlux % checkLWBalance=checkLWBalance                  ! intent(in): flag to check longwave balance
     in_vegNrgFlux % upperBoundTemp=upperBoundTemp                  ! intent(in): temperature of the upper boundary (K) --> NOTE: use air temperature
     in_vegNrgFlux % scalarCanairTempTrial=scalarCanairTempTrial    ! intent(in): trial value of the canopy air space temperature (K)
     in_vegNrgFlux % scalarCanopyTempTrial=scalarCanopyTempTrial    ! intent(in): trial value of canopy temperature (K)
     in_vegNrgFlux % mLayerTempTrial_1=mLayerTempTrial(1)           ! intent(in): trial value of ground temperature (K)
     in_vegNrgFlux % scalarCanopyIceTrial=scalarCanopyIceTrial      ! intent(in): trial value of mass of ice on the vegetation canopy (kg m-2)
     in_vegNrgFlux % scalarCanopyLiqTrial=scalarCanopyLiqTrial      ! intent(in): trial value of mass of liquid water on the vegetation canopy (kg m-2)
     in_vegNrgFlux % dCanLiq_dTcanopy=dCanLiq_dTcanopy              ! intent(in): derivative in canopy liquid storage w.r.t. canopy temperature (kg m-2 K-1)        
    else ! post-processing
     ! intent(out) arguments
     scalarCanopyTranspiration  =out_vegNrgFlux % scalarCanopyTranspiration   ! intent(out): canopy transpiration (kg m-2 s-1)
     scalarCanopyEvaporation    =out_vegNrgFlux % scalarCanopyEvaporation     ! intent(out): canopy evaporation/condensation (kg m-2 s-1)
     scalarGroundEvaporation    =out_vegNrgFlux % scalarGroundEvaporation     ! intent(out): ground evaporation/condensation -- below canopy or non-vegetated (kg m-2 s-1)
     scalarCanairNetNrgFlux     =out_vegNrgFlux % scalarCanairNetNrgFlux      ! intent(out): net energy flux for the canopy air space (W m-2)
     scalarCanopyNetNrgFlux     =out_vegNrgFlux % scalarCanopyNetNrgFlux      ! intent(out): net energy flux for the vegetation canopy (W m-2)
     scalarGroundNetNrgFlux     =out_vegNrgFlux % scalarGroundNetNrgFlux      ! intent(out): net energy flux for the ground surface (W m-2)
     dCanairNetFlux_dCanairTemp =out_vegNrgFlux % dCanairNetFlux_dCanairTemp  ! intent(out): derivative in net canopy air space flux w.r.t. canopy air temperature (W m-2 K-1)
     dCanairNetFlux_dCanopyTemp =out_vegNrgFlux % dCanairNetFlux_dCanopyTemp  ! intent(out): derivative in net canopy air space flux w.r.t. canopy temperature (W m-2 K-1)
     dCanairNetFlux_dGroundTemp =out_vegNrgFlux % dCanairNetFlux_dGroundTemp  ! intent(out): derivative in net canopy air space flux w.r.t. ground temperature (W m-2 K-1)
     dCanopyNetFlux_dCanairTemp =out_vegNrgFlux % dCanopyNetFlux_dCanairTemp  ! intent(out): derivative in net canopy flux w.r.t. canopy air temperature (W m-2 K-1)
     dCanopyNetFlux_dCanopyTemp =out_vegNrgFlux % dCanopyNetFlux_dCanopyTemp  ! intent(out): derivative in net canopy flux w.r.t. canopy temperature (W m-2 K-1)
     dCanopyNetFlux_dGroundTemp =out_vegNrgFlux % dCanopyNetFlux_dGroundTemp  ! intent(out): derivative in net canopy flux w.r.t. ground temperature (W m-2 K-1)
     dGroundNetFlux_dCanairTemp =out_vegNrgFlux % dGroundNetFlux_dCanairTemp  ! intent(out): derivative in net ground flux w.r.t. canopy air temperature (W m-2 K-1)
     dGroundNetFlux_dCanopyTemp =out_vegNrgFlux % dGroundNetFlux_dCanopyTemp  ! intent(out): derivative in net ground flux w.r.t. canopy temperature (W m-2 K-1)
     dGroundNetFlux_dGroundTemp =out_vegNrgFlux % dGroundNetFlux_dGroundTemp  ! intent(out): derivative in net ground flux w.r.t. ground temperature (W m-2 K-1)
     dCanopyEvaporation_dCanWat =out_vegNrgFlux % dCanopyEvaporation_dCanWat  ! intent(out): derivative in canopy evaporation w.r.t. canopy total water content (s-1)
     dCanopyEvaporation_dTCanair=out_vegNrgFlux % dCanopyEvaporation_dTCanair ! intent(out): derivative in canopy evaporation w.r.t. canopy air temperature (kg m-2 s-1 K-1)
     dCanopyEvaporation_dTCanopy=out_vegNrgFlux % dCanopyEvaporation_dTCanopy ! intent(out): derivative in canopy evaporation w.r.t. canopy temperature (kg m-2 s-1 K-1)
     dCanopyEvaporation_dTGround=out_vegNrgFlux % dCanopyEvaporation_dTGround ! intent(out): derivative in canopy evaporation w.r.t. ground temperature (kg m-2 s-1 K-1)
     dGroundEvaporation_dCanWat =out_vegNrgFlux % dGroundEvaporation_dCanWat  ! intent(out): derivative in ground evaporation w.r.t. canopy total water content (s-1)
     dGroundEvaporation_dTCanair=out_vegNrgFlux % dGroundEvaporation_dTCanair ! intent(out): derivative in ground evaporation w.r.t. canopy air temperature (kg m-2 s-1 K-1)
     dGroundEvaporation_dTCanopy=out_vegNrgFlux % dGroundEvaporation_dTCanopy ! intent(out): derivative in ground evaporation w.r.t. canopy temperature (kg m-2 s-1 K-1)
     dGroundEvaporation_dTGround=out_vegNrgFlux % dGroundEvaporation_dTGround ! intent(out): derivative in ground evaporation w.r.t. ground temperature (kg m-2 s-1 K-1)
     dCanopyTrans_dCanWat       =out_vegNrgFlux % dCanopyTrans_dCanWat  ! intent(out): derivative in canopy transpiration w.r.t. canopy total water content (s-1)
     dCanopyTrans_dTCanair      =out_vegNrgFlux % dCanopyTrans_dTCanair ! intent(out): derivative in canopy transpiration w.r.t. canopy air temperature (kg m-2 s-1 K-1)
     dCanopyTrans_dTCanopy      =out_vegNrgFlux % dCanopyTrans_dTCanopy ! intent(out): derivative in canopy transpiration w.r.t. canopy temperature (kg m-2 s-1 K-1)
     dCanopyTrans_dTGround      =out_vegNrgFlux % dCanopyTrans_dTGround ! intent(out): derivative in canopy transpiration w.r.t. ground temperature (kg m-2 s-1 K-1)
     dCanopyNetFlux_dCanWat     =out_vegNrgFlux % dCanopyNetFlux_dCanWat! intent(out): derivative in net canopy fluxes w.r.t. canopy total water content (J kg-1 s-1)
     dGroundNetFlux_dCanWat     =out_vegNrgFlux % dGroundNetFlux_dCanWat! intent(out): derivative in net ground fluxes w.r.t. canopy total water content (J kg-1 s-1)
     err                        =out_vegNrgFlux % err                   ! intent(out): error code
     cmessage                   =out_vegNrgFlux % cmessage              ! intent(out): error message
     ! additional post-processing
     if (err/=0) then; message=trim(message)//trim(cmessage); return; end if  ! check for errors
      ! check fluxes
     if (globalPrintFlag) then
       print*, '**'
       write(*,'(a,1x,10(f30.20))') 'canopyDepth           = ',  canopyDepth
       write(*,'(a,1x,10(f30.20))') 'mLayerDepth(1:2)      = ',  mLayerDepth(1:2)
       write(*,'(a,1x,10(f30.20))') 'scalarCanairTempTrial = ',  scalarCanairTempTrial   ! trial value of the canopy air space temperature (K)
       write(*,'(a,1x,10(f30.20))') 'scalarCanopyTempTrial = ',  scalarCanopyTempTrial   ! trial value of canopy temperature (K)
       write(*,'(a,1x,10(f30.20))') 'mLayerTempTrial(1:2)  = ',  mLayerTempTrial(1:2)    ! trial value of ground temperature (K)
       write(*,'(a,1x,10(f30.20))') 'scalarCanairNetNrgFlux = ', scalarCanairNetNrgFlux
       write(*,'(a,1x,10(f30.20))') 'scalarCanopyNetNrgFlux = ', scalarCanopyNetNrgFlux
       write(*,'(a,1x,10(f30.20))') 'scalarGroundNetNrgFlux = ', scalarGroundNetNrgFlux
       write(*,'(a,1x,10(f30.20))') 'dGroundNetFlux_dGroundTemp = ', dGroundNetFlux_dGroundTemp
     end if ! end if checking fluxes
    end if 
   case(iLookROUTINE%ssdNrgFlux) ! ssdNrgFlux
    if (op==iLookOP%pre) then ! pre-processing
     ! intent(in) arguments
     in_ssdNrgFlux % scalarSolution=scalarSolution .and. .not.firstFluxCall ! intent(in): flag to denote if implementing the scalar solution
     in_ssdNrgFlux % scalarGroundNetNrgFlux=scalarGroundNetNrgFlux          ! intent(in): net energy flux for the ground surface (W m-2)
     in_ssdNrgFlux % iLayerLiqFluxSnow=iLayerLiqFluxSnow                    ! intent(in): liquid flux at the interface of each snow layer (m s-1)
     in_ssdNrgFlux % iLayerLiqFluxSoil=iLayerLiqFluxSoil                    ! intent(in): liquid flux at the interface of each soil layer (m s-1)
     in_ssdNrgFlux % mLayerTempTrial=mLayerTempTrial                        ! intent(in): temperature in each layer at the current iteration (m)
     in_ssdNrgFlux % dThermalC_dWatAbove=dThermalC_dWatAbove                ! intent(in): derivative in the thermal conductivity w.r.t. water state in the layer above
     in_ssdNrgFlux % dThermalC_dWatBelow=dThermalC_dWatBelow                ! intent(in): derivative in the thermal conductivity w.r.t. water state in the layer above
     in_ssdNrgFlux % dThermalC_dTempAbove=dThermalC_dTempAbove              ! intent(in): derivative in the thermal conductivity w.r.t. energy state in the layer above
     in_ssdNrgFlux % dThermalC_dTempBelow=dThermalC_dTempBelow              ! intent(in): derivative in the thermal conductivity w.r.t. energy state in the layer above
     ! intent(inout) arguments
     io_ssdNrgFlux % dGroundNetFlux_dGroundTemp=dGroundNetFlux_dGroundTemp  ! intent(inout): derivative in net ground flux w.r.t. ground temperature (W m-2 K-1)
    else ! post-processing
     ! intent(inout) arguments
     dGroundNetFlux_dGroundTemp=io_ssdNrgFlux % dGroundNetFlux_dGroundTemp  ! intent(inout): derivative in net ground flux w.r.t. ground temperature (W m-2 K-1)
     ! intent(out) arguments
     iLayerNrgFlux      =out_ssdNrgFlux % iLayerNrgFlux                     ! intent(out): energy flux at the layer interfaces (W m-2)
     dNrgFlux_dTempAbove=out_ssdNrgFlux % dNrgFlux_dTempAbove               ! intent(out): derivatives in the flux w.r.t. temperature in the layer above (J m-2 s-1 K-1)
     dNrgFlux_dTempBelow=out_ssdNrgFlux % dNrgFlux_dTempBelow               ! intent(out): derivatives in the flux w.r.t. temperature in the layer below (J m-2 s-1 K-1)
     dNrgFlux_dWatAbove =out_ssdNrgFlux % dNrgFlux_dWatAbove                ! intent(out): derivatives in the flux w.r.t. water state in the layer above (J m-2 s-1 K-1)
     dNrgFlux_dWatBelow =out_ssdNrgFlux % dNrgFlux_dWatBelow                ! intent(out): derivatives in the flux w.r.t. water state in the layer below (J m-2 s-1 K-1)
     err                =out_ssdNrgFlux % err                               ! intent(out): error code
     cmessage           =out_ssdNrgFlux % cmessage                          ! intent(out): error message
     ! error control
     if (err/=0) then; message=trim(message)//trim(cmessage); return; end if
     ! calculate net energy fluxes for each snow and soil layer (J m-3 s-1)
     do iLayer=1,nLayers
       mLayerNrgFlux(iLayer) = -(iLayerNrgFlux(iLayer) - iLayerNrgFlux(iLayer-1))/mLayerDepth(iLayer)
       if (globalPrintFlag) then
         if (iLayer < 10) write(*,'(a,1x,i4,1x,10(f25.15,1x))') 'iLayer, iLayerNrgFlux(iLayer-1:iLayer), mLayerNrgFlux(iLayer)   = ', iLayer, iLayerNrgFlux(iLayer-1:iLayer), mLayerNrgFlux(iLayer)
       end if
     end do
    end if
   case(iLookROUTINE%vegLiqFlux) ! vegLiqFlux
    if (op==iLookOP%pre) then ! pre-processing
     ! intent(in) arguments
     in_vegLiqFlux % computeVegFlux      =computeVegFlux        ! intent(in): flag to denote if computing energy flux over vegetation
     in_vegLiqFlux % scalarCanopyLiqTrial=scalarCanopyLiqTrial  ! intent(in): trial mass of liquid water on the vegetation canopy at the current iteration (kg m-2)
     in_vegLiqFlux % scalarRainfall      =scalarRainfall        ! intent(in): rainfall rate (kg m-2 s-1)
    else ! post-processing
     ! intent(out) arguments
     scalarThroughfallRain       =out_vegLiqFlux % scalarThroughfallRain       ! intent(out): rain that reaches the ground without ever touching the canopy (kg m-2 s-1)
     scalarCanopyLiqDrainage     =out_vegLiqFlux % scalarCanopyLiqDrainage     ! intent(out): drainage of liquid water from the vegetation canopy (kg m-2 s-1)
     scalarThroughfallRainDeriv  =out_vegLiqFlux % scalarThroughfallRainDeriv  ! intent(out): derivative in throughfall w.r.t. canopy liquid water (s-1)
     scalarCanopyLiqDrainageDeriv=out_vegLiqFlux % scalarCanopyLiqDrainageDeriv! intent(out): derivative in canopy drainage w.r.t. canopy liquid water (s-1)
     err                         =out_vegLiqFlux % err                         ! intent(out): error code
     cmessage                    =out_vegLiqFlux % cmessage                    ! intent(out): error control
     ! error control
     if (err/=0) then; message=trim(message)//trim(cmessage); return; end if
     ! calculate the net liquid water flux for the vegetation canopy
     scalarCanopyNetLiqFlux = scalarRainfall + scalarCanopyEvaporation - scalarThroughfallRain - scalarCanopyLiqDrainage
     ! calculate the total derivative in the downward liquid flux
     scalarCanopyLiqDeriv   = scalarThroughfallRainDeriv + scalarCanopyLiqDrainageDeriv
     ! test
     if (globalPrintFlag) then
       print*, '**'
       print*, 'scalarRainfall          = ', scalarRainfall
       print*, 'scalarThroughfallRain   = ', scalarThroughfallRain
       print*, 'scalarCanopyEvaporation = ', scalarCanopyEvaporation
       print*, 'scalarCanopyLiqDrainage = ', scalarCanopyLiqDrainage
       print*, 'scalarCanopyNetLiqFlux  = ', scalarCanopyNetLiqFlux
       print*, 'scalarCanopyLiqTrial    = ', scalarCanopyLiqTrial
     end if
    end if
   case(iLookROUTINE%snowLiqFlx) ! snowLiqFlx
    if (op==iLookOP%pre) then ! pre-processing
     ! intent(in) arguments
     in_snowLiqFlx % nSnow                  =nSnow                          ! intent(in): number of snow layers
     in_snowLiqFlx % firstFluxCall          =firstFluxCall                  ! intent(in): the first flux call (compute variables that are constant over the iterations)
     in_snowLiqFlx % scalarSolution         =(scalarSolution .and. .not.firstFluxCall) ! intent(in): flag to indicate the scalar solution
     in_snowLiqFlx % scalarThroughfallRain  =scalarThroughfallRain          ! intent(in): rain that reaches the snow surface without ever touching vegetation (kg m-2 s-1)
     in_snowLiqFlx % scalarCanopyLiqDrainage=scalarCanopyLiqDrainage        ! intent(in): liquid drainage from the vegetation canopy (kg m-2 s-1)
     in_snowLiqFlx % mLayerVolFracLiqTrial  =mLayerVolFracLiqTrial(1:nSnow) ! intent(in): trial value of volumetric fraction of liquid water at the current iteration (-)
     ! intent(inout) arguments
     io_snowLiqFlx % iLayerLiqFluxSnow      =iLayerLiqFluxSnow       ! intent(inout): vertical liquid water flux at layer interfaces (m s-1)
     io_snowLiqFlx % iLayerLiqFluxSnowDeriv =iLayerLiqFluxSnowDeriv  ! intent(inout): derivative in vertical liquid water flux at layer interfaces (m s-1)
    else ! post-processing
     ! intent(inout) arguments
     iLayerLiqFluxSnow     =io_snowLiqFlx % iLayerLiqFluxSnow        ! intent(inout): vertical liquid water flux at layer interfaces (m s-1)
     iLayerLiqFluxSnowDeriv=io_snowLiqFlx % iLayerLiqFluxSnowDeriv   ! intent(inout): derivative in vertical liquid water flux at layer interfaces (m s-1)
     ! intent(out) arguments
     err     =out_snowLiqFlx % err                                   ! intent(out):   error code
     cmessage=out_snowLiqFlx % cmessage                              ! intent(out):   error message
     ! error control
     if (err/=0) then; message=trim(message)//trim(cmessage); return; end if
     ! define forcing for the soil domain
     scalarRainPlusMelt = iLayerLiqFluxSnow(nSnow)          ! drainage from the base of the snowpack
     ! calculate net liquid water fluxes for each snow layer (s-1)
     do iLayer=1,nSnow
       mLayerLiqFluxSnow(iLayer) = -(iLayerLiqFluxSnow(iLayer) - iLayerLiqFluxSnow(iLayer-1))/mLayerDepth(iLayer)
     end do
     ! compute drainage from the soil zone (needed for mass balance checks)
     scalarSnowDrainage = iLayerLiqFluxSnow(nSnow)
     ! save bottom layer of snow derivatives
     above_soilLiqFluxDeriv = iLayerLiqFluxSnowDeriv(nSnow) ! derivative in vertical liquid water flux at bottom snow layer interface
     above_soildLiq_dTk     = mLayerdTheta_dTk(nSnow)       ! derivative in volumetric liquid water content in bottom snow layer w.r.t. temperature
     above_soilFracLiq      = mLayerFracLiqSnow(nSnow)      ! fraction of liquid water in bottom snow layer (-)
    end if
   case(iLookROUTINE%soilLiqFlx) ! soilLiqFlx
    if (op==iLookOP%pre) then ! pre-processing
     ! intent(in) arguments
     in_soilLiqFlx % nSoil         =nSoil                                         ! intent(in): number of soil layers
     in_soilLiqFlx % firstSplitOper=firstSplitOper                                ! intent(in): flag indicating first flux call in a splitting operation
     in_soilLiqFlx % scalarSolution=(scalarSolution .and. .not.firstFluxCall)     ! intent(in): flag to indicate the scalar solution
     in_soilLiqFlx % deriv_desired =.true.                                        ! intent(in): flag indicating if derivatives are desired
     in_soilLiqFlx % mLayerTempTrial=mLayerTempTrial(nSnow+1:nLayers)             ! intent(in): trial temperature at the current iteration (K)
     in_soilLiqFlx % mLayerMatricHeadTrial   =mLayerMatricHeadTrial(1:nSoil)      ! intent(in): matric potential (m)
     in_soilLiqFlx % mLayerMatricHeadLiqTrial=mLayerMatricHeadLiqTrial(1:nSoil)   ! intent(in): liquid water matric potential (m)
     in_soilLiqFlx % mLayerVolFracLiqTrial=mLayerVolFracLiqTrial(nSnow+1:nLayers) ! intent(in): volumetric fraction of liquid water (-)
     in_soilLiqFlx % mLayerVolFracIceTrial=mLayerVolFracIceTrial(nSnow+1:nLayers) ! intent(in): volumetric fraction of ice (-)
     in_soilLiqFlx % mLayerdTheta_dTk=mLayerdTheta_dTk(nSnow+1:nLayers)           ! intent(in): derivative in volumetric liquid water content w.r.t. temperature (K-1)
     in_soilLiqFlx % dPsiLiq_dTemp=dPsiLiq_dTemp(1:nSoil)                         ! intent(in): derivative in liquid water matric potential w.r.t. temperature (m K-1)
     in_soilLiqFlx % dCanopyTrans_dCanWat  =dCanopyTrans_dCanWat     ! intent(in): derivative in canopy transpiration w.r.t. canopy total water content (s-1)
     in_soilLiqFlx % dCanopyTrans_dTCanair =dCanopyTrans_dTCanair    ! intent(in): derivative in canopy transpiration w.r.t. canopy air temperature (kg m-2 s-1 K-1)
     in_soilLiqFlx % dCanopyTrans_dTCanopy =dCanopyTrans_dTCanopy    ! intent(in): derivative in canopy transpiration w.r.t. canopy temperature (kg m-2 s-1 K-1)
     in_soilLiqFlx % dCanopyTrans_dTGround =dCanopyTrans_dTGround    ! intent(in): derivative in canopy transpiration w.r.t. ground temperature (kg m-2 s-1 K-1)
     in_soilLiqFlx % above_soilLiqFluxDeriv=above_soilLiqFluxDeriv   ! intent(in): derivative in layer above soil (canopy or snow) liquid flux w.r.t. liquid water
     in_soilLiqFlx % above_soildLiq_dTk    =above_soildLiq_dTk       ! intent(in): derivative of layer above soil (canopy or snow) liquid flux w.r.t. temperature
     in_soilLiqFlx % above_soilFracLiq     =above_soilFracLiq        ! intent(in): fraction of liquid water layer above soil (canopy or snow) (-)
     in_soilLiqFlx % scalarCanopyTranspiration=scalarCanopyTranspiration          ! intent(in): canopy transpiration (kg m-2 s-1)
     in_soilLiqFlx % scalarGroundEvaporation  =scalarGroundEvaporation            ! intent(in): ground evaporation (kg m-2 s-1)
     in_soilLiqFlx % scalarRainPlusMelt       =scalarRainPlusMelt                 ! intent(in): rain plus melt (m s-1)
     ! intent(inout) arguments
     io_soilLiqFlx % scalarMaxInfilRate      =scalarMaxInfilRate       ! intent(inout): maximum infiltration rate (m s-1)
     io_soilLiqFlx % scalarInfilArea         =scalarInfilArea          ! intent(inout): fraction of unfrozen area where water can infiltrate (-)
     io_soilLiqFlx % scalarFrozenArea        =scalarFrozenArea         ! intent(inout): fraction of area that is considered impermeable due to soil ice (-)
     io_soilLiqFlx % scalarSurfaceRunoff     =scalarSurfaceRunoff      ! intent(inout): surface runoff (m s-1)
     io_soilLiqFlx % mLayerdTheta_dPsi       =mLayerdTheta_dPsi        ! intent(inout): derivative in the soil water characteristic w.r.t. psi (m-1)
     io_soilLiqFlx % mLayerdPsi_dTheta       =mLayerdPsi_dTheta        ! intent(inout): derivative in the soil water characteristic w.r.t. theta (m)
     io_soilLiqFlx % dHydCond_dMatric        =dHydCond_dMatric         ! intent(inout): derivative in hydraulic conductivity w.r.t matric head (s-1)
     io_soilLiqFlx % scalarInfiltration      =scalarInfiltration       ! intent(inout): surface infiltration rate (m s-1) -- controls on infiltration only computed for iter==1
     io_soilLiqFlx % iLayerLiqFluxSoil       =iLayerLiqFluxSoil        ! intent(inout): liquid fluxes at layer interfaces (m s-1)
     io_soilLiqFlx % mLayerTranspire         =mLayerTranspire          ! intent(inout): transpiration loss from each soil layer (m s-1)
     io_soilLiqFlx % mLayerHydCond           =mLayerHydCond            ! intent(inout): hydraulic conductivity in each layer (m s-1)
     io_soilLiqFlx % dq_dHydStateAbove       =dq_dHydStateAbove        ! intent(inout): derivatives in the flux w.r.t. matric head in the layer above (s-1)
     io_soilLiqFlx % dq_dHydStateBelow       =dq_dHydStateBelow        ! intent(inout): derivatives in the flux w.r.t. matric head in the layer below (s-1)
     io_soilLiqFlx % dq_dHydStateLayerSurfVec=dq_dHydStateLayerSurfVec ! intent(inout): derivative in surface infiltration w.r.t. hydrology state in above soil snow or canopy and every soil layer  (m s-1 or s-1)
     io_soilLiqFlx % dq_dNrgStateAbove       =dq_dNrgStateAbove        ! intent(inout): derivatives in the flux w.r.t. temperature in the layer above (m s-1 K-1)
     io_soilLiqFlx % dq_dNrgStateBelow       =dq_dNrgStateBelow        ! intent(inout): derivatives in the flux w.r.t. temperature in the layer below (m s-1 K-1)
     io_soilLiqFlx % dq_dNrgStateLayerSurfVec=dq_dNrgStateLayerSurfVec ! intent(inout): derivative in surface infiltration w.r.t. energy state in above soil snow or canopy and every soil layer (m s-1 K-1)
     io_soilLiqFlx % mLayerdTrans_dTCanair   =mLayerdTrans_dTCanair    ! intent(inout): derivatives in the soil layer transpiration flux w.r.t. canopy air temperature
     io_soilLiqFlx % mLayerdTrans_dTCanopy   =mLayerdTrans_dTCanopy    ! intent(inout): derivatives in the soil layer transpiration flux w.r.t. canopy temperature
     io_soilLiqFlx % mLayerdTrans_dTGround   =mLayerdTrans_dTGround    ! intent(inout): derivatives in the soil layer transpiration flux w.r.t. ground temperature
     io_soilLiqFlx % mLayerdTrans_dCanWat    =mLayerdTrans_dCanWat     ! intent(inout): derivatives in the soil layer transpiration flux w.r.t. canopy total water 
    else ! post-processing
     ! intent(inout) arguments
     scalarMaxInfilRate      =io_soilLiqFlx % scalarMaxInfilRate       ! intent(inout): maximum infiltration rate (m s-1)
     scalarInfilArea         =io_soilLiqFlx % scalarInfilArea          ! intent(inout): fraction of unfrozen area where water can infiltrate (-)
     scalarFrozenArea        =io_soilLiqFlx % scalarFrozenArea         ! intent(inout): fraction of area that is considered impermeable due to soil ice (-)
     scalarSurfaceRunoff     =io_soilLiqFlx % scalarSurfaceRunoff      ! intent(inout): surface runoff (m s-1)
     mLayerdTheta_dPsi       =io_soilLiqFlx % mLayerdTheta_dPsi        ! intent(inout): derivative in the soil water characteristic w.r.t. psi (m-1)
     mLayerdPsi_dTheta       =io_soilLiqFlx % mLayerdPsi_dTheta        ! intent(inout): derivative in the soil water characteristic w.r.t. theta (m)
     dHydCond_dMatric        =io_soilLiqFlx % dHydCond_dMatric         ! intent(inout): derivative in hydraulic conductivity w.r.t matric head (s-1)
     scalarInfiltration      =io_soilLiqFlx % scalarInfiltration       ! intent(inout): surface infiltration rate (m s-1) -- controls on infiltration only computed for iter==1
     iLayerLiqFluxSoil       =io_soilLiqFlx % iLayerLiqFluxSoil        ! intent(inout): liquid fluxes at layer interfaces (m s-1)
     mLayerTranspire         =io_soilLiqFlx % mLayerTranspire          ! intent(inout): transpiration loss from each soil layer (m s-1)
     mLayerHydCond           =io_soilLiqFlx % mLayerHydCond            ! intent(inout): hydraulic conductivity in each layer (m s-1)
     dq_dHydStateAbove       =io_soilLiqFlx % dq_dHydStateAbove        ! intent(inout): derivatives in the flux w.r.t. matric head in the layer above (s-1)
     dq_dHydStateBelow       =io_soilLiqFlx % dq_dHydStateBelow        ! intent(inout): derivatives in the flux w.r.t. matric head in the layer below (s-1)
     dq_dHydStateLayerSurfVec=io_soilLiqFlx % dq_dHydStateLayerSurfVec ! intent(inout): derivative in surface infiltration w.r.t. hydrology state in above soil snow or canopy and every soil layer  (m s-1 or s-1)
     dq_dNrgStateAbove       =io_soilLiqFlx % dq_dNrgStateAbove        ! intent(inout): derivatives in the flux w.r.t. temperature in the layer above (m s-1 K-1)
     dq_dNrgStateBelow       =io_soilLiqFlx % dq_dNrgStateBelow        ! intent(inout): derivatives in the flux w.r.t. temperature in the layer below (m s-1 K-1)
     dq_dNrgStateLayerSurfVec=io_soilLiqFlx % dq_dNrgStateLayerSurfVec ! intent(inout): derivative in surface infiltration w.r.t. energy state in above soil snow or canopy and every soil layer (m s-1 K-1)
     mLayerdTrans_dTCanair   =io_soilLiqFlx % mLayerdTrans_dTCanair    ! intent(inout): derivatives in the soil layer transpiration flux w.r.t. canopy air temperature
     mLayerdTrans_dTCanopy   =io_soilLiqFlx % mLayerdTrans_dTCanopy    ! intent(inout): derivatives in the soil layer transpiration flux w.r.t. canopy temperature
     mLayerdTrans_dTGround   =io_soilLiqFlx % mLayerdTrans_dTGround    ! intent(inout): derivatives in the soil layer transpiration flux w.r.t. ground temperature
     mLayerdTrans_dCanWat    =io_soilLiqFlx % mLayerdTrans_dCanWat     ! intent(inout): derivatives in the soil layer transpiration flux w.r.t. canopy total water 
     ! intent(out) arguments
     err                     =out_soilLiqFlx % err                     ! intent(out):   error code
     cmessage                =out_soilLiqFlx % cmessage                ! intent(out):   error message
     ! error control
     if (err/=0) then; message=trim(message)//trim(cmessage); return; end if

     ! calculate net liquid water fluxes for each soil layer (s-1)
     do iLayer=1,nSoil
       mLayerLiqFluxSoil(iLayer) = -(iLayerLiqFluxSoil(iLayer) - iLayerLiqFluxSoil(iLayer-1))/mLayerDepth(iLayer+nSnow)
     end do

     ! calculate the soil control on infiltration
     if (nSnow==0) then
       ! * case of infiltration into soil
       if (scalarMaxInfilRate > scalarRainPlusMelt) then  ! infiltration is not rate-limited
         scalarSoilControl = (1._rkind - scalarFrozenArea)*scalarInfilArea
       else
         scalarSoilControl = 0._rkind  ! (scalarRainPlusMelt exceeds maximum infiltration rate
       end if
     else
       ! * case of infiltration into snow
       scalarSoilControl = 1._rkind
     end if

     ! compute drainage from the soil zone (needed for mass balance checks and in aquifer recharge)
     scalarSoilDrainage = iLayerLiqFluxSoil(nSoil)

     ! expand derivatives to the total water matric potential
     ! NOTE: arrays are offset because computing derivatives in interface fluxes, at the top and bottom of the layer respectively
     if (globalPrintFlag) print*, 'dPsiLiq_dPsi0(1:nSoil) = ', dPsiLiq_dPsi0(1:nSoil)
     dq_dHydStateAbove(1:nSoil)   = dq_dHydStateAbove(1:nSoil)  *dPsiLiq_dPsi0(1:nSoil)
     dq_dHydStateBelow(0:nSoil-1) = dq_dHydStateBelow(0:nSoil-1)*dPsiLiq_dPsi0(1:nSoil)
     dq_dHydStateLayerSurfVec(1:nSoil) = dq_dHydStateLayerSurfVec(1:nSoil)*dPsiLiq_dPsi0(1:nSoil)
    end if
   case(iLookROUTINE%groundwatr) ! groundwatr
    if (op==iLookOP%pre) then ! pre-processing
     ! check the derivative matrix is sized appropriately
     if (size(dBaseflow_dMatric,1)/=nSoil .or. size(dBaseflow_dMatric,2)/=nSoil) then
       message=trim(message)//'expect dBaseflow_dMatric to be nSoil x nSoil'
       err=20; return
     end if
     ! intent(in) arguments
     in_groundwatr % nSnow                    = nSnow                                  ! intent(in):    number of snow layers
     in_groundwatr % nSoil                    = nSoil                                  ! intent(in):    number of soil layers
     in_groundwatr % nLayers                  = nLayers                                ! intent(in):    total number of layers
     in_groundwatr % firstFluxCall            = firstFluxCall                          ! intent(in):    logical flag to compute index of the lowest saturated layer
     in_groundwatr % mLayerdTheta_dPsi        = mLayerdTheta_dPsi                      ! intent(in):    derivative in the soil water characteristic w.r.t. matric head in each layer (m-1)
     in_groundwatr % mLayerMatricHeadLiqTrial = mLayerMatricHeadLiqTrial               ! intent(in):    liquid water matric potential (m)
     in_groundwatr % mLayerVolFracLiqTrial    = mLayerVolFracLiqTrial(nSnow+1:nLayers) ! intent(in):    volumetric fraction of liquid water (-)
     in_groundwatr % mLayerVolFracIceTrial    = mLayerVolFracIceTrial(nSnow+1:nLayers) ! intent(in):    volumetric fraction of ice (-)
     ! intent(inout) arguments
     io_groundwatr % ixSaturation = ixSaturation ! intent(inout): index of lowest saturated layer (NOTE: only computed on the first iteration)
    else ! post-processing
     ! intent(inout) arguments
     ixSaturation = io_groundwatr % ixSaturation ! intent(inout): index of lowest saturated layer (NOTE: only computed on the first iteration)
     ! intent(out) arguments
     mLayerBaseflow    = out_groundwatr % mLayerBaseflow                               ! intent(out):   baseflow from each soil layer (m s-1)
     dBaseflow_dMatric = out_groundwatr % dBaseflow_dMatric                            ! intent(out):   derivative in baseflow w.r.t. matric head (s-1)
     err               = out_groundwatr % err                                          ! intent(out):   error code
     cmessage          = out_groundwatr % cmessage                                     ! intent(out):   error message
     ! error control
     if (err/=0) then; message=trim(message)//trim(cmessage); return; end if
    end if
   case(iLookROUTINE%bigAquifer) ! bigAquifer
    if (op==iLookOP%pre) then ! pre-processing
     ! intent(in) arguments
     in_bigAquifer % scalarAquiferStorageTrial = scalarAquiferStorageTrial ! intent(in): trial value of aquifer storage (m)
     in_bigAquifer % scalarCanopyTranspiration = scalarCanopyTranspiration ! intent(in): canopy transpiration (kg m-2 s-1)
     in_bigAquifer % scalarSoilDrainage        = scalarSoilDrainage        ! intent(in): soil drainage (m s-1)
     in_bigAquifer % dCanopyTrans_dCanWat      = dCanopyTrans_dCanWat      ! intent(in): derivative in canopy transpiration w.r.t. canopy total water content (s-1)
     in_bigAquifer % dCanopyTrans_dTCanair     = dCanopyTrans_dTCanair     ! intent(in): derivative in canopy transpiration w.r.t. canopy air temperature (kg m-2 s-1 K-1)
     in_bigAquifer % dCanopyTrans_dTCanopy     = dCanopyTrans_dTCanopy     ! intent(in): derivative in canopy transpiration w.r.t. canopy temperature (kg m-2 s-1 K-1)
     in_bigAquifer % dCanopyTrans_dTGround     = dCanopyTrans_dTGround     ! intent(in): derivative in canopy transpiration w.r.t. ground temperature (kg m-2 s-1 K-1)
     ! intent(inout) arguments
     io_bigAquifer % dAquiferTrans_dTCanair = dAquiferTrans_dTCanair       ! intent(inout): derivatives in the aquifer transpiration flux w.r.t. canopy air temperature
     io_bigAquifer % dAquiferTrans_dTCanopy = dAquiferTrans_dTCanopy       ! intent(inout): derivatives in the aquifer transpiration flux w.r.t. canopy temperature
     io_bigAquifer % dAquiferTrans_dTGround = dAquiferTrans_dTGround       ! intent(inout): derivatives in the aquifer transpiration flux w.r.t. ground temperature
     io_bigAquifer % dAquiferTrans_dCanWat  = dAquiferTrans_dCanWat        ! intent(inout): derivatives in the aquifer transpiration flux w.r.t. canopy total water
    else ! post-processing
     ! intent(inout) arguments
     dAquiferTrans_dTCanair = io_bigAquifer % dAquiferTrans_dTCanair       ! intent(inout): derivatives in the aquifer transpiration flux w.r.t. canopy air temperature
     dAquiferTrans_dTCanopy = io_bigAquifer % dAquiferTrans_dTCanopy       ! intent(inout): derivatives in the aquifer transpiration flux w.r.t. canopy temperature
     dAquiferTrans_dTGround = io_bigAquifer % dAquiferTrans_dTGround       ! intent(inout): derivatives in the aquifer transpiration flux w.r.t. ground temperature
     dAquiferTrans_dCanWat  = io_bigAquifer % dAquiferTrans_dCanWat        ! intent(inout): derivatives in the aquifer transpiration flux w.r.t. canopy total water
     ! intent(out) arguments
     scalarAquiferTranspire = out_bigAquifer % scalarAquiferTranspire      ! intent(out):   transpiration loss from the aquifer (m s-1)
     scalarAquiferRecharge  = out_bigAquifer % scalarAquiferRecharge       ! intent(out):   recharge to the aquifer (m s-1)
     scalarAquiferBaseflow  = out_bigAquifer % scalarAquiferBaseflow       ! intent(out):   total baseflow from the aquifer (m s-1)
     dBaseflow_dAquifer     = out_bigAquifer % dBaseflow_dAquifer          ! intent(out):   change in baseflow flux w.r.t. aquifer storage (s-1)
     err                    = out_bigAquifer % err                         ! intent(out):   error code
     cmessage               = out_bigAquifer % cmessage                    ! intent(out):   error message
     ! error control
     if (err/=0) then; message=trim(message)//trim(cmessage); return; end if
     ! compute total runoff (overwrite previously calculated value before considering aquifer).
     !   (Note:  SoilDrainage goes into aquifer, not runoff)
     scalarTotalRunoff  = scalarSurfaceRunoff + scalarAquiferBaseflow     
    end if
   case default ! Error control for sub argument (must be s subroutine index that is included in the above case blocks)
    err=20
    cmessage="Error in subTools: invalid sub argument requested."
    message=trim(message)//trim(cmessage); return     
  end select

  end associate ! end associate block

 end subroutine subTools

 subroutine initialize_snowLiqFlx
  call in_snowLiqFlx%initialize(nSnow,firstFluxCall,scalarSolution,mLayerVolFracLiqTrial,flux_data)
  call io_snowLiqFlx%initialize(flux_data,deriv_data)
 end subroutine initialize_snowLiqFlx

 subroutine finalize_snowLiqFlx
  call io_snowLiqFlx%finalize(flux_data,deriv_data)
  call out_snowLiqFlx%finalize(err,cmessage) 
  ! error control
  if (err/=0) then; message=trim(message)//trim(cmessage); return; end if
  associate(&
   scalarRainPlusMelt     => flux_data%var(iLookFLUX%scalarRainPlusMelt)%dat(1),   & ! intent(out): [dp] rain plus melt (m s-1)
   mLayerLiqFluxSnow      => flux_data%var(iLookFLUX%mLayerLiqFluxSnow)%dat,       & ! intent(out): [dp] net liquid water flux for each snow layer (s-1)
   iLayerLiqFluxSnow      => flux_data%var(iLookFLUX%iLayerLiqFluxSnow)%dat,       & ! intent(out): [dp(0:)] vertical liquid water flux at snow layer interfaces (-)
   mLayerDepth            => prog_data%var(iLookPROG%mLayerDepth)%dat,             & ! intent(in): [dp(:)]  depth of each layer in the snow-soil sub-domain (m)
   scalarSnowDrainage     => flux_data%var(iLookFLUX%scalarSnowDrainage)%dat(1),   & ! intent(out): [dp]     drainage from the snow profile (m s-1)
   iLayerLiqFluxSnowDeriv => deriv_data%var(iLookDERIV%iLayerLiqFluxSnowDeriv)%dat,& ! intent(out): [dp(:)] derivative in vertical liquid water flux at layer interfaces
   mLayerdTheta_dTk       => deriv_data%var(iLookDERIV%mLayerdTheta_dTk)%dat,      & ! intent(in):  [dp(:)] derivative of volumetric liquid water content w.r.t. temperature
   mLayerFracLiqSnow      => diag_data%var(iLookDIAG%mLayerFracLiqSnow)%dat)         ! intent(inout): [dp(:)] fraction of liquid water in each snow layer (-)
   ! define forcing for the soil domain
   scalarRainPlusMelt = iLayerLiqFluxSnow(nSnow)          ! drainage from the base of the snowpack
   ! calculate net liquid water fluxes for each snow layer (s-1)
   do iLayer=1,nSnow
     mLayerLiqFluxSnow(iLayer) = -(iLayerLiqFluxSnow(iLayer) - iLayerLiqFluxSnow(iLayer-1))/mLayerDepth(iLayer)
   end do
   ! compute drainage from the soil zone (needed for mass balance checks)
   scalarSnowDrainage = iLayerLiqFluxSnow(nSnow)
   ! save bottom layer of snow derivatives
   above_soilLiqFluxDeriv = iLayerLiqFluxSnowDeriv(nSnow) ! derivative in vertical liquid water flux at bottom snow layer interface
   above_soildLiq_dTk     = mLayerdTheta_dTk(nSnow)       ! derivative in volumetric liquid water content in bottom snow layer w.r.t. temperature
   above_soilFracLiq      = mLayerFracLiqSnow(nSnow)      ! fraction of liquid water in bottom snow layer (-)
  end associate
 end subroutine finalize_snowLiqFlx
end subroutine computFlux

! **********************************************************************************************************
! public subroutine soilCmpres: compute soil compressibility (-) and its derivative w.r.t matric head (m-1)
! **********************************************************************************************************
subroutine soilCmpres(&
                      ! input:
                      dt,                                 & ! intent(in):  length of the time step (seconds)
                      ixRichards,                         & ! intent(in):  choice of option for Richards' equation
                      ixBeg,ixEnd,                        & ! intent(in):  start and end indices defining desired layers
                      mLayerMatricHead,                   & ! intent(in):  matric head at the start of the time step (m)
                      mLayerMatricHeadTrial,              & ! intent(in):  trial value of matric head (m)
                      mLayerVolFracLiqTrial,              & ! intent(in):  trial value for the volumetric liquid water content in each soil layer (-)
                      mLayerVolFracIceTrial,              & ! intent(in):  trial value for the volumetric ice content in each soil layer (-)
                      specificStorage,                    & ! intent(in):  specific storage coefficient (m-1)
                      theta_sat,                          & ! intent(in):  soil porosity (-)
                      ! output:
                      compress,                           & ! intent(out): compressibility of the soil matrix (-), per second
                      dCompress_dPsi,                     & ! intent(out): derivative in compressibility w.r.t. matric head (m-1)
                      err,message)                          ! intent(out): error code and error message
  implicit none
  ! input:
  real(rkind),intent(in)         :: dt                        !  length of the time step (seconds)
  integer(i4b),intent(in)        :: ixRichards                ! choice of option for Richards' equation
  integer(i4b),intent(in)        :: ixBeg,ixEnd               ! start and end indices defining desired layers
  real(rkind),intent(in)         :: mLayerMatricHead(:)       ! matric head at the start of the time step (m)
  real(rkind),intent(in)         :: mLayerMatricHeadTrial(:)  ! trial value for matric head (m)
  real(rkind),intent(in)         :: mLayerVolFracLiqTrial(:)  ! trial value for volumetric fraction of liquid water (-)
  real(rkind),intent(in)         :: mLayerVolFracIceTrial(:)  ! trial value for volumetric fraction of ice (-)
  real(rkind),intent(in)         :: specificStorage           ! specific storage coefficient (m-1)
  real(rkind),intent(in)         :: theta_sat(:)              ! soil porosity (-)
  ! output:
  real(rkind),intent(inout)      :: compress(:)               ! soil compressibility (-)
  real(rkind),intent(inout)      :: dCompress_dPsi(:)         ! derivative in soil compressibility w.r.t. matric head (m-1)
  integer(i4b),intent(out)       :: err                       ! error code
  character(*),intent(out)       :: message                   ! error message
  ! local variables
  integer(i4b)                   :: iLayer                    ! index of soil layer
  ! --------------------------------------------------------------
  ! initialize error control
  err=0; message='soilCmpres/'
  ! (only compute for the mixed form of Richards' equation)
  if (ixRichards==mixdform) then
    do iLayer=1,size(mLayerMatricHead)
      if (iLayer>=ixBeg .and. iLayer<=ixEnd) then
      ! compute the derivative for the compressibility term (m-1), no volume expansion for total water
      dCompress_dPsi(iLayer) = specificStorage*(mLayerVolFracLiqTrial(iLayer) + mLayerVolFracIceTrial(iLayer))/theta_sat(iLayer)
      ! compute the compressibility term (-) per second
      compress(iLayer)       = (mLayerMatricHeadTrial(iLayer) - mLayerMatricHead(iLayer))*dCompress_dPsi(iLayer)/dt
      end if
    end do
  else
    compress(:)       = 0._rkind
    dCompress_dPsi(:) = 0._rkind
  end if
end subroutine soilCmpres

! **********************************************************************************************************
! public subroutine soilCmpres: compute soil compressibility (-) and its derivative w.r.t matric head (m-1)
! **********************************************************************************************************
subroutine soilCmpresPrime(&
                          ! input:
                          ixRichards,                         & ! intent(in):  choice of option for Richards' equation
                          ixBeg,ixEnd,                        & ! intent(in):  start and end indices defining desired layers
                          mLayerMatricHeadPrime,              & ! intent(in):  matric head at the start of the time step (m)
                          mLayerVolFracLiqTrial,              & ! intent(in):  trial value for the volumetric liquid water content in each soil layer (-)
                          mLayerVolFracIceTrial,              & ! intent(in):  trial value for the volumetric ice content in each soil layer (-)
                          specificStorage,                    & ! intent(in):  specific storage coefficient (m-1)
                          theta_sat,                          & ! intent(in):  soil porosity (-)
                          ! output:
                          compress,                           & ! intent(out): compressibility of the soil matrix (-)
                          dCompress_dPsi,                     & ! intent(out): derivative in compressibility w.r.t. matric head (m-1)
                          err,message)                          ! intent(out): error code and error message
  implicit none
  ! input:
  integer(i4b),intent(in)           :: ixRichards               ! choice of option for Richards' equation
  integer(i4b),intent(in)           :: ixBeg,ixEnd              ! start and end indices defining desired layers
  real(rkind),intent(in)            :: mLayerMatricHeadPrime(:) ! matric head at the start of the time step (m)
  real(rkind),intent(in)            :: mLayerVolFracLiqTrial(:) ! trial value for volumetric fraction of liquid water (-)
  real(rkind),intent(in)            :: mLayerVolFracIceTrial(:) ! trial value for volumetric fraction of ice (-)
  real(rkind),intent(in)            :: specificStorage          ! specific storage coefficient (m-1)
  real(rkind),intent(in)            :: theta_sat(:)             ! soil porosity (-)
  ! output:
  real(rkind),intent(inout)         :: compress(:)              ! soil compressibility (-)
  real(rkind),intent(inout)         :: dCompress_dPsi(:)        ! derivative in soil compressibility w.r.t. matric head (m-1)
  integer(i4b),intent(out)          :: err                      ! error code
  character(*),intent(out)          :: message                  ! error message
  ! local variables
  integer(i4b)                      :: iLayer                   ! index of soil layer
  ! --------------------------------------------------------------
  ! initialize error control
  err=0; message='soilCmpresPrime/'
  ! (only compute for the mixed form of Richards' equation)
  if (ixRichards==mixdform) then
    do iLayer=1,size(mLayerMatricHeadPrime)
      if (iLayer>=ixBeg .and. iLayer<=ixEnd) then
          ! compute the derivative for the compressibility term (m-1), no volume expansion for total water
          dCompress_dPsi(iLayer) = specificStorage*(mLayerVolFracLiqTrial(iLayer) + mLayerVolFracIceTrial(iLayer))/theta_sat(iLayer)
          ! compute the compressibility term (-) instantaneously
          compress(iLayer)       =   mLayerMatricHeadPrime(iLayer) * dCompress_dPsi(iLayer)
      end if
    end do
  else
    compress(:)       = 0._rkind
    dCompress_dPsi(:) = 0._rkind
  end if
end subroutine soilCmpresPrime

end module computFlux_module

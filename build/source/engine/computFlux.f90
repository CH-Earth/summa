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

! provide access to the derived types and classes used to define data structures and class objects
USE data_types,only:&
                    var_i,              & ! data vector (i4b)
                    var_d,              & ! data vector (rkind)
                    var_ilength,        & ! data vector with variable length dimension (i4b)
                    var_dlength,        & ! data vector with variable length dimension (rkind)
                    model_options,      & ! defines the model decisions
                    in_type_vegNrgFlux,out_type_vegNrgFlux,                   & ! classes for vegNrgFlux call
                    in_type_ssdNrgFlux,io_type_ssdNrgFlux,out_type_ssdNrgFlux,& ! classes for ssdNrgFlux call
                    in_type_vegLiqFlux,out_type_vegLiqFlux,                   & ! classes for vegLiqFlux call
                    in_type_snowLiqFlx,io_type_snowLiqFlx,out_type_snowLiqFlx,& ! classes for snowLiqFlx call                
                    in_type_soilLiqFlx,io_type_soilLiqFlx,out_type_soilLiqFlx,& ! classes for soilLiqFlx call
                    in_type_groundwatr,io_type_groundwatr,out_type_groundwatr,& ! classes for groundwatr call
                    in_type_bigAquifer,io_type_bigAquifer,out_type_bigAquifer   ! classes for bigAquifer call

! indices that define elements of the data structures
USE var_lookup,only:iLookDECISIONS  ! named variables for elements of the decision structure
USE var_lookup,only:iLookPARAM      ! named variables for structure elements
USE var_lookup,only:iLookFORCE      ! named variables for structure elements
USE var_lookup,only:iLookPROG       ! named variables for structure elements
USE var_lookup,only:iLookINDEX      ! named variables for structure elements
USE var_lookup,only:iLookDIAG       ! named variables for structure elements
USE var_lookup,only:iLookFLUX       ! named variables for structure elements
USE var_lookup,only:iLookDERIV      ! named variables for structure elements

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
  ! -------------------------------------------------------------------------------------------------------------------------
  ! * dummy variables
  ! -------------------------------------------------------------------------------------------------------------------------
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
  ! -------------------------------------------------------------------------------------------------------------------------
  ! * local variables
  ! -------------------------------------------------------------------------------------------------------------------------
  integer(i4b)                       :: local_ixGroundwater         ! local index for groundwater representation
  integer(i4b)                       :: iLayer                      ! index of model layers
  logical(lgt)                       :: doVegNrgFlux                ! flag to compute the energy flux over vegetation
  real(rkind),dimension(nSoil)       :: dHydCond_dMatric            ! derivative in hydraulic conductivity w.r.t matric head (s-1)
  character(LEN=256)                 :: cmessage                    ! error message of downwind routine
  real(rkind)                        :: above_soilLiqFluxDeriv      ! derivative in layer above soil (canopy or snow) liquid flux w.r.t. liquid water
  real(rkind)                        :: above_soildLiq_dTk          ! derivative of layer above soil (canopy or snow) liquid flux w.r.t. temperature
  real(rkind)                        :: above_soilFracLiq           ! fraction of liquid water layer above soil (canopy or snow) (-)
  ! ---------------------- classes for flux subroutine arguments (classes defined in data_types module) ----------------------
  !      ** intent(in) arguments **       ||       ** intent(inout) arguments **        ||      ** intent(out) arguments **
  type(in_type_vegNrgFlux) :: in_vegNrgFlux;                                            type(out_type_vegNrgFlux) :: out_vegNrgFlux ! vegNrgFlux arguments
  type(in_type_ssdNrgFlux) :: in_ssdNrgFlux; type(io_type_ssdNrgFlux) :: io_ssdNrgFlux; type(out_type_ssdNrgFlux) :: out_ssdNrgFlux ! ssdNrgFlux arguments
  type(in_type_vegLiqFlux) :: in_vegLiqFlux;                                            type(out_type_vegLiqFlux) :: out_vegLiqFlux ! vegLiqFlux arguments
  type(in_type_snowLiqFlx) :: in_snowLiqFlx; type(io_type_snowLiqFlx) :: io_snowLiqFlx; type(out_type_snowLiqFlx) :: out_snowLiqFlx ! snowLiqFlx arguments
  type(in_type_soilLiqFlx) :: in_soilLiqFlx; type(io_type_soilLiqFlx) :: io_soilLiqFlx; type(out_type_soilLiqFlx) :: out_soilLiqFlx ! soilLiqFlx arguments
  type(in_type_groundwatr) :: in_groundwatr; type(io_type_groundwatr) :: io_groundwatr; type(out_type_groundwatr) :: out_groundwatr ! groundwatr arguments
  type(in_type_bigAquifer) :: in_bigAquifer; type(io_type_bigAquifer) :: io_bigAquifer; type(out_type_bigAquifer) :: out_bigAquifer ! bigAquifer arguments
  ! -------------------------------------------------------------------------------------------------------------------------
  ! initialize error control
  err=0; message='computFlux/'

  call initialize_computFlux ! Preliminary operations to start routine

  ! *** CALCULATE ENERGY FLUXES OVER VEGETATION ***
  associate(&
    ixCasNrg => indx_data%var(iLookINDEX%ixCasNrg)%dat(1), & ! intent(in): [i4b] index of canopy air space energy state variable
    ixVegNrg => indx_data%var(iLookINDEX%ixVegNrg)%dat(1), & ! intent(in): [i4b] index of canopy energy state variable
    ixTopNrg => indx_data%var(iLookINDEX%ixTopNrg)%dat(1)  ) ! intent(in): [i4b] index of upper-most energy state in the snow+soil subdomain
    ! identify the need to calculate the energy flux over vegetation
    doVegNrgFlux = (ixCasNrg/=integerMissing .or. ixVegNrg/=integerMissing .or. ixTopNrg/=integerMissing)
    if (doVegNrgFlux) then ! if necessary, calculate the energy fluxes over vegetation
      call initialize_vegNrgFlux
      call vegNrgFlux(in_vegNrgFlux,type_data,forc_data,mpar_data,indx_data,prog_data,diag_data,flux_data,bvar_data,model_decisions,out_vegNrgFlux)
      call finalize_vegNrgFlux
    end if
  end associate

  ! *** CALCULATE ENERGY FLUXES THROUGH THE SNOW-SOIL DOMAIN ***
  associate(nSnowSoilNrg => indx_data%var(iLookINDEX%nSnowSoilNrg)%dat(1)) ! intent(in): [i4b] number of energy state variables in the snow+soil domain
    if (nSnowSoilNrg>0) then ! if necessary, calculate energy fluxes at layer interfaces through the snow and soil domain
      call initialize_ssdNrgFlux
      call ssdNrgFlux(in_ssdNrgFlux,mpar_data,indx_data,prog_data,diag_data,flux_data,io_ssdNrgFlux,out_ssdNrgFlux)
      call finalize_ssdNrgFlux
    end if
  end associate

  ! *** CALCULATE THE LIQUID FLUX THROUGH VEGETATION ***
  associate(ixVegHyd => indx_data%var(iLookINDEX%ixVegHyd)%dat(1)) ! intent(in): [i4b] index of canopy hydrology state variable (mass)
    if (ixVegHyd/=integerMissing) then ! if necessary, calculate liquid water fluxes through vegetation
      call initialize_vegLiqFlux
      call vegLiqFlux(in_vegLiqFlux,mpar_data,diag_data,out_vegLiqFlux)
      call finalize_vegLiqFlux
    end if
  end associate

  ! *** CALCULATE THE LIQUID FLUX THROUGH SNOW ***
  associate(nSnowOnlyHyd => indx_data%var(iLookINDEX%nSnowOnlyHyd)%dat(1)) ! intent(in): [i4b] number of hydrology variables in the snow domain
    if (nSnowOnlyHyd>0) then ! if necessary, compute liquid fluxes through snow
      call initialize_snowLiqFlx
      call snowLiqFlx(in_snowLiqFlx,indx_data,mpar_data,prog_data,diag_data,io_snowLiqFlx,out_snowLiqFlx)
      call finalize_snowLiqFlx
    else
      call soilForcingNoSnow ! define forcing for the soil domain for the case of no snow layers
    end if
  end associate

  ! *** CALCULATE THE LIQUID FLUX THROUGH SOIL ***
  associate(nSoilOnlyHyd => indx_data%var(iLookINDEX%nSoilOnlyHyd)%dat(1)) ! intent(in): [i4b] number of hydrology variables in the soil domain
    if (nSoilOnlyHyd>0) then ! if necessary, calculate the liquid flux through soil
      call initialize_soilLiqFlx
      call soilLiqFlx(in_soilLiqFlx,mpar_data,indx_data,prog_data,diag_data,flux_data,io_soilLiqFlx,out_soilLiqFlx)
      call finalize_soilLiqFlx
    end if 
  end associate

  ! *** CALCULATE THE GROUNDWATER FLOW ***
  associate(nSoilOnlyHyd => indx_data%var(iLookINDEX%nSoilOnlyHyd)%dat(1)) ! intent(in): [i4b] number of hydrology variables in the soil domain
    if (nSoilOnlyHyd>0) then ! check if computing soil hydrology
      if (local_ixGroundwater/=qbaseTopmodel) then ! set baseflow fluxes to zero if the topmodel baseflow routine is not used
        call zeroBaseflowFluxes
      else ! compute the baseflow flux for topmodel-ish shallow groundwater
        call initialize_groundwatr
        call groundwatr(in_groundwatr,attr_data,mpar_data,prog_data,diag_data,flux_data,io_groundwatr,out_groundwatr)
        call finalize_groundwatr
      end if
      call computeBaseflowRunoff ! compute total baseflow from soil and runoff
    end if
  end associate

  ! *** CALCULATE FLUXES FOR THE DEEP AQUIFER ***
  associate(ixAqWat => indx_data%var(iLookINDEX%ixAqWat)%dat(1)) ! intent(in): [i4b] index of water storage in the aquifer
    if (ixAqWat/=integerMissing) then ! check if computing aquifer fluxes
      if (local_ixGroundwater==bigBucket) then ! compute fluxes for the big bucket
        call initialize_bigAquifer
        call bigAquifer(in_bigAquifer,mpar_data,diag_data,io_bigAquifer,out_bigAquifer)
        call finalize_bigAquifer
      else ! if no aquifer, then fluxes are zero
        call zeroAquiferFluxes
      end if ! end check aquifer model decision
    end if  ! if computing aquifer fluxes
  end associate

  call finalize_computFlux ! final operations to prep for end of routine

contains

 ! **** Subroutines that handle the absence of model features ****
 subroutine soilForcingNoSnow
  ! define forcing for the soil domain for the case of no snow layers
  ! NOTE: in case where nSnowOnlyHyd==0 AND snow layers exist, then scalarRainPlusMelt is taken from the previous flux evaluation
  associate(&
   scalarRainPlusMelt           => flux_data%var(iLookFLUX%scalarRainPlusMelt)%dat(1),      & ! intent(out): [dp] rain plus melt (m s-1)
   scalarThroughfallRain        => flux_data%var(iLookFLUX%scalarThroughfallRain)%dat(1),   & ! intent(out): [dp] rain that reaches the ground without ever touching the canopy (kg m-2 s-1)
   scalarCanopyLiqDrainage      => flux_data%var(iLookFLUX%scalarCanopyLiqDrainage)%dat(1), & ! intent(out): [dp] drainage of liquid water from the vegetation canopy (kg m-2 s-1)
   ixVegHyd                     => indx_data%var(iLookINDEX%ixVegHyd)%dat(1),               & ! intent(in): [i4b]    index of canopy hydrology state variable (mass)
   scalarCanopyLiqDeriv         => deriv_data%var(iLookDERIV%scalarCanopyLiqDeriv)%dat(1),  & ! intent(out): [dp] derivative in (throughfall + drainage) w.r.t. canopy liquid water
   dCanLiq_dTcanopy             => deriv_data%var(iLookDERIV%dCanLiq_dTcanopy)%dat(1),      & ! intent(out): [dp] derivative of canopy liquid storage w.r.t. temperature
   scalarFracLiqVeg             => diag_data%var(iLookDIAG%scalarFracLiqVeg)%dat(1),        & ! intent(inout): [dp] fraction of liquid water on vegetation (-)
   iLayerLiqFluxSnowDeriv       => deriv_data%var(iLookDERIV%iLayerLiqFluxSnowDeriv)%dat,   & ! intent(out): [dp(:)] derivative in vertical liquid water flux at layer interfaces
   mLayerdTheta_dTk             => deriv_data%var(iLookDERIV%mLayerdTheta_dTk)%dat,         & ! intent(in):  [dp(:)] derivative of volumetric liquid water content w.r.t. temperature
   mLayerFracLiqSnow            => diag_data%var(iLookDIAG%mLayerFracLiqSnow)%dat           ) ! intent(inout): [dp(:)] fraction of liquid water in each snow layer (-)
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
   end if ! snow layers or not
  end associate
 end subroutine soilForcingNoSnow

 subroutine zeroBaseflowFluxes
  ! set baseflow fluxes to zero if the topmodel baseflow routine is not used
  associate(&
   scalarExfiltration           => flux_data%var(iLookFLUX%scalarExfiltration)%dat(1), & ! intent(out): [dp] exfiltration from the soil profile (m s-1)
   mLayerColumnOutflow          => flux_data%var(iLookFLUX%mLayerColumnOutflow)%dat,   & ! intent(out): [dp(:)] column outflow from each soil layer (m3 s-1)
   mLayerBaseflow               => flux_data%var(iLookFLUX%mLayerBaseflow)%dat         ) ! intent(out): [dp(:)] baseflow from each soil layer (m s-1)
   ! diagnostic variables in the data structures
   scalarExfiltration     = 0._rkind  ! exfiltration from the soil profile (m s-1)
   mLayerColumnOutflow(:) = 0._rkind  ! column outflow from each soil layer (m3 s-1)
   ! variables needed for the numerical solution
   mLayerBaseflow(:)      = 0._rkind  ! baseflow from each soil layer (m s-1)
  end associate
 end subroutine zeroBaseflowFluxes

 subroutine computeBaseflowRunoff
  ! compute total baseflow from the soil zone (needed for mass balance checks) and total runoff
  ! (Note: scalarSoilBaseflow is zero if topmodel is not used)
  ! (Note: scalarSoilBaseflow may need to re-envisioned in topmodel formulation if parts of it flow into neighboring soil rather than exfiltrate)
  associate(&
   scalarSoilBaseflow           => flux_data%var(iLookFLUX%scalarSoilBaseflow)%dat(1),  & ! intent(out): [dp] total baseflow from the soil profile (m s-1)
   mLayerBaseflow               => flux_data%var(iLookFLUX%mLayerBaseflow)%dat,         & ! intent(out): [dp(:)] baseflow from each soil layer (m s-1)
   scalarTotalRunoff            => flux_data%var(iLookFLUX%scalarTotalRunoff)%dat(1),   & ! intent(out): [dp] total runoff (m s-1)
   scalarSurfaceRunoff          => flux_data%var(iLookFLUX%scalarSurfaceRunoff)%dat(1), & ! intent(out): [dp] surface runoff (m s-1)
   scalarSoilDrainage           => flux_data%var(iLookFLUX%scalarSoilDrainage)%dat(1)   ) ! intent(out): [dp] drainage from the soil profile (m s-1)
   scalarSoilBaseflow = sum(mLayerBaseflow)                                               ! baseflow from the soil zone 
   scalarTotalRunoff  = scalarSurfaceRunoff + scalarSoilDrainage + scalarSoilBaseflow     ! total runoff
  end associate
 end subroutine computeBaseflowRunoff  

 subroutine zeroAquiferFluxes
  ! set aquifer fluxes to zero (if no aquifer exists)
  associate(&
   scalarAquiferTranspire       => flux_data%var(iLookFLUX%scalarAquiferTranspire)%dat(1), & ! intent(out): [dp] transpiration loss from the aquifer (m s-1
   scalarAquiferRecharge        => flux_data%var(iLookFLUX%scalarAquiferRecharge)%dat(1),  & ! intent(out): [dp] recharge to the aquifer (m s-1)
   scalarAquiferBaseflow        => flux_data%var(iLookFLUX%scalarAquiferBaseflow)%dat(1),  & ! intent(out): [dp] total baseflow from the aquifer (m s-1)
   dBaseflow_dAquifer           => deriv_data%var(iLookDERIV%dBaseflow_dAquifer)%dat(1)    ) ! intent(out): [dp(:)] derivative in baseflow flux w.r.t. aquifer storage (s-1)
   scalarAquiferTranspire = 0._rkind  ! transpiration loss from the aquifer (m s-1)
   scalarAquiferRecharge  = 0._rkind  ! recharge to the aquifer (m s-1)
   scalarAquiferBaseflow  = 0._rkind  ! total baseflow from the aquifer (m s-1)
   dBaseflow_dAquifer     = 0._rkind  ! change in baseflow flux w.r.t. aquifer storage (s-1)
  end associate
 end subroutine zeroAquiferFluxes

 ! **** Subroutines for starting/ending operations of computFlux ****
 subroutine initialize_computFlux
  ! operations to prep for the start of computFlux
  associate(&
   numFluxCalls                 => diag_data%var(iLookDIAG%numFluxCalls)%dat(1),         & ! intent(out): [dp] number of flux calls (-)
   ixSpatialGroundwater         => model_decisions(iLookDECISIONS%spatial_gw)%iDecision, & ! intent(in): [i4b] spatial representation of groundwater (local-column or single-basin)
   ixGroundwater                => model_decisions(iLookDECISIONS%groundwatr)%iDecision, & ! intent(in): [i4b] groundwater parameterization
   iLayerLiqFluxSnow            => flux_data%var(iLookFLUX%iLayerLiqFluxSnow)%dat,       & ! intent(out): [dp(0:)] vertical liquid water flux at snow layer interfaces (-)
   iLayerLiqFluxSoil            => flux_data%var(iLookFLUX%iLayerLiqFluxSoil)%dat        ) ! intent(out): [dp(0:)] vertical liquid water flux at soil layer interfaces (-)

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
  end associate
 end subroutine initialize_computFlux

 subroutine finalize_computFlux
  ! operations to prep for the end of computFlux
  associate(&
   ixCasNrg                     => indx_data%var(iLookINDEX%ixCasNrg)%dat(1),              & ! intent(in): [i4b] index of canopy air space energy state variable
   ixVegNrg                     => indx_data%var(iLookINDEX%ixVegNrg)%dat(1),              & ! intent(in): [i4b] index of canopy energy state variable
   ixVegHyd                     => indx_data%var(iLookINDEX%ixVegHyd)%dat(1),              & ! intent(in): [i4b] index of canopy hydrology state variable (mass)
   scalarCanairNetNrgFlux       => flux_data%var(iLookFLUX%scalarCanairNetNrgFlux)%dat(1), & ! intent(out): [dp] net energy flux for the canopy air space (W m-2)
   scalarCanopyNetNrgFlux       => flux_data%var(iLookFLUX%scalarCanopyNetNrgFlux)%dat(1), & ! intent(out): [dp] net energy flux for the vegetation canopy       (W m-2)
   scalarCanopyNetLiqFlux       => flux_data%var(iLookFLUX%scalarCanopyNetLiqFlux)%dat(1), & ! intent(out): [dp] net liquid water flux for the vegetation canopy (kg m-2 s-1)
   canopyDepth                  => diag_data%var(iLookDIAG%scalarCanopyDepth)%dat(1),      & ! intent(in): [dp] canopy depth (m)
   nSnowSoilNrg                 => indx_data%var(iLookINDEX%nSnowSoilNrg)%dat(1),          & ! intent(in): [i4b] number of energy state variables in the snow+soil domain
   ixSnowSoilNrg                => indx_data%var(iLookINDEX%ixSnowSoilNrg)%dat,            & ! intent(in): [i4b(:)] indices for energy states in the snow+soil subdomain
   mLayerNrgFlux                => flux_data%var(iLookFLUX%mLayerNrgFlux)%dat              ) ! intent(out): [dp] net energy flux for each layer within the snow+soil domain (J m-3 s-1)
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
  end associate

  associate(&
   ixAqWat                      => indx_data%var(iLookINDEX%ixAqWat)%dat(1),               & ! intent(in): [i4b] index of water storage in the aquifer
   ixSnowSoilHyd                => indx_data%var(iLookINDEX%ixSnowSoilHyd)%dat,            & ! intent(in): [i4b(:)] indices for hydrology states in the snow+soil subdomain
   nSnowSoilHyd                 => indx_data%var(iLookINDEX%nSnowSoilHyd)%dat(1),          & ! intent(in): [i4b] number of hydrology variables in the snow+soil domain
   layerType                    => indx_data%var(iLookINDEX%layerType)%dat,                & ! intent(in): [i4b(:)] type of layer (iname_soil or iname_snow)
   mLayerLiqFluxSnow            => flux_data%var(iLookFLUX%mLayerLiqFluxSnow)%dat,         & ! intent(out): [dp] net liquid water flux for each snow layer (s-1)
   mLayerLiqFluxSoil            => flux_data%var(iLookFLUX%mLayerLiqFluxSoil)%dat,         & ! intent(out): [dp] net liquid water flux for each soil layer (s-1)
   scalarAquiferTranspire       => flux_data%var(iLookFLUX%scalarAquiferTranspire)%dat(1), & ! intent(out): [dp] transpiration loss from the aquifer (m s-1
   scalarAquiferRecharge        => flux_data%var(iLookFLUX%scalarAquiferRecharge)%dat(1),  & ! intent(out): [dp] recharge to the aquifer (m s-1)
   scalarAquiferBaseflow        => flux_data%var(iLookFLUX%scalarAquiferBaseflow)%dat(1)   ) ! intent(out): [dp] total baseflow from the aquifer (m s-1)
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
  end associate

   firstFluxCall=.false. ! set the first flux call to false
 end subroutine finalize_computFlux

 ! ----------------------- Initialize and Finalize procedures for the flux routines -----------------------
 ! **** vegNrgFlux ****
 subroutine initialize_vegNrgFlux
  associate(&
   dCanLiq_dTcanopy             => deriv_data%var(iLookDERIV%dCanLiq_dTcanopy)%dat(1), & ! intent(out): [dp] derivative of canopy liquid storage w.r.t. temperature
   dTheta_dTkCanopy             => deriv_data%var(iLookDERIV%dTheta_dTkCanopy)%dat(1), & ! intent(in):  [dp] derivative of volumetric liquid water content w.r.t. temperature
   canopyDepth                  => diag_data%var(iLookDIAG%scalarCanopyDepth)%dat(1)   ) ! intent(in): [dp]  canopy depth (m)

   dCanLiq_dTcanopy = dTheta_dTkCanopy*iden_water*canopyDepth     ! derivative in canopy liquid storage w.r.t. canopy temperature (kg m-2 K-1)
  end associate
  call in_vegNrgFlux % initialize(firstSubStep,firstFluxCall,computeVegFlux,checkLWBalance,&
                                  scalarCanairTempTrial,scalarCanopyTempTrial,mLayerTempTrial,scalarCanopyIceTrial,&
                                  scalarCanopyLiqTrial,forc_data,deriv_data)
 end subroutine initialize_vegNrgFlux

 subroutine finalize_vegNrgFlux
  call out_vegNrgFlux%finalize(flux_data,deriv_data,err,cmessage)
  associate(&
   canopyDepth                  => diag_data%var(iLookDIAG%scalarCanopyDepth)%dat(1),      & ! intent(in): [dp   ]  canopy depth (m)
   mLayerDepth                  => prog_data%var(iLookPROG%mLayerDepth)%dat,               & ! intent(in): [dp(:)]  depth of each layer in the snow-soil sub-domain (m)
   scalarCanairNetNrgFlux       => flux_data%var(iLookFLUX%scalarCanairNetNrgFlux)%dat(1), & ! intent(out): [dp] net energy flux for the canopy air space  (W m-2)
   scalarCanopyNetNrgFlux       => flux_data%var(iLookFLUX%scalarCanopyNetNrgFlux)%dat(1), & ! intent(out): [dp] net energy flux for the vegetation canopy (W m-2)
   scalarGroundNetNrgFlux       => flux_data%var(iLookFLUX%scalarGroundNetNrgFlux)%dat(1), & ! intent(out): [dp] net energy flux for the ground surface    (W m-2)
   dGroundNetFlux_dGroundTemp   => deriv_data%var(iLookDERIV%dGroundNetFlux_dGroundTemp)%dat(1) ) ! intent(out): [dp] derivative in net ground flux w.r.t. ground temperature
   ! error control
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
  end associate
 end subroutine finalize_vegNrgFlux
 ! **** end vegNrgFlux ****

 ! **** ssdNrgFlux ****
 subroutine initialize_ssdNrgFlux
  call in_ssdNrgFlux%initialize(scalarSolution,firstFluxCall,mLayerTempTrial,flux_data,deriv_data)
  call io_ssdNrgFlux%initialize(deriv_data)
 end subroutine initialize_ssdNrgFlux

 subroutine finalize_ssdNrgFlux
  call io_ssdNrgFlux%finalize(deriv_data)
  call out_ssdNrgFlux%finalize(flux_data,deriv_data,err,cmessage)
  associate(&
   mLayerNrgFlux                => flux_data%var(iLookFLUX%mLayerNrgFlux)%dat, & ! intent(out): [dp] net energy flux for each layer within the snow+soil domain (J m-3 s-1)
   iLayerNrgFlux                => flux_data%var(iLookFLUX%iLayerNrgFlux)%dat, & ! intent(out): [dp(0:)] vertical energy flux at the interface of snow and soil layers
   mLayerDepth                  => prog_data%var(iLookPROG%mLayerDepth)%dat    ) ! intent(in): [dp(:)]  depth of each layer in the snow-soil sub-domain (m)
   ! error control
   if (err/=0) then; message=trim(message)//trim(cmessage); return; end if
   ! calculate net energy fluxes for each snow and soil layer (J m-3 s-1)
   do iLayer=1,nLayers
     mLayerNrgFlux(iLayer) = -(iLayerNrgFlux(iLayer) - iLayerNrgFlux(iLayer-1))/mLayerDepth(iLayer)
     if (globalPrintFlag) then
       if (iLayer < 10) write(*,'(a,1x,i4,1x,10(f25.15,1x))') 'iLayer, iLayerNrgFlux(iLayer-1:iLayer), mLayerNrgFlux(iLayer)   = ', iLayer, iLayerNrgFlux(iLayer-1:iLayer), mLayerNrgFlux(iLayer)
     end if
   end do
  end associate
 end subroutine finalize_ssdNrgFlux
 ! **** end ssdNrgFlux ****

 ! **** vegLiqFlux ****
 subroutine initialize_vegLiqFlux
  call in_vegLiqFlux%initialize(computeVegFlux,scalarCanopyLiqTrial,flux_data)
 end subroutine initialize_vegLiqFlux
 
 subroutine finalize_vegLiqFlux
  call out_vegLiqFlux%finalize(flux_data,deriv_data,err,cmessage)
  associate( &
   scalarThroughfallRain        => flux_data%var(iLookFLUX%scalarThroughfallRain)%dat(1),         & ! intent(out): [dp] rain that reaches the ground without ever touching the canopy (kg m-2 s-1)
   scalarCanopyLiqDrainage      => flux_data%var(iLookFLUX%scalarCanopyLiqDrainage)%dat(1),       & ! intent(out): [dp] drainage of liquid water from the vegetation canopy (kg m-2 s-1)
   scalarThroughfallRainDeriv   => deriv_data%var(iLookDERIV%scalarThroughfallRainDeriv  )%dat(1),& ! intent(out): [dp] derivative in throughfall w.r.t. canopy liquid water
   scalarCanopyLiqDrainageDeriv => deriv_data%var(iLookDERIV%scalarCanopyLiqDrainageDeriv)%dat(1),& ! intent(out): [dp] derivative in canopy drainage w.r.t. canopy liquid water
   scalarCanopyNetLiqFlux       => flux_data%var(iLookFLUX%scalarCanopyNetLiqFlux)%dat(1),        & ! intent(out): [dp] net liquid water flux for the vegetation canopy (kg m-2 s-1)
   scalarRainfall               => flux_data%var(iLookFLUX%scalarRainfall)%dat(1),                & ! intent(in):  [dp] rainfall rate (kg m-2 s-1)
   scalarCanopyEvaporation      => flux_data%var(iLookFLUX%scalarCanopyEvaporation)%dat(1),       & ! intent(out): [dp] canopy evaporation/condensation (kg m-2 s-1)
   scalarCanopyLiqDeriv         => deriv_data%var(iLookDERIV%scalarCanopyLiqDeriv        )%dat(1) ) ! intent(out): [dp] derivative in (throughfall + drainage) w.r.t. canopy liquid water
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
  end associate
 end subroutine finalize_vegLiqFlux
 ! **** end vegLiqFlux ****

 ! **** snowLiqFlx ****
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
 ! **** end snowLiqFlx ****

 ! **** soilLiqFlx ****
 subroutine initialize_soilLiqFlx
  call in_soilLiqFlx%initialize(nsnow,nSoil,nlayers,firstSplitOper,scalarSolution,firstFluxCall,&
                                mLayerTempTrial,mLayerMatricHeadTrial,mLayerMatricHeadLiqTrial,mLayerVolFracLiqTrial,mLayerVolFracIceTrial,&
                                above_soilLiqFluxDeriv,above_soildLiq_dTk,above_soilFracLiq,flux_data,deriv_data)
  call io_soilLiqFlx%initialize(nsoil,dHydCond_dMatric,flux_data,diag_data,deriv_data)
 end subroutine initialize_soilLiqFlx

 subroutine finalize_soilLiqFlx
  call io_soilLiqFlx%finalize(nsoil,dHydCond_dMatric,flux_data,diag_data,deriv_data)
  call out_soilLiqFlx%finalize(err,cmessage)
  ! error control
  if (err/=0) then; message=trim(message)//trim(cmessage); return; end if

  associate(&
   mLayerLiqFluxSoil            => flux_data%var(iLookFLUX%mLayerLiqFluxSoil)%dat,     & ! intent(out): [dp] net liquid water flux for each soil layer (s-1)
   iLayerLiqFluxSoil            => flux_data%var(iLookFLUX%iLayerLiqFluxSoil)%dat,     & ! intent(out): [dp(0:)] vertical liquid water flux at soil layer interfaces (-)
   mLayerDepth                  => prog_data%var(iLookPROG%mLayerDepth)%dat,           & ! intent(in): [dp(:)]  depth of each layer in the snow-soil sub-domain (m)
   scalarMaxInfilRate           => flux_data%var(iLookFLUX%scalarMaxInfilRate)%dat(1), & ! intent(out): [dp] maximum infiltration rate (m s-1)
   scalarRainPlusMelt           => flux_data%var(iLookFLUX%scalarRainPlusMelt)%dat(1), & ! intent(out): [dp] rain plus melt (m s-1)
   scalarSoilControl            => diag_data%var(iLookDIAG%scalarSoilControl )%dat(1), & ! intent(out): [dp] soil control on infiltration, zero or one
   scalarInfilArea              => diag_data%var(iLookDIAG%scalarInfilArea   )%dat(1), & ! intent(out): [dp] fraction of unfrozen area where water can infiltrate (-)
   scalarFrozenArea             => diag_data%var(iLookDIAG%scalarFrozenArea  )%dat(1), & ! intent(out): [dp] fraction of area that is considered impermeable due to soil ice (-)
   scalarSoilDrainage           => flux_data%var(iLookFLUX%scalarSoilDrainage)%dat(1)  ) ! intent(out): [dp]     drainage from the soil profile (m s-1)

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
  end associate

  associate(&
   dq_dHydStateAbove            => deriv_data%var(iLookDERIV%dq_dHydStateAbove)%dat,        & ! intent(out): [dp(:)] change in flux at layer interfaces w.r.t. states in the layer above
   dq_dHydStateBelow            => deriv_data%var(iLookDERIV%dq_dHydStateBelow)%dat,        & ! intent(out): [dp(:)] change in flux at layer interfaces w.r.t. states in the layer below
   dq_dHydStateLayerSurfVec     => deriv_data%var(iLookDERIV%dq_dHydStateLayerSurfVec)%dat, & ! intent(out): [dp(:)] change in the flux in soil surface interface w.r.t. state variables in layers
   dPsiLiq_dPsi0                => deriv_data%var(iLookDERIV%dPsiLiq_dPsi0   )%dat          ) ! intent(in):  [dp(:)] derivative in liquid water matric pot w.r.t. the total water matric pot (-)
   ! expand derivatives to the total water matric potential
   ! NOTE: arrays are offset because computing derivatives in interface fluxes, at the top and bottom of the layer respectively
   if (globalPrintFlag) print*, 'dPsiLiq_dPsi0(1:nSoil) = ', dPsiLiq_dPsi0(1:nSoil)
   dq_dHydStateAbove(1:nSoil)   = dq_dHydStateAbove(1:nSoil)  *dPsiLiq_dPsi0(1:nSoil)
   dq_dHydStateBelow(0:nSoil-1) = dq_dHydStateBelow(0:nSoil-1)*dPsiLiq_dPsi0(1:nSoil)
   dq_dHydStateLayerSurfVec(1:nSoil) = dq_dHydStateLayerSurfVec(1:nSoil)*dPsiLiq_dPsi0(1:nSoil)
  end associate
 end subroutine finalize_soilLiqFlx
 ! **** end soilLiqFlx ****

 ! **** groundwatr ****
 subroutine initialize_groundwatr
  ! check the derivative matrix is sized appropriately
  if (size(dBaseflow_dMatric,1)/=nSoil .or. size(dBaseflow_dMatric,2)/=nSoil) then
    message=trim(message)//'expect dBaseflow_dMatric to be nSoil x nSoil'
    err=20; return
  end if
  call in_groundwatr%initialize(nSnow,nSoil,nLayers,firstFluxCall,mLayerMatricHeadLiqTrial,mLayerVolFracLiqTrial,mLayerVolFracIceTrial,deriv_data)
  call io_groundwatr%initialize(ixSaturation)
 end subroutine initialize_groundwatr

 subroutine finalize_groundwatr
  call io_groundwatr%finalize(ixSaturation)
  call out_groundwatr%finalize(dBaseflow_dMatric,flux_data,err,cmessage)
  ! error control
  if (err/=0) then; message=trim(message)//trim(cmessage); return; end if
 end subroutine finalize_groundwatr
 ! **** end groundwatr ****

 ! **** bigAquifer ****
 subroutine initialize_bigAquifer
  call in_bigAquifer%initialize(scalarAquiferStorageTrial,flux_data,deriv_data)
  call io_bigAquifer%initialize(deriv_data)
 end subroutine initialize_bigAquifer

 subroutine finalize_bigAquifer
  call io_bigAquifer%finalize(deriv_data)
  call out_bigAquifer%finalize(flux_data,deriv_data,err,cmessage)
  associate(&
   scalarTotalRunoff            => flux_data%var(iLookFLUX%scalarTotalRunoff)%dat(1)               ,&  ! intent(out): [dp] total runoff (m s-1)
   scalarSurfaceRunoff          => flux_data%var(iLookFLUX%scalarSurfaceRunoff)%dat(1)             ,&  ! intent(out): [dp] surface runoff (m s-1)
   scalarAquiferBaseflow        => flux_data%var(iLookFLUX%scalarAquiferBaseflow)%dat(1) )             ! intent(out): [dp] total baseflow from the aquifer (m s-1)
   ! error control
   if (err/=0) then; message=trim(message)//trim(cmessage); return; end if
   ! compute total runoff (overwrite previously calculated value before considering aquifer).
   !   (Note:  SoilDrainage goes into aquifer, not runoff)
   scalarTotalRunoff  = scalarSurfaceRunoff + scalarAquiferBaseflow     
  end associate
 end subroutine finalize_bigAquifer
 ! **** end bigAquifer ****

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

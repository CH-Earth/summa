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
USE nr_type

! provide access to the derived types and classes used to define data structures and class objects
USE data_types,only:&
                    var_i,              & ! data vector (i4b)
                    var_d,              & ! data vector (rkind)
                    var_ilength,        & ! data vector with variable length dimension (i4b)
                    var_dlength,        & ! data vector with variable length dimension (rkind)
                    model_options,      & ! defines the model decisions
                    in_type_vegNrgFlux,out_type_vegNrgFlux,                   & ! classes for vegNrgFlux call
                    in_type_snowLakeSoilGlceNrgFlux,io_type_snowLakeSoilGlceNrgFlux,out_type_snowLakeSoilGlceNrgFlux,& ! classes for snowLakeSoilGlceNrgFlux call
                    in_type_vegLiqFlux,out_type_vegLiqFlux,                   & ! classes for vegLiqFlux call
                    in_type_snowLakeGlceLiqFlux,io_type_snowLakeGlceLiqFlux,out_type_snowLakeGlceLiqFlux,& ! classes for snowLakeGlceLiqFlux call                
                    in_type_soilLiqFlux,io_type_soilLiqFlux,out_type_soilLiqFlux,& ! classes for soilLiqFlux call
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
USE globalData,only:iname_glce      ! named variables for glacier ice
USE globalData,only:iname_lake      ! named variables for lake

! constants
USE multiconst,only:&
                    Tfreeze,      & ! temperature at freezing              (K)
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

! look-up values for the choice of boundary conditions for hydrology
USE mDecisions_module,only:       &
 prescribedHead,                  & ! prescribed head
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
                      nLake,                    & ! intent(in):    number of lake layers
                      nSoil,                    & ! intent(in):    number of soil layers
                      nGlce,                    & ! intent(in):    number of glacier ice layers
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
                      mLayerTempTrial,          & ! intent(in):    trial value for the temperature of each layer (K)
                      mLayerMatricHeadLiqTrial, & ! intent(in):    trial value for the liquid water matric potential in each soil layer (m)
                      mLayerMatricHeadTrial,    & ! intent(in):    trial vector of total water matric potential (m)
                      scalarAquiferStorageTrial,& ! intent(in):    trial value of storage of water in the aquifer (m)
                      ! input: diagnostic variables defining the liquid water and ice content
                      scalarCanopyLiqTrial,     & ! intent(in):    trial value for the liquid water on the vegetation canopy (kg m-2)
                      scalarCanopyIceTrial,     & ! intent(in):    trial value for the ice on the vegetation canopy (kg m-2)
                      mLayerVolFracLiqTrial,    & ! intent(in):    trial value for the volumetric liquid water content in each layer (-)
                      mLayerVolFracIceTrial,    & ! intent(in):    trial value for the volumetric ice in each layer (-)
                      ! input: data structures
                      model_decisions,          & ! intent(in):    model decisions
                      type_data,                & ! intent(in):    type of vegetation and soil
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
                      dBaseflow_dWat,           & ! intent(out):   derivative in baseflow w.r.t. soil water characteristic
                      dBaseflow_dTk,            & ! intent(out):   derivative in baseflow w.r.t. temperature (m s-1 K-1)
                      fluxVec,                  & ! intent(out):   flux vector (mixed units)
                      ! output: error control
                      err,message)                ! intent(out):   error code and error message
  ! provide access to flux subroutines
  USE vegNrgFlux_module,only:vegNrgFlux                           ! compute energy fluxes over vegetation
  USE snowLakeSoilGlceNrgFlux_module,only:snowLakeSoilGlceNrgFlux ! compute energy fluxes throughout the layers
  USE vegLiqFlux_module,only:vegLiqFlux                           ! compute liquid water fluxes through vegetation
  USE snowLakeGlceLiqFlux_module,only:snowLakeGlceLiqFlux         ! compute liquid water fluxes through non-soil layers
  USE soilLiqFlux_module,only:soilLiqFlux                           ! compute liquid water fluxes through soil
  USE groundwatr_module,only:groundwatr                           ! compute the baseflow flux
  USE bigAquifer_module,only:bigAquifer                           ! compute fluxes for the big aquifer
  implicit none
  ! -------------------------------------------------------------------------------------------------------------------------
  ! * dummy variables
  ! -------------------------------------------------------------------------------------------------------------------------
  ! input-output: control
  integer(i4b),intent(in)            :: nSnow                       ! number of snow layers
  integer(i4b),intent(in)            :: nLake                       ! number of lake layers
  integer(i4b),intent(in)            :: nSoil                       ! number of soil layers
  integer(i4b),intent(in)            :: nGlce                       ! number of glacier ice layers
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
  real(rkind),intent(out)            :: dBaseflow_dWat(:,:)         ! derivative in baseflow w.r.t. soil water characteristic
  real(rkind),intent(out)            :: dBaseflow_dTk(:,:)          ! derivative in baseflow w.r.t. temperature (m s-1 K-1)
  real(rkind),intent(out)            :: fluxVec(:)                  ! model flux vector (mixed units)
  ! output: error control
  integer(i4b),intent(out)           :: err                         ! error code
  character(*),intent(out)           :: message                     ! error message
  ! -------------------------------------------------------------------------------------------------------------------------
  ! * local variables
  ! -------------------------------------------------------------------------------------------------------------------------
  logical(lgt)                       :: lake_frozen                 ! flag that lake is frozen
  integer(i4b)                       :: nLake_frz                   ! number of frozen lake layers
  integer(i4b)                       :: local_ixGroundwater         ! local index for groundwater representation
  integer(i4b)                       :: iLayer,nStart               ! index control of model layers
  logical(lgt)                       :: doVegNrgFlux                ! flag to compute the energy flux over vegetation
  real(rkind),dimension(nSoil)       :: dHydCond_dMatric            ! derivative in hydraulic conductivity w.r.t matric head (s-1)
  character(LEN=256)                 :: cmessage                    ! error message of downwind routine
  real(rkind)                        :: surface_flux                ! surface flux (m s-1) into snow or ice
  real(rkind)                        :: bottom_flux                 ! bottom flux (m s-1) out of snow or ice
  ! ---------------------- classes for flux subroutine arguments (classes defined in data_types module) ----------------------
  !      ** intent(in) arguments **       ||       ** intent(inout) arguments **        ||      ** intent(out) arguments **
  type(in_type_vegNrgFlux) :: in_vegNrgFlux;                                            type(out_type_vegNrgFlux) :: out_vegNrgFlux ! vegNrgFlux arguments
  type(in_type_snowLakeSoilGlceNrgFlux) :: in_snowLakeSoilGlceNrgFlux; type(io_type_snowLakeSoilGlceNrgFlux) :: io_snowLakeSoilGlceNrgFlux; type(out_type_snowLakeSoilGlceNrgFlux) :: out_snowLakeSoilGlceNrgFlux ! snowLakeSoilGlceNrgFlux arguments
  type(in_type_vegLiqFlux) :: in_vegLiqFlux;                                            type(out_type_vegLiqFlux) :: out_vegLiqFlux ! vegLiqFlux arguments
  type(in_type_snowLakeGlceLiqFlux) :: in_snowLakeGlceLiqFlux; type(io_type_snowLakeGlceLiqFlux) :: io_snowLakeGlceLiqFlux; type(out_type_snowLakeGlceLiqFlux) :: out_snowLakeGlceLiqFlux ! snowLakeGlceLiqFlux arguments
  type(in_type_soilLiqFlux) :: in_soilLiqFlux; type(io_type_soilLiqFlux) :: io_soilLiqFlux; type(out_type_soilLiqFlux) :: out_soilLiqFlux ! soilLiqFlux arguments
  type(in_type_groundwatr) :: in_groundwatr; type(io_type_groundwatr) :: io_groundwatr; type(out_type_groundwatr) :: out_groundwatr ! groundwatr arguments
  type(in_type_bigAquifer) :: in_bigAquifer; type(io_type_bigAquifer) :: io_bigAquifer; type(out_type_bigAquifer) :: out_bigAquifer ! bigAquifer arguments
  ! -------------------------------------------------------------------------------------------------------------------------
  ! initialize error control
  err=0; message='computFlux/'

  call initialize_computFlux; if(err/=0)then; return; endif ! Preliminary operations to start routine

  ! *** CALCULATE ENERGY FLUXES OVER VEGETATION ***
  associate(&
    ixCasNrg => indx_data%var(iLookINDEX%ixCasNrg)%dat(1), & ! intent(in): [i4b] index of canopy air space energy state variable
    ixVegNrg => indx_data%var(iLookINDEX%ixVegNrg)%dat(1), & ! intent(in): [i4b] index of canopy energy state variable
    ixTopNrg => indx_data%var(iLookINDEX%ixTopNrg)%dat(1)  ) ! intent(in): [i4b] index of upper-most energy state in the snow+lake+soil+glce subdomain
    ! identify the need to calculate the energy flux over vegetation
    doVegNrgFlux = (ixCasNrg/=integerMissing .or. ixVegNrg/=integerMissing .or. ixTopNrg/=integerMissing)
    if (doVegNrgFlux) then ! if necessary, calculate the energy fluxes over vegetation
      call initialize_vegNrgFlux
      call vegNrgFlux(in_vegNrgFlux,type_data,forc_data,mpar_data,indx_data,prog_data,diag_data,flux_data,bvar_data,model_decisions,out_vegNrgFlux)
      call finalize_vegNrgFlux; if(err/=0)then; return; endif
    end if
  end associate

  ! *** CALCULATE ENERGY FLUXES THROUGH THE LAYER DOMAIN ***
  associate(nSnLaSoGlNrg => indx_data%var(iLookINDEX%nSnLaSoGlNrg)%dat(1)) ! intent(in): [i4b] number of energy state variables in the layer domains
    if (nSnLaSoGlNrg>0) then ! if necessary, calculate energy fluxes at layer interfaces through the snow and soil domain
      call initialize_snowLakeSoilGlceNrgFlux
      call snowLakeSoilGlceNrgFlux(in_snowLakeSoilGlceNrgFlux,mpar_data,indx_data,prog_data,diag_data,flux_data,io_snowLakeSoilGlceNrgFlux,out_snowLakeSoilGlceNrgFlux)
      call finalize_snowLakeSoilGlceNrgFlux; if(err/=0)then; return; endif
    end if
  end associate

  ! *** CALCULATE THE LIQUID FLUX THROUGH VEGETATION ***
  associate(ixVegHyd => indx_data%var(iLookINDEX%ixVegHyd)%dat(1)) ! intent(in): [i4b] index of canopy hydrology state variable (mass)
    if (ixVegHyd/=integerMissing) then ! if necessary, calculate liquid water fluxes through vegetation
      call initialize_vegLiqFlux
      call vegLiqFlux(in_vegLiqFlux,mpar_data,diag_data,out_vegLiqFlux)
      call finalize_vegLiqFlux; if(err/=0)then; return; endif
    end if
  end associate

  ! *** CALCULATE THE LIQUID FLUX FROM GLACIER ICE, impermeable ice so fluxes go upwards
  ! NOTE: in a domain, there is no aquifer or veetation if there are glce layers
  associate(nGlceOnlyHyd => indx_data%var(iLookINDEX%nGlceOnlyHyd)%dat(1)) ! intent(in): [i4b] number of hydrology variables in the glacier ice
    if (nGlceOnlyHyd>0) then ! if necessary, calculate the liquid flux through glacier ice
      call initialize_glceLiqFlux
      call snowLakeGlceLiqFlux(in_snowLakeGlceLiqFlux,mpar_data,indx_data,prog_data,diag_data,io_snowLakeGlceLiqFlux,out_snowLakeGlceLiqFlux)
      call finalize_glceLiqFlux; if(err/=0)then; return; endif
    else
      call zeroGlacierFluxes ! set glacier ice fluxes to zero if there are no glacier ice layers
    end if 
  end associate

  ! *** CALCULATE THE LIQUID FLUX THROUGH LAKE ***
  associate(nLakeOnlyHyd => indx_data%var(iLookINDEX%nLakeOnlyHyd)%dat(1)) ! intent(in): [i4b] number of hydrology variables in the lake
    if (nLakeOnlyHyd>0) then ! if necessary, compute liquid fluxes through lake
      !call initialize_lakeLiqFlux ! only liquid flux is out top (possibly into ice) or out bottom of lake
      !call lakeSpillFrzFlux(in_lakeSpillFrzFlux,mpar_data,indx_data,prog_data,diag_data,io_snowLakeGlceLiqFlux,out_snowLakeGlceLiqFlux)
      !call finalize_lakeLiqFlux
      !if(lake_frozen)then ! NOTE: now possible to have snow layers
      !  call initialize_frzlakeLiqFlux
      !  call snowLakeGlceLiqFlux(in_snowLakeGlceLiqFlux,mpar_data,indx_data,prog_data,diag_data,io_snowLakeGlceLiqFlux,out_snowLakeGlceLiqFlux)
      !  call finalize_frzlakeLiqFlux
      !endif
      print*, 'Lake liquid fluxes are not yet implemented'; stop
    else
      call forcingNoLake ! define forcing for the domain beneath for the case of no lake layers
    end if
  end associate

  ! *** CALCULATE THE LIQUID FLUX THROUGH SNOW ***
  associate(nSnowOnlyHyd => indx_data%var(iLookINDEX%nSnowOnlyHyd)%dat(1)) ! intent(in): [i4b] number of hydrology variables in the snow
    if (nSnowOnlyHyd>0) then ! if necessary, compute liquid fluxes through snow
      call initialize_snowLiqFlux
      call snowLakeGlceLiqFlux(in_snowLakeGlceLiqFlux,mpar_data,indx_data,prog_data,diag_data,io_snowLakeGlceLiqFlux,out_snowLakeGlceLiqFlux)
      call finalize_snowLiqFlux; if(err/=0)then; return; endif
    else
      call forcingNoSnow ! define forcing for the domain beneath for the case of no snow layers
    end if
  end associate

  ! *** CALCULATE THE LIQUID FLUX THROUGH SOIL ***
  associate(nSoilOnlyHyd => indx_data%var(iLookINDEX%nSoilOnlyHyd)%dat(1)) ! intent(in): [i4b] number of hydrology variables in the soil
    if (nSoilOnlyHyd>0) then ! if necessary, calculate the liquid flux through soil
      call initialize_soilLiqFlux
      call soilLiqFlux(in_soilLiqFlux,mpar_data,indx_data,prog_data,diag_data,flux_data,io_soilLiqFlux,out_soilLiqFlux)
      call finalize_soilLiqFlux; if(err/=0)then; return; endif
    else
      call forcingNoSoil ! define forcing for the domain beneath for the case of no soil layers
    end if 
  end associate

  ! *** CALCULATE THE SHALLOW GROUNDWATER FLOW OR DEBRIS LATERAL FLOW ***
  associate(nSoilOnlyHyd => indx_data%var(iLookINDEX%nSoilOnlyHyd)%dat(1)) ! intent(in): [i4b] number of hydrology variables in the soil domain
    if (nSoilOnlyHyd>0) then ! check if computing soil hydrology
      if (local_ixGroundwater/=qbaseTopmodel .and. nGlce==0) then ! set baseflow fluxes to zero if the topmodel baseflow routine is not used
        call zeroBaseflowFluxes
      else ! compute the baseflow flux for topmodel-ish shallow groundwater or lateral flow for glacier debris
        call initialize_groundwatr; if(err/=0)then; return; endif
        call groundwatr(in_groundwatr,mpar_data,prog_data,flux_data,io_groundwatr,out_groundwatr)
        call finalize_groundwatr;   if(err/=0)then; return; endif
      end if
      call computBaseflowRunoff ! compute total baseflow from soil and runoff
    end if
  end associate

  ! *** CALCULATE FLUXES FOR THE DEEP AQUIFER ***
  associate(ixAqWat => indx_data%var(iLookINDEX%ixAqWat)%dat(1)) ! intent(in): [i4b] index of water storage in the aquifer
    if (ixAqWat/=integerMissing) then ! check if computing aquifer fluxes
      if (local_ixGroundwater==bigBucket) then ! compute fluxes for the big bucket
        call initialize_bigAquifer
        call bigAquifer(in_bigAquifer,mpar_data,diag_data,io_bigAquifer,out_bigAquifer)
        call finalize_bigAquifer; if(err/=0)then; return; endif
      else ! if no deep aquifer, then fluxes are zero
        call zeroAquiferFluxes
      end if ! end check aquifer model decision
    end if  ! if computing aquifer fluxes
  end associate

  call finalize_computFlux; if(err/=0)then; return; endif ! final operations to prep for end of routine

contains

 ! **** Subroutines that handle the absence of model features ****
 subroutine zeroGlacierFluxes
  ! set glacier ice fluxes to zero if no glacier ice layers
  associate(&
   scalarThroughfallRain     => flux_data%var(iLookFLUX%scalarThroughfallRain)%dat(1),       & ! intent(in):  [dp] rain that reaches the ground without ever touching the canopy (kg m-2 s-1)
   scalarCanopyLiqDrainage   => flux_data%var(iLookFLUX%scalarCanopyLiqDrainage)%dat(1),     & ! intent(in):  [dp] drainage of liquid water from the vegetation canopy (kg m-2 s-1)
   scalarSurfaceIceMelt      => flux_data%var(iLookFLUX%scalarSurfaceIceMelt)%dat(1),        & ! intent(out): [dp] liquid water flux at the top of the glacier ice layer (m s-1)
   scalarSurfaceIceMeltDeriv => deriv_data%var(iLookDERIV%scalarSurfaceIceMeltDeriv)%dat(1), & ! intent(out): [dp] derivative in liquid water flux at the top of the glacier ice layer (s-1)
   scalarGlceMelt            => flux_data%var(iLookFLUX%scalarGlceMelt)%dat(1),              & ! intent(out): [dp] glacier ice melt (m s-1)
   scalarGlacierMelt         => flux_data%var(iLookFLUX%scalarGlacierMelt)%dat(1)            ) ! intent(out): [dp] glacier ice melt plus snow and debris drainage (m s-1)
   if(nGlce==0)then ! no glacier ice layers
    scalarGlceMelt    = 0._rkind ! glacier ice melt (m s-1)
    scalarGlacierMelt = 0._rkind ! glacier rain + snow+debris melt (m s-1)
    scalarSurfaceIceMelt = 0._rkind
    scalarSurfaceIceMeltDeriv = 0._rkind
   else ! glacier ice layers, take from previous flux calculation
    scalarGlacierMelt = (scalarThroughfallRain + scalarCanopyLiqDrainage)/iden_water  & ! liquid flux from the canopy (m s-1)
                        + drainageMeltPond/iden_water & ! melt of the snow without a layer (m s-1)
                        - scalarGlceMelt ! (= above layer positive flux (0) - upward flux) save for glacier melt flow calculations, may be overwritten with addition of above domain fluxes
   endif
  end associate
 end subroutine zeroGlacierFluxes

 subroutine forcingNoLake
  ! define forcing for the beneath domains for the case of no lake layers
  associate(&
   scalarLakeDrainage        => flux_data%var(iLookFLUX%scalarLakeDrainage)%dat(1),          & ! intent(out): [dp] drainage from the lake profile (m s-1)
   scalarRainPlusMelt        => flux_data%var(iLookFLUX%scalarRainPlusMelt)%dat(1),          & ! intent(out): [dp] rain plus melt plus lake drainage (m s-1)
   scalarSurfaceRunoff       => flux_data%var(iLookFLUX%scalarSurfaceRunoff)%dat(1),         & ! intent(in):  [dp] surface runoff (m s-1)
   scalarSurfaceIceMelt      => flux_data%var(iLookFLUX%scalarSurfaceIceMelt)%dat(1),        & ! intent(out): [dp] liquid water flux at the top of the glacier ice layer (m s-1)
   scalarSurfaceIceMeltDeriv => deriv_data%var(iLookDERIV%scalarSurfaceIceMeltDeriv)%dat(1), & ! intent(out): [dp] derivative in liquid water flux at the top of the glacier ice layer (s-1)
   scalarGlceMelt            => flux_data%var(iLookFLUX%scalarGlceMelt)%dat(1),              & ! intent(out): [dp] glacier ice melt (m s-1)
   scalarGlacierMelt         => flux_data%var(iLookFLUX%scalarGlacierMelt)%dat(1)            ) ! intent(out): [dp] glacier ice melt plus snow and debris drainage (m s-1)
   if(nLake==0)then ! no lake layers, will have above_* variables from snow calculations (or forcingNoSnow if no snow)
    scalarLakeDrainage = 0._rkind ! drainage from the above layer (m s-1)
    scalarSurfaceRunoff = scalarRainPlusMelt ! all rain plus melt becomes surface runoff
    if(nGlce==0)then
     scalarSurfaceIceMelt = 0._rkind
     scalarSurfaceIceMeltDeriv = 0._rkind
    endif
   else ! lake layers, take from previous flux calculation
    scalarRainPlusMelt = scalarLakeDrainage                    ! drainage from the base of the lake
    if(nGlce>0) scalarGlacierMelt = scalarLakeDrainage - scalarGlceMelt + scalarSurfaceRunoff ! save for glacier melt flow calculations, may be overwritten with addition of below domain fluxes
   end if ! lake layers or not
  end associate
 end subroutine forcingNoLake

 subroutine forcingNoSnow
  ! define forcing for the beneath domains for the case of no snow layers
  ! NOTE: in case where nSnowOnlyHyd==0 AND snow layers exist, then scalarRainPlusMelt is taken from the previous flux evaluation
  associate(&
   scalarRainPlusMelt           => flux_data%var(iLookFLUX%scalarRainPlusMelt)%dat(1),      & ! intent(out): [dp] rain plus melt plus lake drainage (m s-1)
   scalarThroughfallRain        => flux_data%var(iLookFLUX%scalarThroughfallRain)%dat(1),   & ! intent(in):  [dp] rain that reaches the ground without ever touching the canopy (kg m-2 s-1)
   scalarCanopyLiqDrainage      => flux_data%var(iLookFLUX%scalarCanopyLiqDrainage)%dat(1), & ! intent(in):  [dp] drainage of liquid water from the vegetation canopy (kg m-2 s-1)
   scalarSnowDrainage           => flux_data%var(iLookFLUX%scalarSnowDrainage)%dat(1),      & ! intent(out): [dp] drainage from the snow layers (m s-1)
   scalarGlceMelt               => flux_data%var(iLookFLUX%scalarGlceMelt)%dat(1),          & ! intent(out): [dp] glacier ice melt (m s-1)
   scalarGlacierMelt            => flux_data%var(iLookFLUX%scalarGlacierMelt)%dat(1)        ) ! intent(out): [dp] glacier ice melt plus snow and debris drainage (m s-1)
   if (nSnow==0) then ! no snow layers
    scalarRainPlusMelt = (scalarThroughfallRain + scalarCanopyLiqDrainage)/iden_water  & ! liquid flux from the canopy (m s-1)
                        + drainageMeltPond/iden_water  ! melt of the snow without a layer (m s-1)
    scalarSnowDrainage = 0._rkind ! drainage from the snow layers (m s-1)
   else ! snow layers, take from previous flux calculation
    if(nGlce>0) scalarGlacierMelt = scalarSnowDrainage - scalarGlceMelt ! save for glacier melt flow calculations, may be overwritten with addition of above domain fluxes
   end if ! snow layers or not
  end associate
 end subroutine forcingNoSnow

 subroutine forcingNoSoil
  ! define forcing for the beneath domains for the case of no soil layers
  associate(&
   scalarRainPlusMelt           => flux_data%var(iLookFLUX%scalarRainPlusMelt)%dat(1),      & ! intent(out): [dp] rain plus melt plus lake drainage (m s-1)
   scalarSurfaceRunoff          => flux_data%var(iLookFLUX%scalarSurfaceRunoff)%dat(1),     & ! intent(in):  [dp] surface runoff (m s-1)
   scalarGlceMelt               => flux_data%var(iLookFLUX%scalarGlceMelt)%dat(1),          & ! intent(out): [dp] glacier ice melt (m s-1)
   scalarGlacierMelt            => flux_data%var(iLookFLUX%scalarGlacierMelt)%dat(1),       & ! intent(out): [dp] glacier ice melt plus snow and debris drainage (m s-1)
   scalarSoilControl            => diag_data%var(iLookDIAG%scalarSoilControl )%dat(1),      & ! intent(out): [dp] soil control on infiltration for derivative, zero or one
   scalarSoilDrainage           => flux_data%var(iLookFLUX%scalarSoilDrainage)%dat(1),      & ! intent(out): [dp] drainage from the soil profile (m s-1)
   iLayerLiqFluxSoil            => flux_data%var(iLookFLUX%iLayerLiqFluxSoil)%dat           ) ! intent(out): [dp(0:)] vertical liquid water flux at soil layer interfaces (-)
   if(nSoil==0)then ! no soil layers
    scalarSoilDrainage = 0._rkind ! drainage from the above layer (m s-1)
    scalarSoilControl  = 1._rkind ! infiltration not controlled by soil
    scalarSurfaceRunoff = scalarRainPlusMelt ! all rain plus melt becomes surface runoff
   else ! soil layers, take from previous flux calculation
    if(nGlce>0) scalarGlacierMelt = scalarSoilDrainage - scalarGlceMelt + scalarSurfaceRunoff ! save for glacier melt flow calculations, may be overwritten with addition of below domain fluxes
   end if
  end associate
 end subroutine forcingNoSoil

 subroutine zeroBaseflowFluxes
  ! set baseflow fluxes to zero if the topmodel baseflow routine is not used
  associate(&
   scalarSoilBaseflow           => flux_data%var(iLookFLUX%scalarSoilBaseflow)%dat(1), & ! intent(out): [dp] total baseflow from the soil profile (m s-1)
   scalarExfiltration           => flux_data%var(iLookFLUX%scalarExfiltration)%dat(1), & ! intent(out): [dp] exfiltration from the soil profile (m s-1)
   mLayerColumnOutflow          => flux_data%var(iLookFLUX%mLayerColumnOutflow)%dat,   & ! intent(out): [dp(:)] column outflow from each soil layer (m3 s-1)
   mLayerBaseflow               => flux_data%var(iLookFLUX%mLayerBaseflow)%dat         ) ! intent(out): [dp(:)] baseflow from each soil layer (m s-1)
   ! diagnostic variables in the data structures
   scalarExfiltration     = 0._rkind  ! exfiltration from the soil profile (m s-1)
   mLayerColumnOutflow(:) = 0._rkind  ! column outflow from each soil layer (m3 s-1)
   ! flux variables in the data structures
   mLayerBaseflow(:)      = 0._rkind  ! baseflow from each soil layer (m s-1)
   scalarSoilBaseflow     = 0._rkind  ! total baseflow from the soil profile (m s-1)
  end associate
 end subroutine zeroBaseflowFluxes

 subroutine computBaseflowRunoff
  ! compute total baseflow from the soil zone (needed for mass balance checks) and total runoff, before aquifer fluxes
  ! (Note: scalarSoilBaseflow is nonzero only if topmodel or have glacier debris layers)
  associate(&
   scalarSoilBaseflow           => flux_data%var(iLookFLUX%scalarSoilBaseflow)%dat(1),  & ! intent(out): [dp] total baseflow from the soil profile (m s-1)
   mLayerBaseflow               => flux_data%var(iLookFLUX%mLayerBaseflow)%dat,         & ! intent(out): [dp(:)] baseflow from each soil layer (m s-1)
   scalarTotalRunoff            => flux_data%var(iLookFLUX%scalarTotalRunoff)%dat(1),   & ! intent(out): [dp] total runoff (m s-1)
   scalarSurfaceRunoff          => flux_data%var(iLookFLUX%scalarSurfaceRunoff)%dat(1), & ! intent(out): [dp] surface runoff (m s-1)
   scalarSoilDrainage           => flux_data%var(iLookFLUX%scalarSoilDrainage)%dat(1),  & ! intent(out): [dp] drainage from the soil profile (m s-1)
   scalarGlceMelt               => flux_data%var(iLookFLUX%scalarGlceMelt)%dat(1),      & ! intent(out): [dp] glacier ice melt (m s-1)
   scalarGlacierMelt            => flux_data%var(iLookFLUX%scalarGlacierMelt)%dat(1)    ) ! intent(out): [dp] glacier ice melt plus snow and soil drainage (m s-1)
   ! baseflow from the soil zone
   scalarSoilBaseflow = sum(mLayerBaseflow) 
   ! compute total runoff
   scalarTotalRunoff = scalarSurfaceRunoff + scalarSoilBaseflow + scalarSoilDrainage - scalarGlceMelt ! total runoff (m s-1)
   if(nGlce>0) scalarGlacierMelt = scalarTotalRunoff ! if got here then needs debris outflow which is total runoff
  end associate
 end subroutine computBaseflowRunoff  

 subroutine zeroAquiferFluxes
  ! set aquifer fluxes to zero (if no aquifer exists)
  associate(&
   scalarAquiferTranspire      => flux_data%var(iLookFLUX%scalarAquiferTranspire)%dat(1), & ! intent(out): [dp] transpiration loss from the aquifer (m s-1
   scalarAquiferRecharge       => flux_data%var(iLookFLUX%scalarAquiferRecharge)%dat(1),  & ! intent(out): [dp] recharge to the aquifer (m s-1)
   scalarAquiferBaseflow       => flux_data%var(iLookFLUX%scalarAquiferBaseflow)%dat(1),  & ! intent(out): [dp] total baseflow from the aquifer (m s-1)
   dBaseflow_dAquifer          => deriv_data%var(iLookDERIV%dBaseflow_dAquifer)%dat(1)    ) ! intent(out): [dp(:)] derivative in baseflow flux w.r.t. aquifer storage (s-1)
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
   ixSpatialGroundwater         => model_decisions(iLookDECISIONS%spatial_gw)%iDecision, & ! intent(in):  [i4b] spatial representation of groundwater (local-column or single-basin)
   ixGroundwater                => model_decisions(iLookDECISIONS%groundwatr)%iDecision, & ! intent(in):  [i4b] groundwater parameterization
   iLayerLiqFluxSnLaGl          => flux_data%var(iLookFLUX%iLayerLiqFluxSnLaGl)%dat,     & ! intent(out): [dp(0:)] vertical liquid water flux at snow layer interfaces (-)
   iLayerLiqFluxSoil            => flux_data%var(iLookFLUX%iLayerLiqFluxSoil)%dat        ) ! intent(out): [dp(0:)] vertical liquid water flux at soil layer interfaces (-)

   numFluxCalls = numFluxCalls+1 ! increment the number of flux calls
   lake_frozen = .false.
   if(nLake>0 .and. mLayerTempTrial(nSnow+1)<=Tfreeze) lake_frozen= .true.

   ! modify the groundwater representation for this single-column implementation
   select case(ixSpatialGroundwater)
     case(singleBasin); local_ixGroundwater = noExplicit    ! force no explicit representation of groundwater at the local scale
     case(localColumn); local_ixGroundwater = ixGroundwater ! go with the specified decision
     case default; err=20; message=trim(message)//'unable to identify spatial representation of groundwater'; return
   end select ! end modify the groundwater representation for this single-column implementation

   ! initialize liquid water fluxes throughout the layer domains
   ! NOTE: used in the energy routines, which is called before the hydrology routines
   if (firstFluxCall) then
     iLayerLiqFluxSnLaGl(0:nLayers) = 0._rkind
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
   scalarCanopyNetNrgFlux       => flux_data%var(iLookFLUX%scalarCanopyNetNrgFlux)%dat(1), & ! intent(out): [dp] net energy flux for the vegetation canopy (W m-2)
   scalarCanopyNetLiqFlux       => flux_data%var(iLookFLUX%scalarCanopyNetLiqFlux)%dat(1), & ! intent(out): [dp] net liquid water flux for the vegetation canopy (kg m-2 s-1)
   canopyDepth                  => diag_data%var(iLookDIAG%scalarCanopyDepth)%dat(1),      & ! intent(in): [dp] canopy depth (m)
   nSnLaSoGlNrg                 => indx_data%var(iLookINDEX%nSnLaSoGlNrg)%dat(1),          & ! intent(in): [i4b] number of energy state variables in the layer domains
   ixSnLaSoGlNrg                => indx_data%var(iLookINDEX%ixSnLaSoGlNrg)%dat,            & ! intent(in): [i4b(:)] indices for energy states in the layer domains
   mLayerNrgFlux                => flux_data%var(iLookFLUX%mLayerNrgFlux)%dat              ) ! intent(out): [dp] net energy flux for each layer within the layer domains (J m-3 s-1)
   ! *** WRAP UP ***
   ! define model flux vector for the vegetation sub-domain
   if (ixCasNrg/=integerMissing) fluxVec(ixCasNrg) = scalarCanairNetNrgFlux/canopyDepth
   if (ixVegNrg/=integerMissing) fluxVec(ixVegNrg) = scalarCanopyNetNrgFlux/canopyDepth
   if (ixVegHyd/=integerMissing) fluxVec(ixVegHyd) = scalarCanopyNetLiqFlux   ! NOTE: solid fluxes are handled separately
   if (nSnLaSoGlNrg>0) then ! if necessary, populate the flux vector for energy
     do concurrent (iLayer=1:nLayers,ixSnLaSoGlNrg(iLayer)/=integerMissing)   ! loop through non-missing energy state variables in the layer domains
       fluxVec( ixSnLaSoGlNrg(iLayer) ) = mLayerNrgFlux(iLayer)
     end do
   end if
  end associate

  associate(&
   ixAqWat                      => indx_data%var(iLookINDEX%ixAqWat)%dat(1),               & ! intent(in): [i4b] index of water storage in the aquifer
   ixSnLaSoGlHyd                => indx_data%var(iLookINDEX%ixSnLaSoGlHyd)%dat,            & ! intent(in): [i4b(:)] indices for hydrology states in the layer domains
   nSnLaSoGlHyd                 => indx_data%var(iLookINDEX%nSnLaSoGlHyd)%dat(1),          & ! intent(in): [i4b] number of hydrology variables in the layer domains
   layerType                    => indx_data%var(iLookINDEX%layerType)%dat,                & ! intent(in): [i4b(:)] type of layer (iname_*)
   mLayerLiqFluxSnLaGl          => flux_data%var(iLookFLUX%mLayerLiqFluxSnLaGl)%dat,       & ! intent(out): [dp] net liquid water flux for each non-soil layer (s-1)
   mLayerLiqFluxSoil            => flux_data%var(iLookFLUX%mLayerLiqFluxSoil)%dat,         & ! intent(out): [dp] net liquid water flux for each soil layer (s-1)
   scalarAquiferTranspire       => flux_data%var(iLookFLUX%scalarAquiferTranspire)%dat(1), & ! intent(out): [dp] transpiration loss from the aquifer (m s-1)
   scalarAquiferRecharge        => flux_data%var(iLookFLUX%scalarAquiferRecharge)%dat(1),  & ! intent(out): [dp] recharge to the aquifer (m s-1)
   scalarAquiferBaseflow        => flux_data%var(iLookFLUX%scalarAquiferBaseflow)%dat(1)   ) ! intent(out): [dp] total baseflow from the aquifer (m s-1)
   ! populate the flux vector for hydrology
   ! NOTE: ixVolFracWat  and ixVolFracLiq can also include states in the soil domain, hence enable primary variable switching
   if (nSnLaSoGlHyd>0) then  ! check if any hydrology states exist
     do iLayer=1,nLayers     ! loop through non-missing energy state variables in the layer domains
       if (ixSnLaSoGlHyd(iLayer)/=integerMissing) then   ! check if a given hydrology state exists
         select case(layerType(iLayer))
           case(iname_snow,iname_lake,iname_glce); fluxVec(ixSnLaSoGlHyd(iLayer)) = mLayerLiqFluxSnLaGl(iLayer)
           case(iname_soil);                       fluxVec(ixSnLaSoGlHyd(iLayer)) = mLayerLiqFluxSoil(iLayer-nSnow-nLake)
           case default; err=20; message=trim(message)//'expect layerType to be iname_snow, iname_lake, iname_soil, or iname_glce'; return
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
   canopyDepth                  => diag_data%var(iLookDIAG%scalarCanopyDepth)%dat(1)   ) ! intent(in):  [dp]  canopy depth (m)

   dCanLiq_dTcanopy = dTheta_dTkCanopy*iden_water*canopyDepth     ! derivative in canopy liquid storage w.r.t. canopy temperature (kg m-2 K-1)
  end associate
  call in_vegNrgFlux % initialize(firstSubStep,firstFluxCall,computeVegFlux,checkLWBalance,&
                                  scalarCanairTempTrial,scalarCanopyTempTrial,mLayerTempTrial,scalarCanopyIceTrial,&
                                  scalarCanopyLiqTrial,forc_data,deriv_data)
 end subroutine initialize_vegNrgFlux

 subroutine finalize_vegNrgFlux
  call out_vegNrgFlux%finalize(flux_data,deriv_data,err,cmessage)
  ! error control
  if (err/=0) then; message=trim(message)//trim(cmessage); return; end if  ! check for errors
 end subroutine finalize_vegNrgFlux
 ! **** end vegNrgFlux ****

 ! **** snowLakeSoilGlceNrgFlux ****
 subroutine initialize_snowLakeSoilGlceNrgFlux
  call in_snowLakeSoilGlceNrgFlux%initialize(scalarSolution,firstFluxCall,mLayerTempTrial,flux_data,deriv_data)
  call io_snowLakeSoilGlceNrgFlux%initialize(deriv_data)
 end subroutine initialize_snowLakeSoilGlceNrgFlux

 subroutine finalize_snowLakeSoilGlceNrgFlux
  call io_snowLakeSoilGlceNrgFlux%finalize(deriv_data)
  call out_snowLakeSoilGlceNrgFlux%finalize(flux_data,deriv_data,err,cmessage)
  ! error control
  if (err/=0) then; message=trim(message)//trim(cmessage); return; end if
  associate(&
   mLayerNrgFlux                => flux_data%var(iLookFLUX%mLayerNrgFlux)%dat, & ! intent(out): [dp] net energy flux for each layer within the layer domains (J m-3 s-1)
   iLayerNrgFlux                => flux_data%var(iLookFLUX%iLayerNrgFlux)%dat, & ! intent(in):  [dp(0:)] vertical energy flux at the interface of layers
   mLayerDepth                  => prog_data%var(iLookPROG%mLayerDepth)%dat    ) ! intent(in):  [dp(:)]  depth of each layer in the layer domains (m)
   ! calculate net energy fluxes for each layer (J m-3 s-1)
   do iLayer=1,nLayers
     mLayerNrgFlux(iLayer) = -(iLayerNrgFlux(iLayer) - iLayerNrgFlux(iLayer-1))/mLayerDepth(iLayer)
   end do
  end associate
 end subroutine finalize_snowLakeSoilGlceNrgFlux
 ! **** end snowLakeSoilGlceNrgFlux ****

 ! **** vegLiqFlux ****
 subroutine initialize_vegLiqFlux
  call in_vegLiqFlux%initialize(computeVegFlux,scalarCanopyLiqTrial,flux_data)
 end subroutine initialize_vegLiqFlux
 
 subroutine finalize_vegLiqFlux
  call out_vegLiqFlux%finalize(flux_data,deriv_data,err,cmessage)
  ! error control
  if (err/=0) then; message=trim(message)//trim(cmessage); return; end if
  associate( &
   scalarThroughfallRain       => flux_data%var(iLookFLUX%scalarThroughfallRain)%dat(1),         & ! intent(in):  [dp] rain that reaches the ground without ever touching the canopy (kg m-2 s-1)
   scalarCanopyLiqDrainage     => flux_data%var(iLookFLUX%scalarCanopyLiqDrainage)%dat(1),       & ! intent(in):  [dp] drainage of liquid water from the vegetation canopy (kg m-2 s-1)
   scalarThroughfallRainDeriv  => deriv_data%var(iLookDERIV%scalarThroughfallRainDeriv  )%dat(1),& ! intent(in):  [dp] derivative in throughfall w.r.t. canopy liquid water
   scalarCanopyLiqDrainageDeriv=> deriv_data%var(iLookDERIV%scalarCanopyLiqDrainageDeriv)%dat(1),& ! intent(out): [dp] derivative in canopy drainage w.r.t. canopy liquid water
   scalarCanopyNetLiqFlux      => flux_data%var(iLookFLUX%scalarCanopyNetLiqFlux)%dat(1),        & ! intent(out): [dp] net liquid water flux for the vegetation canopy (kg m-2 s-1)
   scalarRainfall              => flux_data%var(iLookFLUX%scalarRainfall)%dat(1),                & ! intent(in):  [dp] rainfall rate (kg m-2 s-1)
   scalarCanopyEvaporation     => flux_data%var(iLookFLUX%scalarCanopyEvaporation)%dat(1),       & ! intent(in):  [dp] canopy evaporation/condensation (kg m-2 s-1)
   scalarCanopyLiqDeriv        => deriv_data%var(iLookDERIV%scalarCanopyLiqDeriv        )%dat(1) ) ! intent(in):  [dp] derivative in (throughfall + drainage) w.r.t. canopy liquid water
   ! calculate the net liquid water flux for the vegetation canopy
   scalarCanopyNetLiqFlux = scalarRainfall + scalarCanopyEvaporation - scalarThroughfallRain - scalarCanopyLiqDrainage
   ! calculate the total derivative in the downward liquid flux
   scalarCanopyLiqDeriv   = scalarThroughfallRainDeriv + scalarCanopyLiqDrainageDeriv
  end associate
 end subroutine finalize_vegLiqFlux
 ! **** end vegLiqFlux ****

  ! **** glceLiqFlux ****
 subroutine initialize_glceLiqFlux
  associate(&
   noThetaChange        => indx_data%var(iLookINDEX%noThetaChange)%dat(1)      ) ! intent(in): [int] number of layers with no change in total water content (bottom layers)
   surface_flux = 0._rkind ! no surface flux for glacier ice layers since impermeable
   bottom_flux = 0._rkind ! no bottom flux for glacier ice layers
   nStart = nSnow + nLake + nSoil
   call in_snowLakeGlceLiqFlux%initialize(nGlce-noThetaChange,nStart,.true.,.false.,surface_flux,bottom_flux,firstFluxCall,scalarSolution,mLayerVolFracLiqTrial)
   call io_snowLakeGlceLiqFlux%initialize(flux_data,deriv_data) ! only compute liquid water fluxes for top layers
  end associate
 end subroutine initialize_glceLiqFlux

 subroutine finalize_glceLiqFlux
  nStart = nSnow + nLake + nSoil
  call io_snowLakeGlceLiqFlux%finalize(flux_data,deriv_data) ! only compute liquid water fluxes for top layers
  call out_snowLakeGlceLiqFlux%finalize(err,cmessage) 
  ! error control
  if (err/=0) then; message=trim(message)//trim(cmessage); return; end if
  associate(&
   noThetaChange               => indx_data%var(iLookINDEX%noThetaChange)%dat(1),              & ! intent(in):    [int] number of layers with no change in total water content (bottom layers)
   mLayerLiqFluxSnLaGl         => flux_data%var(iLookFLUX%mLayerLiqFluxSnLaGl)%dat,            & ! intent(out):   [dp] net liquid water flux for each snow lake glce layer (s-1)
   iLayerLiqFluxSnLaGl         => flux_data%var(iLookFLUX%iLayerLiqFluxSnLaGl)%dat,            & ! intent(in):    [dp(0:)] vertical liquid water flux at snow lake glce layer interfaces (-)
   iLayerLiqFluxSnLaGlDeriv    => deriv_data%var(iLookDERIV%iLayerLiqFluxSnLaGlDeriv)%dat,     & ! intent(inout): [dp(:)] derivative in vertical liquid water flux at layer interfaces
   scalarSurfaceIceMelt        => flux_data%var(iLookFLUX%scalarSurfaceIceMelt)%dat(1),        & ! intent(out):   [dp] liquid water flux at the top of the glacier ice layer (m s-1)
   scalarSurfaceIceMeltDeriv   => deriv_data%var(iLookDERIV%scalarSurfaceIceMeltDeriv)%dat(1), & ! intent(out):   [dp] derivative in liquid water flux at the top of the glacier ice layer (s-1)
   mLayerDepth                 => prog_data%var(iLookPROG%mLayerDepth)%dat,                    & ! intent(in):    [dp(:)]  depth of each layer (m)
   scalarThroughfallRain       => flux_data%var(iLookFLUX%scalarThroughfallRain)%dat(1),       & ! intent(in):    [dp] rain that reaches the ground without ever touching the canopy (kg m-2 s-1)
   scalarCanopyLiqDrainage     => flux_data%var(iLookFLUX%scalarCanopyLiqDrainage)%dat(1),     & ! intent(in):    [dp] drainage of liquid water from the vegetation canopy (kg m-2 s-1)
   scalarGlceMelt              => flux_data%var(iLookFLUX%scalarGlceMelt)%dat(1),              & ! intent(out):   [dp] glacier ice melt (m s-1)
   scalarGlacierMelt           => flux_data%var(iLookFLUX%scalarGlacierMelt)%dat(1)            ) ! intent(out):   [dp] glacier ice melt plus snow and soil drainage (m s-1)
   ! iLayer surface will be later written over by snow domain or lake domain derivative (if no soil), save here
   scalarSurfaceIceMelt = iLayerLiqFluxSnLaGl(nStart) ! won't be used as this if there is soil
   scalarSurfaceIceMeltDeriv = iLayerLiqFluxSnLaGlDeriv(nStart) ! won't be used as this if there is soil
   ! calculate net liquid water fluxes for top layers of glacier ice layer only (s-1)
   do iLayer=1,nGlce-noThetaChange
     mLayerLiqFluxSnLaGl(iLayer+nStart) = -(iLayerLiqFluxSnLaGl(iLayer+nStart) - iLayerLiqFluxSnLaGl(iLayer-1+nStart))/mLayerDepth(iLayer+nStart)
   end do
   if(noThetaChange>0)then ! no water flux in lower glacier ice layers
     do iLayer=nGlce-noThetaChange+1,nGlce
       iLayerLiqFluxSnLaGl(iLayer+nStart) = 0._rkind
       iLayerLiqFluxSnLaGlDeriv(iLayer+nStart) = 0._rkind
       mLayerLiqFluxSnLaGl(iLayer+nStart) = 0._rkind
     end do
   end if
   ! compute melt from the glacier ice zone (all melt goes to top of glacier ice), make positive since will be negative as upward flux
   scalarGlceMelt = iLayerLiqFluxSnLaGl(nStart) ! glacier ice melt is the liquid water flux at the top of the glacier ice layer
   scalarGlacierMelt = (scalarThroughfallRain + scalarCanopyLiqDrainage)/iden_water  & ! liquid flux from the canopy (m s-1), since no canopy on glacier, this is just the rain
                       + drainageMeltPond/iden_water & ! melt of the snow without a layer (m s-1)
                       - scalarGlceMelt ! (= above layer positive flux (0) - upward flux) save for glacier melt flow calculations, may be overwritten with addition of above domain fluxes
  end associate
 end subroutine finalize_glceLiqFlux
 ! **** end glceLiqFlux ****

    ! **** frozen lakeLiqFlux ****
 subroutine initialize_frzlakeLiqFlux
  associate(&
   noThetaChange               => indx_data%var(iLookINDEX%noThetaChange)%dat(1),    & ! intent(in): [int] number of layers with no change in total water content (bottom layers)
   scalarGlceMelt              => flux_data%var(iLookFLUX%scalarGlceMelt)%dat(1)     ) ! intent(in): [dp]  glacier ice melt (m s-1)
   surface_flux = 0._rkind ! no surface flux for glacier ice layers since impermeable
   bottom_flux = 0._rkind ! no bottom flux for frozen lake layers
   nStart = nSnow
   call in_snowLakeGlceLiqFlux%initialize(nLake-noThetaChange,nStart,nGlce>0,.false.,surface_flux,bottom_flux,firstFluxCall,scalarSolution,mLayerVolFracLiqTrial)
   call io_snowLakeGlceLiqFlux%initialize(flux_data,deriv_data)
  end associate
 end subroutine initialize_frzlakeLiqFlux

 subroutine finalize_frzlakeLiqFlux
  nStart = nSnow
  call io_snowLakeGlceLiqFlux%finalize(flux_data,deriv_data)
  call out_snowLakeGlceLiqFlux%finalize(err,cmessage) 
  ! error control
  if (err/=0) then; message=trim(message)//trim(cmessage); return; end if
  associate(&
   noThetaChange               => indx_data%var(iLookINDEX%noThetaChange)%dat(1),              & ! intent(in): [int] number of layers with no change in total water content (bottom layers)
   scalarRainPlusMelt          => flux_data%var(iLookFLUX%scalarRainPlusMelt)%dat(1),          & ! intent(out): [dp] rain plus melt plus lake drainage (m s-1)
   mLayerLiqFluxSnLaGl         => flux_data%var(iLookFLUX%mLayerLiqFluxSnLaGl)%dat,            & ! intent(out): [dp] net liquid water flux for each snow layer (s-1)
   iLayerLiqFluxSnLaGl         => flux_data%var(iLookFLUX%iLayerLiqFluxSnLaGl)%dat,            & ! intent(in):  [dp(0:)] vertical liquid water flux at snow layer interfaces (-)
   iLayerLiqFluxSnLaGlDeriv    => deriv_data%var(iLookDERIV%iLayerLiqFluxSnLaGlDeriv)%dat,     & ! intent(out): [dp(:)] derivative in vertical liquid water flux at layer interfaces
   scalarSurfaceIceMelt        => flux_data%var(iLookFLUX%scalarSurfaceIceMelt)%dat(1),        & ! intent(out): [dp] liquid water flux at the top of the glacier ice layer (m s-1)
   scalarSurfaceIceMeltDeriv   => deriv_data%var(iLookDERIV%scalarSurfaceIceMeltDeriv)%dat(1), & ! intent(out): [dp] derivative in liquid water flux at the top of the glacier ice layer (s-1)
   mLayerDepth                 => prog_data%var(iLookPROG%mLayerDepth)%dat,                    & ! intent(in):  [dp(:)]  depth of each layer (m)
   scalarSurfaceRunoff         => flux_data%var(iLookFLUX%scalarSurfaceRunoff)%dat(1),         & ! intent(in):  [dp] surface runoff (m s-1)
   scalarLakeDrainage          => flux_data%var(iLookFLUX%scalarLakeDrainage)%dat(1),          & ! intent(out): [dp] drainage from the lake profile (m s-1)
   scalarGlceMelt              => flux_data%var(iLookFLUX%scalarGlceMelt)%dat(1),              & ! intent(in):  [dp]  glacier ice melt (m s-1)
   scalarGlacierMelt           => flux_data%var(iLookFLUX%scalarGlacierMelt)%dat(1)            ) ! intent(out): [dp] glacier ice melt plus snow and soil drainage (m s-1)
   ! iLayer surface will be later written over by snow domain, save here
   scalarSurfaceIceMelt = iLayerLiqFluxSnLaGl(nStart)
   scalarSurfaceIceMeltDeriv = iLayerLiqFluxSnLaGlDeriv(nStart)
   ! define forcing for the beneath domain
   scalarRainPlusMelt = iLayerLiqFluxSnLaGl(nLake+nStart) ! drainage from the base of the lake
   ! calculate net liquid water fluxes for each ice layer (s-1)
   do iLayer=1,nLake-noThetaChange
     mLayerLiqFluxSnLaGl(iLayer+nStart) = -(iLayerLiqFluxSnLaGl(iLayer+nStart) - iLayerLiqFluxSnLaGl(iLayer-1+nStart))/mLayerDepth(iLayer+nStart)
   end do
   scalarLakeDrainage = 0._rkind ! no drainage from solid frozen lake
   if(nGlce>0) scalarGlacierMelt = scalarLakeDrainage + scalarSurfaceRunoff - scalarGlceMelt ! save for glacier melt flow calculations, may be overwritten with addition of below domain fluxes
  end associate
 end subroutine finalize_frzlakeLiqFlux

 ! **** unfrozen lakeLiqFlux ****
 subroutine initialize_lakeLiqFlux
  associate(&
   noThetaChange               => indx_data%var(iLookINDEX%noThetaChange)%dat(1),   & ! intent(in): [int] number of layers with no change in total water content (bottom layers)
   scalarSnowfall              => flux_data%var(iLookFLUX%scalarSnowfall)%dat(1),   & ! intent(in): [dp] computed snowfall rate (kg m-2 s-1)
   scalarRainfall              => flux_data%var(iLookFLUX%scalarRainfall)%dat(1),   & ! intent(in): [dp] computed rainfall rate (kg m-2 s-1)
   scalarGlceMelt              => flux_data%var(iLookFLUX%scalarGlceMelt)%dat(1)    ) ! intent(in): [dp]  glacier ice melt (m s-1)
   ! This should include all liquid and all solid fluxes
   surface_flux = (scalarSnowfall + scalarRainfall)/iden_water 
   !surface_flux = (scalarThroughfallRain + scalarCanopyLiqDrainage)/iden_water  & ! liquid flux from the canopy (m s-1), since no canopy on lake, this is just the rain + snow which melts immediately
   !                    + drainageMeltPond/iden_water
   bottom_flux = 0._rkind ! bottom flux for lake layers (m s-1)
   !if(nSoil==0) bottom_flux = scalarGlceMelt ! leave this here in case want to couple with glacier ice melt, or soil, or aquifer
   nLake_frz = 0
   if(lake_frozen) nLake_frz = nLake-noThetaChange
   nStart = nSnow + nLake_frz
   call in_snowLakeGlceLiqFlux%initialize(nLake-nLake_frz,nStart,nGlce>0,.false.,surface_flux,bottom_flux,firstFluxCall,scalarSolution,mLayerVolFracLiqTrial)
   call io_snowLakeGlceLiqFlux%initialize(flux_data,deriv_data)
  end associate
 end subroutine initialize_lakeLiqFlux

 subroutine finalize_lakeLiqFlux
  nStart = nSnow
  call io_snowLakeGlceLiqFlux%finalize(flux_data,deriv_data)
  call out_snowLakeGlceLiqFlux%finalize(err,cmessage) 
  ! error control
  if (err/=0) then; message=trim(message)//trim(cmessage); return; end if
  associate(&
   scalarRainPlusMelt          => flux_data%var(iLookFLUX%scalarRainPlusMelt)%dat(1),          & ! intent(out): [dp] rain plus melt plus lake drainage (m s-1)
   mLayerLiqFluxSnLaGl         => flux_data%var(iLookFLUX%mLayerLiqFluxSnLaGl)%dat,            & ! intent(out): [dp] net liquid water flux for each snow layer (s-1)
   iLayerLiqFluxSnLaGl         => flux_data%var(iLookFLUX%iLayerLiqFluxSnLaGl)%dat,            & ! intent(in):  [dp(0:)] vertical liquid water flux at snow layer interfaces (-)
   iLayerLiqFluxSnLaGlDeriv    => deriv_data%var(iLookDERIV%iLayerLiqFluxSnLaGlDeriv)%dat,     & ! intent(out): [dp(:)] derivative in vertical liquid water flux at layer interfaces
   scalarSurfaceIceMelt        => flux_data%var(iLookFLUX%scalarSurfaceIceMelt)%dat(1),        & ! intent(out): [dp] liquid water flux at the top of the glacier ice layer (m s-1)
   scalarSurfaceIceMeltDeriv   => deriv_data%var(iLookDERIV%scalarSurfaceIceMeltDeriv)%dat(1), & ! intent(out): [dp] derivative in liquid water flux at the top of the glacier ice layer (s-1)
   mLayerDepth                 => prog_data%var(iLookPROG%mLayerDepth)%dat,                    & ! intent(in):  [dp(:)]  depth of each layer (m)
   scalarSurfaceRunoff         => flux_data%var(iLookFLUX%scalarSurfaceRunoff)%dat(1),         & ! intent(in):  [dp] surface runoff (m s-1)
   scalarLakeDrainage          => flux_data%var(iLookFLUX%scalarLakeDrainage)%dat(1),          & ! intent(out): [dp] drainage from the lake profile (m s-1)
   scalarGlceMelt              => flux_data%var(iLookFLUX%scalarGlceMelt)%dat(1),              & ! intent(in):  [dp]  glacier ice melt (m s-1)
   scalarGlacierMelt           => flux_data%var(iLookFLUX%scalarGlacierMelt)%dat(1)            ) ! intent(out): [dp] glacier ice melt plus snow and soil drainage (m s-1)
   ! define forcing for the beneath domain
   scalarRainPlusMelt = iLayerLiqFluxSnLaGl(nLake+nStart) ! drainage from the base of the lake
   ! calculate net liquid water fluxes for each snow layer (s-1)
   do iLayer=1,nLake
     mLayerLiqFluxSnLaGl(iLayer+nStart) = -(iLayerLiqFluxSnLaGl(iLayer+nStart) - iLayerLiqFluxSnLaGl(iLayer-1+nStart))/mLayerDepth(iLayer+nStart)
   end do
   ! compute drainage from the lake zone (needed for mass balance checks)
   scalarLakeDrainage = iLayerLiqFluxSnLaGl(nLake+nStart)
   ! ice melt is 0 if not frozen
   if(.not.lake_frozen)then
     scalarSurfaceIceMelt = 0._rkind
     scalarSurfaceIceMeltDeriv = 0._rkind
   endif
   if(nGlce>0) scalarGlacierMelt = scalarLakeDrainage + scalarSurfaceRunoff - scalarGlceMelt ! save for glacier melt flow calculations, may be overwritten with addition of below domain fluxes
  end associate
 end subroutine finalize_lakeLiqFlux

 ! **** snowLiqFlux ****
 subroutine initialize_snowLiqFlux
  associate(&
   scalarThroughfallRain       => flux_data%var(iLookFLUX%scalarThroughfallRain)%dat(1),   & ! intent(in): [dp] rain that reaches the ground without ever touching the canopy (kg m-2 s-1)
   scalarCanopyLiqDrainage     => flux_data%var(iLookFLUX%scalarCanopyLiqDrainage)%dat(1), & ! intent(in): [dp] drainage of liquid water from the vegetation canopy (kg m-2 s-1)
   scalarGlceMelt              => flux_data%var(iLookFLUX%scalarGlceMelt)%dat(1)           ) ! intent(in): [dp]  glacier ice melt (m s-1)
   surface_flux = (scalarThroughfallRain + scalarCanopyLiqDrainage)/iden_water
   bottom_flux = 0._rkind ! bottom flux for snow layers (m s-1)
   !if (nLake==0 .and. nSoil==0) bottom_flux = scalarGlceMelt ! leave this here in case want to couple with glacier ice melt for slush layer, will change derivatives
   nStart = 0
   call in_snowLakeGlceLiqFlux%initialize(nSnow,nStart,nGlce>0,.true.,surface_flux,bottom_flux,firstFluxCall,scalarSolution,mLayerVolFracLiqTrial)
   call io_snowLakeGlceLiqFlux%initialize(flux_data,deriv_data)
  end associate
 end subroutine initialize_snowLiqFlux

 subroutine finalize_snowLiqFlux
  nStart = 0
  call io_snowLakeGlceLiqFlux%finalize(flux_data,deriv_data)
  call out_snowLakeGlceLiqFlux%finalize(err,cmessage) 
  ! error control
  if (err/=0) then; message=trim(message)//trim(cmessage); return; end if
  associate(&
   scalarRainPlusMelt          => flux_data%var(iLookFLUX%scalarRainPlusMelt)%dat(1), & ! intent(out): [dp] rain plus melt plus lake drainage (m s-1)
   mLayerLiqFluxSnLaGl         => flux_data%var(iLookFLUX%mLayerLiqFluxSnLaGl)%dat,   & ! intent(out): [dp] net liquid water flux for each snow layer (s-1)
   iLayerLiqFluxSnLaGl         => flux_data%var(iLookFLUX%iLayerLiqFluxSnLaGl)%dat,   & ! intent(in):  [dp(0:)] vertical liquid water flux at snow layer interfaces (-)
   mLayerDepth                 => prog_data%var(iLookPROG%mLayerDepth)%dat,           & ! intent(in):  [dp(:)]  depth of each layer (m)
   scalarSnowDrainage          => flux_data%var(iLookFLUX%scalarSnowDrainage)%dat(1), & ! intent(out): [dp] drainage from the snow profile (m s-1)
   scalarGlceMelt              => flux_data%var(iLookFLUX%scalarGlceMelt)%dat(1),     & ! intent(in):  [dp]  glacier ice melt (m s-1)
   scalarGlacierMelt           => flux_data%var(iLookFLUX%scalarGlacierMelt)%dat(1)   ) ! intent(out): [dp] glacier ice melt plus snow and soil drainage (m s-1)
   ! define forcing for the beneath domain
   scalarRainPlusMelt = iLayerLiqFluxSnLaGl(nSnow+nStart) ! drainage from the base of the snowpack
   ! calculate net liquid water fluxes for each snow layer (s-1)
   do iLayer=1,nSnow
     mLayerLiqFluxSnLaGl(iLayer+nStart) = -(iLayerLiqFluxSnLaGl(iLayer+nStart) - iLayerLiqFluxSnLaGl(iLayer-1+nStart))/mLayerDepth(iLayer+nStart)
   end do
   ! compute drainage from the snow zone (needed for mass balance checks)
   scalarSnowDrainage = iLayerLiqFluxSnLaGl(nSnow+nStart)
   if(nGlce>0) scalarGlacierMelt = scalarSnowDrainage - scalarGlceMelt ! save for glacier melt flow calculations, may be overwritten with addition of below domain fluxes
  end associate
 end subroutine finalize_snowLiqFlux
 ! **** end snowLiqFlux ****

 ! **** soilLiqFlux ****
 subroutine initialize_soilLiqFlux
  call in_soilLiqFlux%initialize(nSnow,nLake,nSoil,firstSplitOper,scalarSolution,firstFluxCall,scalarAquiferStorageTrial,&
                                mLayerTempTrial,mLayerMatricHeadTrial,mLayerMatricHeadLiqTrial,mLayerVolFracLiqTrial,mLayerVolFracIceTrial,&
                                flux_data,deriv_data)
  call io_soilLiqFlux%initialize(nSoil,dHydCond_dMatric,flux_data,diag_data,deriv_data)
 end subroutine initialize_soilLiqFlux

 subroutine finalize_soilLiqFlux
  nStart = nSnow + nLake
  call io_soilLiqFlux%finalize(nSoil,dHydCond_dMatric,flux_data,diag_data,deriv_data)
  call out_soilLiqFlux%finalize(err,cmessage)
  ! error control
  if (err/=0) then; message=trim(message)//trim(cmessage); return; end if
  associate(&
   mLayerLiqFluxSoil           => flux_data%var(iLookFLUX%mLayerLiqFluxSoil)%dat,       & ! intent(out):   [dp] net liquid water flux for each soil layer (s-1)
   iLayerLiqFluxSoil           => flux_data%var(iLookFLUX%iLayerLiqFluxSoil)%dat,       & ! intent(in):    [dp(0:)] vertical liquid water flux at soil layer interfaces (-)
   mLayerLiqFluxSnLaGl         => flux_data%var(iLookFLUX%mLayerLiqFluxSnLaGl)%dat,     & ! intent(inout): [dp] net liquid water flux for each layer, 0 in soil (s-1)
   iLayerLiqFluxSnLaGl         => flux_data%var(iLookFLUX%iLayerLiqFluxSnLaGl)%dat,     & ! intent(inout): [dp(0:)] vertical liquid water flux at layer interfaces, 0 in soil (-) 
   mLayerDepth                 => prog_data%var(iLookPROG%mLayerDepth)%dat,             & ! intent(in):    [dp(:)]  depth of each layer in the sub-domain (m)
   scalarSurfaceRunoff         => flux_data%var(iLookFLUX%scalarSurfaceRunoff)%dat(1),  & ! intent(in):    [dp] surface runoff (m s-1)
   scalarSoilDrainage          => flux_data%var(iLookFLUX%scalarSoilDrainage)%dat(1),   & ! intent(in):    [dp] drainage from the soil profile (m s-1)
   scalarGlceMelt              => flux_data%var(iLookFLUX%scalarGlceMelt)%dat(1),       & ! intent(in):  [dp]  glacier ice melt (m s-1)
   scalarGlacierMelt           => flux_data%var(iLookFLUX%scalarGlacierMelt)%dat(1)     ) ! intent(out):   [dp] glacier ice melt plus snow and soil drainage (m s-1)
   ! calculate net liquid water fluxes for each soil layer (s-1)
   if (nStart==0) iLayerLiqFluxSnLaGl(0) = 0._rkind ! then 0 layer is top of soil, iLayerLiqFluxSnLaGl does not exist in soil
   do iLayer=1,nSoil
     if(iLayer/=nSoil) iLayerLiqFluxSnLaGl(iLayer+nStart) = realMissing ! iLayerLiqFluxSnLaGl does not exist in soil but could exist at the bottom of the soil domain
     mLayerLiqFluxSnLaGl(iLayer+nStart) = realMissing ! iLayerLiqFluxSnLaGl does not exist in soil
     mLayerLiqFluxSoil(iLayer) = -(iLayerLiqFluxSoil(iLayer) - iLayerLiqFluxSoil(iLayer-1))/mLayerDepth(iLayer+nStart)
   end do
   if(nGlce==0) iLayerLiqFluxSnLaGl(nSoil+nStart) = realMissing ! if nothing below the soil domain, then does not exist
   ! compute drainage from the soil zone (needed for mass balance checks and in aquifer recharge)
   scalarSoilDrainage = iLayerLiqFluxSoil(nSoil)
   if(nGlce>0) scalarGlacierMelt = scalarSoilDrainage + scalarSurfaceRunoff - scalarGlceMelt ! save for glacier melt flow calculations, may be overwritten with addition of below domain fluxes
  end associate

  associate(&
   ixBcUpper                   => model_decisions(iLookDECISIONS%bcUpprSoiH)%iDecision,    & ! index defining the type of boundary conditions
   dq_dHydStateAbove           => deriv_data%var(iLookDERIV%dq_dHydStateAbove)%dat,        & ! intent(out): [dp(:)] change in flux at layer interfaces w.r.t. states in the layer above
   dq_dHydStateBelow           => deriv_data%var(iLookDERIV%dq_dHydStateBelow)%dat,        & ! intent(out): [dp(:)] change in flux at layer interfaces w.r.t. states in the layer below
   dq_dHydStateLayerSurfVec    => deriv_data%var(iLookDERIV%dq_dHydStateLayerSurfVec)%dat, & ! intent(out): [dp(:)] change in the flux in soil surface interface w.r.t. state variables in layers
   dq_dNrgStateAbove           => deriv_data%var(iLookDERIV%dq_dNrgStateAbove)%dat,        & ! intent(out): [dp(:)] change in flux at layer interfaces w.r.t. states in the layer above
   dq_dNrgStateBelow           => deriv_data%var(iLookDERIV%dq_dNrgStateBelow)%dat,        & ! intent(out): [dp(:)] change in flux at layer interfaces w.r.t. states in the layer below
   dq_dNrgStateLayerSurfVec    => deriv_data%var(iLookDERIV%dq_dNrgStateLayerSurfVec)%dat, & ! intent(out): [dp(:)] change in the flux in soil surface interface w.r.t. state variables in layers
   mLayerdTheta_dTk            => deriv_data%var(iLookDERIV%mLayerdTheta_dTk)%dat(nStart+1:nStart+nSoil), & ! intent(in): [dp(:)]  derivatives in volumetric liquid water content w.r.t. temperature
   iLayerLiqFluxSnLaGlDeriv    => deriv_data%var(iLookDERIV%iLayerLiqFluxSnLaGlDeriv)%dat, & ! intent(inout): [dp(:)]  derivative in vertical liquid water flux at layer interfaces
   dPsiLiq_dPsi0               => deriv_data%var(iLookDERIV%dPsiLiq_dPsi0)%dat             ) ! intent(in):    [dp(:)] derivative in liquid water matric pot w.r.t. the total water matric pot (-)
   ! expand derivatives to the total water matric potential
   ! NOTE: arrays are offset because computing derivatives in interface fluxes, at the top and bottom of the layer respectively
   dq_dHydStateAbove(1:nSoil)   = dq_dHydStateAbove(1:nSoil)  *dPsiLiq_dPsi0(1:nSoil)
   dq_dHydStateBelow(0:nSoil-1) = dq_dHydStateBelow(0:nSoil-1)*dPsiLiq_dPsi0(1:nSoil)
   if (ixBcUpper==prescribedHead) dq_dHydStateLayerSurfVec(1) = dq_dHydStateLayerSurfVec(1)*dPsiLiq_dPsi0(1)
  ! iLayerLiqFluxSnLaGlDeriv does not exist in soil but could exist at the bottom of the soil domain
   do iLayer=1,nSoil-1
     iLayerLiqFluxSnLaGlDeriv(iLayer+nStart) = realMissing
   end do
   if(nGlce==0) iLayerLiqFluxSnLaGlDeriv(nSoil+nStart) = realMissing ! if nothing below the soil domain, then does not exist
  end associate
 end subroutine finalize_soilLiqFlux
 ! **** end soilLiqFlux ****

 ! **** groundwatr ****
 subroutine initialize_groundwatr
  ! check the derivative matrix is sized appropriately
  if (size(dBaseflow_dWat,1)/=nSoil .or. size(dBaseflow_dWat,2)/=nSoil .or. size(dBaseflow_dTk,1)/=nSoil .or. size(dBaseflow_dTk,2)/=nSoil) then
    message=trim(message)//'expect dBaseflow_dWat and dBaseflow_dTk to be nSoil x nSoil'
    err=20; return
  end if
  call in_groundwatr%initialize(nSnow,nLake,nSoil,nGlce,firstFluxCall,mLayerVolFracLiqTrial,mLayerVolFracIceTrial,deriv_data)
  call io_groundwatr%initialize(ixSaturation)
 end subroutine initialize_groundwatr

 subroutine finalize_groundwatr
  call io_groundwatr%finalize(ixSaturation)
  call out_groundwatr%finalize(dBaseflow_dWat,dBaseflow_dTk,flux_data,err,cmessage)
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
  ! error control
  if (err/=0) then; message=trim(message)//trim(cmessage); return; end if
  associate(&
   scalarTotalRunoff           => flux_data%var(iLookFLUX%scalarTotalRunoff)%dat(1)     ,&  ! intent(out): [dp] total runoff (m s-1)
   scalarSurfaceRunoff         => flux_data%var(iLookFLUX%scalarSurfaceRunoff)%dat(1)   ,&  ! intent(out): [dp] surface runoff (m s-1)
   scalarAquiferBaseflow       => flux_data%var(iLookFLUX%scalarAquiferBaseflow)%dat(1)  )  ! intent(out): [dp] total baseflow from the aquifer (m s-1)
   ! compute total runoff (overwrite previously calculated value before considering aquifer)
   scalarTotalRunoff = scalarSurfaceRunoff + scalarAquiferBaseflow     
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
                      ixTop,ixBot,                        & ! intent(in):  top and bottom defining desired layers
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
  real(rkind),intent(in)         :: dt                        ! length of the time step (seconds)
  integer(i4b),intent(in)        :: ixTop,ixBot               ! top and bottom defining desired layers
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
  do iLayer=1,size(mLayerMatricHead)
    if (iLayer>=ixTop .and. iLayer<=ixBot) then
      ! compute the derivative for the compressibility term (m-1), no volume expansion for total water
      dCompress_dPsi(iLayer) = specificStorage*(mLayerVolFracLiqTrial(iLayer) + mLayerVolFracIceTrial(iLayer))/theta_sat(iLayer)
      ! compute the compressibility term (-) per second
      compress(iLayer) = (mLayerMatricHeadTrial(iLayer) - mLayerMatricHead(iLayer))*dCompress_dPsi(iLayer)/dt
    end if
  end do
end subroutine soilCmpres

! **********************************************************************************************************
! public subroutine soilCmpres: compute soil compressibility (-) and its derivative w.r.t matric head (m-1)
! **********************************************************************************************************
subroutine soilCmpresPrime(&
                          ! input:
                          ixTop,ixBot,                        & ! intent(in):  top and bottom defining desired layers
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
  integer(i4b),intent(in)           :: ixTop,ixBot              ! top and bottom defining desired layers
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
  do iLayer=1,size(mLayerMatricHeadPrime)
    if (iLayer>=ixTop .and. iLayer<=ixBot) then
      ! compute the derivative for the compressibility term (m-1), no volume expansion for total water
      dCompress_dPsi(iLayer) = specificStorage*(mLayerVolFracLiqTrial(iLayer) + mLayerVolFracIceTrial(iLayer))/theta_sat(iLayer)
      ! compute the compressibility term (-) instantaneously
      compress(iLayer) = mLayerMatricHeadPrime(iLayer) * dCompress_dPsi(iLayer)
    end if
  end do
end subroutine soilCmpresPrime

end module computFlux_module

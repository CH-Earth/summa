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

module computJacobWithPrime_module

! data types
USE nr_type

! derived types to define the data structures
USE data_types,only:&
                    var_ilength,  & ! data vector with variable length dimension (i4b)
                    var_dlength,  & ! data vector with variable length dimension (rkind)
                    model_options   ! defines the model decisions

! named variables for structure elements
USE var_lookup,only:iLookDECISIONS  ! named variables for elements of the decision structure
USE var_lookup,only:iLookPARAM      ! named variables for structure elements
USE var_lookup,only:iLookPROG       ! named variables for structure elements
USE var_lookup,only:iLookINDEX      ! named variables for structure elements
USE var_lookup,only:iLookDIAG       ! named variables for structure elements
USE var_lookup,only:iLookDERIV      ! named variables for structure elements

! access the global print flag
USE globalData,only:globalPrintFlag

! access missing values
USE globalData,only:integerMissing  ! missing integer
USE globalData,only:realMissing     ! missing real number

! named variables to describe the state variable type
USE globalData,only:iname_watLayer   ! named variable defining the total water state variable for layers
USE globalData,only:maxVolIceContent ! snow maximum volumetric ice content to store water (-)

! access named variables to describe the form and structure of the matrices used in the numerical solver
USE globalData,only: kl             ! number of sub-diagonal bands, assume kl>=4
USE globalData,only: ixDiag         ! index for the diagonal band
USE globalData,only: nBands         ! length of the leading dimension of the band diagonal matrix
USE globalData,only: ixFullMatrix   ! named variable for the full Jacobian matrix
USE globalData,only: ixBandMatrix   ! named variable for the band diagonal matrix
USE globalData,only: iJac1          ! first layer of the Jacobian to print
USE globalData,only: iJac2          ! last layer of the Jacobian to print

! constants
USE multiconst,only:&
                    LH_fus,       & ! latent heat of fusion                (J kg-1)
                    iden_water      ! intrinsic density of liquid water    (kg m-3)

! look-up values for the choice of groundwater parameterization
USE mDecisions_module,only:       &
 qbaseTopmodel,                   & ! TOPMODEL-ish baseflow parameterization
 bigBucket,                       & ! a big bucket (lumped aquifer model)
 noExplicit                         ! no explicit groundwater parameterization

! look-up values for the choice of variable in energy equations (BE residual or IDA state variable)
USE mDecisions_module,only:       &
 closedForm,                      & ! use temperature with closed form heat capacity
 enthalpyForm,                    & ! use enthalpy with soil temperature-enthalpy lookup tables
 enthalpyFormAN                     ! use enthalpy with soil temperature-enthalpy analytical solution

implicit none
private
public::computJacobWithPrime
public::computJacob4ida
contains


! **********************************************************************************************************
! public subroutine computJacobWithPrime: compute the Jacobian matrix
! **********************************************************************************************************
subroutine computJacobWithPrime(&
                      ! input: model control
                      cj,                         & ! intent(in):    this scalar changes whenever the step size or method order changes
                      dt,                         & ! intent(in):    length of the time step (seconds)
                      nSnow,                      & ! intent(in):    number of snow layers
                      nLake,                      & ! intent(in):    number of lake layers
                      nSoil,                      & ! intent(in):    number of soil layers
                      nGlce,                      & ! intent(in):    number of glacier ice layers
                      nLayers,                    & ! intent(in):    total number of layers
                      computeVegFlux,             & ! intent(in):    flag to indicate if we need to compute fluxes over vegetation
                      computeBaseflow,            & ! intent(in):    flag to indicate if we need to compute baseflow
                      ixMatrix,                   & ! intent(in):    form of the Jacobian matrix
                      specificStorage,            & ! intent(in):    specific storage coefficient (m-1)
                      theta_sat,                  & ! intent(in):    soil porosity (-)
                      enthalpyStateVec,           & ! intent(in):    flag if enthalpy is state variable
                      ! input: data structures
                      indx_data,                  & ! intent(in):    index data
                      prog_data,                  & ! intent(in):    model prognostic variables for a local HRU
                      diag_data,                  & ! intent(in):    model diagnostic variables for a local HRU
                      deriv_data,                 & ! intent(in):    derivatives in model fluxes w.r.t. relevant state variables
                      dBaseflow_dWat,             & ! intent(in):    derivative in baseflow w.r.t. soil water characteristic
                      dBaseflow_dTk,              & ! intent(in):    derivative in baseflow w.r.t. temperature (m s-1 K-1)
                      ! input: state variables
                      mLayerTempPrime,            & ! intent(in):    vector of derivative value for layer temperature (K)
                      mLayerMatricHeadPrime,      & ! intent(in):    vector of derivative value for layer matric head
                      mLayerVolFracWatPrime,      & ! intent(in):    vector of derivative value for layer water volume fraction
                      scalarCanopyTempPrime,      & ! intent(in):    derivative value for temperature of the vegetation canopy (K)
                      scalarCanopyWatPrime,       & ! intent(in):    derivative value for water content of the vegetation canopy
                      ! input-output: Jacobian and its diagonal
                      dMat0,                      & ! intent(in):    diagonal of the Jacobian matrix excluding fluxes, not depending on the state vector
                      aJac,                       & ! intent(out):   Jacobian matrix
                      ! output: error control
                      err,message)                  ! intent(out):   error code and error message
  ! -----------------------------------------------------------------------------------------------------------------
  ! provide access to subroutines
  use computJacob_module,only:fluxJacAdd
  use computJacob_module,only:ixInd
  ! -----------------------------------------------------------------------------------------------------------------
  implicit none
  ! input: model control
  real(rkind),intent(in)               :: cj                         ! this scalar changes whenever the step size or method order changes
  real(rkind),intent(in)               :: dt                         ! length of the time step (seconds)
  integer(i4b),intent(in)              :: nSnow                      ! number of snow layers
  integer(i4b),intent(in)              :: nLake                      ! number of lake layers
  integer(i4b),intent(in)              :: nSoil                      ! number of soil layers
  integer(i4b),intent(in)              :: nGlce                      ! number of glacier ice layers
  integer(i4b),intent(in)              :: nLayers                    ! total number of layers in the layer domains
  logical(lgt),intent(in)              :: computeVegFlux             ! flag to indicate if computing fluxes over vegetation
  logical(lgt),intent(in)              :: computeBaseflow            ! flag to indicate if computing baseflow
  integer(i4b),intent(in)              :: ixMatrix                   ! form of the Jacobian matrix
  real(rkind),intent(in)               :: specificStorage            ! specific storage coefficient (m-1)
  real(rkind),intent(in)               :: theta_sat(:)               ! soil porosity (-)
  logical(lgt),intent(in)              :: enthalpyStateVec           ! flag if enthalpy is state variable
  ! input: data structures
  type(var_ilength),intent(in)         :: indx_data                  ! indices defining model states and layers
  type(var_dlength),intent(in)         :: prog_data                  ! prognostic variables for a local HRU
  type(var_dlength),intent(in)         :: diag_data                  ! diagnostic variables for a local HRU
  type(var_dlength),intent(in)         :: deriv_data                 ! derivatives in model fluxes w.r.t. relevant state variables
  real(rkind),intent(in)               :: dBaseflow_dWat(:,:)        ! derivative in baseflow w.r.t. soil water characteristic
  real(rkind),intent(in)               :: dBaseflow_dTk(:,:)         ! derivative in baseflow w.r.t. temperature (m s-1 K-1)
  ! input: state variables
  real(rkind),intent(in)               :: mLayerTempPrime(:)         ! vector of derivative value for layer temperature
  real(rkind),intent(in)               :: mLayerMatricHeadPrime(:)   ! vector of derivative value for layer matric head
  real(rkind),intent(in)               :: mLayerVolFracWatPrime(:)   ! vector of derivative value for layer water volume fraction
  real(rkind),intent(in)               :: scalarCanopyTempPrime      ! derivative value for temperature of the vegetation canopy (K)
  real(rkind),intent(in)               :: scalarCanopyWatPrime       ! derivative value for water content of the vegetation canopy
  ! input-output: Jacobian and its diagonal
  real(rkind),intent(in)               :: dMat0(:)                   ! diagonal of the Jacobian matrix excluding fluxes, not depending on the state vector
  real(rkind),intent(out)              :: aJac(:,:)                  ! Jacobian matrix
  ! output variables
  integer(i4b),intent(out)             :: err                        ! error code
  character(*),intent(out)             :: message                    ! error message
  ! --------------------------------------------------------------
  ! * local variables
  ! --------------------------------------------------------------
  real(rkind),allocatable              :: dMat(:)                    ! diagonal of the Jacobian matrix excluding fluxes, depending on the state vector
  ! indices of model state variables
  integer(i4b)                         :: nrgState                   ! energy state variable
  integer(i4b)                         :: watState                   ! hydrology state variable
  integer(i4b)                         :: nState                     ! number of state variables
  integer(i4b),allocatable             :: nrgRows(:)                 ! indices of rows for energy column
  integer(i4b),allocatable             :: watRows(:)                 ! indices of rows for hydrology column
  ! indices of model layers
  integer(i4b)                         :: iLayer                     ! index of model layer
  integer(i4b)                         :: jLayer                     ! index of model layer within the full state vector (hydrology)
  integer(i4b)                         :: qLayer                     ! indices of snow+lake+glce layers
  ! conversion factors
  real(rkind)                          :: LH_fu0                     ! latent heat of fusion, modified to be 0 if using enthalpy formulation and not using
  real(rkind)                          :: maxVolIceContent_use       ! maximum volumetric ice content depending if snow or firn
  character(LEN=256)                   :: cmessage                   ! error message of downwind routine
  logical(lgt)                         :: full                       ! flag to indicate if the matrix is full (true) or banded (false)
  ! --------------------------------------------------------------
  ! associate variables from data structures
  associate(&
    ! indices of model state variables
    ixCasNrg                     => indx_data%var(iLookINDEX%ixCasNrg)%dat(1)                  ,& ! intent(in): [i4b]    index of canopy air space energy state variable
    ixVegNrg                     => indx_data%var(iLookINDEX%ixVegNrg)%dat(1)                  ,& ! intent(in): [i4b]    index of canopy energy state variable
    ixVegHyd                     => indx_data%var(iLookINDEX%ixVegHyd)%dat(1)                  ,& ! intent(in): [i4b]    index of canopy hydrology state variable (mass)
    ixAqWat                      => indx_data%var(iLookINDEX%ixAqWat)%dat(1)                   ,& ! intent(in): [i4b]    index of water storage in the aquifer
    ! vector of energy indices for the layer domains
    ! NOTE: states not in the subset are equal to integerMissing
    ixSnLaSoGlNrg                => indx_data%var(iLookINDEX%ixSnLaSoGlNrg)%dat                ,& ! intent(in): [i4b(:)] index in the state subset for energy state variables in the layer domains
    ixSoilOnlyNrg                => indx_data%var(iLookINDEX%ixSoilOnlyNrg)%dat                ,& ! intent(in): [i4b(:)] index in the state subset for energy state variables in the soil domain
    ! vector of hydrology indices for the layer domains
    ! NOTE: states not in the subset are equal to integerMissing
    noThetaChange                => indx_data%var(iLookINDEX%noThetaChange)%dat(1)             ,& ! intent(in): [i4b]    number of layers with no change in total water content (bottom layers)  
    ixSnLaSoGlHyd                => indx_data%var(iLookINDEX%ixSnLaSoGlHyd)%dat                ,& ! intent(in): [i4b(:)] index in the state subset for hydrology state variables in the layer domains
    ixSoilOnlyHyd                => indx_data%var(iLookINDEX%ixSoilOnlyHyd)%dat                ,& ! intent(in): [i4b(:)] index in the state subset for hydrology state variables in the soil domain
    ! number of state variables of a specific type
    nSnLaSoGlNrg                 => indx_data%var(iLookINDEX%nSnLaSoGlNrg)%dat(1)              ,& ! intent(in): [i4b]    number of energy state variables in the layer domains
    nSnowOnlyNrg                 => indx_data%var(iLookINDEX%nSnowOnlyNrg)%dat(1)              ,& ! intent(in): [i4b]    number of energy state variables in the snow domain
    nSoilOnlyNrg                 => indx_data%var(iLookINDEX%nSoilOnlyNrg)%dat(1)              ,& ! intent(in): [i4b]    number of energy state variables in the soil domain
    nLakeOnlyNrg                 => indx_data%var(iLookINDEX%nLakeOnlyNrg)%dat(1)              ,& ! intent(in): [i4b]    number of energy state variables in the lake domain
    nGlceOnlyNrg                 => indx_data%var(iLookINDEX%nGlceOnlyNrg)%dat(1)              ,& ! intent(in): [i4b]    number of energy state variables in the glacier ice domain
    nSnLaSoGlHyd                 => indx_data%var(iLookINDEX%nSnLaSoGlHyd)%dat(1)              ,& ! intent(in): [i4b]    number of hydrology variables in the layer domains
    nSnowOnlyHyd                 => indx_data%var(iLookINDEX%nSnowOnlyHyd)%dat(1)              ,& ! intent(in): [i4b]    number of hydrology variables in the snow domain
    nSoilOnlyHyd                 => indx_data%var(iLookINDEX%nSoilOnlyHyd)%dat(1)              ,& ! intent(in): [i4b]    number of hydrology variables in the soil domain
    nLakeOnlyHyd                 => indx_data%var(iLookINDEX%nLakeOnlyHyd)%dat(1)              ,& ! intent(in): [i4b]    number of hydrology variables in the lake domain
    nGlceOnlyHyd                 => indx_data%var(iLookINDEX%nGlceOnlyHyd)%dat(1)              ,& ! intent(in): [i4b]    number of hydrology variables in the glacier ice domain
    ! derivatives in net vegetation energy fluxes w.r.t. relevant state variables
    dCanopyNetFlux_dCanWat       => deriv_data%var(iLookDERIV%dCanopyNetFlux_dCanWat)%dat(1)   ,& ! intent(in): [dp]     derivative in net canopy fluxes w.r.t. canopy total water content
    ! derivatives in canopy water w.r.t canopy temperature
    dTheta_dTkCanopy             => deriv_data%var(iLookDERIV%dTheta_dTkCanopy)%dat(1)         ,& ! intent(in): [dp]     derivative in volumetric liquid water content w.r.t. temperature
    d2Theta_dTkCanopy2           => deriv_data%var(iLookDERIV%d2Theta_dTkCanopy2)%dat(1)       ,& ! intent(in): [dp]     second derivative of volumetric liquid water content w.r.t. temperature
    dFracLiqVeg_dTkCanopy        => deriv_data%var(iLookDERIV%dFracLiqVeg_dTkCanopy)%dat(1)    ,& ! intent(in): [dp]     derivative in fraction of (throughfall + drainage)  w.r.t. temperature
    ! derivatives in energy fluxes at the interface of layers w.r.t. water state in layers above and below
    dNrgFlux_dWatAbove           => deriv_data%var(iLookDERIV%dNrgFlux_dWatAbove)%dat          ,& ! intent(in): [dp(:)]  derivatives in the flux w.r.t. water state in the layer above
    dNrgFlux_dWatBelow           => deriv_data%var(iLookDERIV%dNrgFlux_dWatBelow)%dat          ,& ! intent(in): [dp(:)]  derivatives in the flux w.r.t. water state in the layer below
    ! derivative in liquid water fluxes for the soil domain w.r.t hydrology state variables
    dVolTot_dPsi0                => deriv_data%var(iLookDERIV%dVolTot_dPsi0)%dat               ,& ! intent(in): [dp(:)]  derivatives in total water content w.r.t. total water matric potential
    d2VolTot_dPsi02              => deriv_data%var(iLookDERIV%d2VolTot_dPsi02)%dat             ,& ! intent(in): [dp(:)]  second derivatives in total water content w.r.t. total water matric potential
    dCompress_dPsi               => deriv_data%var(iLookDERIV%dCompress_dPsi)%dat              ,& ! intent(in): [dp(:)]  derivatives in compressibility w.r.t matric head
    ! derivative in liquid water fluxes for the layer domains w.r.t temperature
    dFracLiqWat_dTk              => deriv_data%var(iLookDERIV%dFracLiqWat_dTk)%dat             ,& ! intent(in): [dp(:)]  derivatives in fraction of liquid w.r.t. temperature
    mLayerdTheta_dTk             => deriv_data%var(iLookDERIV%mLayerdTheta_dTk)%dat            ,& ! intent(in): [dp(:)]  derivatives in volumetric liquid water content w.r.t. temperature
    mLayerd2Theta_dTk2           => deriv_data%var(iLookDERIV%mLayerd2Theta_dTk2)%dat          ,& ! intent(in): [dp(:)]  second derivative in volumetric liquid water content w.r.t. temperature
    ! derivative in bulk heat capacity w.r.t. relevant state variables
    dVolHtCapBulk_dPsi0          => deriv_data%var(iLookDERIV%dVolHtCapBulk_dPsi0)%dat         ,& ! intent(in): [dp(:)]  derivatives in bulk heat capacity w.r.t. matric potential
    dVolHtCapBulk_dTheta         => deriv_data%var(iLookDERIV%dVolHtCapBulk_dTheta)%dat        ,& ! intent(in): [dp(:)]  derivatives in bulk heat capacity w.r.t. volumetric water content
    dVolHtCapBulk_dCanWat        => deriv_data%var(iLookDERIV%dVolHtCapBulk_dCanWat)%dat(1)    ,& ! intent(in): [dp   ]  derivative in bulk heat capacity w.r.t. volumetric water content
    dVolHtCapBulk_dTk            => deriv_data%var(iLookDERIV%dVolHtCapBulk_dTk)%dat           ,& ! intent(in): [dp(:)]  derivatives in bulk heat capacity w.r.t. temperature
    dVolHtCapBulk_dTkCanopy      => deriv_data%var(iLookDERIV%dVolHtCapBulk_dTkCanopy)%dat(1)  ,& ! intent(in): [dp   ]  derivative in bulk heat capacity w.r.t. temperature
    ! derivative in Cm w.r.t. relevant state variables
    dCm_dPsi0                    => deriv_data%var(iLookDERIV%dCm_dPsi0)%dat                   ,& ! intent(in): [dp(:)]  derivatives in heat capacity w.r.t. matric potential (J kg-1)
    dCm_dTk                      => deriv_data%var(iLookDERIV%dCm_dTk)%dat                     ,& ! intent(in): [dp(:)]  derivatives in heat capacity w.r.t. temperature (J kg-1 K-2)
    dCm_dTkCanopy                => deriv_data%var(iLookDERIV%dCm_dTkCanopy)%dat(1)            ,& ! intent(in): [dp   ]  derivative in heat capacity w.r.t. canopy temperature (J kg-1 K-2)
    ! derivatives of temperature if enthalpy is the state variable
    dCanairTemp_dEnthalpy        => deriv_data%var(iLookDERIV%dCanairTemp_dEnthalpy)%dat(1)    ,& ! intent(in): [dp]     derivative in canopy air temperature w.r.t. enthalpy
    dCanopyTemp_dEnthalpy        => deriv_data%var(iLookDERIV%dCanopyTemp_dEnthalpy)%dat(1)    ,& ! intent(in): [dp]     derivative in canopy temperature w.r.t. enthalpy 
    dTemp_dEnthalpy              => deriv_data%var(iLookDERIV%dTemp_dEnthalpy)%dat             ,& ! intent(in): [dp(:)]  derivatives in temperature w.r.t. enthalpy
    dCanopyTemp_dCanWat          => deriv_data%var(iLookDERIV%dCanopyTemp_dCanWat)%dat(1)      ,& ! intent(in): [dp]     derivative in canopy temperature w.r.t. volumetric water content
    dTemp_dTheta                 => deriv_data%var(iLookDERIV%dTemp_dTheta)%dat                ,& ! intent(in): [dp(:)]  derivatives in temperature w.r.t. volumetric water content
    dTemp_dPsi0                  => deriv_data%var(iLookDERIV%dTemp_dPsi0)%dat                 ,& ! intent(in): [dp(:)]  derivatives in temperature w.r.t. total water matric potential
    ! diagnostic variables
    scalarFracLiqVeg             => diag_data%var(iLookDIAG%scalarFracLiqVeg)%dat(1)           ,& ! intent(in): [dp]     fraction of liquid water on vegetation (-)
    scalarBulkVolHeatCapVeg      => diag_data%var(iLookDIAG%scalarBulkVolHeatCapVeg)%dat(1)    ,& ! intent(in): [dp]     bulk volumetric heat capacity of vegetation (J m-3 K-1)
    scalarCanopyCm               => diag_data%var(iLookDIAG%scalarCanopyCm)%dat(1)             ,& ! intent(in): [dp]     Cm for canopy vegetation (J kg-1)
    mLayerFracLiq                => diag_data%var(iLookDIAG%mLayerFracLiq)%dat                 ,& ! intent(in): [dp(:)]  fraction of liquid water in each snow, lake, or glce layer (-)
    mLayerVolHtCapBulk           => diag_data%var(iLookDIAG%mLayerVolHtCapBulk)%dat            ,& ! intent(in): [dp(:)]  bulk volumetric heat capacity in each layer (J m-3 K-1)
    mLayerCm                     => diag_data%var(iLookDIAG%mLayerCm)%dat                      ,& ! intent(in): [dp(:)]  Cm for each layer (J m-3)
    ! canopy and layer depth
    canopyDepth                  => diag_data%var(iLookDIAG%scalarCanopyDepth)%dat(1)          ,& ! intent(in): [dp   ]  canopy depth (m)
    mLayerDepth                  => prog_data%var(iLookPROG%mLayerDepth)%dat                    & ! intent(in): [dp(:)]  depth of each layer in the sub-domain (m)
    ) ! making association with data in structures
    ! --------------------------------------------------------------
    ! initialize error control
    err=0; message='computJacobWithPrime/'

    ! *********************************************************************************************************************************************************
    ! * PART 0: PRELIMINARIES (INITIALIZE JACOBIAN AND COMPUTE TIME-VARIABLE DIAGONAL TERMS)
    ! *********************************************************************************************************************************************************
    ! get the number of state variables
    nState = size(dMat0)

    ! initialize the Jacobian and diagonal
    ! NOTE: this needs to be done every time, since Jacobian matrix is modified in the solver and dMat is modified below
    aJac(:,:) = 0._rkind  ! analytical Jacobian matrix
    allocate(dMat(nState)) 
    dMat = dMat0 * cj ! dMat0(ixCasNrg) = Cp_air*iden_air and dMat0(Wat states) = 1.0

    if(computeVegFlux)then
      ! compute terms in the Jacobian for vegetation (excluding fluxes)
      if(ixVegNrg/=integerMissing)&
          dMat(ixVegNrg) = ( scalarBulkVolHeatCapVeg + LH_fus*iden_water*dTheta_dTkCanopy ) * cj &
                          + dVolHtCapBulk_dTkCanopy * scalarCanopyTempPrime &
                          + dCm_dTkCanopy * scalarCanopyWatPrime/canopyDepth &
                          + LH_fus*iden_water * scalarCanopyTempPrime * d2Theta_dTkCanopy2 &
                          + LH_fus * dFracLiqVeg_dTkCanopy * scalarCanopyWatPrime/canopyDepth
    endif

    ! comput terms for the Jacobian for the layer domain (excluding fluxes)
    do iLayer=1,nLayers
      if(ixSnLaSoGlNrg(iLayer)/=integerMissing)&
          dMat(ixSnLaSoGlNrg(iLayer)) = ( mLayerVolHtCapBulk(iLayer) + LH_fus*iden_water*mLayerdTheta_dTk(iLayer) ) * cj &
                                       + dVolHtCapBulk_dTk(iLayer) * mLayerTempPrime(iLayer) &
                                       + dCm_dTk(iLayer) * mLayerVolFracWatPrime(iLayer) &
                                       + LH_fus*iden_water * mLayerTempPrime(iLayer) * mLayerd2Theta_dTk2(iLayer) &
                                       + LH_fus*iden_water * dFracLiqWat_dTk(iLayer) * mLayerVolFracWatPrime(iLayer)
    end do

    ! compute terms for the Jacobian for the soil domain (excluding fluxes)
    do iLayer=1,nSoil
      if(ixSoilOnlyHyd(iLayer)/=integerMissing)then ! writes over dMat(ixSoilOnlyHyd(iLayer) = 1.0 * cj above
        dMat(ixSoilOnlyHyd(iLayer)) = ( dVolTot_dPsi0(iLayer) + dCompress_dPsi(iLayer) ) * cj + d2VolTot_dPsi02(iLayer) * mLayerMatricHeadPrime(iLayer)
        dMat(ixSoilOnlyHyd(iLayer)) = dMat(ixSoilOnlyHyd(iLayer)) + specificStorage * dVolTot_dPsi0(iLayer) * mLayerMatricHeadPrime(iLayer)/theta_sat(iLayer)
      endif
    end do

    ! if using enthalpy as a state variable, zero out usual RHS terms and add them end of the iteration loop 
    ! NOTE: other terms on RHS that are not fluxes are zeroed out by not computing heat capacity and Cm and their derivatives
    if(enthalpyStateVec)then 
      if(ixCasNrg/=integerMissing) dMat(ixCasNrg) = 0._rkind
      if(ixVegNrg/=integerMissing) dMat(ixVegNrg) = 0._rkind
      do iLayer=1,nLayers
        if(ixSnLaSoGlNrg(iLayer)/=integerMissing) dMat(ixSnLaSoGlNrg(iLayer)) = 0._rkind
      end do
      LH_fu0 = 0._rkind ! set to 0 to not use RHS terms
    else
      LH_fu0 = LH_fus ! use regular value
    endif

    ! compute the maximum volumetric ice content for the layer domains
    if(nGlce>0)then ! snow can be firn
      maxVolIceContent_use = min(maxVolIceContent+0.15,0.85_rkind) ! firn maximum volumetric ice content to store water (-)
    else ! snow
      maxVolIceContent_use = maxVolIceContent ! snow maximum volumetric ice content to store water (-)
    endif

    ! define the form of the matrix
    select case(ixMatrix)
     case(ixBandMatrix)
       ! check
       if(size(aJac,1)/=nBands .or. size(aJac,2)/=size(dMat))then
         message=trim(message)//'unexpected shape of the Jacobian matrix: expect aJac(nBands,nState)'
         err=20; return
       endif
       full = .false.
     case(ixFullMatrix)
       ! check
       if(size(aJac,1)/=size(dMat) .or. size(aJac,2)/=size(dMat))then
         message=trim(message)//'unexpected shape of the Jacobian matrix: expect aJac(nState,nState)'
         err=20; return
       endif
       full = .true.
     case default; err=20; message=trim(message)//'unable to identify option for the type of matrix'; return
    end select

    ! *********************************************************************************************************************************************************
    ! * PART 1: COMPUTE CROSS-DERIVATIVE JACOBIAN TERMS
    ! *********************************************************************************************************************************************************
    ! -----
    ! * cross derivatives in the vegetation...
    ! ---------------------------------------------
    if(computeVegFlux)then ! (derivatives only defined when vegetation protrudes over the surface)
      if(ixVegHyd/=integerMissing .and. ixVegNrg/=integerMissing)&
          ! NOTE: dIce/dLiq = (1 - scalarFracLiqVeg); dIce*LH_fu0/canopyDepth = J m-3; dLiq = kg m-2
          aJac(ixInd(full,ixVegNrg,ixVegHyd),ixVegHyd) = (-1._rkind + scalarFracLiqVeg)*LH_fu0/canopyDepth * cj &
                                                   + dVolHtCapBulk_dCanWat * scalarCanopyTempPrime + scalarCanopyCm/canopyDepth * cj &
                                                   - (dt/canopyDepth) * dCanopyNetFlux_dCanWat &
                                                   + LH_fu0 * scalarCanopyTempPrime * dFracLiqVeg_dTkCanopy/canopyDepth
    endif  ! if there is a need to compute energy fluxes within vegetation

    ! -----
    ! * cross derivatives in the snow, lake, glce domains...
    ! ----------------------------------------
    if((nSnowOnlyHyd>0 .and. nSnowOnlyNrg>0) .or. (nLakeOnlyHyd>0 .and. nLakeOnlyNrg>0) .or. (nGlceOnlyHyd>0 .and. nGlceOnlyNrg>0))then
      do qLayer=1,nSnow+nLake+nGlce-noThetaChange ! loop through layers in the snow, lake, glce domains

        if(qLayer<=nSnow+nLake)then
          jLayer = qLayer
          iLayer = qLayer
        else
          jLayer = qLayer + nLake + nSoil
          iLayer = qLayer - nSnow - nLake
        endif
        ! - check that the layer is desired
        if(ixSnLaSoGlNrg(jLayer)==integerMissing) cycle
        ! (define the energy state)
        nrgState = ixSnLaSoGlNrg(jLayer)       ! index within the full state vector
        ! - define state indices for the current layer
        watState = ixSnLaSoGlHyd(jLayer)   ! hydrology state index within the state subset

        if(watState/=integerMissing)then       ! (water state for the current layer is within the state subset)
          ! - include derivatives of energy fluxes w.r.t water fluxes for current layer
          aJac(ixInd(full,nrgState,watState),watState) = (-1._rkind + mLayerFracLiq(jLayer))*LH_fu0*iden_water * cj &
                                      + dVolHtCapBulk_dTheta(jLayer) * mLayerTempPrime(jLayer) + mLayerCm(jLayer) * cj &
                                      + (dt/mLayerDepth(jLayer))*(-dNrgFlux_dWatBelow(jLayer-1) + dNrgFlux_dWatAbove(jLayer)) &
                                      + LH_fu0*iden_water * mLayerTempPrime(jLayer) * dFracLiqWat_dTk(jLayer)    ! (dF/dLiq)
        endif ! (if the water state for the current layer is within the state subset)

      end do ! (looping through snow, lake, glce layers)
    endif ! (if there are state variables for both water and energy in the snow, lake, glce domains)

    ! -----
    ! * cross derivatives in the soil domain...
    ! ----------------------------------------
    if(nSoilOnlyHyd>0 .and. nSoilOnlyNrg>0)then
      do iLayer=1,nSoilOnlyNrg

        ! - check that the soil layer is desired
        if(ixSoilOnlyNrg(iLayer)==integerMissing) cycle
        ! - define indices of the soil layers
        jLayer   = iLayer+nSnow+nLake            ! index of layer in the layer system
        ! - define the energy state variable
        nrgState = ixSoilOnlyNrg(iLayer)         ! index within the full state vector
        ! - define index of hydrology state variable within the state subset
        watState = ixSoilOnlyHyd(iLayer)

        ! only compute derivatives if the water state for the current layer is within the state subset
        if(watState/=integerMissing)then
          ! - include derivatives in energy fluxes w.r.t. with respect to water for current layer
          aJac(ixInd(full,nrgState,watState),watState) = (dt/mLayerDepth(jLayer))*(-dNrgFlux_dWatBelow(jLayer-1) + dNrgFlux_dWatAbove(jLayer))
          aJac(ixInd(full,nrgState,watState),watState) = mLayerCm(jLayer) * dVolTot_dPsi0(iLayer) * cj + dVolHtCapBulk_dPsi0(iLayer) * mLayerTempPrime(jLayer) &
                                                        + dCm_dPsi0(iLayer) * mLayerVolFracWatPrime(jLayer) + mLayerCm(jLayer) * d2VolTot_dPsi02(iLayer) * mLayerMatricHeadPrime(iLayer) &
                                                        + aJac(ixInd(full,nrgState,watState),watState) 
          if(mLayerdTheta_dTk(jLayer) > tiny(1.0_rkind))&  ! ice is present
             aJac(ixInd(full,nrgState,watState),watState) = -LH_fu0*iden_water * dVolTot_dPsi0(iLayer) * cj &
                                                           - LH_fu0*iden_water * mLayerMatricHeadPrime(iLayer) * d2VolTot_dPsi02(iLayer) + aJac(ixInd(full,nrgState,watState),watState) ! dNrg/dMat (J m-3 m-1) -- dMat changes volumetric water, and hence ice content
        endif ! (if the water state for the current layer is within the state subset)

      end do ! (looping through energy states in the soil domain)
    endif ! (if there are state variables for both water and energy in the soil domain)

    ! *********************************************************************************************************************************************************
    ! * PART 2: COMPUTE FLUX JACOBIAN TERMS 
    ! *********************************************************************************************************************************************************
    call fluxJacAdd(full,maxVolIceContent_use,dt,nSnow,nLake,nSoil,nGlce,nLayers,computeVegFlux,computeBaseflow,&
                    indx_data,prog_data,diag_data,deriv_data,dBaseflow_dWat,dBaseflow_dTk,dMat,aJac,err,cmessage)
    if(err/=0)then; message=trim(message)//trim(cmessage); return; endif
    deallocate(dMat)

    ! *********************************************************************************************************************************************************
    ! * PART 3: CLEAN UP JACOBIAN (IF USING ENTHALPY AS A STATE VARIABLE) AND PRINT (IF DESIRED)
    ! *********************************************************************************************************************************************************
    ! * if desired, modify to use enthalpy as a state variable instead of temperature 
    ! NOTE, dMat(Nrg states) was used as 0 and now 1._rkind * cj is added instead 
    ! ----------------------------------------
    if(enthalpyStateVec)then 

      if(full) then
        allocate(watRows(nState),nrgRows(nState)) ! all rows are used
        do jLayer=1,nState
          watRows(jLayer) = jLayer
          nrgRows(jLayer) = jLayer
        end do
      else
        allocate(watRows(nBands-1),nrgRows(nBands-1)) ! only the bands are used
        do jLayer=1,nBands-1 ! water row nBand would add a 0 value from the nrg matrix since it's out of range
          watRows(jLayer) = jLayer
          nrgRows(jLayer) = jLayer + 1
        enddo
      endif

      if(ixCasNrg/=integerMissing)then
        aJac(:,ixCasNrg) = aJac(:,ixCasNrg) * dCanairTemp_dEnthalpy
        aJac(ixInd(full,ixCasNrg,ixCasNrg),ixCasNrg) = aJac(ixInd(full,ixCasNrg,ixCasNrg),ixCasNrg) + 1._rkind * cj
      endif
      
      if(ixVegNrg/=integerMissing)then
        if(ixVegHyd/=integerMissing) aJac(watRows,ixVegHyd) = aJac(watRows,ixVegHyd) + aJac(nrgRows,ixVegNrg) * dCanopyTemp_dCanWat
        aJac(:,ixVegNrg) = aJac(:,ixVegNrg) * dCanopyTemp_dEnthalpy
        aJac(ixInd(full,ixVegNrg,ixVegNrg),ixVegNrg) = aJac(ixInd(full,ixVegNrg,ixVegNrg),ixVegNrg) + 1._rkind * cj
      endif
      
      if(nSnLaSoGlNrg>0)then
        do iLayer=1,nLayers
          nrgState = ixSnLaSoGlNrg(iLayer)       
          if(nrgState==integerMissing) cycle
          watState = ixSnLaSoGlHyd(iLayer)
          if(watState/=integerMissing)then 
            if(iLayer<=nSnow+nLake .or. iLayer>nSnow+nLake+nSoil)then
              aJac(watRows,watState) = aJac(watRows,watState) + aJac(nrgRows,nrgState) * dTemp_dTheta(iLayer)
            else
              aJac(watRows,watState) = aJac(watRows,watState) + aJac(nrgRows,nrgState) * dTemp_dPsi0(iLayer-nSnow-nLake)
            endif
          endif
          aJac(:,nrgState) = aJac(:,nrgState) * dTemp_dEnthalpy(iLayer)
          aJac(ixInd(full,nrgState,nrgState),nrgState) = aJac(ixInd(full,nrgState,nrgState),nrgState) + 1._rkind * cj
        enddo
      endif
    else
      allocate(watRows(0),nrgRows(0)) ! dummy allocation to avoid compiler warning
    endif
    deallocate(watRows,nrgRows)
    
    ! print the Jacobian
    if(globalPrintFlag .or. any(isNan(aJac)))then
      if(full) then
        print*, '** full analytical Jacobian:'
        write(*,'(a4,1x,100(i12,1x))') 'xCol', (iLayer, iLayer=min(iJac1,nState),min(iJac2,nState))
        do iLayer=min(iJac1,nState),min(iJac2,nState)
          write(*,'(i4,1x,100(e12.5,1x))') iLayer, aJac(min(iJac1,nState):min(iJac2,nState),iLayer)
        end do
      else
        print*, '** banded analytical Jacobian:'
        write(*,'(a4,1x,100(i17,1x))') 'xCol', (iLayer, iLayer=min(iJac1,nState),min(iJac2,nState))
        do iLayer=kl+1,nBands
          write(*,'(i4,1x,100(e17.10,1x))') iLayer, (aJac(iLayer,jLayer),jLayer=min(iJac1,nState),min(iJac2,nState))
        end do
      endif
    endif
    if(any(isNan(aJac)))then; message=trim(message)//'NaN in Jacobian';err=20; return; endif

  ! end association to variables in the data structures
  end associate

end subroutine computJacobWithPrime

! **********************************************************************************************************
! public function computJacob4ida: the interface to compute the Jacobian matrix dF/dy + c dF/dy' for IDA solver
! **********************************************************************************************************
! Return values:
!    0 = success,
!    1 = recoverable error,
!   -1 = non-recoverable error
! ----------------------------------------------------------------
integer(c_int) function computJacob4ida(t, cj, sunvec_y, sunvec_yp, sunvec_r, &
                    sunmat_J, user_data, sunvec_temp1, sunvec_temp2, sunvec_temp3) &
                    result(ierr) bind(C,name='computJacob4ida')

  !======= Inclusions ===========
  use, intrinsic :: iso_c_binding
  use fsundials_core_mod
  use fnvector_serial_mod
  use fsunmatrix_band_mod
  use fsunmatrix_dense_mod
  use type4ida

  !======= Declarations =========
  implicit none

  ! calling variables
  real(rkind), value            :: t              ! current time
  real(rkind), value            :: cj             ! step size scaling factor
  type(N_Vector)                :: sunvec_y       ! solution N_Vector
  type(N_Vector)                :: sunvec_yp      ! derivative N_Vector
  type(N_Vector)                :: sunvec_r       ! residual N_Vector
  type(SUNMatrix)               :: sunmat_J       ! Jacobian SUNMatrix
  type(c_ptr), value            :: user_data      ! user-defined data
  type(N_Vector)                :: sunvec_temp1   ! temporary N_Vector
  type(N_Vector)                :: sunvec_temp2   ! temporary N_Vector
  type(N_Vector)                :: sunvec_temp3   ! temporary N_Vector

  ! pointers to data in SUNDIALS vectors
  real(rkind), pointer          :: Jac(:,:)       ! Jacobian matrix
  type(data4ida), pointer       :: eqns_data      ! equations data
  ! ----------------------------------------------------------------

  ! get equations data from user-defined data
  call c_f_pointer(user_data, eqns_data)

  ! get data arrays from SUNDIALS vectors
  if(eqns_data%ixMatrix==ixBandMatrix) Jac(1:nBands, 1:eqns_data%nState) => FSUNBandMatrix_Data(sunmat_J)
  if(eqns_data%ixMatrix==ixFullMatrix) Jac(1:eqns_data%nState, 1:eqns_data%nState) => FSUNDenseMatrix_Data(sunmat_J)

  ! compute the analytical Jacobian matrix
  ! NOTE: The derivatives were computed in the previous call to computFlux
  !       This occurred either at the call to eval8summaWithPrime at the start of systemSolv
  !        or in the call to eval8summaWithPrime in the previous iteration
  call computJacobWithPrime(&
                ! input: model control
                cj,                                       & ! intent(in):    this scalar changes whenever the step size or method order changes
                1._qp,                                    & ! intent(in):    length of the time step (seconds)
                eqns_data%nSnow,                          & ! intent(in):    number of snow layers
                eqns_data%nLake,                          & ! intent(in):    number of lake layers
                eqns_data%nSoil,                          & ! intent(in):    number of soil layers
                eqns_data%nGlce,                          & ! intent(in):    number of glacier ice layers
                eqns_data%nLayers,                        & ! intent(in):    total number of layers
                eqns_data%computeVegFlux,                 & ! intent(in):    flag to indicate if we need to compute fluxes over vegetation
                eqns_data%model_decisions(iLookDECISIONS%groundwatr)%iDecision==qbaseTopmodel, & ! intent(in): flag to indicate if we need to compute baseflow
                eqns_data%ixMatrix,                                                            & ! intent(in): form of the Jacobian matrix
                eqns_data%mpar_data%var(iLookPARAM%specificStorage)%dat(1),                    & ! intent(in): specific storage coefficient (m-1)
                eqns_data%mpar_data%var(iLookPARAM%theta_sat)%dat,                             & ! intent(in): soil porosity (-)
                eqns_data%model_decisions(iLookDECISIONS%nrgConserv)%iDecision.ne.closedForm,  & ! intent(in): flag if enthalpy is state variable
                ! input: data structures
                eqns_data%indx_data,                      & ! intent(in):    index data
                eqns_data%prog_data,                      & ! intent(in):    model prognostic variables for a local HRU
                eqns_data%diag_data,                      & ! intent(in):    model diagnostic variables for a local HRU
                eqns_data%deriv_data,                     & ! intent(in):    derivatives in model fluxes w.r.t. relevant state variables
                eqns_data%dBaseflow_dWat,                 & ! intent(in):    derivative in baseflow w.r.t. soil water characteristic
                eqns_data%dBaseflow_dTk,                  & ! intent(in):    derivative in baseflow w.r.t. temperature (m s-1 K-1)
                ! input: state variables
                eqns_data%mLayerTempPrime,                & ! intent(in):    derivative value for temperature of each layer (K)
                eqns_data%mLayerMatricHeadPrime,          & ! intent(in):    derivative value for matric head of each layer (m)
                eqns_data%mLayerVolFracWatPrime,          & ! intent(in):    derivative value for volumetric total water content of each layer (-)
                eqns_data%scalarCanopyTempPrime,          & ! intent(in):    derivative value for temperature of the vegetation canopy (K)
                eqns_data%scalarCanopyWatPrime,           & ! intent(in):    derivative value for total water content of the vegetation canopy (kg m-2)
                ! input-output: Jacobian and its diagonal
                eqns_data%dMat,                           & ! intent(in):    diagonal of the Jacobian matrix excluding fluxes, not depending on the state vector
                Jac,                                      & ! intent(out):   Jacobian matrix
                ! output: error control
                eqns_data%err,eqns_data%message)            ! intent(out):   error code and error message
  if(eqns_data%err > 0)then; eqns_data%message=trim(eqns_data%message); ierr=-1; return; endif
  if(eqns_data%err < 0)then; eqns_data%message=trim(eqns_data%message); ierr=1; return; endif

  ! return success
  ierr = 0
  return

end function computJacob4ida

end module computJacobWithPrime_module

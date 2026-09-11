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

module groundwatr_module

! data types
USE nr_type

! model constants
USE multiconst,only:iden_water   ! density of water (kg m-3)

! derived types to define the data structures
USE data_types,only:&
                    var_dlength,        & ! data vector with variable length dimension (rkind)
                    in_type_groundwatr, & ! intent(in) arguments for groundwatr call
                    io_type_groundwatr, & ! intent(inout) arguments for groundwatr call
                    out_type_groundwatr   ! intent(out) arguments for groundwatr call

! named variables defining elements in the data structures
USE var_lookup,only:iLookPROG    ! named variables for structure elements
USE var_lookup,only:iLookDIAG    ! named variables for structure elements
USE var_lookup,only:iLookFLUX    ! named variables for structure elements
USE var_lookup,only:iLookPARAM   ! named variables for structure elements

! model decision structures
USE globalData,only:model_decisions        ! model decision structure
USE var_lookup,only:iLookDECISIONS         ! named variables for elements of the decision structure

! look-up values for the choice of hydraulic conductivity profile
USE mDecisions_module,only: &
 constant,                  & ! constant hydraulic conductivity with depth
 powerLaw_profile,          & ! power-law profile
 expLaw_profile               ! exponential profile

! privacy
implicit none
private
public :: groundwatr
contains

! ************************************************************************************************
! public subroutine groundwatr: compute the groundwater sink term in Richards' equation
! ************************************************************************************************
!
! Method
! ------
!
! Here we assume that water available for shallow groundwater flow includes is all water above
! "field capacity" below the depth zCrit, where zCrit is defined as the lowest point in the soil
! profile where the volumetric liquid water content is less than field capacity.
!
! We further assume that transmssivity (m2 s-1) for each layer is defined assuming that the water
! available for saturated flow is located at the bottom of the soil profile. Specifically:
!  trTotal(iLayer) = tran0*xTrans(zActive(iLayer))
!  trSoil(iLayer)  = trTotal(iLayer) - trTotal(iLayer+1)
! where zActive(iLayer) is the effective water table thickness for all layers up to and including
! the current layer (working from the bottom to the top).
!
! Transmissivity is the vertical integral of the hydraulic conductivity profile, so xTrans follows
! whichever profile satHydCond used, for saturated thickness s and total soil depth D:
!  powerLaw_profile: tran0 = kAnisotropic*K_0*D/zScale_TOPMODEL, xTrans = (s/D)**zScale_TOPMODEL
!  expLaw_profile:   tran0 = kAnisotropic*K_0/f_hydCond,         xTrans = exp(-f*(D-s)) - exp(-f*D)
! Both give xTrans=0 at s=0. The exponential form keeps conductivity finite at the base of the soil,
! which matters for glacier debris where melt enters through the bottom boundary.
!
! The outflow from each layer is then (m3 s-1)
!  mLayerOutflow(iLayer) = trSoil(iLayer)*tan_slope*contourLength
! where contourLength is the width of a hillslope (m) parallel to a stream
!
! ************************************************************************************************
subroutine groundwatr(&
                      ! input: model control, state variables, and diagnostic variables
                      in_groundwatr,                          & ! intent(in): model control, state variables, and diagnostic variables
                      ! input/output: data structures
                      mpar_data,                              & ! intent(in):    model parameters
                      prog_data,                              & ! intent(in):    model prognostic variables for a local HRU
                      flux_data,                              & ! intent(inout): model fluxes for a local HRU
                      ! input-output: baseflow
                      io_groundwatr,                          & ! intent(inout): index of lowest saturated layer (NOTE: only computed on the first iteration)
                      ! output: baseflow and error control
                      out_groundwatr)                           ! intent(out):   baseflow and error control
  ! ---------------------------------------------------------------------------------------
  implicit none
  ! input: model control, state variables, and diagnostic variables
  type(in_type_groundwatr),intent(in)    :: in_groundwatr     ! model control, state variables, and diagnostic variables
  ! input-output: data structures
  type(var_dlength),intent(in)           :: mpar_data         ! model parameters
  type(var_dlength),intent(in)           :: prog_data         ! prognostic variables for a local HRU
  type(var_dlength),intent(inout)        :: flux_data         ! model fluxes for a local HRU
  ! input-output: baseflow
  type(io_type_groundwatr),intent(inout) :: io_groundwatr     ! index of lowest saturated layer (NOTE: only computed on the first iteration)
  ! output: baseflow and error control
  type(out_type_groundwatr),intent(out)  :: out_groundwatr    ! baseflow and error control
  ! general local variables
  integer(i4b)                           :: iLayer            ! index of soil layer
  real(rkind)                            :: fieldCapacity_use ! field capacity used to determine the "active" portion of the soil profile
  character(len=256)                     :: cmessage          ! error message
  ! ***************************************************************************************
  ! associate variables in data structures
  allocate(out_groundwatr % mLayerBaseflow(in_groundwatr%nSoil),out_groundwatr % dBaseflow_dWat(in_groundwatr%nSoil,in_groundwatr%nSoil),&
           out_groundwatr % dBaseflow_dTk(in_groundwatr%nSoil,in_groundwatr%nSoil)) ! allocate intent(out) data structure components
  associate(&
    ! input: model control
    nSnow               => in_groundwatr % nSnow,                              & ! intent(in):    [i4b] number of snow layers
    nLake               => in_groundwatr % nLake,                              & ! intent(in):    [i4b] number of lake layers
    nSoil               => in_groundwatr % nSoil,                              & ! intent(in):    [i4b] number of soil layers
    nGlce               => in_groundwatr % nGlce,                              & ! intent(in):    [i4b] number of glacier ice layers
    getSatDepth         => in_groundwatr % firstFluxCall,                      & ! intent(in):    [lgt] logical flag to compute index of the lowest saturated layer
    ! input: diagnostic variables
    mLayerVolFracLiq    => in_groundwatr % mLayerVolFracLiqTrial,              & ! intent(in):    [dp] volumetric fraction of liquid water (-)
    mLayerVolFracIce    => in_groundwatr % mLayerVolFracIceTrial,              & ! intent(in):    [dp] volumetric fraction of ice (-)
    ! input: derivatives
    dVolTot_dPsi0       => in_groundwatr % dVolTot_dPsi0,                      & ! intent(in):    [dp] derivative in total volumetric water content w.r.t. matric head (m-1)
    mLayerdTheta_dPsi   => in_groundwatr % mLayerdTheta_dPsi,                  & ! intent(in):    [dp] derivative in liquid water content w.r.t. matric potential (m-1)
    mLayerdTheta_dTk    => in_groundwatr % mLayerdTheta_dTk,                   & ! intent(in):    [dp] derivative in volumetric liquid water content w.r.t. temperature (K-1)
    ! input: baseflow parameters
    fieldCapacity       => mpar_data%var(iLookPARAM%fieldCapacity)%dat(1),     & ! intent(in):    [dp] field capacity (-)
    theta_sat           => mpar_data%var(iLookPARAM%theta_sat)%dat,            & ! intent(in):    [dp] soil porosity (-)
    ! input-output: baseflow
    ixSaturation        => io_groundwatr % ixSaturation,                       & ! intent(inout): [i4b] index of lowest saturated layer (NOTE: only computed on the first iteration)
    ! output: diagnostic variables
    scalarExfiltration  => flux_data%var(iLookFLUX%scalarExfiltration)%dat(1), & ! intent(out):   [dp]    exfiltration from the soil profile (m s-1)
    mLayerColumnOutflow => flux_data%var(iLookFLUX%mLayerColumnOutflow)%dat,   & ! intent(out):   [dp(:)] column outflow from each soil layer (m3 s-1)
    ! output: baseflow
    mLayerBaseflow      => out_groundwatr % mLayerBaseflow,                    & ! intent(out):   [dp(:)]   baseflow from each soil layer (m s-1)
    dBaseflow_dWat      => out_groundwatr % dBaseflow_dWat,                    & ! intent(out):   [dp(:,:)] derivative in baseflow w.r.t. soil water characteristic
    dBaseflow_dTk       => out_groundwatr % dBaseflow_dTk,                     & ! intent(out):   [dp(:,:)] derivative in baseflow w.r.t. temperature (m s-1 K-1)
    ! output: error control
    err                 => out_groundwatr % err,                               & ! intent(out):   [i4b]       error code
    message             => out_groundwatr % cmessage                           & ! intent(out):   [character] error message
    )  ! end association to variables in data structures

    ! initialize error control
    err=0; message='groundwatr/'

    ! ************************************************************************************************
    ! (1) compute the "active" portion of the soil profile
    ! ************************************************************************************************
    fieldCapacity_use = fieldCapacity
    if(nGlce>0) fieldCapacity_use = 0._rkind ! if glacier ice layers are present, set field capacity to zero (i.e. all water is "active" for flow)

    ! get index of the layer closest to surface that is more than field capacity (NOTE: only compute on the first flux call)
    if (getSatDepth) then
      ixSaturation = nSoil+1  ! unsaturated profile when ixSaturation>nSoil
      do iLayer=nSoil,1,-1  ! start at the lowest soil layer and work upwards to the top layer
        if (mLayerVolFracLiq(iLayer) > fieldCapacity_use) then; ixSaturation = iLayer  ! index of saturated layer -- keeps getting over-written as move upwards
        else; exit; end if                                                         ! only consider saturated layer at the bottom of the soil profile
      end do  ! end looping through soil layers
    end if

    ! check for an early return (no layers are "active")
    if (ixSaturation > nSoil) then
      scalarExfiltration     = 0._rkind   ! exfiltration from the soil profile (m s-1)
      mLayerColumnOutflow(:) = 0._rkind   ! column outflow from each soil layer (m3 s-1)
      mLayerBaseflow(:)      = 0._rkind   ! baseflow from each soil layer (m s-1)
      dBaseflow_dWat(:,:)    = 0._rkind   ! derivative in baseflow w.r.t. soil water characteristic
      dBaseflow_dTk(:,:)     = 0._rkind   ! derivative in baseflow w.r.t. temperature (m s-1 K-1)
      return
    end if  ! if some layers are saturated

    ! ************************************************************************************************
    ! (2) compute the baseflow flux and its derivative w.r.t matric head and temperature
    ! ************************************************************************************************

    ! use private subroutine to compute baseflow (for multiple calls for numerical Jacobian)
    call computBaseflow(&
                          ! input: control and state variables
                          nSnow,                   & ! intent(in):    number of snow layers
                          nLake,                   & ! intent(in):    number of lake layers
                          nSoil,                   & ! intent(in):    number of soil layers
                          nGlce,                   & ! intent(in):    number of glacier ice layers
                          ixSaturation,            & ! intent(in):    index of upper-most "saturated" layer
                          mLayerVolFracLiq,        & ! intent(in):    volumetric fraction of liquid water in each soil layer (-)
                          mLayerVolFracIce,        & ! intent(in):    volumetric fraction of ice in each soil layer (-)
                          ! input/output: data structures
                          mpar_data,               & ! intent(in):    model parameters
                          prog_data,               & ! intent(in):    model prognostic variables for a local HRU
                          flux_data,               & ! intent(inout): model fluxes for a local HRU
                          ! derivatives
                          dVolTot_dPsi0,           & ! intent(in):    derivative in total volumetric water content w.r.t. matric head (m-1)
                          mLayerdTheta_dPsi,       & ! intent(in):    derivative in liquid water content w.r.t. matric potential (m-1)
                          mLayerdTheta_dTk,        & ! intent(in):    derivative in volumetric liquid water content w.r.t. temperature (K-1)
                          ! output: fluxes and derivatives
                          mLayerBaseflow,          & ! intent(out):   baseflow flux in each soil layer (m s-1)
                          dBaseflow_dWat,          & ! intent(out):   derivative in baseflow w.r.t. soil water characteristic
                          dBaseflow_dTk,           & ! intent(out):   derivative in baseflow w.r.t. temperature (m s-1 K-1)
                          err, cmessage)             ! intent(out):   error control
   if(err/=0)then; err=20; message=trim(message)//trim(cmessage); return; endif
   ! end association to variables in data structures
  end associate

end subroutine groundwatr


! ***********************************************************************************************************************
! * private subroutine computBaseflow: compute the baseflow flux and its derivative w.r.t. volumetric liquid water content
! ***********************************************************************************************************************
subroutine computBaseflow(&
                          ! input: control and state variables
                          nSnow,                         & ! intent(in):    number of snow layers
                          nLake,                         & ! intent(in):    number of lake layers
                          nSoil,                         & ! intent(in):    number of soil layers
                          nGlce,                         & ! intent(in):    number of glacier ice layers
                          ixSaturation,                  & ! intent(in):    index of upper-most "saturated" layer
                          mLayerVolFracLiq,              & ! intent(in):    volumetric fraction of liquid water in each soil layer (-)
                          mLayerVolFracIce,              & ! intent(in):    volumetric fraction of ice in each soil layer (-)
                          ! input/output: data structures
                          mpar_data,                     & ! intent(in):    model parameters
                          prog_data,                     & ! intent(in):    model prognostic variables for a local HRU
                          flux_data,                     & ! intent(inout): model fluxes for a local HRU
                          ! derivatives
                          dVolTot_dPsi0,                 & ! intent(in):    derivative in total volumetric water content w.r.t. matric head (m-1)
                          mLayerdTheta_dPsi,             & ! intent(in):    derivative in liquid water content w.r.t. matric potential (m-1)
                          mLayerdTheta_dTk,              & ! intent(in):    derivative in volumetric liquid water content w.r.t. temperature (K-1)
                          ! output: fluxes and derivatives
                          mLayerBaseflow,                & ! intent(out):   baseflow flux in each soil layer (m s-1)
                          dBaseflow_dWat,                & ! intent(out):   derivative in baseflow w.r.t. soil water characteristic
                          dBaseflow_dTk,                 & ! intent(out):   derivative in baseflow w.r.t. temperature (m s-1 K-1)
                          ! error handling
                          err, message)                    ! intent(out):   error control
  implicit none
  ! ---------------------------------------------------------------------------------------
  ! * dummy variables
  ! ---------------------------------------------------------------------------------------
  ! input: control and state variables
  integer(i4b),intent(in)          :: nSnow                   ! number of snow layers
  integer(i4b),intent(in)          :: nLake                   ! number of lake layers
  integer(i4b),intent(in)          :: nSoil                   ! number of soil layers
  integer(i4b),intent(in)          :: nGlce                   ! number of glacier ice layers
  integer(i4b),intent(in)          :: ixSaturation            ! index of upper-most "saturated" layer
  real(rkind),intent(in)           :: mLayerVolFracLiq(:)     ! volumetric fraction of liquid water (-)
  real(rkind),intent(in)           :: mLayerVolFracIce(:)     ! volumetric fraction of ice (-)
  ! input/output: data structures
  type(var_dlength),intent(in)     :: mpar_data               ! model parameters
  type(var_dlength),intent(in)     :: prog_data               ! prognostic variables for a local HRU
  type(var_dlength),intent(inout)  :: flux_data               ! model fluxes for a local HRU
  ! derivatives
  real(rkind),intent(in)           :: dVolTot_dPsi0(:)        ! derivative in total volumetric water content w.r.t. matric head (m-1)
  real(rkind),intent(in)           :: mLayerdTheta_dPsi(:)    ! derivative in liquid water content w.r.t. matric potential (m-1)
  real(rkind),intent(in)           :: mLayerdTheta_dTk(:)     ! derivative in volumetric liquid water content w.r.t. temperature (K-1)
  ! output: baseflow
  real(rkind),intent(out)          :: mLayerBaseflow(:)       ! baseflow from each soil layer (m s-1)
  real(rkind),intent(out)          :: dBaseflow_dWat(:,:)     ! derivative in baseflow w.r.t. soil water characteristic
  real(rkind),intent(out)          :: dBaseflow_dTk(:,:)      ! derivative in baseflow w.r.t. temperature (m s-1 K-1)
  ! error handling
  integer(i4b), intent(out)        :: err                     ! error code
  character(*), intent(out)        :: message                 ! error message
  ! ---------------------------------------------------------------------------------------
  ! * local variables
  ! ---------------------------------------------------------------------------------------
  ! general local variables
  integer(i4b)                       :: iLayer,jLayer         ! index of model layer
  ! local variables for the exfiltration
  real(rkind)                        :: totalColumnInflow     ! total column inflow (m s-1)
  real(rkind)                        :: totalColumnOutflow    ! total column outflow (m s-1)
  real(rkind)                        :: availStorage          ! available storage (m)
  real(rkind),parameter              :: xMinEval=0.002_rkind  ! minimum value to evaluate the exfiltration function (m)
  real(rkind),parameter              :: xCenter=0.001_rkind   ! center of the exfiltration function (m)
  real(rkind),parameter              :: xWidth=0.0001_rkind   ! width of the exfiltration function (m)
  real(rkind)                        :: expF,logF             ! logistic smoothing function (-)
  real(rkind)                        :: fieldCapacity_use     ! field capacity used to determine the "active" portion of the soil profile
  real(rkind)                        :: kAnisotropic_use     ! anisotropy factor used to compute transmissivity
  ! local variables for the lateral flux among soil columns
  real(rkind)                        :: activePorosity        ! "active" porosity associated with storage above a threshold (-)
  real(rkind)                        :: drainableWater        ! drainable water in each layer (m)
  real(rkind)                        :: tran0                 ! maximum transmissivity (m2 s-1)
  real(rkind)                        :: surfaceHydCond_use    ! macropore conductivity at the soil surface (m s-1)
  real(rkind)                        :: refDepth              ! reference depth for the power-law profile scaling (m)
  real(rkind)                        :: cDepth                ! compacted depth, limited to the soil column (m)
  real(rkind)                        :: scaleFacC             ! power-law scale factor at the compacted depth (-)
  real(rkind)                        :: wtDepth               ! depth to the water table (m)
  integer(i4b)                       :: ix_hc_profile         ! index for the choice of the hydraulic conductivity profile
  real(rkind),dimension(nSoil)       :: xTrans                ! dimensionless transmissivity profile, trTotal = tran0*xTrans (-)
  real(rkind),dimension(nSoil)       :: zActive               ! water table thickness associated with storage below and including the given layer (m)
  real(rkind),dimension(nSoil)       :: trTotal               ! total transmissivity associated with total water table depth zActive (m2 s-1)
  real(rkind),dimension(nSoil)       :: trSoil                ! transmissivity of water in a given layer (m2 s-1)
  ! local variables for the derivatives
  real(rkind),dimension(nSoil,nSoil) :: dBaseflow_dVolLiq     ! derivative in baseflow w.r.t. volumetric liquid water content (s-1)
  real(rkind)                        :: qbTotal               ! total baseflow (m s-1)
  real(rkind)                        :: length2area           ! ratio of hillslope width to hillslope area (m m-2)
  real(rkind),dimension(nSoil)       :: depth2capacity        ! ratio of layer depth to total subsurface storage capacity (-)
  real(rkind),dimension(nSoil)       :: dXdS                  ! change in dimensionless flux w.r.t. change in dimensionless storage (-)
  real(rkind),dimension(nSoil)       :: dLogFunc_dWat         ! derivative in the logistic function w.r.t. soil water characteristic
  real(rkind),dimension(nSoil)       :: dExfiltrate_dWat      ! derivative in exfiltration w.r.t. soil water characteristic
  real(rkind),dimension(nSoil)       :: dExfiltrate_dTk       ! derivative in exfiltration w.r.t. temperature (K-1)
  ! ---------------------------------------------------------------------------------------
  ! * association to data in structures
  ! ---------------------------------------------------------------------------------------
  associate(&
    ! input: coordinate variables
    soilDepth               => prog_data%var(iLookPROG%iLayerHeight)%dat(nSnow+nLake+nSoil),             & ! intent(in):  [dp]    total soil depth (m)
    mLayerDepth             => prog_data%var(iLookPROG%mLayerDepth)%dat(nSnow+nLake+1:nSnow+nLake+nSoil),& ! intent(in):  [dp(:)] depth of each soil layer (m)
    ! input: diagnostic variables
    surfaceHydCond          => flux_data%var(iLookFLUX%mLayerSatHydCondMP)%dat(1),       & ! intent(in):  [dp]    macropore conductivity at the first soil layer midpoint (m s-1)
    mLayerSatHydCond        => flux_data%var(iLookFLUX%mLayerSatHydCond)%dat,            & ! intent(in):  [dp(:)] micropore conductivity at the mid-point of each layer (m s-1)
    iLayerSatHydCond        => flux_data%var(iLookFLUX%iLayerSatHydCond)%dat,            & ! intent(in):  [dp(:)] micropore conductivity at layer interfaces, index 0 is the soil surface (m s-1)
    mLayerColumnInflow      => flux_data%var(iLookFLUX%mLayerColumnInflow)%dat,          & ! intent(in):  [dp(:)] inflow into each soil layer (m3/s)
    ! input: local attributes
    area                    => prog_data%var(iLookPROG%DOMarea)%dat(1),                  & ! intent(in):  [dp]    Domain area in HRU (m2)
    tan_slope               => prog_data%var(iLookPROG%DOMtan_slope)%dat(1),             & ! intent(in):  [dp]    tan water table slope, taken as tan local ground surface slope (-)
    contourLength           => prog_data%var(iLookPROG%DOMcontourLength)%dat(1),         & ! intent(in):  [dp]    length of contour at downslope edge of HRU (m)
    ! input: baseflow parameters
    zScale_TOPMODEL         => mpar_data%var(iLookPARAM%zScale_TOPMODEL)%dat(1),         & ! intent(in):  [dp]    TOPMODEL exponent (-)
    f_hydCond               => mpar_data%var(iLookPARAM%f_hydCond)%dat(1),               & ! intent(in):  [dp]    decay rate of hydraulic conductivity with depth (m-1)
    compactedDepth          => mpar_data%var(iLookPARAM%compactedDepth)%dat(1),          & ! intent(in):  [dp]    depth where k_soil reaches the compacted value (m)
    kAnisotropic            => mpar_data%var(iLookPARAM%kAnisotropic)%dat(1),            & ! intent(in):  [dp]    anisotropy factor for lateral hydraulic conductivity (-)
    fieldCapacity           => mpar_data%var(iLookPARAM%fieldCapacity)%dat(1),           & ! intent(in):  [dp]    field capacity (-)
    theta_sat               => mpar_data%var(iLookPARAM%theta_sat)%dat,                  & ! intent(in):  [dp(:)] soil porosity (-)
    ! output: diagnostic variables
    scalarExfiltration      => flux_data%var(iLookFLUX%scalarExfiltration)%dat(1),       & ! intent(out): [dp]    exfiltration from the soil profile (m s-1)
    mLayerColumnOutflow     => flux_data%var(iLookFLUX%mLayerColumnOutflow)%dat          & ! intent(out): [dp(:)] column outflow from each soil layer (m3 s-1)
    )  ! end association to variables in data structures
    ! -----------------------------------------------------------------------------------------------------------------------------------------
    ! initialize error control
    err=0; message="computBaseflow/"
    ! ***********************************************************************************************************************
    ! (1) compute the baseflow flux in each soil layer
    ! ***********************************************************************************************************************
    fieldCapacity_use = fieldCapacity
    kAnisotropic_use = kAnisotropic
    if(nGlce>0)then
      fieldCapacity_use = 0._rkind ! if glacier ice layers are present, set field capacity to zero (i.e. all water is "active" for flow)
      kAnisotropic_use = kAnisotropic*10._rkind ! if glacier ice layers are present, increase anisotropy factor to reflect higher hydraulic conductivity in glacier debris
    end if

    ! the transmissivity profile is the vertical integral of the hydraulic conductivity profile, so it must
    ! match the profile satHydCond used to build the conductivity itself
    ix_hc_profile = model_decisions(iLookDECISIONS%hc_profile)%iDecision
    if(nGlce>0) ix_hc_profile = expLaw_profile ! must match the override in satHydCond

    ! the transmissivity integrals below are written in terms of the conductivity at the soil surface, but
    ! mLayerSatHydCondMP(1) is the macropore value at the first layer midpoint, already scaled by the depth
    ! profile. satHydCond applies that same scaling to the micropore conductivity, and evaluates it at both
    ! the midpoint and the surface interface, so their ratio recovers the surface value for any profile.
    surfaceHydCond_use = surfaceHydCond*iLayerSatHydCond(0)/mLayerSatHydCond(1)

    ! compute the water table thickness (m) in each layer, working from the bottom of the profile up
    do iLayer=nSoil,ixSaturation,-1  ! loop through "active" soil layers, from lowest to highest
      ! define drainable water in each layer (m)
      activePorosity = theta_sat(iLayer) - fieldCapacity_use ! "active" porosity (-)
      drainableWater = mLayerDepth(iLayer)*(max(0._rkind,mLayerVolFracLiq(iLayer) - fieldCapacity_use))/activePorosity
      if (iLayer==nSoil) then
        zActive(iLayer) = drainableWater                        ! water table thickness associated with storage in a given layer (m)
      else
        zActive(iLayer) = zActive(iLayer+1) + drainableWater
      end if
    end do  ! end looping through soil layers

    ! set un-used portions of the water table thickness to zero, both profiles give xTrans=0 there
    if (ixSaturation>1) zActive(1:ixSaturation-1) = 0._rkind

    ! compute the maximum transmissivity (m2 s-1) and the dimensionless transmissivity profile xTrans(zActive),
    ! along with dXdS = d(xTrans)/d(zActive/soilDepth) used to build the derivative matrix below
    ! NOTE: tran0 can be done as a pre-processing step
    select case(ix_hc_profile)

      ! K(z) = K_0*exp(-f*z) integrated from the water table up to the base of the soil gives
      !  T(s) = (K_0/f)*[exp(-f*(D-s)) - exp(-f*D)], for saturated thickness s and soil depth D
      ! NOTE: written so that no exponential ever takes a positive argument
      case(expLaw_profile)
        tran0 = kAnisotropic_use*surfaceHydCond_use/f_hydCond
        xTrans(1:nSoil) = exp(-f_hydCond*(soilDepth - zActive(1:nSoil))) - exp(-f_hydCond*soilDepth)
        dXdS(1:nSoil)   = soilDepth*f_hydCond*exp(-f_hydCond*(soilDepth - zActive(1:nSoil)))

      ! power-law transmissivity, the vertical integral of the same floored profile satHydCond builds,
      !  K(z) = K_0*(1 - min(z,zc)/R)**(zScale_TOPMODEL-1), for zc = compactedDepth and R = refDepth.
      ! For saturated thickness s, water table depth d = D-s, and Kc the scale factor at zc, that integrates to
      !  s <= D-zc:  T = K_0*Kc*s                                                  (water table below the scaling zone)
      !  s >  D-zc:  T = K_0*Kc*(D-zc) + K_0*(R/n)*[(1-d/R)**n - (1-zc/R)**n]
      ! NOTE: zc < R always, so 1-min(z,zc)/R is bounded away from zero and this never evaluates 0**(n-1). That is why
      !       satHydCond floors the profile in the first place, so the floor has to be carried through here too
      case(powerLaw_profile)
        refDepth  = soilDepth
        if (soilDepth < compactedDepth) refDepth = compactedDepth + 1._rkind ! as in satHydCond
        cDepth    = min(compactedDepth, soilDepth)
        scaleFacC = (1._rkind - cDepth/refDepth)**(zScale_TOPMODEL - 1._rkind)
        tran0     = kAnisotropic_use*surfaceHydCond_use*soilDepth
        do iLayer=1,nSoil
          if (zActive(iLayer) <= soilDepth - cDepth) then ! water table at or below the compacted depth, uniform conductivity
            xTrans(iLayer) = scaleFacC*zActive(iLayer)/soilDepth
            dXdS(iLayer)   = scaleFacC
          else                                            ! water table up inside the scaling zone
            wtDepth        = soilDepth - zActive(iLayer)
            xTrans(iLayer) = scaleFacC*(soilDepth - cDepth)/soilDepth                                    &
                             + ( refDepth/(zScale_TOPMODEL*soilDepth) )                                  &
                               *( (1._rkind - wtDepth/refDepth)**zScale_TOPMODEL                         &
                                - (1._rkind -  cDepth/refDepth)**zScale_TOPMODEL )
            dXdS(iLayer)   = (1._rkind - wtDepth/refDepth)**(zScale_TOPMODEL - 1._rkind)
          end if
        end do

      ! uniform conductivity with depth, so transmissivity is simply linear in the saturated thickness
      ! NOTE: grouped here only for completeness, mDecisions does not allow constant with qbaseTopmodel
      case(constant)
        tran0 = kAnisotropic_use*surfaceHydCond_use*soilDepth
        xTrans(1:nSoil) = zActive(1:nSoil)/soilDepth
        dXdS(1:nSoil)   = 1._rkind

      case default
        message=trim(message)//"unknown hydraulic conductivity profile for the baseflow transmissivity"
        err=20; return

    end select

    ! compute the transmissivity of each layer (m2 s-1)
    trTotal(1:nSoil) = tran0*xTrans(1:nSoil)                    ! total transmissivity for total depth zActive (m2 s-1)
    trSoil(nSoil)    = trTotal(nSoil)                           ! transmissivity of water in a given layer (m2 s-1)
    do iLayer=nSoil-1,1,-1
      trSoil(iLayer) = trTotal(iLayer) - trTotal(iLayer+1)
    end do

    ! set un-used portions of the vectors to zero
    if (ixSaturation>1) trSoil(1:ixSaturation-1) = 0._rkind

    ! compute the outflow from each layer (m3 s-1)
    mLayerColumnOutflow(1:nSoil) = trSoil(1:nSoil)*tan_slope*contourLength

    ! compute total column inflow and total column outflow (m s-1)
    totalColumnInflow  = sum(mLayerColumnInflow(1:nSoil))/area
    totalColumnOutflow = sum(mLayerColumnOutflow(1:nSoil))/area

    ! compute the available storage (m)
    availStorage = sum(mLayerDepth(1:nSoil)*(theta_sat(1:nSoil) - (mLayerVolFracLiq(1:nSoil)+mLayerVolFracIce(1:nSoil))))

    ! compute the smoothing function (-)
    if (availStorage < xMinEval) then
      ! compute the logistic function
      expF = exp((availStorage - xCenter)/xWidth)
      logF = 1._rkind / (1._rkind + expF)
      ! compute the derivative in the logistic function w.r.t. volumetric water content in each soil layer, NOTE dLogFunc_dTemp = 0
      dLogFunc_dWat(1:nSoil) = mLayerDepth(1:nSoil)*(expF/xWidth)/(1._rkind + expF)**2_i4b * dVolTot_dPsi0(1:nSoil)
    else
      logF             = 0._rkind
      dLogFunc_dWat(:) = 0._rkind
    end if

    ! compute the exfiltration (m s-1)
    if (totalColumnInflow > totalColumnOutflow .and. logF > tiny(1._rkind)) then
      scalarExfiltration = logF*(totalColumnInflow - totalColumnOutflow)  ! m s-1
    else
      scalarExfiltration = 0._rkind
    end if

    ! compute the baseflow in each layer (m s-1)
    mLayerBaseflow(1:nSoil) = (mLayerColumnOutflow(1:nSoil) - mLayerColumnInflow(1:nSoil))/area

    ! compute the total baseflow
    qbTotal = sum(mLayerBaseflow)

    ! add exfiltration to the baseflow flux at the top layer
    mLayerBaseflow(1)      = mLayerBaseflow(1) + scalarExfiltration
    mLayerColumnOutflow(1) = mLayerColumnOutflow(1) + scalarExfiltration*area

    ! ***********************************************************************************************************************
    ! (2) compute the derivative in the baseflow flux w.r.t. volumetric liquid water content (m s-1)
    ! ***********************************************************************************************************************

    ! initialize the derivative matrix
    dBaseflow_dVolLiq(:,:) = 0._rkind
    dBaseflow_dWat(:,:) = 0._rkind
    dBaseflow_dTk(:,:) = 0._rkind

    ! compute ratio of hillslope width to hillslope area (m m-2)
    length2area = tan_slope*contourLength/area

    ! compute the ratio of layer depth to maximum water holding capacity (-)
    depth2capacity(1:nSoil) = mLayerDepth(1:nSoil)/(theta_sat(1:nSoil) - fieldCapacity_use)/soilDepth
    do iLayer=1,nSoil
      if (mLayerVolFracLiq(iLayer) <= fieldCapacity_use) depth2capacity(iLayer) = 0._rkind
    end do

    ! loop through soil layers
    do iLayer=1,nSoil
      ! compute diagonal terms (s-1)
      dBaseflow_dVolLiq(iLayer,iLayer) = tran0*dXdS(iLayer)*depth2capacity(iLayer)*length2area
      dBaseflow_dWat(iLayer,iLayer) = dBaseflow_dVolLiq(iLayer,iLayer)*mLayerdTheta_dPsi(iLayer)
      dBaseflow_dTk(iLayer,iLayer) = dBaseflow_dVolLiq(iLayer,iLayer)*mLayerdTheta_dTk(iLayer)
      ! compute off-diagonal terms
      do jLayer=iLayer+1,nSoil  ! only dependent on layers below
        dBaseflow_dVolLiq(iLayer,jLayer) = tran0*(dXdS(iLayer) - dXdS(iLayer+1))*depth2capacity(jLayer)*length2area
        dBaseflow_dWat(iLayer,jLayer) = dBaseflow_dVolLiq(iLayer,jLayer)*mLayerdTheta_dPsi(jLayer)
        dBaseflow_dTk(iLayer,jLayer) = dBaseflow_dVolLiq(iLayer,jLayer)*mLayerdTheta_dTk(jLayer)
      end do  ! end looping through soil layers
    end do  ! end looping through soil layers

    ! trSoil is clamped to zero above the saturated zone, so the outflow of those layers does not respond to
    ! the state at all and their derivative rows must be zero to match. Without this the rows depend on the
    ! shape of xTrans near zActive=0, which differs between profiles (and is non-zero for the exponential,
    ! and for the power law whenever zScale_TOPMODEL=1).
    ! NOTE: this is done before the exfiltration derivative, which is a real flux and still belongs in row 1
    if (ixSaturation>1) then
      dBaseflow_dVolLiq(1:ixSaturation-1,:) = 0._rkind
      dBaseflow_dWat(1:ixSaturation-1,:)    = 0._rkind
      dBaseflow_dTk(1:ixSaturation-1,:)     = 0._rkind
    end if

    ! compute the derivative in the exfiltration flux and add to the baseflow derivative matrix
    if (totalColumnInflow > totalColumnOutflow .and. logF > tiny(1._rkind)) then
      do iLayer=1,nSoil
        dExfiltrate_dWat(iLayer) = -sum(dBaseflow_dWat(1:nSoil,iLayer))*logF - dLogFunc_dWat(iLayer)*qbTotal
        dExfiltrate_dTk(iLayer) = -sum(dBaseflow_dTk(1:nSoil,iLayer))*logF
      end do  ! end looping through soil layers
      dBaseflow_dWat(1,1:nSoil) = dBaseflow_dWat(1,1:nSoil) + dExfiltrate_dWat(1:nSoil)
      dBaseflow_dTk(1,1:nSoil) = dBaseflow_dTk(1,1:nSoil) + dExfiltrate_dTk(1:nSoil)
    end if

  end associate ! end association to data in structures

end subroutine computBaseflow

end module groundwatr_module

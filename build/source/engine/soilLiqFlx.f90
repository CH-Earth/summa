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

module soilLiqFlx_module
! -----------------------------------------------------------------------------------------------------------

! data types
USE nrtype
USE data_types,only:var_d                  ! x%var(:)       (rkind)
USE data_types,only:var_ilength            ! x%var(:)%dat   (i4b)
USE data_types,only:var_dlength            ! x%var(:)%dat   (rkind)

! missing values
USE globalData,only:integerMissing  ! missing integer
USE globalData,only:realMissing     ! missing real number

! physical constants
USE multiconst,only:&
                    LH_fus,  & ! latent heat of fusion         (J kg-1)
                    LH_vap,  & ! latent heat of vaporization   (J kg-1)
                    LH_sub,  & ! latent heat of sublimation    (J kg-1)
                    gravity, & ! gravitational acceleteration  (m s-2)
                    Tfreeze, & ! freezing point of pure water  (K)
                    iden_air,& ! intrinsic density of air      (kg m-3)
                    iden_ice,& ! intrinsic density of ice      (kg m-3)
                    iden_water ! intrinsic density of water    (kg m-3)

! named variables
USE var_lookup,only:iLookPROG              ! named variables for structure elements
USE var_lookup,only:iLookDIAG              ! named variables for structure elements
USE var_lookup,only:iLookFLUX              ! named variables for structure elements
USE var_lookup,only:iLookPARAM             ! named variables for structure elements
USE var_lookup,only:iLookINDEX             ! named variables for structure elements

! model decisions
USE globalData,only:model_decisions        ! model decision structure
USE var_lookup,only:iLookDECISIONS         ! named variables for elements of the decision structure

! provide access to look-up values for model decisions
USE mDecisions_module,only:  &
  ! look-up values for method used to compute derivative
  numerical,                  & ! numerical solution
  analytical,                 & ! analytical solution
  ! look-up values for the form of Richards' equation
  moisture,                   & ! moisture-based form of Richards' equation
  mixdform,                   & ! mixed form of Richards' equation
  ! look-up values for the type of hydraulic conductivity profile
  constant,                   & ! constant hydraulic conductivity with depth
  powerLaw_profile,           & ! power-law profile
  ! look-up values for the choice of groundwater parameterization
  qbaseTopmodel,              & ! TOPMODEL-ish baseflow parameterization
  bigBucket,                  & ! a big bucket (lumped aquifer model)
  noExplicit,                 & ! no explicit groundwater parameterization
  ! look-up values for the choice of boundary conditions for hydrology
  prescribedHead,             & ! prescribed head (volumetric liquid water content for mixed form of Richards' eqn)
  funcBottomHead,             & ! function of matric head in the lower-most layer
  freeDrainage,               & ! free drainage
  liquidFlux,                 & ! liquid water flux
  zeroFlux                      ! zero flux

! -----------------------------------------------------------------------------------------------------------
implicit none
private
public::soilLiqFlx
! constant parameters
real(rkind),parameter     :: verySmall=1.e-12_rkind       ! a very small number (used to avoid divide by zero)
real(rkind),parameter     :: dx=1.e-8_rkind               ! finite difference increment
contains


! ***************************************************************************************************************
! public subroutine soilLiqFlx: compute liquid water fluxes and their derivatives
! ***************************************************************************************************************
subroutine soilLiqFlx(&
                      ! input: model control
                      nSoil,                        & ! intent(in): number of soil layers
                      firstSplitOper,               & ! intent(in): flag to compute infiltration, if firstSplitOper
                      scalarSolution,               & ! intent(in): flag to indicate the scalar solution
                      deriv_desired,                & ! intent(in): flag indicating if derivatives are desired
                      ! input: trial state variables
                      mLayerTempTrial,              & ! intent(in): temperature (K)
                      mLayerMatricHeadTrial,        & ! intent(in): matric head (m)
                      mLayerMatricHeadLiqTrial,     & ! intent(in): liquid matric head (m)
                      mLayerVolFracLiqTrial,        & ! intent(in): volumetric fraction of liquid water (-)
                      mLayerVolFracIceTrial,        & ! intent(in): volumetric fraction of ice (-)
                      ! input: pre-computed derivatives
                      mLayerdTheta_dTk,             & ! intent(in): derivative in volumetric liquid water content w.r.t. temperature (K-1)
                      dPsiLiq_dTemp,                & ! intent(in): derivative in liquid water matric potential w.r.t. temperature (m K-1)
                      dCanopyTrans_dCanWat,         & ! intent(in): derivative in canopy transpiration w.r.t. canopy total water content (s-1)
                      dCanopyTrans_dTCanair,        & ! intent(in): derivative in canopy transpiration w.r.t. canopy air temperature (kg m-2 s-1 K-1)
                      dCanopyTrans_dTCanopy,        & ! intent(in): derivative in canopy transpiration w.r.t. canopy temperature (kg m-2 s-1 K-1)
                      dCanopyTrans_dTGround,        & ! intent(in): derivative in canopy transpiration w.r.t. ground temperature (kg m-2 s-1 K-1)
                      above_soilLiqFluxDeriv,       & ! intent(in): derivative in layer above soil (canopy or snow) liquid flux w.r.t. liquid water
                      above_soildLiq_dTk,           & ! intent(in): derivative of layer above soil (canopy or snow) liquid flux w.r.t. temperature
                      above_soilFracLiq,            & ! intent(in): fraction of liquid water layer above soil (canopy or snow) (-)
                      ! input: fluxes
                      scalarCanopyTranspiration,    & ! intent(in): canopy transpiration (kg m-2 s-1)
                      scalarGroundEvaporation,      & ! intent(in): ground evaporation (kg m-2 s-1)
                      scalarRainPlusMelt,           & ! intent(in): rain plus melt (m s-1)
                      ! input-output: data structures
                      mpar_data,                    & ! intent(in):    model parameters
                      indx_data,                    & ! intent(in):    model indices
                      prog_data,                    & ! intent(in):    model prognostic variables for a local HRU
                      diag_data,                    & ! intent(inout): model diagnostic variables for a local HRU
                      flux_data,                    & ! intent(inout): model fluxes for a local HRU
                      ! output: diagnostic variables for surface runoff
                      xMaxInfilRate,                & ! intent(inout): maximum infiltration rate (m s-1)
                      scalarInfilArea,              & ! intent(inout): fraction of unfrozen area where water can infiltrate (-)
                      scalarFrozenArea,             & ! intent(inout): fraction of area that is considered impermeable due to soil ice (-)
                      scalarSurfaceRunoff,          & ! intent(out): surface runoff (m s-1)
                      ! output: diagnostic variables for model layers
                      mLayerdTheta_dPsi,            & ! intent(out): derivative in the soil water characteristic w.r.t. psi (m-1)
                      mLayerdPsi_dTheta,            & ! intent(out): derivative in the soil water characteristic w.r.t. theta (m)
                      dHydCond_dMatric,             & ! intent(out): derivative in hydraulic conductivity w.r.t matric head (s-1)
                      ! output: fluxes
                      scalarSurfaceInfiltration,    & ! intent(out): surface infiltration rate (m s-1)
                      iLayerLiqFluxSoil,            & ! intent(out): liquid fluxes at layer interfaces (m s-1)
                      mLayerTranspire,              & ! intent(out): transpiration loss from each soil layer (m s-1)
                      mLayerHydCond,                & ! intent(out): hydraulic conductivity in each soil layer (m s-1)
                      ! output: derivatives in fluxes w.r.t. hydrology state variables -- matric head or volumetric lquid water -- in the layer above and layer below (m s-1 or s-1)
                      dq_dHydStateAbove,            & ! intent(out): derivatives in the flux w.r.t. volumetric liquid water content in the layer above (m s-1)
                      dq_dHydStateBelow,            & ! intent(out): derivatives in the flux w.r.t. volumetric liquid water content in the layer below (m s-1)
                      dq_dHydStateLayerSurfVec,     & ! intent(out): derivative in surface infiltration w.r.t. hydrology state in above soil snow or canopy and every soil layer  (m s-1 or s-1)
                      ! output: derivatives in fluxes w.r.t. energy state variables -- now just temperature -- in the layer above and layer below (m s-1 K-1)
                      dq_dNrgStateAbove,            & ! intent(out): derivatives in the flux w.r.t. temperature in the layer above (m s-1 K-1)
                      dq_dNrgStateBelow,            & ! intent(out): derivatives in the flux w.r.t. temperature in the layer below (m s-1 K-1)
                      dq_dNrgStateLayerSurfVec,     & ! intent(out): derivative in surface infiltration w.r.t. energy state in above soil snow or canopy and every soil layer (m s-1 K-1)
                      ! output: derivatives in transpiration w.r.t. canopy state variables
                      mLayerdTrans_dTCanair,        & ! intent(out): derivatives in the soil layer transpiration flux w.r.t. canopy air temperature
                      mLayerdTrans_dTCanopy,        & ! intent(out): derivatives in the soil layer transpiration flux w.r.t. canopy temperature
                      mLayerdTrans_dTGround,        & ! intent(out): derivatives in the soil layer transpiration flux w.r.t. ground temperature
                      mLayerdTrans_dCanWat,         & ! intent(out): derivatives in the soil layer transpiration flux w.r.t. canopy total water
                      ! output: error control
                      err,message)                    ! intent(out): error control
  ! utility modules
  USE soil_utils_module,only:volFracLiq      ! compute volumetric fraction of liquid water
  USE soil_utils_module,only:matricHead      ! compute matric head (m)
  USE soil_utils_module,only:dTheta_dPsi     ! compute derivative of the soil moisture characteristic w.r.t. psi (m-1)
  USE soil_utils_module,only:dPsi_dTheta     ! compute derivative of the soil moisture characteristic w.r.t. theta (m)
  USE soil_utils_module,only:hydCond_psi     ! compute hydraulic conductivity as a function of matric head
  USE soil_utils_module,only:hydCond_liq     ! compute hydraulic conductivity as a function of volumetric liquid water content
  USE soil_utils_module,only:hydCondMP_liq   ! compute hydraulic conductivity of macropores as a function of volumetric liquid water content
  ! -------------------------------------------------------------------------------------------------------------------------------------------------
  implicit none
  ! input: model control
  integer(i4b),intent(in)             :: nSoil                         ! number of soil layers
  logical(lgt),intent(in)             :: firstSplitOper                ! flag to compute infiltration
  logical(lgt),intent(in)             :: scalarSolution                ! flag to denote if implementing the scalar solution
  logical(lgt),intent(in)             :: deriv_desired                 ! flag indicating if derivatives are desired
  ! input: trial model state variables
  real(rkind),intent(in)              :: mLayerTempTrial(:)            ! temperature in each layer at the current iteration (m)
  real(rkind),intent(in)              :: mLayerMatricHeadTrial(:)      ! matric head in each layer at the current iteration (m)
  real(rkind),intent(in)              :: mLayerMatricHeadLiqTrial(:)   ! liquid matric head in each layer at the current iteration (m)
  real(rkind),intent(in)              :: mLayerVolFracLiqTrial(:)      ! volumetric fraction of liquid water at the current iteration (-)
  real(rkind),intent(in)              :: mLayerVolFracIceTrial(:)      ! volumetric fraction of ice at the current iteration (-)
  ! input: pre-computed derivatves
  real(rkind),intent(in)              :: mLayerdTheta_dTk(:)           ! derivative in volumetric liquid water content w.r.t. temperature (K-1)
  real(rkind),intent(in)              :: dPsiLiq_dTemp(:)              ! derivative in liquid water matric potential w.r.t. temperature (m K-1)
  real(rkind),intent(in)              :: dCanopyTrans_dCanWat          ! derivative in canopy transpiration w.r.t. canopy total water content (s-1)
  real(rkind),intent(in)              :: dCanopyTrans_dTCanair         ! derivative in canopy transpiration w.r.t. canopy air temperature (kg m-2 s-1 K-1)
  real(rkind),intent(in)              :: dCanopyTrans_dTCanopy         ! derivative in canopy transpiration w.r.t. canopy temperature (kg m-2 s-1 K-1)
  real(rkind),intent(in)              :: dCanopyTrans_dTGround         ! derivative in canopy transpiration w.r.t. ground temperature (kg m-2 s-1 K-1)
  real(rkind),intent(in)              :: above_soilLiqFluxDeriv        ! derivative in layer above soil (canopy or snow) liquid flux w.r.t. liquid water
  real(rkind),intent(in)              :: above_soildLiq_dTk            ! derivative of layer above soil (canopy or snow) liquid flux w.r.t. temperature
  real(rkind),intent(in)              :: above_soilFracLiq             ! fraction of liquid water layer above soil (canopy or snow) (-)
  ! input: model fluxes
  real(rkind),intent(in)              :: scalarCanopyTranspiration     ! canopy transpiration (kg m-2 s-1)
  real(rkind),intent(in)              :: scalarGroundEvaporation       ! ground evaporation (kg m-2 s-1)
  real(rkind),intent(in)              :: scalarRainPlusMelt            ! rain plus melt (m s-1)
  ! input-output: data structures
  type(var_dlength),intent(in)        :: mpar_data                     ! model parameters
  type(var_ilength),intent(in)        :: indx_data                     ! state vector geometry
  type(var_dlength),intent(in)        :: prog_data                     ! prognostic variables for a local HRU
  type(var_dlength),intent(inout)     :: diag_data                     ! diagnostic variables for a local HRU
  type(var_dlength),intent(inout)     :: flux_data                     ! model fluxes for a local HRU
  ! output: diagnostic variables for surface runoff
  real(rkind),intent(inout)           :: xMaxInfilRate                 ! maximum infiltration rate (m s-1)
  real(rkind),intent(inout)           :: scalarInfilArea               ! fraction of unfrozen area where water can infiltrate (-)
  real(rkind),intent(inout)           :: scalarFrozenArea              ! fraction of area that is considered impermeable due to soil ice (-)
  real(rkind),intent(inout)           :: scalarSurfaceRunoff           ! surface runoff (m s-1)
  ! output: diagnostic variables for each layer
  real(rkind),intent(inout)           :: mLayerdTheta_dPsi(:)          ! derivative in the soil water characteristic w.r.t. psi (m-1)
  real(rkind),intent(inout)           :: mLayerdPsi_dTheta(:)          ! derivative in the soil water characteristic w.r.t. theta (m)
  real(rkind),intent(inout)           :: dHydCond_dMatric(:)           ! derivative in hydraulic conductivity w.r.t matric head (s-1)
  ! output: liquid fluxes
  real(rkind),intent(inout)           :: scalarSurfaceInfiltration     ! surface infiltration rate (m s-1)
  real(rkind),intent(inout)           :: iLayerLiqFluxSoil(0:)         ! liquid flux at soil layer interfaces (m s-1)
  real(rkind),intent(inout)           :: mLayerTranspire(:)            ! transpiration loss from each soil layer (m s-1)
  real(rkind),intent(inout)           :: mLayerHydCond(:)              ! hydraulic conductivity in each soil layer (m s-1)
  ! output: derivatives in fluxes w.r.t. state variables in the layer above and layer below (m s-1)
  real(rkind),intent(inout)           :: dq_dHydStateAbove(0:)         ! derivative in the flux in layer interfaces w.r.t. state variables in the layer above
  real(rkind),intent(inout)           :: dq_dHydStateBelow(0:)         ! derivative in the flux in layer interfaces w.r.t. state variables in the layer below
  real(rkind),intent(inout)           :: dq_dHydStateLayerSurfVec(0:)  ! derivative in surface infiltration w.r.t. hydrology state in above soil snow or canopy and every soil layer  (m s-1 or s-1)
  ! output: derivatives in fluxes w.r.t. energy state variables -- now just temperature -- in the layer above and layer below (m s-1 K-1)
  real(rkind),intent(inout)           :: dq_dNrgStateAbove(0:)         ! derivatives in the flux w.r.t. temperature in the layer above (m s-1 K-1)
  real(rkind),intent(inout)           :: dq_dNrgStateBelow(0:)         ! derivatives in the flux w.r.t. temperature in the layer below (m s-1 K-1)
  real(rkind),intent(inout)           :: dq_dNrgStateLayerSurfVec(0:)  ! derivative in surface infiltration w.r.t. temperature in above soil snow or canopy and every soil layer  (m s-1 or s-1)
  ! output: derivatives in transpiration w.r.t. canopy state variables
  real(rkind),intent(inout)           :: mLayerdTrans_dTCanair(:)      ! derivatives in the soil layer transpiration flux w.r.t. canopy air temperature
  real(rkind),intent(inout)           :: mLayerdTrans_dTCanopy(:)      ! derivatives in the soil layer transpiration flux w.r.t. canopy temperature
  real(rkind),intent(inout)           :: mLayerdTrans_dTGround(:)      ! derivatives in the soil layer transpiration flux w.r.t. ground temperature
  real(rkind),intent(inout)           :: mLayerdTrans_dCanWat(:)       ! derivatives in the soil layer transpiration flux w.r.t. canopy total water
  ! output: error control
  integer(i4b),intent(out)            :: err                           ! error code
  character(*),intent(out)            :: message                       ! error message
  ! -----------------------------------------------------------------------------------------------------------------------------------------------------
  ! local variables: general
  character(LEN=256)                  :: cmessage                     ! error message of downwind routine
  integer(i4b)                        :: ibeg,iend                    ! start and end indices of the soil layers in concatanated snow-soil vector
  integer(i4b)                        :: iLayer,iSoil                 ! index of soil layer
  integer(i4b)                        :: ixLayerDesired(1)            ! layer desired (scalar solution)
  integer(i4b)                        :: ixTop                        ! top layer in subroutine call
  integer(i4b)                        :: ixBot                        ! bottom layer in subroutine call
  ! transpiration sink term
  real(rkind),dimension(nSoil)        :: mLayerTranspireFrac          ! fraction of transpiration allocated to each soil layer (-)
  ! diagnostic variables
  real(rkind),dimension(nSoil)        :: iceImpedeFac                 ! ice impedence factor at layer mid-points (-)
  real(rkind),dimension(nSoil)        :: mLayerDiffuse                ! diffusivity at layer mid-point (m2 s-1)
  real(rkind),dimension(nSoil)        :: dHydCond_dVolLiq             ! derivative in hydraulic conductivity w.r.t volumetric liquid water content (m s-1)
  real(rkind),dimension(nSoil)        :: dDiffuse_dVolLiq             ! derivative in hydraulic diffusivity w.r.t volumetric liquid water content (m2 s-1)
  real(rkind),dimension(nSoil)        :: dHydCond_dTemp               ! derivative in hydraulic conductivity w.r.t temperature (m s-1 K-1)
  real(rkind),dimension(0:nSoil)      :: iLayerHydCond                ! hydraulic conductivity at layer interface (m s-1)
  real(rkind),dimension(0:nSoil)      :: iLayerDiffuse                ! diffusivity at layer interface (m2 s-1)
  ! compute surface flux
  integer(i4b)                        :: nRoots                       ! number of soil layers with roots
  integer(i4b)                        :: ixIce                        ! index of the lowest soil layer that contains ice
  real(rkind),dimension(0:nSoil)      :: iLayerHeight                 ! height of the layer interfaces (m)
  ! compute fluxes and derivatives at layer interfaces
   real(rkind)                         :: scalardPsi_dTheta            ! derivative in soil water characteristix, used for perturbations when computing numerical derivatives
  ! -------------------------------------------------------------------------------------------------------------------------------------------------
  ! initialize error control
  err=0; message='soilLiqFlx/'

  ! get indices for the data structures
  ibeg = indx_data%var(iLookINDEX%nSnow)%dat(1) + 1
  iend = indx_data%var(iLookINDEX%nSnow)%dat(1) + indx_data%var(iLookINDEX%nSoil)%dat(1)

  ! get a copy of iLayerHeight
  ! NOTE: performance hit, though cannot define the shape (0:) with the associate construct
  iLayerHeight(0:nSoil) = prog_data%var(iLookPROG%iLayerHeight)%dat(ibeg-1:iend)  ! height of the layer interfaces (m)

  ! make association between local variables and the information in the data structures
  associate(&
    ! input: model control
    ixRichards             => model_decisions(iLookDECISIONS%f_Richards)%iDecision,   & ! intent(in): index of the form of Richards' equation
    ixBcUpperSoilHydrology => model_decisions(iLookDECISIONS%bcUpprSoiH)%iDecision,   & ! intent(in): index of the upper boundary conditions for soil hydrology
    ixBcLowerSoilHydrology => model_decisions(iLookDECISIONS%bcLowrSoiH)%iDecision,   & ! intent(in): index of the lower boundary conditions for soil hydrology
    ! input: model indices
    ixMatricHead           => indx_data%var(iLookINDEX%ixMatricHead)%dat,             & ! intent(in): indices of soil layers where matric head is the state variable
    ixSoilOnlyHyd          => indx_data%var(iLookINDEX%ixSoilOnlyHyd)%dat,            & ! intent(in): index in the state subset for hydrology state variables in the soil domain
    ! input: model coordinate variables -- NOTE: use of ibeg and iend
    mLayerDepth            => prog_data%var(iLookPROG%mLayerDepth)%dat(ibeg:iend),    & ! intent(in): depth of the layer (m)
    mLayerHeight           => prog_data%var(iLookPROG%mLayerHeight)%dat(ibeg:iend),   & ! intent(in): height of the layer mid-point (m)
    ! input: upper boundary conditions
    upperBoundHead         => mpar_data%var(iLookPARAM%upperBoundHead)%dat(1),        & ! intent(in): upper boundary condition for matric head (m)
    upperBoundTheta        => mpar_data%var(iLookPARAM%upperBoundTheta)%dat(1),       & ! intent(in): upper boundary condition for volumetric liquid water content (-)
    ! input: lower boundary conditions
    lowerBoundHead         => mpar_data%var(iLookPARAM%lowerBoundHead)%dat(1),        & ! intent(in): lower boundary condition for matric head (m)
    lowerBoundTheta        => mpar_data%var(iLookPARAM%lowerBoundTheta)%dat(1),       & ! intent(in): lower boundary condition for volumetric liquid water content (-)
    ! input: vertically variable soil parameters
    vGn_m                  => diag_data%var(iLookDIAG%scalarVGn_m)%dat,               & ! intent(in): van Genutchen "m" parameter (-)
    vGn_n                  => mpar_data%var(iLookPARAM%vGn_n)%dat,                    & ! intent(in): van Genutchen "n" parameter (-)
    vGn_alpha              => mpar_data%var(iLookPARAM%vGn_alpha)%dat,                & ! intent(in): van Genutchen "alpha" parameter (m-1)
    theta_sat              => mpar_data%var(iLookPARAM%theta_sat)%dat,                & ! intent(in): soil porosity (-)
    theta_res              => mpar_data%var(iLookPARAM%theta_res)%dat,                & ! intent(in): soil residual volumetric water content (-)
    ! input: vertically constant soil parameters
    wettingFrontSuction    => mpar_data%var(iLookPARAM%wettingFrontSuction)%dat(1),   & ! intent(in): Green-Ampt wetting front suction (m)
    rootingDepth           => mpar_data%var(iLookPARAM%rootingDepth)%dat(1),          & ! intent(in): rooting depth (m)
    kAnisotropic           => mpar_data%var(iLookPARAM%kAnisotropic)%dat(1),          & ! intent(in): anisotropy factor for lateral hydraulic conductivity (-)
    zScale_TOPMODEL        => mpar_data%var(iLookPARAM%zScale_TOPMODEL)%dat(1),       & ! intent(in): TOPMODEL scaling factor (m)
    qSurfScale             => mpar_data%var(iLookPARAM%qSurfScale)%dat(1),            & ! intent(in): scaling factor in the surface runoff parameterization (-)
    f_impede               => mpar_data%var(iLookPARAM%f_impede)%dat(1),              & ! intent(in): ice impedence factor (-)
    soilIceScale           => mpar_data%var(iLookPARAM%soilIceScale)%dat(1),          & ! intent(in): scaling factor for depth of soil ice, used to get frozen fraction (m)
    soilIceCV              => mpar_data%var(iLookPARAM%soilIceCV)%dat(1),             & ! intent(in): CV of depth of soil ice, used to get frozen fraction (-)
    theta_mp               => mpar_data%var(iLookPARAM%theta_mp)%dat(1),              & ! intent(in): volumetric liquid water content when macropore flow begins (-)
    mpExp                  => mpar_data%var(iLookPARAM%mpExp)%dat(1),                 & ! intent(in): empirical exponent in macropore flow equation (-)
    ! input: saturated hydraulic conductivity
    mLayerSatHydCondMP     => flux_data%var(iLookFLUX%mLayerSatHydCondMP)%dat,        & ! intent(in): saturated hydraulic conductivity of macropores at the mid-point of each layer (m s-1)
    mLayerSatHydCond       => flux_data%var(iLookFLUX%mLayerSatHydCond)%dat,          & ! intent(in): saturated hydraulic conductivity at the mid-point of each layer (m s-1)
    iLayerSatHydCond       => flux_data%var(iLookFLUX%iLayerSatHydCond)%dat,          & ! intent(in): saturated hydraulic conductivity at the interface of each layer (m s-1)
    ! input: factors limiting transpiration (from vegFlux routine)
    mLayerRootDensity      => diag_data%var(iLookDIAG%mLayerRootDensity)%dat,         & ! intent(in): root density in each layer (-)
    scalarTranspireLim     => diag_data%var(iLookDIAG%scalarTranspireLim)%dat(1),     & ! intent(in): weighted average of the transpiration limiting factor (-)
    mLayerTranspireLim     => diag_data%var(iLookDIAG%mLayerTranspireLim)%dat         & ! intent(in): transpiration limiting factor in each layer (-)
    )  ! associating local variables with the information in the data structures

    ! -------------------------------------------------------------------------------------------------------------------------------------------------
    ! preliminaries
    ! -------------------------------------------------------------------------------------------------------------------------------------------------

    ! get the indices for the soil layers
    if(scalarSolution)then
      ixLayerDesired = pack(ixMatricHead, ixSoilOnlyHyd/=integerMissing)
      ixTop = ixLayerDesired(1)
      ixBot = ixLayerDesired(1)
    else
      ixTop = 1
      ixBot = nSoil
    endif

    ! identify the number of layers that contain roots
    nRoots = count(iLayerHeight(0:nSoil-1) < rootingDepth-verySmall)
    if(nRoots==0)then
      message=trim(message)//'no layers with roots'
      err=20; return
    endif

    ! identify lowest soil layer with ice
    ! NOTE: cannot use count because there may be an unfrozen wedge
    ixIce = 0  ! initialize the index of the ice layer (0 means no ice in the soil profile)
    do iLayer=1,nSoil ! (loop through soil layers)
      if(mLayerVolFracIceTrial(iLayer) > verySmall) ixIce = iLayer
    end do

    ! -------------------------------------------------------------------------------------------------------------------------------------------------
    ! compute the transpiration sink term
    ! -------------------------------------------------------------------------------------------------------------------------------------------------

    ! check the need to compute transpiration (NOTE: intent=inout)
    if( .not. (scalarSolution .and. ixTop>1) )then
      
      ! compute the fraction of transpiration loss from each soil layer
      if(scalarTranspireLim > tiny(scalarTranspireLim))then ! (transpiration may be non-zero even if the soil moisture limiting factor is zero)
        mLayerTranspireFrac(:) = mLayerRootDensity(:)*mLayerTranspireLim(:)/scalarTranspireLim
      else ! (possible for there to be non-zero conductance and therefore transpiration in this case)
        mLayerTranspireFrac(:) = mLayerRootDensity(:) / sum(mLayerRootDensity)
      end if
      ! check fractions sum to one
      if(abs(sum(mLayerTranspireFrac) - 1._rkind) > verySmall)then
        message=trim(message)//'fraction transpiration in soil layers does not sum to one'
        err=20; return
      endif

      ! compute transpiration loss from each soil layer (kg m-2 s-1 --> m s-1)
      mLayerTranspire(:) = mLayerTranspireFrac(:)*scalarCanopyTranspiration/iden_water
      ! derivatives in transpiration w.r.t. canopy state variables
      mLayerdTrans_dCanWat(:)  = mLayerTranspireFrac(:)*dCanopyTrans_dCanWat /iden_water
      mLayerdTrans_dTCanair(:) = mLayerTranspireFrac(:)*dCanopyTrans_dTCanair/iden_water
      mLayerdTrans_dTCanopy(:) = mLayerTranspireFrac(:)*dCanopyTrans_dTCanopy/iden_water
      mLayerdTrans_dTGround(:) = mLayerTranspireFrac(:)*dCanopyTrans_dTGround/iden_water

      ! special case of prescribed head -- no transpiration
      if(ixBcUpperSoilHydrology==prescribedHead) then
        mLayerTranspire(:)      = 0._rkind
        ! derivatives in transpiration w.r.t. canopy state variables
        mLayerdTrans_dCanWat(:) = 0._rkind
        mLayerdTrans_dTCanair(:)= 0._rkind
        mLayerdTrans_dTCanopy(:)= 0._rkind
        mLayerdTrans_dTGround(:)= 0._rkind
      endif

    endif  ! if need to compute transpiration

    ! -------------------------------------------------------------------------------------------------------------------------------------------------
    ! compute diagnostic variables at the nodes throughout the soil profile
    ! -------------------------------------------------------------------------------------------------------------------------------------------------
    do iSoil=ixTop,min(ixBot+1,nSoil) ! (loop through soil layers)
      call diagv_node(&
                      ! input: model control
                      deriv_desired,                   & ! intent(in): flag indicating if derivatives are desired
                      ixRichards,                      & ! intent(in): index defining the option for Richards' equation (moisture or mixdform)
                      ! input: state variables
                      mLayerTempTrial(iSoil),          & ! intent(in): temperature (K)
                      mLayerMatricHeadLiqTrial(iSoil), & ! intent(in): liquid  matric head in each layer (m)
                      mLayerVolFracLiqTrial(iSoil),    & ! intent(in): volumetric liquid water content in each soil layer (-)
                      mLayerVolFracIceTrial(iSoil),    & ! intent(in): volumetric ice content in each soil layer (-)
                      ! input: pre-computed deriavatives
                      mLayerdTheta_dTk(iSoil),         & ! intent(in): derivative in volumetric liquid water content w.r.t. temperature (K-1)
                      dPsiLiq_dTemp(iSoil),            & ! intent(in): derivative in liquid water matric potential w.r.t. temperature (m K-1)
                      ! input: soil parameters
                      vGn_alpha(iSoil),                & ! intent(in): van Genutchen "alpha" parameter (m-1)
                      vGn_n(iSoil),                    & ! intent(in): van Genutchen "n" parameter (-)
                      vGn_m(iSoil),                    & ! intent(in): van Genutchen "m" parameter (-)
                      mpExp,                           & ! intent(in): empirical exponent in macropore flow equation (-)
                      theta_sat(iSoil),                & ! intent(in): soil porosity (-)
                      theta_res(iSoil),                & ! intent(in): soil residual volumetric water content (-)
                      theta_mp,                        & ! intent(in): volumetric liquid water content when macropore flow begins (-)
                      f_impede,                        & ! intent(in): ice impedence factor (-)
                      ! input: saturated hydraulic conductivity
                      mLayerSatHydCond(iSoil),         & ! intent(in): saturated hydraulic conductivity at the mid-point of each layer (m s-1)
                      mLayerSatHydCondMP(iSoil),       & ! intent(in): saturated hydraulic conductivity of macropores at the mid-point of each layer (m s-1)
                      ! output: derivative in the soil water characteristic
                      mLayerdPsi_dTheta(iSoil),        & ! intent(out): derivative in the soil water characteristic
                      mLayerdTheta_dPsi(iSoil),        & ! intent(out): derivative in the soil water characteristic
                      ! output: transmittance
                      mLayerHydCond(iSoil),            & ! intent(out): hydraulic conductivity at layer mid-points (m s-1)
                      mLayerDiffuse(iSoil),            & ! intent(out): diffusivity at layer mid-points (m2 s-1)
                      iceImpedeFac(iSoil),             & ! intent(out): ice impedence factor in each layer (-)
                      ! output: transmittance derivatives
                      dHydCond_dVolLiq(iSoil),         & ! intent(out): derivative in hydraulic conductivity w.r.t volumetric liquid water content (m s-1)
                      dDiffuse_dVolLiq(iSoil),         & ! intent(out): derivative in hydraulic diffusivity w.r.t volumetric liquid water content (m2 s-1)
                      dHydCond_dMatric(iSoil),         & ! intent(out): derivative in hydraulic conductivity w.r.t matric head (m s-1)
                      dHydCond_dTemp(iSoil),           & ! intent(out): derivative in hydraulic conductivity w.r.t temperature (m s-1 K-1)
                      ! output: error control
                      err,cmessage)                      ! intent(out): error control
      if(err/=0)then; message=trim(message)//trim(cmessage); return; end if
    end do  ! (looping through soil layers)

    ! -------------------------------------------------------------------------------------------------------------------------------------------------
    ! compute infiltraton at the surface and its derivative w.r.t. mass in the upper soil layer
    ! -------------------------------------------------------------------------------------------------------------------------------------------------
  
    ! set derivative w.r.t. state above to zero (does not exist)
    dq_dHydStateAbove(0) = 0._rkind
    dq_dNrgStateAbove(0) = 0._rkind

    ! compute surface flux and its derivative...
    call surfaceFlx(&
                    ! input: model control
                    firstSplitOper,                     & ! intent(in): flag indicating if desire to compute infiltration
                    deriv_desired,                      & ! intent(in): flag indicating if derivatives are desired
                    ixRichards,                         & ! intent(in): index defining the form of Richards' equation (moisture or mixdform)
                    ixBcUpperSoilHydrology,             & ! intent(in): index defining the type of boundary conditions (neumann or diriclet)
                    nRoots,                             & ! intent(in): number of layers that contain roots
                    ixIce,                              & ! intent(in): index of lowest ice layer
                    nSoil,                              & ! intent(in): number of soil layers
                    ! input: state variables
                    mLayerTempTrial,                    & ! intent(in): temperature (K)
                    mLayerMatricHeadLiqTrial(1),        & ! intent(in): liquid matric head in the upper-most soil layer (m)
                    mLayerMatricHeadTrial,              & ! intent(in): matric head in each soil layer (m)
                    mLayerVolFracLiqTrial(1),           & ! intent(in): volumetric liquid water content the upper-most soil layer (-)
                    mLayerVolFracLiqTrial,              & ! intent(in): volumetric liquid water content in each soil layer (-)
                    mLayerVolFracIceTrial,              & ! intent(in): volumetric ice content in each soil layer (-)
                    ! input: pre-computed deriavatives
                    mLayerdTheta_dTk,                   & ! intent(in): derivative in volumetric liquid water content w.r.t. temperature (K-1)
                    mLayerdTheta_dPsi,                  & ! intent(in): derivative in the soil water characteristic w.r.t. psi (m-1)
                    mLayerdPsi_dTheta,                  & ! intent(in): derivative in the soil water characteristic w.r.t. theta (m)
                    above_soilLiqFluxDeriv,             & ! intent(in): derivative in layer above soil (canopy or snow) liquid flux w.r.t. liquid water
                    above_soildLiq_dTk,                 & ! intent(in): derivative of layer above soil (canopy or snow) liquid flux w.r.t. temperature
                    above_soilFracLiq,                  & ! intent(in): fraction of liquid water layer above soil (canopy or snow) (-)
                    ! input: depth of upper-most soil layer (m)
                    mLayerDepth,                        & ! intent(in): depth of each soil layer (m)
                    iLayerHeight,                       & ! intent(in): height at the interface of each layer (m)
                    ! input: boundary conditions
                    upperBoundHead,                     & ! intent(in): upper boundary condition (m)
                    upperBoundTheta,                    & ! intent(in): upper boundary condition (-)
                    ! input: flux at the upper boundary
                    scalarRainPlusMelt,                 & ! intent(in): rain plus melt (m s-1)
                    ! input: transmittance
                    iLayerSatHydCond(0),                & ! intent(in): saturated hydraulic conductivity at the surface (m s-1)
                    dHydCond_dTemp(1),                  & ! intent(in): derivative in hydraulic conductivity w.r.t temperature (m s-1 K-1)
                    iceImpedeFac(1),                    & ! intent(in): ice impedence factor in the upper-most soil layer (-)
                    ! input: soil parameters
                    vGn_alpha(1),                       & ! intent(in): van Genutchen "alpha" parameter (m-1)
                    vGn_n(1),                           & ! intent(in): van Genutchen "n" parameter (-)
                    vGn_m(1),                           & ! intent(in): van Genutchen "m" parameter (-)
                    theta_sat(1),                       & ! intent(in): soil porosity (-)
                    theta_res(1),                       & ! intent(in): soil residual volumetric water content (-)
                    qSurfScale,                         & ! intent(in): scaling factor in the surface runoff parameterization (-)
                    zScale_TOPMODEL,                    & ! intent(in): scaling factor used to describe decrease in hydraulic conductivity with depth (m)
                    rootingDepth,                       & ! intent(in): rooting depth (m)
                    wettingFrontSuction,                & ! intent(in): Green-Ampt wetting front suction (m)
                    soilIceScale,                       & ! intent(in): soil ice scaling factor in Gamma distribution used to define frozen area (m)
                    soilIceCV,                          & ! intent(in): soil ice CV in Gamma distribution used to define frozen area (-)
                    ! input-output: hydraulic conductivity and diffusivity at the surface
                    iLayerHydCond(0),                   & ! intent(inout): hydraulic conductivity at the surface (m s-1)
                    iLayerDiffuse(0),                   & ! intent(inout): hydraulic diffusivity at the surface (m2 s-1)
                    ! input-output: fluxes at layer interfaces and surface runoff
                    xMaxInfilRate,                      & ! intent(inout): maximum infiltration rate (m s-1)
                    scalarInfilArea,                    & ! intent(inout): fraction of unfrozen area where water can infiltrate (-)
                    scalarFrozenArea,                   & ! intent(inout): fraction of area that is considered impermeable due to soil ice (-)
                    scalarSurfaceRunoff,                & ! intent(out): surface runoff (m s-1)
                    scalarSurfaceInfiltration,          & ! intent(out): surface infiltration (m s-1)
                    ! input-output: deriavtives in surface infiltration w.r.t. volumetric liquid water (m s-1) and matric head (s-1) in the upper-most soil layer
                    dq_dHydStateLayerSurfVec,           & ! intent(inout): derivative in surface infiltration w.r.t. hydrology state in above soil snow or canopy and every soil layer  (m s-1 or s-1)
                    dq_dNrgStateLayerSurfVec,           & ! intent(inout): derivative in surface infiltration w.r.t. energy state in above soil snow or canopy and every soil layer (m s-1 K-1)
                    ! output: error control
                    err,cmessage)                         ! intent(out): error control
    if(err/=0)then; message=trim(message)//trim(cmessage); return; end if

    ! include base soil evaporation as the upper boundary flux
    iLayerLiqFluxSoil(0) = scalarGroundEvaporation/iden_water + scalarSurfaceInfiltration

    dq_dHydStateBelow(0) = 0._rkind ! contribution will be in dq_dHydStateLayerSurfVec(1)
    dq_dNrgStateBelow(0) = 0._rkind ! contribution will be in dq_dNrgStateLayerSurfVec(1)

    ! -------------------------------------------------------------------------------------------------------------------------------------------------
    ! * compute fluxes and derivatives at layer interfaces...
    ! -------------------------------------------------------------------------------------------------------------------------------------------------

    ! computing flux at the bottom of the layer
    do iLayer=ixTop,min(ixBot,nSoil-1)
      call iLayerFlux(&
                      ! input: model control
                      deriv_desired,                      & ! intent(in): flag indicating if derivatives are desired
                      ixRichards,                         & ! intent(in): index defining the form of Richards' equation (moisture or mixdform)
                      ! input: state variables (adjacent layers)
                      mLayerMatricHeadLiqTrial(iLayer:iLayer+1), & ! intent(in): liquid matric head at the soil nodes (m)
                      mLayerVolFracLiqTrial(iLayer:iLayer+1),    & ! intent(in): volumetric liquid water content at the soil nodes (-)
                      ! input: model coordinate variables (adjacent layers)
                      mLayerHeight(iLayer:iLayer+1),      & ! intent(in): height of the soil nodes (m)
                      ! input: temperature derivatives
                      dPsiLiq_dTemp(iLayer:iLayer+1),     & ! intent(in): derivative in liquid water matric potential w.r.t. temperature (m K-1)
                      dHydCond_dTemp(iLayer:iLayer+1),    & ! intent(in): derivative in hydraulic conductivity w.r.t temperature (m s-1 K-1)
                      ! input: transmittance (adjacent layers)
                      mLayerHydCond(iLayer:iLayer+1),     & ! intent(in): hydraulic conductivity at the soil nodes (m s-1)
                      mLayerDiffuse(iLayer:iLayer+1),     & ! intent(in): hydraulic diffusivity at the soil nodes (m2 s-1)
                      ! input: transmittance derivatives (adjacent layers)
                      dHydCond_dVolLiq(iLayer:iLayer+1),  & ! intent(in): change in hydraulic conductivity w.r.t. change in volumetric liquid water content (m s-1)
                      dDiffuse_dVolLiq(iLayer:iLayer+1),  & ! intent(in): change in hydraulic diffusivity w.r.t. change in volumetric liquid water content (m2 s-1)
                      dHydCond_dMatric(iLayer:iLayer+1),  & ! intent(in): change in hydraulic conductivity w.r.t. change in matric head (s-1)
                      ! output: tranmsmittance at the layer interface (scalars)
                      iLayerHydCond(iLayer),              & ! intent(out): hydraulic conductivity at the interface between layers (m s-1)
                      iLayerDiffuse(iLayer),              & ! intent(out): hydraulic diffusivity at the interface between layers (m2 s-1)
                      ! output: vertical flux at the layer interface (scalars)
                      iLayerLiqFluxSoil(iLayer),          & ! intent(out): vertical flux of liquid water at the layer interface (m s-1)
                      ! output: derivatives in fluxes w.r.t. state variables -- matric head or volumetric lquid water -- in the layer above and layer below (m s-1 or s-1)
                      dq_dHydStateAbove(iLayer),          & ! intent(out): derivatives in the flux w.r.t. matric head or volumetric lquid water in the layer above (m s-1 or s-1)
                      dq_dHydStateBelow(iLayer),          & ! intent(out): derivatives in the flux w.r.t. matric head or volumetric lquid water in the layer below (m s-1 or s-1)
                      ! output: derivatives in fluxes w.r.t. energy state variables -- now just temperature -- in the layer above and layer below (m s-1 K-1)
                      dq_dNrgStateAbove(iLayer),          & ! intent(out): derivatives in the flux w.r.t. temperature in the layer above (m s-1 K-1)
                      dq_dNrgStateBelow(iLayer),          & ! intent(out): derivatives in the flux w.r.t. temperature in the layer below (m s-1 K-1)
                      ! output: error control
                      err,cmessage)                         ! intent(out): error control
      if(err/=0)then; message=trim(message)//trim(cmessage); return; end if
    end do  ! (looping through soil layers)

    ! -------------------------------------------------------------------------------------------------------------------------------------------------
    ! * compute drainage flux from the bottom of the soil profile, and its derivative
    ! -------------------------------------------------------------------------------------------------------------------------------------------------

    ! define the need to compute drainage
    if( .not. (scalarSolution .and. ixTop<nSoil) )then
      ! compute drainage flux and its derivative...
      call qDrainFlux(&
                      ! input: model control
                      deriv_desired,                      & ! intent(in): flag indicating if derivatives are desired
                      ixRichards,                      & ! intent(in): index defining the form of Richards' equation (moisture or mixdform)
                      ixBcLowerSoilHydrology,          & ! intent(in): index defining the type of boundary conditions
                      ! input: state variables
                      mLayerMatricHeadLiqTrial(nSoil), & ! intent(in): liquid matric head in the lowest unsaturated node (m)
                      mLayerVolFracLiqTrial(nSoil),    & ! intent(in): volumetric liquid water content the lowest unsaturated node (-)
                      ! input: model coordinate variables
                      mLayerDepth(nSoil),              & ! intent(in): depth of the lowest unsaturated soil layer (m)
                      mLayerHeight(nSoil),             & ! intent(in): height of the lowest unsaturated soil node (m)
                      ! input: boundary conditions
                      lowerBoundHead,                  & ! intent(in): lower boundary condition (m)
                      lowerBoundTheta,                 & ! intent(in): lower boundary condition (-)
                      ! input: derivative in the soil water characteristic
                      mLayerdPsi_dTheta(nSoil),        & ! intent(in): derivative in the soil water characteristic
                      ! input: transmittance
                      iLayerSatHydCond(0),             & ! intent(in): saturated hydraulic conductivity at the surface (m s-1)
                      iLayerSatHydCond(nSoil),         & ! intent(in): saturated hydraulic conductivity at the bottom of the unsaturated zone (m s-1)
                      mLayerHydCond(nSoil),            & ! intent(in): hydraulic conductivity at the node itself (m s-1)
                      iceImpedeFac(nSoil),             & ! intent(in): ice impedence factor in the lower-most soil layer (-)
                      ! input: transmittance derivatives
                      dHydCond_dVolLiq(nSoil),         & ! intent(in): derivative in hydraulic conductivity w.r.t. volumetric liquid water content (m s-1)
                      dHydCond_dMatric(nSoil),         & ! intent(in): derivative in hydraulic conductivity w.r.t. matric head (s-1)
                      dHydCond_dTemp(nSoil),           & ! intent(in): derivative in hydraulic conductivity w.r.t temperature (m s-1 K-1)
                      ! input: soil parameters
                      vGn_alpha(nSoil),                & ! intent(in): van Genutchen "alpha" parameter (m-1)
                      vGn_n(nSoil),                    & ! intent(in): van Genutchen "n" parameter (-)
                      vGn_m(nSoil),                    & ! intent(in): van Genutchen "m" parameter (-)
                      theta_sat(nSoil),                & ! intent(in): soil porosity (-)
                      theta_res(nSoil),                & ! intent(in): soil residual volumetric water content (-)
                      kAnisotropic,                    & ! intent(in): anisotropy factor for lateral hydraulic conductivity (-)
                      zScale_TOPMODEL,                 & ! intent(in): TOPMODEL scaling factor (m)
                      ! output: hydraulic conductivity and diffusivity at the surface
                      iLayerHydCond(nSoil),            & ! intent(out): hydraulic conductivity at the bottom of the unsatuarted zone (m s-1)
                      iLayerDiffuse(nSoil),            & ! intent(out): hydraulic diffusivity at the bottom of the unsatuarted zone (m2 s-1)
                      ! output: drainage flux
                      iLayerLiqFluxSoil(nSoil),        & ! intent(out): drainage flux (m s-1)
                      ! output: derivatives in drainage flux
                      dq_dHydStateAbove(nSoil),        & ! intent(out): change in drainage flux w.r.t. change in hydrology state in lowest unsaturated node (m s-1 or s-1)
                      dq_dNrgStateAbove(nSoil),        & ! intent(out): change in drainage flux w.r.t. change in energy state in lowest unsaturated node (m s-1 or s-1)
                      ! output: error control
                      err,cmessage)                ! intent(out): error control
      if(err/=0)then; message=trim(message)//trim(cmessage); return; end if

      ! no dependence on the aquifer for drainage
      dq_dHydStateBelow(nSoil) = 0._rkind  ! keep this here in case we want to couple some day....
      dq_dNrgStateBelow(nSoil) = 0._rkind  ! keep this here in case we want to couple some day....

    endif  ! if computing drainage

  end associate

end subroutine soilLiqFlx

! ***************************************************************************************************************
! private subroutine diagv_node: compute transmittance and derivatives for model nodes
! ***************************************************************************************************************
subroutine diagv_node(&
                      ! input: model control
                      deriv_desired,         & ! intent(in): flag indicating if derivatives are desired
                      ixRichards,            & ! intent(in): index defining the option for Richards' equation (moisture or mixdform)
                      ! input: state variables
                      scalarTempTrial,       & ! intent(in): temperature (K)
                      scalarMatricHeadLiqTrial, & ! intent(in): liquid matric head in a given layer (m)
                      scalarVolFracLiqTrial, & ! intent(in): volumetric liquid water content in a given soil layer (-)
                      scalarVolFracIceTrial, & ! intent(in): volumetric ice content in a given soil layer (-)
                      ! input: pre-computed deriavatives
                      dTheta_dTk,            & ! intent(in): derivative in volumetric liquid water content w.r.t. temperature (K-1)
                      dPsiLiq_dTemp,         & ! intent(in): derivative in liquid water matric potential w.r.t. temperature (m K-1)
                      ! input: soil parameters
                      vGn_alpha,             & ! intent(in): van Genutchen "alpha" parameter (m-1)
                      vGn_n,                 & ! intent(in): van Genutchen "n" parameter (-)
                      vGn_m,                 & ! intent(in): van Genutchen "m" parameter (-)
                      mpExp,                 & ! intent(in): empirical exponent in macropore flow equation (-)
                      theta_sat,             & ! intent(in): soil porosity (-)
                      theta_res,             & ! intent(in): soil residual volumetric water content (-)
                      theta_mp,              & ! intent(in): volumetric liquid water content when macropore flow begins (-)
                      f_impede,              & ! intent(in): ice impedence factor (-)
                      ! input: saturated hydraulic conductivity
                      scalarSatHydCond,      & ! intent(in): saturated hydraulic conductivity at the mid-point of a given layer (m s-1)
                      scalarSatHydCondMP,    & ! intent(in): saturated hydraulic conductivity of macropores at the mid-point of a given layer (m s-1)
                      ! output: derivative in the soil water characteristic
                      scalardPsi_dTheta,     & ! derivative in the soil water characteristic
                      scalardTheta_dPsi,     & ! derivative in the soil water characteristic
                      ! output: transmittance
                      scalarHydCond,         & ! intent(out): hydraulic conductivity at layer mid-points (m s-1)
                      scalarDiffuse,         & ! intent(out): diffusivity at layer mid-points (m2 s-1)
                      iceImpedeFac,          & ! intent(out): ice impedence factor in each layer (-)
                      ! output: transmittance derivatives
                      dHydCond_dVolLiq,      & ! intent(out): derivative in hydraulic conductivity w.r.t volumetric liquid water content (m s-1)
                      dDiffuse_dVolLiq,      & ! intent(out): derivative in hydraulic diffusivity w.r.t volumetric liquid water content (m2 s-1)
                      dHydCond_dMatric,      & ! intent(out): derivative in hydraulic conductivity w.r.t matric head (m s-1)
                      dHydCond_dTemp,        & ! intent(out): derivative in hydraulic conductivity w.r.t temperature (m s-1 K-1)
                      ! output: error control
                      err,message)             ! intent(out): error control
  USE soil_utils_module,only:iceImpede           ! compute the ice impedence factor
  USE soil_utils_module,only:volFracLiq          ! compute volumetric fraction of liquid water as a function of matric head
  USE soil_utils_module,only:matricHead          ! compute matric head (m)
  USE soil_utils_module,only:hydCond_psi         ! compute hydraulic conductivity as a function of matric head
  USE soil_utils_module,only:hydCond_liq         ! compute hydraulic conductivity as a function of volumetric liquid water content
  USE soil_utils_module,only:hydCondMP_liq       ! compute hydraulic conductivity of macropores as a function of volumetric liquid water content
  USE soil_utils_module,only:dTheta_dPsi         ! compute derivative of the soil moisture characteristic w.r.t. psi (m-1)
  USE soil_utils_module,only:dPsi_dTheta         ! compute derivative of the soil moisture characteristic w.r.t. theta (m)
  USE soil_utils_module,only:dPsi_dTheta2        ! compute derivative in dPsi_dTheta (m)
  USE soil_utils_module,only:dHydCond_dLiq       ! compute derivative in hydraulic conductivity w.r.t. volumetric liquid water content (m s-1)
  USE soil_utils_module,only:dHydCond_dPsi       ! compute derivative in hydraulic conductivity w.r.t. matric head (s-1)
  USE soil_utils_module,only:dIceImpede_dTemp    ! compute the derivative in the ice impedance factor w.r.t. temperature (K-1)
  ! compute hydraulic transmittance and derivatives for all layers
  implicit none
  ! input: model control
  logical(lgt),intent(in)       :: deriv_desired             ! flag indicating if derivatives are desired
  integer(i4b),intent(in)       :: ixRichards                ! index defining the option for Richards' equation (moisture or mixdform)
  ! input: state and diagnostic variables
  real(rkind),intent(in)           :: scalarTempTrial           ! temperature in each layer (K)
  real(rkind),intent(in)           :: scalarMatricHeadLiqTrial  ! liquid matric head in each layer (m)
  real(rkind),intent(in)           :: scalarVolFracLiqTrial     ! volumetric fraction of liquid water in a given layer (-)
  real(rkind),intent(in)           :: scalarVolFracIceTrial     ! volumetric fraction of ice in a given layer (-)
  ! input: pre-computed deriavatives
  real(rkind),intent(in)           :: dTheta_dTk                ! derivative in volumetric liquid water content w.r.t. temperature (K-1)
  real(rkind),intent(in)           :: dPsiLiq_dTemp             ! derivative in liquid water matric potential w.r.t. temperature (m K-1)
  ! input: soil parameters
  real(rkind),intent(in)           :: vGn_alpha                 ! van Genutchen "alpha" parameter (m-1)
  real(rkind),intent(in)           :: vGn_n                     ! van Genutchen "n" parameter (-)
  real(rkind),intent(in)           :: vGn_m                     ! van Genutchen "m" parameter (-)
  real(rkind),intent(in)           :: mpExp                     ! empirical exponent in macropore flow equation (-)
  real(rkind),intent(in)           :: theta_sat                 ! soil porosity (-)
  real(rkind),intent(in)           :: theta_res                 ! soil residual volumetric water content (-)
  real(rkind),intent(in)           :: theta_mp                  ! volumetric liquid water content when macropore flow begins (-)
  real(rkind),intent(in)           :: f_impede                  ! ice impedence factor (-)
  ! input: saturated hydraulic conductivity
  real(rkind),intent(in)           :: scalarSatHydCond          ! saturated hydraulic conductivity at the mid-point of a given layer (m s-1)
  real(rkind),intent(in)           :: scalarSatHydCondMP        ! saturated hydraulic conductivity of macropores at the mid-point of a given layer (m s-1)
  ! output: derivative in the soil water characteristic
  real(rkind),intent(out)          :: scalardPsi_dTheta         ! derivative in the soil water characteristic
  real(rkind),intent(out)          :: scalardTheta_dPsi         ! derivative in the soil water characteristic
  ! output: transmittance
  real(rkind),intent(out)          :: scalarHydCond             ! hydraulic conductivity at layer mid-points (m s-1)
  real(rkind),intent(out)          :: scalarDiffuse             ! diffusivity at layer mid-points (m2 s-1)
  real(rkind),intent(out)          :: iceImpedeFac              ! ice impedence factor in each layer (-)
  ! output: transmittance derivatives
  real(rkind),intent(out)          :: dHydCond_dVolLiq          ! derivative in hydraulic conductivity w.r.t volumetric liquid water content (m s-1)
  real(rkind),intent(out)          :: dDiffuse_dVolLiq          ! derivative in hydraulic diffusivity w.r.t volumetric liquid water content (m2 s-1)
  real(rkind),intent(out)          :: dHydCond_dMatric          ! derivative in hydraulic conductivity w.r.t matric head (s-1)
  real(rkind),intent(out)          :: dHydCond_dTemp            ! derivative in hydraulic conductivity w.r.t temperature (m s-1 K-1)
  ! output: error control
  integer(i4b),intent(out)      :: err                       ! error code
  character(*),intent(out)      :: message                   ! error message
  ! local variables
  real(rkind)                      :: localVolFracLiq           ! local volumetric fraction of liquid water
  real(rkind)                      :: scalarHydCondMP           ! hydraulic conductivity of macropores at layer mid-points (m s-1)
  real(rkind)                      :: dIceImpede_dT             ! derivative in ice impedance factor w.r.t. temperature (K-1)
  real(rkind)                      :: dHydCondMacro_dVolLiq     ! derivative in hydraulic conductivity of macropores w.r.t volumetric liquid water content (m s-1)
  real(rkind)                      :: dHydCondMacro_dMatric     ! derivative in hydraulic conductivity of macropores w.r.t matric head (s-1)
  real(rkind)                      :: dHydCondMicro_dMatric     ! derivative in hydraulic conductivity of micropores w.r.t matric head (s-1)
  real(rkind)                      :: dHydCondMicro_dTemp       ! derivative in hydraulic conductivity of micropores w.r.t temperature (m s-1 K-1)
  real(rkind)                      :: dPsi_dTheta2a             ! derivative in dPsi_dTheta (analytical)
  real(rkind)                      :: dIceImpede_dLiq           ! derivative in ice impedence factor w.r.t. volumetric liquid water content (-)
  real(rkind)                      :: hydCond_noIce             ! hydraulic conductivity in the absence of ice (m s-1)
  real(rkind)                      :: dK_dLiq__noIce            ! derivative in hydraulic conductivity w.r.t volumetric liquid water content, in the absence of ice (m s-1)
  real(rkind)                      :: dK_dPsi__noIce            ! derivative in hydraulic conductivity w.r.t matric head, in the absence of ice (s-1)
  real(rkind)                      :: relSatMP                  ! relative saturation of macropores (-)
  ! initialize error control
  err=0; message="diagv_node/"

  ! *****
  ! compute the derivative in the soil water characteristic
  select case(ixRichards)
    case(moisture)
      scalardPsi_dTheta = dPsi_dTheta(scalarVolFracLiqTrial,vGn_alpha,theta_res,theta_sat,vGn_n,vGn_m)
      scalardTheta_dPsi = realMissing  ! (deliberately cause problems if this is ever used)
    case(mixdform)
      scalardTheta_dPsi = dTheta_dPsi(scalarMatricHeadLiqTrial,vGn_alpha,theta_res,theta_sat,vGn_n,vGn_m)
      scalardPsi_dTheta = dPsi_dTheta(scalarVolFracLiqTrial,vGn_alpha,theta_res,theta_sat,vGn_n,vGn_m)
    case default; err=10; message=trim(message)//"unknown form of Richards' equation"; return
  end select

  ! *****
  ! compute hydraulic conductivity and its derivative in each soil layer

  ! compute the ice impedence factor and its derivative w.r.t. volumetric liquid water content (-)
  call iceImpede(scalarVolFracIceTrial,f_impede, &  ! input
                  iceImpedeFac,dIceImpede_dLiq)      ! output

  select case(ixRichards)
    ! ***** moisture-based form of Richards' equation
    case(moisture)
      ! haven't included macropores yet
      err=20; message=trim(message)//'still need to include macropores for the moisture-based form of Richards eqn'; return
      ! compute the hydraulic conductivity (m s-1) and diffusivity (m2 s-1) for a given layer
      hydCond_noIce = hydCond_liq(scalarVolFracLiqTrial,scalarSatHydCond,theta_res,theta_sat,vGn_m)
      scalarHydCond = hydCond_noIce*iceImpedeFac
      scalarDiffuse = scalardPsi_dTheta * scalarHydCond
      ! compute derivative in hydraulic conductivity (m s-1) and hydraulic diffusivity (m2 s-1)
      if(deriv_desired)then
        if(scalarVolFracIceTrial > epsilon(iceImpedeFac))then
          dK_dLiq__noIce   = dHydCond_dLiq(scalarVolFracLiqTrial,scalarSatHydCond,theta_res,theta_sat,vGn_m,.true.)  ! [.true. = analytical]
          dHydCond_dVolLiq = hydCond_noIce*dIceImpede_dLiq + dK_dLiq__noIce*iceImpedeFac
        else
          dHydCond_dVolLiq = dHydCond_dLiq(scalarVolFracLiqTrial,scalarSatHydCond,theta_res,theta_sat,vGn_m,.true.)
        end if
          dPsi_dTheta2a    = dPsi_dTheta2(scalarVolFracLiqTrial,vGn_alpha,theta_res,theta_sat,vGn_n,vGn_m,.true.)   ! [.true. = analytical] compute derivative in dPsi_dTheta (m)
          dDiffuse_dVolLiq = dHydCond_dVolLiq*scalardPsi_dTheta + scalarHydCond*dPsi_dTheta2a
          dHydCond_dMatric = realMissing ! not used, so cause problems
      end if

    ! ***** mixed form of Richards' equation -- just compute hydraulic condictivity
    case(mixdform)
      ! compute the hydraulic conductivity (m s-1) and diffusivity (m2 s-1) for a given layer
      hydCond_noIce = hydCond_psi(scalarMatricHeadLiqTrial,scalarSatHydCond,vGn_alpha,vGn_n,vGn_m)
      scalarDiffuse = realMissing ! not used, so cause problems
      ! compute the hydraulic conductivity of macropores (m s-1)
      localVolFracLiq = volFracLiq(scalarMatricHeadLiqTrial,vGn_alpha,theta_res,theta_sat,vGn_n,vGn_m)
      scalarHydCondMP = hydCondMP_liq(localVolFracLiq,theta_sat,theta_mp,mpExp,scalarSatHydCondMP,scalarSatHydCond)
      scalarHydCond   = hydCond_noIce*iceImpedeFac + scalarHydCondMP

      ! compute derivative in hydraulic conductivity (m s-1)
      if(deriv_desired)then 
        ! (compute derivative for macropores)
        if(localVolFracLiq > theta_mp)then
          relSatMP              = (localVolFracLiq - theta_mp)/(theta_sat - theta_mp)
          dHydCondMacro_dVolLiq = ((scalarSatHydCondMP - scalarSatHydCond)/(theta_sat - theta_mp))*mpExp*(relSatMP**(mpExp - 1._rkind))
          dHydCondMacro_dMatric = scalardTheta_dPsi*dHydCondMacro_dVolLiq
        else
          dHydCondMacro_dVolLiq = 0._rkind
          dHydCondMacro_dMatric = 0._rkind
        end if
        ! (compute derivatives for micropores)
        if(scalarVolFracIceTrial > verySmall)then
          dK_dPsi__noIce        = dHydCond_dPsi(scalarMatricHeadLiqTrial,scalarSatHydCond,vGn_alpha,vGn_n,vGn_m,.true.)  ! analytical
          dHydCondMicro_dTemp   = dPsiLiq_dTemp*dK_dPsi__noIce  ! m s-1 K-1
          dHydCondMicro_dMatric = hydCond_noIce*dIceImpede_dLiq*scalardTheta_dPsi + dK_dPsi__noIce*iceImpedeFac
        else
          dHydCondMicro_dTemp   = 0._rkind
          dHydCondMicro_dMatric = dHydCond_dPsi(scalarMatricHeadLiqTrial,scalarSatHydCond,vGn_alpha,vGn_n,vGn_m,.true.)
        end if
        ! (combine derivatives)
        dHydCond_dMatric = dHydCondMicro_dMatric + dHydCondMacro_dMatric

        ! (compute analytical derivative for change in ice impedance factor w.r.t. temperature)
        call dIceImpede_dTemp(scalarVolFracIceTrial, & ! intent(in): trial value of volumetric ice content (-)
                              dTheta_dTk,            & ! intent(in): derivative in volumetric liquid water content w.r.t. temperature (K-1)
                              f_impede,              & ! intent(in): ice impedance parameter (-)
                              dIceImpede_dT          ) ! intent(out): derivative in ice impedance factor w.r.t. temperature (K-1)
        ! (compute derivative in hydraulic conductivity w.r.t. temperature)
        dHydCond_dTemp = hydCond_noIce*dIceImpede_dT + dHydCondMicro_dTemp*iceImpedeFac
        ! (set values that are not used to missing)
        dHydCond_dVolLiq = realMissing ! not used, so cause problems
        dDiffuse_dVolLiq = realMissing ! not used, so cause problems
      end if
    case default; err=10; message=trim(message)//"unknown form of Richards' equation"; return
  end select ! select form of Richards' equation

  ! if derivatives are not desired, then set values to missing
  if(.not.deriv_desired)then
    dHydCond_dVolLiq   = realMissing ! not used, so cause problems
    dDiffuse_dVolLiq   = realMissing ! not used, so cause problems
    dHydCond_dMatric   = realMissing ! not used, so cause problems
  end if

end subroutine diagv_node


! ***************************************************************************************************************
! private subroutine surfaceFlx: compute the surface flux and its derivative
! ***************************************************************************************************************
subroutine surfaceFlx(&
                      ! input: model control
                      firstSplitOper,            & ! intent(in): flag indicating if desire to compute infiltration
                      deriv_desired,             & ! intent(in): flag indicating if derivatives are desired
                      ixRichards,                & ! intent(in): index defining the form of Richards' equation (moisture or mixdform)
                      bc_upper,                  & ! intent(in): index defining the type of boundary conditions (neumann or diriclet)
                      nRoots,                    & ! intent(in): number of layers that contain roots
                      ixIce,                     & ! intent(in): index of lowest ice layer
                      nSoil,                     & ! intent(in): number of soil layers
                      ! input: state variables
                      mLayerTemp,                & ! intent(in): temperature (K)
                      scalarMatricHeadLiq,       & ! intent(in): liquid matric head in the upper-most soil layer (m)
                      mLayerMatricHead,          & ! intent(in): matric head in each soil layer (m)
                      scalarVolFracLiq,          & ! intent(in): volumetric liquid water content in the upper-most soil layer (-)
                      mLayerVolFracLiq,          & ! intent(in): volumetric liquid water content in each soil layer (-)
                      mLayerVolFracIce,          & ! intent(in): volumetric ice content in each soil layer (-)
                      ! input: pre-computed derivatives
                      dTheta_dTk,                & ! intent(in): derivative in volumetric liquid water content w.r.t. temperature (K-1)
                      dTheta_dPsi,               & ! intent(in): derivative in the soil water characteristic w.r.t. psi (m-1)
                      mLayerdPsi_dTheta,         & ! intent(in): derivative in the soil water characteristic w.r.t. theta (m)
                      above_soilLiqFluxDeriv,    & ! intent(in): derivative in layer above soil (canopy or snow) liquid flux w.r.t. liquid water
                      above_soildLiq_dTk,        & ! intent(in): derivative of layer above soil (canopy or snow) liquid flux w.r.t. temperature
                      above_soilFracLiq,         & ! intent(in): fraction of liquid water layer above soil (canopy or snow) (-)
                      ! input: depth of upper-most soil layer (m)
                      mLayerDepth,               & ! intent(in): depth of each soil layer (m)
                      iLayerHeight,              & ! intent(in): height at the interface of each layer (m)
                      ! input: boundary conditions
                      upperBoundHead,            & ! intent(in): upper boundary condition (m)
                      upperBoundTheta,           & ! intent(in): upper boundary condition (-)
                      ! input: flux at the upper boundary
                      scalarRainPlusMelt,        & ! intent(in): rain plus melt (m s-1)
                      ! input: transmittance
                      surfaceSatHydCond,         & ! intent(in): saturated hydraulic conductivity at the surface (m s-1)
                      dHydCond_dTemp,            & ! intent(in): derivative in hydraulic conductivity w.r.t temperature (m s-1 K-1)
                      iceImpedeFac,              & ! intent(in): ice impedence factor in the upper-most soil layer (-)
                      ! input: soil parameters
                      vGn_alpha,                 & ! intent(in): van Genutchen "alpha" parameter (m-1)
                      vGn_n,                     & ! intent(in): van Genutchen "n" parameter (-)
                      vGn_m,                     & ! intent(in): van Genutchen "m" parameter (-)
                      theta_sat,                 & ! intent(in): soil porosity (-)
                      theta_res,                 & ! intent(in): soil residual volumetric water content (-)
                      qSurfScale,                & ! intent(in): scaling factor in the surface runoff parameterization (-)
                      zScale_TOPMODEL,           & ! intent(in): scaling factor used to describe decrease in hydraulic conductivity with depth (m)
                      rootingDepth,              & ! intent(in): rooting depth (m)
                      wettingFrontSuction,       & ! intent(in): Green-Ampt wetting front suction (m)
                      soilIceScale,              & ! intent(in): soil ice scaling factor in Gamma distribution used to define frozen area (m)
                      soilIceCV,                 & ! intent(in): soil ice CV in Gamma distribution used to define frozen area (-)
                      ! input-output: hydraulic conductivity and diffusivity at the surface
                      surfaceHydCond,            & ! intent(inout): hydraulic conductivity at the surface (m s-1)
                      surfaceDiffuse,            & ! intent(inout): hydraulic diffusivity at the surface (m2 s-1)
                      ! input-output: fluxes at layer interfaces and surface runoff
                      xMaxInfilRate,             & ! intent(inout): maximum infiltration rate (m s-1)
                      scalarInfilArea,           & ! intent(inout): fraction of unfrozen area where water can infiltrate (-)
                      scalarFrozenArea,          & ! intent(inout): fraction of area that is considered impermeable due to soil ice (-)
                      scalarSurfaceRunoff,       & ! intent(out): surface runoff (m s-1)
                      scalarSurfaceInfiltration, & ! intent(out): surface infiltration (m s-1)
                      ! input-output: deriavtives in surface infiltration w.r.t. volumetric liquid water (m s-1) and matric head (s-1) in the upper-most soil layer
                      dq_dHydStateVec,           & ! intent(inout): derivative in surface infiltration w.r.t. hydrology state in above soil snow or canopy and every soil layer (m s-1 or s-1)
                      dq_dNrgStateVec,           & ! intent(inout): derivative in surface infiltration w.r.t. energy state in above soil snow or canopy and every soil layer (m s-1 K-1)
                      ! output: error control
                      err,message)                 ! intent(out): error control
  USE soil_utils_module,only:volFracLiq      ! compute volumetric fraction of liquid water as a function of matric head (-)
  USE soil_utils_module,only:hydCond_psi     ! compute hydraulic conductivity as a function of matric head (m s-1)
  USE soil_utils_module,only:hydCond_liq     ! compute hydraulic conductivity as a function of volumetric liquid water content (m s-1)
  USE soil_utils_module,only:dPsi_dTheta     ! compute derivative of the soil moisture characteristic w.r.t. theta (m)
  USE soil_utils_module,only:crit_soilT      ! compute critical temperature below which ice exists
  USE soil_utils_module,only:gammp           ! compute the cumulative probabilty based on the Gamma distribution
  ! compute infiltraton at the surface and its derivative w.r.t. mass in the upper soil layer
  implicit none
  ! -----------------------------------------------------------------------------------------------------------------------------
  ! input: model control
  logical(lgt),intent(in)          :: firstSplitOper            ! flag indicating if desire to compute infiltration
  logical(lgt),intent(in)          :: deriv_desired             ! flag to indicate if derivatives are desired
  integer(i4b),intent(in)          :: bc_upper                  ! index defining the type of boundary conditions
  integer(i4b),intent(in)          :: ixRichards                ! index defining the option for Richards' equation (moisture or mixdform)
  integer(i4b),intent(in)          :: nRoots                    ! number of layers that contain roots
  integer(i4b),intent(in)          :: ixIce                     ! index of lowest ice layer
  integer(i4b),intent(in)          :: nSoil                     ! number of soil layers
  ! input: state and diagnostic variables
  real(rkind),intent(in)           :: mLayerTemp(:)             ! temperature (K)
  real(rkind),intent(in)           :: scalarMatricHeadLiq       ! liquid matric head in the upper-most soil layer (m)
  real(rkind),intent(in)           :: mLayerMatricHead(:)       ! matric head in each soil layer (m)
  real(rkind),intent(in)           :: scalarVolFracLiq          ! volumetric liquid water content in the upper-most soil layer (-)
  real(rkind),intent(in)           :: mLayerVolFracLiq(:)       ! volumetric liquid water content in each soil layer (-)
  real(rkind),intent(in)           :: mLayerVolFracIce(:)       ! volumetric ice content in each soil layer (-)
  ! input: pre-computed derivatives, note all of these would need to be recomputed if wanted a numerical derivative
  real(rkind),intent(in)           :: dTheta_dTk(:)             ! derivative in volumetric liquid water content w.r.t. temperature (K-1)
  real(rkind),intent(in)           :: dTheta_dPsi(:)            ! derivative in the soil water characteristic w.r.t. psi (m-1)
  real(rkind),intent(in)           :: mLayerdPsi_dTheta(:)      ! derivative in the soil water characteristic w.r.t. theta (m)
  real(rkind),intent(in)           :: above_soilLiqFluxDeriv    ! derivative in layer above soil (canopy or snow) liquid flux w.r.t. liquid water
  real(rkind),intent(in)           :: above_soildLiq_dTk        ! derivative of layer above soil (canopy or snow) liquid flux w.r.t. temperature
  real(rkind),intent(in)           :: above_soilFracLiq         ! fraction of liquid water layer above soil (canopy or snow) (-)
  ! input: depth of upper-most soil layer (m)
  real(rkind),intent(in)           :: mLayerDepth(:)            ! depth of upper-most soil layer (m)
  real(rkind),intent(in)           :: iLayerHeight(0:)          ! height at the interface of each layer (m)
  ! input: diriclet boundary conditions
  real(rkind),intent(in)           :: upperBoundHead            ! upper boundary condition for matric head (m)
  real(rkind),intent(in)           :: upperBoundTheta           ! upper boundary condition for volumetric liquid water content (-)
  ! input: flux at the upper boundary
  real(rkind),intent(in)           :: scalarRainPlusMelt        ! rain plus melt, used as input to the soil zone before computing surface runoff (m s-1)
  ! input: transmittance
  real(rkind),intent(in)           :: surfaceSatHydCond         ! saturated hydraulic conductivity at the surface (m s-1)
  real(rkind),intent(in)           :: dHydCond_dTemp            ! derivative in hydraulic conductivity w.r.t temperature (m s-1 K-1)
  real(rkind),intent(in)           :: iceImpedeFac              ! ice impedence factor in the upper-most soil layer (-)
  ! input: soil parameters
  real(rkind),intent(in)           :: vGn_alpha                 ! van Genutchen "alpha" parameter (m-1)
  real(rkind),intent(in)           :: vGn_n                     ! van Genutchen "n" parameter (-)
  real(rkind),intent(in)           :: vGn_m                     ! van Genutchen "m" parameter (-)
  real(rkind),intent(in)           :: theta_sat                 ! soil porosity (-)
  real(rkind),intent(in)           :: theta_res                 ! soil residual volumetric water content (-)
  real(rkind),intent(in)           :: qSurfScale                ! scaling factor in the surface runoff parameterization (-)
  real(rkind),intent(in)           :: zScale_TOPMODEL           ! scaling factor used to describe decrease in hydraulic conductivity with depth (m)
  real(rkind),intent(in)           :: rootingDepth              ! rooting depth (m)
  real(rkind),intent(in)           :: wettingFrontSuction       ! Green-Ampt wetting front suction (m)
  real(rkind),intent(in)           :: soilIceScale              ! soil ice scaling factor in Gamma distribution used to define frozen area (m)
  real(rkind),intent(in)           :: soilIceCV                 ! soil ice CV in Gamma distribution used to define frozen area (-)
  ! -----------------------------------------------------------------------------------------------------------------------------
  ! input-output: hydraulic conductivity and diffusivity at the surface
  ! NOTE: intent(inout) because infiltration may only be computed for the first iteration
  real(rkind),intent(inout)        :: surfaceHydCond            ! hydraulic conductivity (m s-1)
  real(rkind),intent(inout)        :: surfaceDiffuse            ! hydraulic diffusivity at the surface (m
  ! output: surface runoff and infiltration flux (m s-1)
  real(rkind),intent(inout)        :: xMaxInfilRate             ! maximum infiltration rate (m s-1)
  real(rkind),intent(inout)        :: scalarInfilArea           ! fraction of unfrozen area where water can infiltrate (-)
  real(rkind),intent(inout)        :: scalarFrozenArea          ! fraction of area that is considered impermeable due to soil ice (-)
  real(rkind),intent(out)          :: scalarSurfaceRunoff       ! surface runoff (m s-1)
  real(rkind),intent(out)          :: scalarSurfaceInfiltration ! surface infiltration (m s-1)
  ! output: derivatives in surface infiltration w.r.t. states in above soil snow or canopy and every soil layer
  real(rkind),intent(out)          :: dq_dHydStateVec(0:)       ! derivative in surface infiltration w.r.t. hydrology state in above soil snow or canopy and every soil layer (m s-1 or s-1)
  real(rkind),intent(out)          :: dq_dNrgStateVec(0:)       ! derivative in surface infiltration w.r.t. energy state in above soil snow or canopy and every soil layer  (m s-1 K-1)
  ! output: error control
  integer(i4b),intent(out)         :: err                       ! error code
  character(*),intent(out)         :: message                   ! error message
  ! -----------------------------------------------------------------------------------------------------------------------------
  ! local variables
  ! general
  integer(i4b)                     :: iLayer                    ! index of soil layer
  real(rkind)                      :: Tcrit                     ! temperature where all water is unfrozen (K)
  real(rkind)                      :: fpart1,fpart2             ! different parts of a function
  real(rkind)                      :: dpart1(1:nSoil),dpart2(1:nSoil)     ! derivatives for different parts of a function
  real(rkind)                      :: dfracCap(1:nSoil),dfInfRaw(1:nSoil)  ! derivatives for different parts of a function
  ! head boundary condition
  real(rkind)                      :: cFlux                     ! capillary flux (m s-1)
  real(rkind)                      :: dNum                      ! numerical derivative
  ! simplified Green-Ampt infiltration
  real(rkind)                      :: rootZoneLiq               ! depth of liquid water in the root zone (m)
  real(rkind)                      :: rootZoneIce               ! depth of ice in the root zone (m)
  real(rkind)                      :: availCapacity             ! available storage capacity in the root zone (m)
  real(rkind)                      :: depthWettingFront         ! depth to the wetting front (m)
  real(rkind)                      :: hydCondWettingFront       ! hydraulic conductivity at the wetting front (m s-1)
  ! saturated area associated with variable storage capacity
  real(rkind)                      :: fracCap                   ! fraction of pore space filled with liquid water and ice (-)
  real(rkind)                      :: fInfRaw                   ! infiltrating area before imposing solution constraints (-)
  real(rkind),parameter            :: maxFracCap=0.995_rkind       ! maximum fraction capacity -- used to avoid numerical problems associated with an enormous derivative
  real(rkind),parameter            :: scaleFactor=0.000001_rkind   ! scale factor for the smoothing function (-)
  real(rkind),parameter            :: qSurfScaleMax=1000._rkind    ! maximum surface runoff scaling factor (-)
  ! fraction of impermeable area associated with frozen ground
  real(rkind)                      :: alpha                     ! shape parameter in the Gamma distribution
  real(rkind)                      :: xLimg                     ! upper limit of the integral
  ! derivatives
  real(rkind)                      :: dVolFracLiq_dWat(1:nSoil)  ! derivative in vol fraction of liquid w.r.t. water state variable in root layers
  real(rkind)                      :: dVolFracIce_dWat(1:nSoil)  ! derivative in vol fraction of ice w.r.t. water state variable in root layers
  real(rkind)                      :: dVolFracLiq_dTk(1:nSoil)   ! derivative in vol fraction of liquid w.r.t. temperature in root layers
  real(rkind)                      :: dVolFracIce_dTk(1:nSoil)   ! derivative in vol fraction of ice w.r.t. temperature in root layers
  real(rkind)                      :: dRootZoneLiq_dWat(1:nSoil)       ! derivative in vol fraction of scalar root zone liquid w.r.t. water state variable in root layers
  real(rkind)                      :: dRootZoneIce_dWat(1:nSoil)       ! derivative in vol fraction of scalar root zone ice w.r.t. water state variable in root layers
  real(rkind)                      :: dRootZoneLiq_dTk(1:nSoil)        ! derivative in vol fraction of scalar root zone liquid w.r.t. temperature in root layers
  real(rkind)                      :: dRootZoneIce_dTk(1:nSoil)        ! derivative in vol fraction of scalar root zone ice w.r.t. temperature in root layers
  real(rkind)                      :: dDepthWettingFront_dWat(1:nSoil) ! derivative in scalar depth of wetting front w.r.t. water state variable in root layers
  real(rkind)                      :: dDepthWettingFront_dTk(1:nSoil)  ! derivative in scalar depth of wetting front w.r.t. temperature in root layers
  real(rkind)                      :: dxMaxInfilRate_dWat(1:nSoil)     ! derivative in scalar max infiltration rate w.r.t. water state variable in root layers
  real(rkind)                      :: dxMaxInfilRate_dTk(1:nSoil)      ! derivative in scalar max infiltration rate w.r.t. temperature in root layers
  real(rkind)                      :: dInfilArea_dWat(0:nSoil)         ! derivative in scalar infiltration rate w.r.t. water state variable in canopy or snow and root layers
  real(rkind)                      :: dInfilArea_dTk(0:nSoil)          ! derivative in scalar infiltration rate w.r.t. temperature in canopy or snow and root layers
  real(rkind)                      :: dFrozenArea_dWat(0:nSoil)        ! derivative in scalar frozen area w.r.t. water state variable in canopy or snow and root layers
  real(rkind)                      :: dFrozenArea_dTk(0:nSoil)         ! derivative in scalar frozen area w.r.t. temperature in canopy or snow and root layers
  real(rkind)                      :: dInfilRate_dWat(0:nSoil)         ! derivative in scalar infiltration rate w.r.t. water state variable in canopy or snow and root layers
  real(rkind)                      :: dInfilRate_dTk(0:nSoil)          ! derivative in scalar infiltration rate w.r.t. temperature in canopy or snow and root layers

  ! initialize error control
  err=0; message="surfaceFlx/"

  ! initialize derivatives
  dq_dHydStateVec(:) = 0._rkind
  dq_dNrgStateVec(:) = 0._rkind

  ! *****
  ! compute the surface flux and its derivative
  select case(bc_upper)

    ! *****
    ! head condition
    case(prescribedHead)
      ! surface runoff iz zero for the head condition
      scalarSurfaceRunoff = 0._rkind
      ! compute transmission and the capillary flux
      select case(ixRichards)  ! (form of Richards' equation)
        case(moisture)
          ! compute the hydraulic conductivity and diffusivity at the boundary
          surfaceHydCond = hydCond_liq(upperBoundTheta,surfaceSatHydCond,theta_res,theta_sat,vGn_m) * iceImpedeFac
          surfaceDiffuse = dPsi_dTheta(upperBoundTheta,vGn_alpha,theta_res,theta_sat,vGn_n,vGn_m) * surfaceHydCond
          ! compute the capillary flux
          cflux = -surfaceDiffuse*(scalarVolFracLiq - upperBoundTheta) / (mLayerDepth(1)*0.5_rkind)
        case(mixdform)
          ! compute the hydraulic conductivity and diffusivity at the boundary
          surfaceHydCond = hydCond_psi(upperBoundHead,surfaceSatHydCond,vGn_alpha,vGn_n,vGn_m) * iceImpedeFac
          surfaceDiffuse = realMissing
          ! compute the capillary flux
          cflux = -surfaceHydCond*(scalarMatricHeadLiq - upperBoundHead) / (mLayerDepth(1)*0.5_rkind)
        case default; err=10; message=trim(message)//"unknown form of Richards' equation"; return
      end select  ! (form of Richards' eqn)
      ! compute the total flux
      scalarSurfaceInfiltration = cflux + surfaceHydCond
      ! compute the derivative
      if(deriv_desired)then
        ! compute the hydrology derivative at the surface
        select case(ixRichards)  ! (form of Richards' equation)
          case(moisture); dq_dHydStateVec(1) = -surfaceDiffuse/(mLayerDepth(1)/2._rkind)
          case(mixdform); dq_dHydStateVec(1) = -surfaceHydCond/(mLayerDepth(1)/2._rkind)
          case default; err=10; message=trim(message)//"unknown form of Richards' equation"; return
        end select
        ! compute the energy derivative at the surface
        dq_dNrgStateVec(1) = -(dHydCond_dTemp/2._rkind)*(scalarMatricHeadLiq - upperBoundHead)/(mLayerDepth(1)*0.5_rkind) + dHydCond_dTemp/2._rkind
      else
        dNum         = 0._rkind
      end if

    ! *****
    ! flux condition
    case(liquidFlux)
      ! force infiltration to be constant over the iterations
      if(firstSplitOper)then
        ! process root layers only liquid and ice derivatives
        dVolFracLiq_dWat(:) = 0._rkind
        dVolFracIce_dWat(:) = 0._rkind
        dVolFracLiq_dTk(:)  = 0._rkind
        dVolFracIce_dTk(:)  = 0._rkind
        if(deriv_desired .and. nRoots > 0)then
          select case(ixRichards)  ! form of Richards' equation
            case(moisture)
              dVolFracLiq_dWat(:) = 1._rkind
              dVolFracIce_dWat(:) = mLayerdPsi_dTheta(:) - 1._rkind
            case(mixdform)
              do iLayer=1,nRoots
                Tcrit = crit_soilT( mLayerMatricHead(iLayer) )
                if(mLayerTemp(iLayer) < Tcrit)then
                  dVolFracLiq_dWat(iLayer) = 0._rkind
                  dVolFracIce_dWat(iLayer) = dTheta_dPsi(iLayer)
                else
                  dVolFracLiq_dWat(iLayer) = dTheta_dPsi(iLayer)
                  dVolFracIce_dWat(iLayer) = 0._rkind
                endif
              enddo
          end select ! (form of Richards' equation)
          dVolFracLiq_dTk(:) = dTheta_dTk(:) !already zeroed out if not below critical temperature
          dVolFracIce_dTk(:) = -dVolFracLiq_dTk(:) !often can and will simplify one of these terms out
        endif

        ! define the storage in the root zone (m) and derivatives
        rootZoneLiq = 0._rkind
        rootZoneIce = 0._rkind
        dRootZoneLiq_dWat(:) = 0._rkind
        dRootZoneIce_dWat(:) = 0._rkind
        dRootZoneLiq_dTk(:)  = 0._rkind
        dRootZoneIce_dTk(:)  = 0._rkind

        ! process layers where the roots extend to the bottom of the layer
        if(nRoots > 1)then
          do iLayer=1,nRoots-1
            rootZoneLiq = rootZoneLiq + mLayerVolFracLiq(iLayer)*mLayerDepth(iLayer)
            rootZoneIce = rootZoneIce + mLayerVolFracIce(iLayer)*mLayerDepth(iLayer)
            dRootZoneLiq_dWat(iLayer) = dVolFracLiq_dWat(iLayer)*mLayerDepth(iLayer)
            dRootZoneIce_dWat(iLayer) = dVolFracIce_dWat(iLayer)*mLayerDepth(iLayer)
            dRootZoneLiq_dTk(iLayer)  = dVolFracLiq_dTk(iLayer) *mLayerDepth(iLayer)
            dRootZoneIce_dTk(iLayer)  = dVolFracIce_dTk(iLayer) *mLayerDepth(iLayer)
          end do
        end if
        ! process layers where the roots end in the current layer
        rootZoneLiq = rootZoneLiq + mLayerVolFracLiq(nRoots)*(rootingDepth - iLayerHeight(nRoots-1))
        rootZoneIce = rootZoneIce + mLayerVolFracIce(nRoots)*(rootingDepth - iLayerHeight(nRoots-1))
        dRootZoneLiq_dWat(nRoots) = dVolFracLiq_dWat(nRoots)*(rootingDepth - iLayerHeight(nRoots-1))
        dRootZoneIce_dWat(nRoots) = dVolFracIce_dWat(nRoots)*(rootingDepth - iLayerHeight(nRoots-1))
        dRootZoneLiq_dTk(nRoots)  = dVolFracLiq_dTk(nRoots)* (rootingDepth - iLayerHeight(nRoots-1))
        dRootZoneIce_dTk(nRoots)  = dVolFracIce_dTk(nRoots)* (rootingDepth - iLayerHeight(nRoots-1))

        ! define available capacity to hold water (m)
        availCapacity = theta_sat*rootingDepth - rootZoneIce
        if(rootZoneLiq > availCapacity+verySmall)then
          message=trim(message)//'liquid water in the root zone exceeds capacity'
          err=20; return
        end if

        ! define the depth to the wetting front (m) and derivatives
        depthWettingFront = (rootZoneLiq/availCapacity)*rootingDepth
        dDepthWettingFront_dWat(:)=( dRootZoneLiq_dWat(:)*rootingDepth + dRootZoneIce_dWat(:)*depthWettingFront )/availCapacity
        dDepthWettingFront_dTk(:) =( dRootZoneLiq_dTk(:)*rootingDepth  + dRootZoneIce_dTk(:)*depthWettingFront  )/availCapacity

        ! define the hydraulic conductivity at depth=depthWettingFront (m s-1)
        hydCondWettingFront =  surfaceSatHydCond * ( (1._rkind - depthWettingFront/sum(mLayerDepth))**(zScale_TOPMODEL - 1._rkind) )

        ! define the maximum infiltration rate (m s-1) and derivatives
        xMaxInfilRate = hydCondWettingFront*( (wettingFrontSuction + depthWettingFront)/depthWettingFront )  ! maximum infiltration rate (m s-1)
        !write(*,'(a,1x,f9.3,1x,10(e20.10,1x))') 'depthWettingFront, surfaceSatHydCond, hydCondWettingFront, xMaxInfilRate = ', depthWettingFront, surfaceSatHydCond, hydCondWettingFront, xMaxInfilRate
        fPart1    = hydCondWettingFront
        fPart2    = (wettingFrontSuction + depthWettingFront)/depthWettingFront
        dPart1(:) = surfaceSatHydCond*(zScale_TOPMODEL - 1._rkind) * ( (1._rkind - depthWettingFront/sum(mLayerDepth))**(zScale_TOPMODEL - 2._rkind) ) * (-dDepthWettingFront_dWat(:))/sum(mLayerDepth)
        dPart2(:) = -dDepthWettingFront_dWat(:)*wettingFrontSuction / (depthWettingFront**2_i4b)
        dxMaxInfilRate_dWat(:) = fPart1*dpart2(:) + fPart2*dPart1(:)
        dPart1(:) = surfaceSatHydCond*(zScale_TOPMODEL - 1._rkind) * ( (1._rkind - depthWettingFront/sum(mLayerDepth))**(zScale_TOPMODEL - 2._rkind) ) * (-dDepthWettingFront_dTk(:))/sum(mLayerDepth)
        dPart2(:) = -dDepthWettingFront_dTk(:)*wettingFrontSuction / (depthWettingFront**2_i4b)
        dxMaxInfilRate_dTk(:)  = fPart1*dpart2(:) + fPart2*dPart1(:)

        ! define the infiltrating area and derivatives for the non-frozen part of the cell/basin
        if(qSurfScale < qSurfScaleMax)then
          fracCap         = rootZoneLiq/(maxFracCap*availCapacity)                              ! fraction of available root zone filled with water
          fInfRaw         = 1._rkind - exp(-qSurfScale*(1._rkind - fracCap))                          ! infiltrating area -- allowed to violate solution constraints
          scalarInfilArea = min(0.5_rkind*(fInfRaw + sqrt(fInfRaw**2_i4b + scaleFactor)), 1._rkind)   ! infiltrating area -- constrained
          if (0.5_rkind*(fInfRaw + sqrt(fInfRaw**2_i4b + scaleFactor))< 1._rkind) then
            dfracCap(:) = ( dRootZoneLiq_dWat(:)/maxFracCap + dRootZoneIce_dWat(:)*fracCap )/availCapacity
            dfInfRaw(:) = -qSurfScale*dfracCap(:) * exp(-qSurfScale*(1._rkind - fracCap))
            dInfilArea_dWat(1:nSoil) = 0.5_rkind*dfInfRaw(:) * (1._rkind + fInfRaw/sqrt(fInfRaw**2_i4b + scaleFactor))
            dfracCap(:) = ( dRootZoneLiq_dTk(:)/maxFracCap + dRootZoneIce_dTk(:)*fracCap )/availCapacity
            dfInfRaw(:) = -qSurfScale*dfracCap(:) * exp(-qSurfScale*(1._rkind - fracCap))
            dInfilArea_dTk(1:nSoil)  = 0.5_rkind*dfInfRaw(:) * (1._rkind + fInfRaw/sqrt(fInfRaw**2_i4b + scaleFactor))
          else ! scalarInfilArea = 1._rkind
            dInfilArea_dWat(1:nSoil) = 0._rkind
            dInfilArea_dTk(1:nSoil)  = 0._rkind
          endif
        else
          scalarInfilArea = 1._rkind
          dInfilArea_dWat(1:nSoil) = 0._rkind
          dInfilArea_dTk(1:nSoil)  = 0._rkind
        endif
        dInfilArea_dWat(0) = 0._rkind
        dInfilArea_dTk(0)  = 0._rkind

        ! check to ensure we are not infiltrating into a fully saturated column
        if(ixIce<nRoots)then
          if(sum(mLayerVolFracLiq(ixIce+1:nRoots)*mLayerDepth(ixIce+1:nRoots)) > 0.9999_rkind*theta_sat*sum(mLayerDepth(ixIce+1:nRoots))) scalarInfilArea=0._rkind
        endif

        ! define the impermeable area and derivatives due to frozen ground
        if(rootZoneIce > tiny(rootZoneIce))then  ! (avoid divide by zero)
          alpha            = 1._rkind/(soilIceCV**2_i4b)        ! shape parameter in the Gamma distribution
          xLimg            = alpha*soilIceScale/rootZoneIce  ! upper limit of the integral

          !if we use this, we will have a derivative of scalarFrozenArea w.r.t. water and temperature in each layer (through mLayerVolFracIce)
          scalarFrozenArea = 0._rkind
          dFrozenArea_dWat(1:nSoil) = 0._rkind
          dFrozenArea_dTk(1:nSoil)  = 0._rkind
        else
          scalarFrozenArea = 0._rkind
          dFrozenArea_dWat(1:nSoil) = 0._rkind
          dFrozenArea_dTk(1:nSoil)  = 0._rkind
        end if
        dFrozenArea_dWat(0) = 0._rkind
        dFrozenArea_dTk(0)  = 0._rkind


        if (xMaxInfilRate < scalarRainPlusMelt) then ! = dxMaxInfilRate_d, dependent on layers not at surface
          dInfilRate_dWat(0) = 0._rkind
          dInfilRate_dTk(0)  = 0._rkind
          dInfilRate_dWat(1:nSoil) = dxMaxInfilRate_dWat(:)
          dInfilRate_dTk(1:nSoil)  = dxMaxInfilRate_dTk(:)
        else ! = dRainPlusMelt_d, dependent on above layer (canopy or snow) water and temp
          dInfilRate_dWat(0) = above_soilLiqFluxDeriv*above_soilFracLiq
          dInfilRate_dTk(0)  = above_soilLiqFluxDeriv*above_soildLiq_dTk
          dInfilRate_dWat(1:nSoil) = 0._rkind
          dInfilRate_dTk(1:nSoil)  = 0._rkind
        endif

        ! dq w.r.t. infiltration only, scalarRainPlusMelt accounted for in computJacob module
        dq_dHydStateVec(:) = (1._rkind - scalarFrozenArea) * ( dInfilArea_dWat(:)*min(scalarRainPlusMelt,xMaxInfilRate) + scalarInfilArea*dInfilRate_dWat(:) ) +&
                              (-dFrozenArea_dWat(:))*scalarInfilArea*min(scalarRainPlusMelt,xMaxInfilRate)
        dq_dNrgStateVec(:) = (1._rkind - scalarFrozenArea) * ( dInfilArea_dTk(:) *min(scalarRainPlusMelt,xMaxInfilRate) + scalarInfilArea*dInfilRate_dTk(:)  ) +&
                              (-dFrozenArea_dTk(:)) *scalarInfilArea*min(scalarRainPlusMelt,xMaxInfilRate)

      else ! do not compute infiltration after first flux call in a splitting operation
        dq_dHydStateVec(:) = 0._rkind
        dq_dNrgStateVec(:) = 0._rkind

      end if ! (if desire to compute infiltration)

      ! compute infiltration (m s-1), if after first flux call in a splitting operation does not change
      scalarSurfaceInfiltration = (1._rkind - scalarFrozenArea)*scalarInfilArea*min(scalarRainPlusMelt,xMaxInfilRate)

      ! compute surface runoff (m s-1)
      scalarSurfaceRunoff = scalarRainPlusMelt - scalarSurfaceInfiltration

      ! set surface hydraulic conductivity and diffusivity to missing (not used for flux condition)
      surfaceHydCond = realMissing
      surfaceDiffuse = realMissing

    case default; err=20; message=trim(message)//'unknown upper boundary condition for soil hydrology'; return

  end select  ! (type of upper boundary condition)

end subroutine surfaceFlx


! ***************************************************************************************************************
! private subroutine iLayerFlux: compute the fluxes and derivatives at layer interfaces
! ***************************************************************************************************************
subroutine iLayerFlux(&
                      ! input: model control
                      deriv_desired,             & ! intent(in): flag indicating if derivatives are desired
                      ixRichards,                & ! intent(in): index defining the form of Richards' equation (moisture or mixdform)
                      ! input: state variables (adjacent layers)
                      nodeMatricHeadLiqTrial,    & ! intent(in): liquid matric head at the soil nodes (m)
                      nodeVolFracLiqTrial,       & ! intent(in): volumetric liquid water content at the soil nodes (-)
                      ! input: model coordinate variables (adjacent layers)
                      nodeHeight,                & ! intent(in): height of the soil nodes (m)
                      ! input: temperature derivatives
                      dPsiLiq_dTemp,             & ! intent(in): derivative in liquid water matric potential w.r.t. temperature (m K-1)
                      dHydCond_dTemp,            & ! intent(in): derivative in hydraulic conductivity w.r.t temperature (m s-1 K-1)
                      ! input: transmittance (adjacent layers)
                      nodeHydCondTrial,          & ! intent(in): hydraulic conductivity at the soil nodes (m s-1)
                      nodeDiffuseTrial,          & ! intent(in): hydraulic diffusivity at the soil nodes (m2 s-1)
                      ! input: transmittance derivatives (adjacent layers)
                      dHydCond_dVolLiq,          & ! intent(in): derivative in hydraulic conductivity w.r.t. change in volumetric liquid water content (m s-1)
                      dDiffuse_dVolLiq,          & ! intent(in): derivative in hydraulic diffusivity w.r.t. change in volumetric liquid water content (m2 s-1)
                      dHydCond_dMatric,          & ! intent(in): derivative in hydraulic conductivity w.r.t. change in matric head (s-1)
                      ! output: tranmsmittance at the layer interface (scalars)
                      iLayerHydCond,             & ! intent(out): hydraulic conductivity at the interface between layers (m s-1)
                      iLayerDiffuse,             & ! intent(out): hydraulic diffusivity at the interface between layers (m2 s-1)
                      ! output: vertical flux at the layer interface (scalars)
                      iLayerLiqFluxSoil,         & ! intent(out): vertical flux of liquid water at the layer interface (m s-1)
                      ! output: derivatives in fluxes w.r.t. state variables -- matric head or volumetric lquid water -- in the layer above and layer below (m s-1 or s-1)
                      dq_dHydStateAbove,         & ! intent(out): derivatives in the flux w.r.t. matric head or volumetric lquid water in the layer above (m s-1 or s-1)
                      dq_dHydStateBelow,         & ! intent(out): derivatives in the flux w.r.t. matric head or volumetric lquid water in the layer below (m s-1 or s-1)
                      ! output: derivatives in fluxes w.r.t. energy state variables -- now just temperature -- in the layer above and layer below (m s-1 K-1)
                      dq_dNrgStateAbove,         & ! intent(out): derivatives in the flux w.r.t. temperature in the layer above (m s-1 K-1)
                      dq_dNrgStateBelow,         & ! intent(out): derivatives in the flux w.r.t. temperature in the layer below (m s-1 K-1)
                      ! output: error control
                      err,message)                 ! intent(out): error control
  ! ------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  ! input: model control
  logical(lgt),intent(in)       :: deriv_desired               ! flag indicating if derivatives are desired
  integer(i4b),intent(in)       :: ixRichards                  ! index defining the option for Richards' equation (moisture or mixdform)
  ! input: state variables
  real(rkind),intent(in)           :: nodeMatricHeadLiqTrial(:)   ! liquid matric head at the soil nodes (m)
  real(rkind),intent(in)           :: nodeVolFracLiqTrial(:)      ! volumetric fraction of liquid water at the soil nodes (-)
  ! input: model coordinate variables
  real(rkind),intent(in)           :: nodeHeight(:)               ! height at the mid-point of the lower layer (m)
  ! input: temperature derivatives
  real(rkind),intent(in)           :: dPsiLiq_dTemp(:)            ! derivative in liquid water matric potential w.r.t. temperature (m K-1)
  real(rkind),intent(in)           :: dHydCond_dTemp(:)           ! derivative in hydraulic conductivity w.r.t temperature (m s-1 K-1)
  ! input: transmittance
  real(rkind),intent(in)           :: nodeHydCondTrial(:)         ! hydraulic conductivity at layer mid-points (m s-1)
  real(rkind),intent(in)           :: nodeDiffuseTrial(:)         ! diffusivity at layer mid-points (m2 s-1)
  ! input: transmittance derivatives
  real(rkind),intent(in)           :: dHydCond_dVolLiq(:)         ! derivative in hydraulic conductivity w.r.t volumetric liquid water content (m s-1)
  real(rkind),intent(in)           :: dDiffuse_dVolLiq(:)         ! derivative in hydraulic diffusivity w.r.t volumetric liquid water content (m2 s-1)
  real(rkind),intent(in)           :: dHydCond_dMatric(:)         ! derivative in hydraulic conductivity w.r.t matric head (m s-1)
  ! output: tranmsmittance at the layer interface (scalars)
  real(rkind),intent(out)          :: iLayerHydCond               ! hydraulic conductivity at the interface between layers (m s-1)
  real(rkind),intent(out)          :: iLayerDiffuse               ! hydraulic diffusivity at the interface between layers (m2 s-1)
  ! output: vertical flux at the layer interface (scalars)
  real(rkind),intent(out)          :: iLayerLiqFluxSoil           ! vertical flux of liquid water at the layer interface (m s-1)
  ! output: derivatives in fluxes w.r.t. state variables -- matric head or volumetric lquid water -- in the layer above and layer below (m s-1 or s-1)
  real(rkind),intent(out)          :: dq_dHydStateAbove           ! derivatives in the flux w.r.t. matric head or volumetric lquid water in the layer above (m s-1 or s-1)
  real(rkind),intent(out)          :: dq_dHydStateBelow           ! derivatives in the flux w.r.t. matric head or volumetric lquid water in the layer below (m s-1 or s-1)
  ! output: derivatives in fluxes w.r.t. energy state variables -- now just temperature -- in the layer above and layer below (m s-1 K-1)
  real(rkind),intent(out)          :: dq_dNrgStateAbove           ! derivatives in the flux w.r.t. temperature in the layer above (m s-1 K-1)
  real(rkind),intent(out)          :: dq_dNrgStateBelow           ! derivatives in the flux w.r.t. temperature in the layer below (m s-1 K-1)
  ! output: error control
  integer(i4b),intent(out)      :: err                         ! error code
  character(*),intent(out)      :: message                     ! error message
  ! ------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  ! local variables (named variables to provide index of 2-element vectors)
  integer(i4b),parameter        :: ixUpper=1                   ! index of upper node in the 2-element vectors
  integer(i4b),parameter        :: ixLower=2                   ! index of lower node in the 2-element vectors
  logical(lgt),parameter        :: useGeometric=.false.        ! switch between the arithmetic and geometric mean
  ! local variables (Darcy flux)
  real(rkind)                      :: dPsi                        ! spatial difference in matric head (m)
  real(rkind)                      :: dLiq                        ! spatial difference in volumetric liquid water (-)
  real(rkind)                      :: dz                          ! spatial difference in layer mid-points (m)
  real(rkind)                      :: cflux                       ! capillary flux (m s-1)
  ! local variables (derivative in Darcy's flux)
  real(rkind)                      :: dHydCondIface_dVolLiqAbove  ! derivative in hydraulic conductivity at layer interface w.r.t. volumetric liquid water content in layer above
  real(rkind)                      :: dHydCondIface_dVolLiqBelow  ! derivative in hydraulic conductivity at layer interface w.r.t. volumetric liquid water content in layer below
  real(rkind)                      :: dDiffuseIface_dVolLiqAbove  ! derivative in hydraulic diffusivity at layer interface w.r.t. volumetric liquid water content in layer above
  real(rkind)                      :: dDiffuseIface_dVolLiqBelow  ! derivative in hydraulic diffusivity at layer interface w.r.t. volumetric liquid water content in layer below
  real(rkind)                      :: dHydCondIface_dMatricAbove  ! derivative in hydraulic conductivity at layer interface w.r.t. matric head in layer above
  real(rkind)                      :: dHydCondIface_dMatricBelow  ! derivative in hydraulic conductivity at layer interface w.r.t. matric head in layer below
 ! ------------------------------------------------------------------------------------------------------------------------------------------------------------------------
 ! initialize error control
 err=0; message="iLayerFlux/"

  ! *****
  ! compute the vertical flux of liquid water
  ! compute the hydraulic conductivity at the interface
  if(useGeometric)then
    iLayerHydCond   = sqrt(nodeHydCondTrial(ixLower)   * nodeHydCondTrial(ixUpper))
  else
    iLayerHydCond   = (nodeHydCondTrial(ixLower)   + nodeHydCondTrial(ixUpper))*0.5_rkind
  end if
  
  dz = nodeHeight(ixLower) - nodeHeight(ixUpper)
  ! compute the capillary flux
  select case(ixRichards)  ! (form of Richards' equation)
    case(moisture)
    iLayerDiffuse = sqrt(nodeDiffuseTrial(ixLower) * nodeDiffuseTrial(ixUpper))
    dLiq          = nodeVolFracLiqTrial(ixLower) - nodeVolFracLiqTrial(ixUpper)
    cflux         = -iLayerDiffuse * dLiq/dz
    case(mixdform)
    iLayerDiffuse = realMissing
    dPsi          = nodeMatricHeadLiqTrial(ixLower) - nodeMatricHeadLiqTrial(ixUpper)
    cflux         = -iLayerHydCond * dPsi/dz
    case default; err=10; message=trim(message)//"unable to identify option for Richards' equation"; return
  end select
  ! compute the total flux (add gravity flux, positive downwards)
  iLayerLiqFluxSoil = cflux + iLayerHydCond

  ! ** compute the derivatives
  if(deriv_desired)then
    select case(ixRichards)  ! (form of Richards' equation)
      case(moisture)
        ! still need to implement arithmetric mean for the moisture-based form
        if(.not.useGeometric)then
          message=trim(message)//'only currently implemented for geometric mean -- change local flag'
          err=20; return
        end if
        ! derivatives in hydraulic conductivity at the layer interface (m s-1)
        dHydCondIface_dVolLiqAbove = dHydCond_dVolLiq(ixUpper)*nodeHydCondTrial(ixLower) * 0.5_rkind/max(iLayerHydCond,verySmall)
        dHydCondIface_dVolLiqBelow = dHydCond_dVolLiq(ixLower)*nodeHydCondTrial(ixUpper) * 0.5_rkind/max(iLayerHydCond,verySmall)
        ! derivatives in hydraulic diffusivity at the layer interface (m2 s-1)
        dDiffuseIface_dVolLiqAbove = dDiffuse_dVolLiq(ixUpper)*nodeDiffuseTrial(ixLower) * 0.5_rkind/max(iLayerDiffuse,verySmall)
        dDiffuseIface_dVolLiqBelow = dDiffuse_dVolLiq(ixLower)*nodeDiffuseTrial(ixUpper) * 0.5_rkind/max(iLayerDiffuse,verySmall)
        ! derivatives in the flux w.r.t. volumetric liquid water content
        dq_dHydStateAbove = -dDiffuseIface_dVolLiqAbove*dLiq/dz + iLayerDiffuse/dz + dHydCondIface_dVolLiqAbove
        dq_dHydStateBelow = -dDiffuseIface_dVolLiqBelow*dLiq/dz - iLayerDiffuse/dz + dHydCondIface_dVolLiqBelow
      case(mixdform)
        ! derivatives in hydraulic conductivity
        if(useGeometric)then
          dHydCondIface_dMatricAbove = dHydCond_dMatric(ixUpper)*nodeHydCondTrial(ixLower) * 0.5_rkind/max(iLayerHydCond,verySmall)
          dHydCondIface_dMatricBelow = dHydCond_dMatric(ixLower)*nodeHydCondTrial(ixUpper) * 0.5_rkind/max(iLayerHydCond,verySmall)
        else
          dHydCondIface_dMatricAbove = dHydCond_dMatric(ixUpper)/2._rkind
          dHydCondIface_dMatricBelow = dHydCond_dMatric(ixLower)/2._rkind
        end if
        ! derivatives in the flux w.r.t. matric head
        dq_dHydStateAbove = -dHydCondIface_dMatricAbove*dPsi/dz + iLayerHydCond/dz + dHydCondIface_dMatricAbove
        dq_dHydStateBelow = -dHydCondIface_dMatricBelow*dPsi/dz - iLayerHydCond/dz + dHydCondIface_dMatricBelow
        ! derivative in the flux w.r.t. temperature
        dq_dNrgStateAbove = -(dHydCond_dTemp(ixUpper)/2._rkind)*dPsi/dz + iLayerHydCond*dPsiLiq_dTemp(ixUpper)/dz + dHydCond_dTemp(ixUpper)/2._rkind
        dq_dNrgStateBelow = -(dHydCond_dTemp(ixLower)/2._rkind)*dPsi/dz - iLayerHydCond*dPsiLiq_dTemp(ixLower)/dz + dHydCond_dTemp(ixLower)/2._rkind
      case default; err=10; message=trim(message)//"unknown form of Richards' equation"; return
    end select
  else
    dq_dHydStateAbove = realMissing
    dq_dHydStateBelow = realMissing
  end if

end subroutine iLayerFlux


! ***************************************************************************************************************
! private subroutine qDrainFlux: compute the drainage flux from the bottom of the soil profile and its derivative
! ***************************************************************************************************************
subroutine qDrainFlux(&
                      ! input: model control
                      deriv_desired,             & ! intent(in): flag indicating if derivatives are desired
                      ixRichards,                & ! intent(in): index defining the form of Richards' equation (moisture or mixdform)
                      bc_lower,                  & ! intent(in): index defining the type of boundary conditions
                      ! input: state variables
                      nodeMatricHeadLiq,         & ! intent(in): liquid matric head in the lowest unsaturated node (m)
                      nodeVolFracLiq,            & ! intent(in): volumetric liquid water content the lowest unsaturated node (-)
                      ! input: model coordinate variables
                      nodeDepth,                 & ! intent(in): depth of the lowest unsaturated soil layer (m)
                      nodeHeight,                & ! intent(in): height of the lowest unsaturated soil node (m)
                      ! input: boundary conditions
                      lowerBoundHead,            & ! intent(in): lower boundary condition (m)
                      lowerBoundTheta,           & ! intent(in): lower boundary condition (-)
                      ! input: derivative in soil water characteristix
                      node_dPsi_dTheta,         & ! intent(in): derivative of the soil moisture characteristic w.r.t. theta (m)
                      ! input: transmittance
                      surfaceSatHydCond,         & ! intent(in): saturated hydraulic conductivity at the surface (m s-1)
                      bottomSatHydCond,          & ! intent(in): saturated hydraulic conductivity at the bottom of the unsaturated zone (m s-1)
                      nodeHydCond,               & ! intent(in): hydraulic conductivity at the node itself (m s-1)
                      iceImpedeFac,              & ! intent(in): ice impedence factor in the lower-most soil layer (-)
                      ! input: transmittance derivatives
                      dHydCond_dVolLiq,          & ! intent(in): derivative in hydraulic conductivity w.r.t. volumetric liquid water content (m s-1)
                      dHydCond_dMatric,          & ! intent(in): derivative in hydraulic conductivity w.r.t. matric head (s-1)
                      dHydCond_dTemp,            & ! intent(in): derivative in hydraulic conductivity w.r.t temperature (m s-1 K-1)
                      ! input: soil parameters
                      vGn_alpha,                 & ! intent(in): van Genutchen "alpha" parameter (m-1)
                      vGn_n,                     & ! intent(in): van Genutchen "n" parameter (-)
                      vGn_m,                     & ! intent(in): van Genutchen "m" parameter (-)
                      theta_sat,                 & ! intent(in): soil porosity (-)
                      theta_res,                 & ! intent(in): soil residual volumetric water content (-)
                      kAnisotropic,              & ! intent(in): anisotropy factor for lateral hydraulic conductivity (-)
                      zScale_TOPMODEL,           & ! intent(in): TOPMODEL scaling factor (m)
                      ! output: hydraulic conductivity and diffusivity at the surface
                      bottomHydCond,             & ! intent(out): hydraulic conductivity at the bottom of the unsatuarted zone (m s-1)
                      bottomDiffuse,             & ! intent(out): hydraulic diffusivity at the bottom of the unsatuarted zone (m2 s-1)
                      ! output: drainage flux from the bottom of the soil profile
                      scalarDrainage,            & ! intent(out): drainage flux from the bottom of the soil profile (m s-1)
                      ! output: derivatives in drainage flux
                      dq_dHydStateUnsat,         & ! intent(out): change in drainage flux w.r.t. change in hydrology state variable in lowest unsaturated node (m s-1 or s-1)
                      dq_dNrgStateUnsat,         & ! intent(out): change in drainage flux w.r.t. change in energy state variable in lowest unsaturated node (m s-1 K-1)
                      ! output: error control
                      err,message)                 ! intent(out): error control
  USE soil_utils_module,only:volFracLiq      ! compute volumetric fraction of liquid water as a function of matric head (-)
  USE soil_utils_module,only:matricHead      ! compute matric head as a function of volumetric fraction of liquid water (m)
  USE soil_utils_module,only:hydCond_psi     ! compute hydraulic conductivity as a function of matric head (m s-1)
  USE soil_utils_module,only:hydCond_liq     ! compute hydraulic conductivity as a function of volumetric liquid water content (m s-1)
  USE soil_utils_module,only:dPsi_dTheta     ! compute derivative of the soil moisture characteristic w.r.t. theta (m)
  ! compute infiltraton at the surface and its derivative w.r.t. mass in the upper soil layer
  implicit none
  ! -----------------------------------------------------------------------------------------------------------------------------
  ! input: model control
  logical(lgt),intent(in)       :: deriv_desired             ! flag to indicate if derivatives are desired
  integer(i4b),intent(in)       :: ixRichards                ! index defining the option for Richards' equation (moisture or mixdform)
  integer(i4b),intent(in)       :: bc_lower                  ! index defining the type of boundary conditions
  ! input: state and diagnostic variables
  real(rkind),intent(in)           :: nodeMatricHeadLiq         ! liquid matric head in the lowest unsaturated node (m)
  real(rkind),intent(in)           :: nodeVolFracLiq            ! volumetric liquid water content in the lowest unsaturated node (-)
  ! input: model coordinate variables
  real(rkind),intent(in)           :: nodeDepth                 ! depth of the lowest unsaturated soil layer (m)
  real(rkind),intent(in)           :: nodeHeight                ! height of the lowest unsaturated soil node (m)
  ! input: diriclet boundary conditions
  real(rkind),intent(in)           :: lowerBoundHead            ! lower boundary condition for matric head (m)
  real(rkind),intent(in)           :: lowerBoundTheta           ! lower boundary condition for volumetric liquid water content (-)
  ! input: derivative in soil water characteristix
  real(rkind),intent(in)           :: node_dPsi_dTheta         ! derivative of the soil moisture characteristic w.r.t. theta (m)
  ! input: transmittance
  real(rkind),intent(in)           :: surfaceSatHydCond         ! saturated hydraulic conductivity at the surface (m s-1)
  real(rkind),intent(in)           :: bottomSatHydCond          ! saturated hydraulic conductivity at the bottom of the unsaturated zone (m s-1)
  real(rkind),intent(in)           :: nodeHydCond               ! hydraulic conductivity at the node itself (m s-1)
  real(rkind),intent(in)           :: iceImpedeFac              ! ice impedence factor in the upper-most soil layer (-)
  ! input: transmittance derivatives
  real(rkind),intent(in)           :: dHydCond_dVolLiq          ! derivative in hydraulic conductivity w.r.t. volumetric liquid water content (m s-1)
  real(rkind),intent(in)           :: dHydCond_dMatric          ! derivative in hydraulic conductivity w.r.t. matric head (s-1)
  real(rkind),intent(in)           :: dHydCond_dTemp            ! derivative in hydraulic conductivity w.r.t temperature (m s-1 K-1)
  ! input: soil parameters
  real(rkind),intent(in)           :: vGn_alpha                 ! van Genutchen "alpha" parameter (m-1)
  real(rkind),intent(in)           :: vGn_n                     ! van Genutchen "n" parameter (-)
  real(rkind),intent(in)           :: vGn_m                     ! van Genutchen "m" parameter (-)
  real(rkind),intent(in)           :: theta_sat                 ! soil porosity (-)
  real(rkind),intent(in)           :: theta_res                 ! soil residual volumetric water content (-)
  real(rkind),intent(in)           :: kAnisotropic              ! anisotropy factor for lateral hydraulic conductivity (-)
  real(rkind),intent(in)           :: zScale_TOPMODEL           ! scale factor for TOPMODEL-ish baseflow parameterization (m)
  ! -----------------------------------------------------------------------------------------------------------------------------
  ! output: hydraulic conductivity at the bottom of the unsaturated zone
  real(rkind),intent(out)          :: bottomHydCond             ! hydraulic conductivity at the bottom of the unsaturated zone (m s-1)
  real(rkind),intent(out)          :: bottomDiffuse             ! hydraulic diffusivity at the bottom of the unsatuarted zone (m2 s-1)
  ! output: drainage flux from the bottom of the soil profile
  real(rkind),intent(out)          :: scalarDrainage            ! drainage flux from the bottom of the soil profile (m s-1)
  ! output: derivatives in drainage flux
  real(rkind),intent(out)          :: dq_dHydStateUnsat         ! change in drainage flux w.r.t. change in state variable in lowest unsaturated node (m s-1 or s-1)
  real(rkind),intent(out)          :: dq_dNrgStateUnsat         ! change in drainage flux w.r.t. change in energy state variable in lowest unsaturated node (m s-1 K-1)
  ! output: error control
  integer(i4b),intent(out)      :: err                       ! error code
  character(*),intent(out)      :: message                   ! error message
  ! -----------------------------------------------------------------------------------------------------------------------------
  ! local variables
  real(rkind)                      :: zWater                    ! effective water table depth (m)
  real(rkind)                      :: nodePsi                   ! matric head in the lowest unsaturated node (m)
  real(rkind)                      :: cflux                     ! capillary flux (m s-1)
  ! -----------------------------------------------------------------------------------------------------------------------------
  ! initialize error control
  err=0; message="qDrainFlux/"

  ! determine lower boundary condition
  select case(bc_lower)

    case(prescribedHead)
      ! compute flux
      select case(ixRichards)
        case(moisture) ! (moisture-based form of Richards' equation)
          ! compute the hydraulic conductivity and diffusivity at the boundary
          bottomHydCond = hydCond_liq(lowerBoundTheta,bottomSatHydCond,theta_res,theta_sat,vGn_m) * iceImpedeFac
          bottomDiffuse = dPsi_dTheta(lowerBoundTheta,vGn_alpha,theta_res,theta_sat,vGn_n,vGn_m) * bottomHydCond
          ! compute the capillary flux
          cflux = -bottomDiffuse*(lowerBoundTheta - nodeVolFracLiq) / (nodeDepth*0.5_rkind)
        case(mixdform)
          ! compute the hydraulic conductivity and diffusivity at the boundary
          bottomHydCond = hydCond_psi(lowerBoundHead,bottomSatHydCond,vGn_alpha,vGn_n,vGn_m) * iceImpedeFac
          bottomDiffuse = realMissing
          ! compute the capillary flux
          cflux = -bottomHydCond*(lowerBoundHead  - nodeMatricHeadLiq) / (nodeDepth*0.5_rkind)
        case default; err=10; message=trim(message)//"unknown form of Richards' equation"; return
      end select  ! (form of Richards' eqn)
      scalarDrainage = cflux + bottomHydCond

      ! compute derivatives
      if(deriv_desired)then
        ! hydrology derivatives
        select case(ixRichards)  ! (form of Richards' equation)
          case(moisture); dq_dHydStateUnsat = bottomDiffuse/(nodeDepth/2._rkind)
          case(mixdform); dq_dHydStateUnsat = bottomHydCond/(nodeDepth/2._rkind)
          case default; err=10; message=trim(message)//"unknown form of Richards' equation"; return
        end select
        ! energy derivatives
        dq_dNrgStateUnsat = -(dHydCond_dTemp/2._rkind)*(lowerBoundHead  - nodeMatricHeadLiq)/(nodeDepth*0.5_rkind) + dHydCond_dTemp/2._rkind
      else     ! (do not desire derivatives)
        dq_dHydStateUnsat = realMissing
        dq_dNrgStateUnsat = realMissing
      end if

    case(funcBottomHead) !function of matric head in the bottom layer
      ! compute flux
      select case(ixRichards)
        case(moisture); nodePsi = matricHead(nodeVolFracLiq,vGn_alpha,theta_res,theta_sat,vGn_n,vGn_m)
        case(mixdform); nodePsi = nodeMatricHeadLiq
      end select
      zWater = nodeHeight - nodePsi
      scalarDrainage = kAnisotropic*surfaceSatHydCond * exp(-zWater/zScale_TOPMODEL)

      ! compute derivatives
      if(deriv_desired)then
        ! hydrology derivatives
        select case(ixRichards)  ! (form of Richards' equation)
          case(moisture); dq_dHydStateUnsat = kAnisotropic*surfaceSatHydCond * node_dPsi_dTheta*exp(-zWater/zScale_TOPMODEL)/zScale_TOPMODEL
          case(mixdform); dq_dHydStateUnsat = kAnisotropic*surfaceSatHydCond * exp(-zWater/zScale_TOPMODEL)/zScale_TOPMODEL
          case default; err=10; message=trim(message)//"unknown form of Richards' equation"; return
        end select
        ! energy derivatives
        err=20; message=trim(message)//"not yet implemented energy derivatives"; return
      else     ! (do not desire derivatives)
        dq_dHydStateUnsat = realMissing
        dq_dNrgStateUnsat = realMissing
      end if

    case(freeDrainage)
      ! compute flux
      scalarDrainage = nodeHydCond*kAnisotropic

      ! compute derivatives
      if(deriv_desired)then
        ! hydrology derivatives
        select case(ixRichards)  ! (form of Richards' equation)
          case(moisture); dq_dHydStateUnsat = dHydCond_dVolLiq*kAnisotropic
          case(mixdform); dq_dHydStateUnsat = dHydCond_dMatric*kAnisotropic
          case default; err=10; message=trim(message)//"unknown form of Richards' equation"; return
        end select
        ! energy derivatives
        dq_dNrgStateUnsat = dHydCond_dTemp*kAnisotropic
      else     ! (do not desire derivatives)
        dq_dHydStateUnsat = realMissing
        dq_dNrgStateUnsat = realMissing
      end if

    case(zeroFlux)
      scalarDrainage = 0._rkind
      if(deriv_desired)then
        dq_dHydStateUnsat = 0._rkind
        dq_dNrgStateUnsat = 0._rkind
      else
        dq_dHydStateUnsat = realMissing
        dq_dNrgStateUnsat = realMissing
      end if

    case default; err=20; message=trim(message)//'unknown lower boundary condition for soil hydrology'; return

  end select ! (type of boundary condition)

end subroutine qDrainFlux

end module soilLiqFlx_module

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

module ssdNrgFlux_module

! data types
USE nrtype

! data types
USE data_types,only:var_d           ! x%var(:)       (rkind)
USE data_types,only:var_dlength     ! x%var(:)%dat   (rkind)
USE data_types,only:var_ilength     ! x%var(:)%dat   (i4b)

! physical constants
USE multiconst,only:&
                    sb,          & ! Stefan Boltzman constant      (W m-2 K-4)
                    Em_Sno,      & ! emissivity of snow            (-)
                    LH_fus,      & ! latent heat of fusion         (J kg-1)
                    LH_vap,      & ! latent heat of vaporization   (J kg-1)
                    LH_sub,      & ! latent heat of sublimation    (J kg-1)
                    gravity,     & ! gravitational acceleteration  (m s-2)
                    Tfreeze,     & ! freezing point of pure water  (K)
                    iden_air,    & ! intrinsic density of air      (kg m-3)
                    iden_ice,    & ! intrinsic density of ice      (kg m-3)
                    iden_water,  &  ! intrinsic density of water    (kg m-3)
                    ! specific heat
                    Cp_air,      & ! specific heat of air          (J kg-1 K-1)
                    Cp_water,    & ! specific heat of liquid water (J kg-1 K-1)
                    ! thermal conductivity
                    lambda_air,  & ! thermal conductivity of air   (J s-1 m-1)
                    lambda_ice,  & ! thermal conductivity of ice   (J s-1 m-1)
                    lambda_water   ! thermal conductivity of water (J s-1 m-1)


! missing values
USE globalData,only:integerMissing  ! missing integer
USE globalData,only:realMissing     ! missing real number

! named variables for snow and soil
USE globalData,only:iname_snow     ! named variables for snow
USE globalData,only:iname_soil     ! named variables for soil

! named variables
USE var_lookup,only:iLookPROG       ! named variables for structure elements
USE var_lookup,only:iLookDIAG       ! named variables for structure elements
USE var_lookup,only:iLookFLUX       ! named variables for structure elements
USE var_lookup,only:iLookPARAM      ! named variables for structure elements
USE var_lookup,only:iLookINDEX      ! named variables for structure elements

! model decisions
USE globalData,only:model_decisions                         ! model decision structure
USE var_lookup,only:iLookDECISIONS                          ! named variables for elements of the decision structure

! provide access to look-up values for model decisions
USE mDecisions_module,only:      &
 ! look-up values for method used to compute derivative
 numerical,                      & ! numerical solution
 analytical,                     & ! analytical solution
 ! look-up values for choice of boundary conditions for thermodynamics
 prescribedTemp,                 & ! prescribed temperature
 energyFlux,                     & ! energy flux
 zeroFlux,                       & ! zero flux
 ! look-up values for choice of boundary conditions for soil hydrology
 prescribedHead,                 & ! prescribed head
 ! look-up values for choice of thermal conductivity representation for snow
 Yen1965,                        & ! Yen (1965)
 Mellor1977,                     & ! Mellor (1977)
 Jordan1991,                     & ! Jordan (1991)
 Smirnova2000,                   & ! Smirnova et al. (2000)
 ! look-up values for choice of thermal conductivity representation for soil
 funcSoilWet,                    & ! function of soil wetness
 mixConstit,                     & ! mixture of constituents
 hanssonVZJ,                     & ! test case for the mizoguchi lab experiment, Hansson et al. VZJ 2004
 ! look-up values for the form of Richards' equation
 moisture,                       & ! moisture-based form of Richards' equation
 mixdform                          ! mixed form of Richards' equation

! -------------------------------------------------------------------------------------------------
implicit none
private
public::ssdNrgFlux
! global parameters
real(rkind),parameter            :: dx=1.e-10_rkind             ! finite difference increment (K)
contains


! **********************************************************************************************************
! public subroutine ssdNrgFlux: compute energy fluxes and derivatives at layer interfaces
! **********************************************************************************************************
subroutine ssdNrgFlux(&
                      ! input: model control
                      scalarSolution,                     & ! intent(in):    flag to indicate the scalar solution
                      deriv_desired,                & ! intent(in): flag indicating if derivatives are desired
                      ! input: fluxes and derivatives at the upper boundary
                      groundNetFlux,                      & ! intent(in):    total flux at the ground surface (W m-2)
                      dGroundNetFlux_dGroundTemp,         & ! intent(in):    derivative in total ground surface flux w.r.t. ground temperature (W m-2 K-1)
                      ! input: liquid water fluxes
                      iLayerLiqFluxSnow,                  & ! intent(in):    liquid flux at the interface of each snow layer (m s-1)
                      iLayerLiqFluxSoil,                  & ! intent(in):    liquid flux at the interface of each soil layer (m s-1)
                      ! input: trial value of model state variables
                      mLayerTempTrial,                    & ! intent(in):    trial temperature at the current iteration (K)
                      mLayerMatricHeadTrial,              & ! intent(in):    trial matric head at the current iteration(m)
                      mLayerVolFracLiqTrial,              & ! intent(in):    trial volumetric fraction of liquid water at the current iteration(-)
                      mLayerVolFracIceTrial,              & ! intent(in):    trial volumetric fraction of ice water at the current iteration(-)
                      ! input: pre-computed derivatives
                      mLayerdTheta_dTk,                   & ! intent(in):    derivative in volumetric liquid water content w.r.t. temperature (K-1)
                      mLayerFracLiqSnow,                  & ! intent(in):    fraction of liquid water (-)
                      ! input-output: data structures
                      mpar_data,                          & ! intent(in):    model parameters
                      indx_data,                          & ! intent(in):    model indices
                      prog_data,                          & ! intent(in):    model prognostic variables for a local HRU
                      diag_data,                          & ! intent(inout):    model diagnostic variables for a local HRU
                      flux_data,                          & ! intent(inout): model fluxes for a local HRU
                      ! output: fluxes and derivatives at all layer interfaces
                      iLayerNrgFlux,                      & ! intent(out):   energy flux at the layer interfaces (W m-2)
                      dFlux_dTempAbove,                   & ! intent(out):   derivatives in the flux w.r.t. temperature in the layer above (W m-2 K-1)
                      dFlux_dTempBelow,                   & ! intent(out):   derivatives in the flux w.r.t. temperature in the layer below (W m-2 K-1)
                      dFlux_dWatAbove,                    & ! intent(out):   derivatives in the flux w.r.t. water state in the layer above (W m-2 K-1)
                      dFlux_dWatBelow,                    & ! intent(out):   derivatives in the flux w.r.t. water state in the layer below (W m-2 K-1)
                    ! output: error control
                      err,message)                          ! intent(out): error control
  ! utility modules
  USE soil_utils_module,only:volFracLiq     ! compute volumetric fraction of liquid water
  USE soil_utils_module,only:matricHead     ! compute the matric head based on volumetric water content
  USE soil_utils_module,only:crit_soilT     ! compute critical temperature below which ice exists
  USE snow_utils_module,only:fracliquid     ! compute fraction of liquid water at a given temperature
  USE soil_utils_module,only:dTheta_dPsi    ! compute derivative of the soil moisture characteristic w.r.t. psi (m-1)
  USE soil_utils_module,only:dPsi_dTheta    ! compute derivative of the soil moisture characteristic w.r.t. theta (m)

  ! constants
  USE multiconst, only: gravity, &                          ! gravitational acceleration (m s-1)
                        Tfreeze, &                          ! freezing point of water (K)
                        iden_water,iden_ice,&      ! intrinsic density of water and ice (kg m-3)
                        LH_fus                              ! latent heat of fusion (J kg-1)
  ! -------------------------------------------------------------------------------------------------------------------------------------------------
  implicit none
  ! input: model control
  logical(lgt),intent(in)            :: scalarSolution              ! flag to denote if implementing the scalar solution
  logical(lgt),intent(in)            :: deriv_desired               ! flag indicating if derivatives are desired
  ! input: fluxes and derivatives at the upper boundary
  real(rkind),intent(in)             :: groundNetFlux               ! net energy flux for the ground surface (W m-2)
  real(rkind),intent(inout)          :: dGroundNetFlux_dGroundTemp  ! derivative in net ground flux w.r.t. ground temperature (W m-2 K-1)
  ! input: liquid water fluxes
  real(rkind),intent(in)             :: iLayerLiqFluxSnow(0:)       ! liquid flux at the interface of each snow layer (m s-1)
  real(rkind),intent(in)             :: iLayerLiqFluxSoil(0:)       ! liquid flux at the interface of each soil layer (m s-1)
  ! input: trial model state variables
  real(rkind),intent(in)              :: mLayerTempTrial(:)         ! temperature in each layer at the current iteration (m)
  real(rkind),intent(in)              :: mLayerMatricHeadTrial(:)   ! matric head in each layer at the current iteration (m)
  real(rkind),intent(in)              :: mLayerVolFracLiqTrial(:)   ! volumetric fraction of liquid at the current iteration (-)
  real(rkind),intent(in)              :: mLayerVolFracIceTrial(:)   ! volumetric fraction of ice at the current iteration (-)
  ! input: pre-computed derivatives
  real(rkind),intent(in)              :: mLayerdTheta_dTk(:)        ! derivative in volumetric liquid water content w.r.t. temperature (K-1)
  real(rkind),intent(in)              :: mLayerFracLiqSnow(:)       ! fraction of liquid water (-)
  ! input-output: data structures
  type(var_dlength),intent(in)        :: mpar_data                  ! model parameters
  type(var_ilength),intent(in)        :: indx_data                  ! state vector geometry
  type(var_dlength),intent(in)        :: prog_data                  ! prognostic variables for a local HRU
  type(var_dlength),intent(inout)        :: diag_data                  ! diagnostic variables for a local HRU
  type(var_dlength),intent(inout)     :: flux_data                  ! model fluxes for a local HRU
  ! output: fluxes and derivatives at all layer interfaces
  real(rkind),intent(out)             :: iLayerNrgFlux(0:)          ! energy flux at the layer interfaces (W m-2)
  real(rkind),intent(out)             :: dFlux_dTempAbove(0:)       ! derivatives in the flux w.r.t. temperature in the layer above (J m-2 s-1 K-1)
  real(rkind),intent(out)             :: dFlux_dTempBelow(0:)       ! derivatives in the flux w.r.t. temperature in the layer below (J m-2 s-1 K-1)
  real(rkind),intent(out)             :: dFlux_dWatAbove(0:)        ! derivatives in the flux w.r.t. water state in the layer above (J m-2 s-1 K-1)
  real(rkind),intent(out)             :: dFlux_dWatBelow(0:)        ! derivatives in the flux w.r.t. water state in the layer below (J m-2 s-1 K-1)
  ! output: error control
  integer(i4b),intent(out)            :: err                        ! error code
  character(*),intent(out)            :: message                    ! error message
  ! ------------------------------------------------------------------------------------------------------------------------------------------------------
  ! local variables
  character(LEN=256)               :: cmessage                     ! error message of downwind routine
  integer(i4b)                     :: i,j,iLayer,iLayer0           ! index of model layers
  integer(i4b)                     :: ixLayerDesired(1)            ! layer desired (scalar solution)
  integer(i4b)                     :: ixTop                        ! top layer in subroutine call
  integer(i4b)                     :: ixBot                        ! bottom layer in subroutine call
  real(rkind)                      :: qFlux                        ! liquid flux at layer interfaces (m s-1)
  real(rkind)                      :: dz                           ! height difference (m)
  ! additional variables to compute numerical derivatives
  integer(i4b)                     :: nFlux                        ! number of flux calculations required (>1 = numerical derivatives with one-sided finite differences)
  integer(i4b)                     :: itry                         ! index of different flux calculations
  integer(i4b),parameter           :: unperturbed=0                ! named variable to identify the case of unperturbed state variables
  integer(i4b),parameter           :: perturbState=1               ! named variable to identify the case where we perturb the state in the current layer
  integer(i4b),parameter           :: perturbStateTempAbove=2      ! named variable to identify the case where we perturb the state layer above
  integer(i4b),parameter           :: perturbStateTempBelow=3      ! named variable to identify the case where we perturb the state layer below
  integer(i4b),parameter           :: perturbStateWatAbove=4       ! named variable to identify the case where we perturb the state layer above
  integer(i4b),parameter           :: perturbStateWatBelow=5       ! named variable to identify the case where we perturb the state layer below
  integer(i4b)                     :: ixPerturb                    ! index of element in 2-element vector to perturb
  integer(i4b)                     :: ixOriginal                   ! index of perturbed element in the original vector
  real(rkind)                      :: scalarThermCFlux               ! thermal conductivity (W m-1 K-1)
  real(rkind)                      :: scalarThermCFlux_dTempAbove    ! thermal conductivity with perturbation to the temperature state above (W m-1 K-1)
  real(rkind)                      :: scalarThermCFlux_dTempBelow    ! thermal conductivity with perturbation to the temperature state below (W m-1 K-1)
  real(rkind)                      :: scalarThermCFlux_dWatAbove     ! thermal conductivity with perturbation to the water state above
  real(rkind)                      :: scalarThermCFlux_dWatBelow     ! thermal conductivity with perturbation to the water state below
  real(rkind)                      :: flux0,flux1,flux2            ! fluxes used to calculate derivatives (W m-2)
  ! compute fluxes and derivatives at layer interfaces
  integer(i4b),dimension(2)        :: mLayer_ind                   ! indices of above and below layers
  integer(i4b),dimension(2)        :: iLayer_ind                   ! indices of above and below interfaces
  real(rkind)                      :: matricFHead                  ! matric head for frozen soil
  real(rkind)                      :: Tcrit                        ! temperature where all water is unfrozen (K)
  real(rkind)                      :: fLiq                         ! fraction of liquid water (-)
  real(rkind),dimension(2)         :: vectorMatricHeadTrial        ! trial value of matric head (m)
  real(rkind),dimension(2)         :: vectorVolFracLiqTrial        ! trial value of volumetric liquid content (-)
  real(rkind),dimension(2)         :: vectorVolFracIceTrial        ! trial value of volumetric ice content (-)
  real(rkind),dimension(2)         :: vectorTempTrial              ! trial value of temperature (K)
  real(rkind),dimension(2)         :: vectordTheta_dPsi            ! derivative in the soil water characteristic w.r.t. psi (m-1)
  real(rkind),dimension(2)         :: vectordPsi_dTheta            ! derivative in the soil water characteristic w.r.t. theta (m)
  real(rkind),dimension(2)         :: vectorFracLiqSnow            ! fraction of liquid water (-)
  real(rkind),dimension(2)         :: vectortheta_sat              ! layer above and below soil porosity (-)
  real(rkind),dimension(2)         :: vectoriden_soil              ! layer above and below density of soil (kg m-3)
  real(rkind),dimension(2)         :: vectorthCond_soil            ! layer above and below thermal conductivity of soil (W m-1 K-1)
  real(rkind),dimension(2)         :: vectorfrac_sand              ! layer above and below fraction of sand (-)
  real(rkind),dimension(2)         :: vectorfrac_clay              ! layer above and below fraction of clay (-)
  ! recompute the perturbed version of iLayerThermalC, this could be the only version and remove the computThermConduct_module
  real(rkind)                      :: dThermalC_dHydStateAbove     ! derivative in the thermal conductivity w.r.t. water state in the layer above
  real(rkind)                      :: dThermalC_dHydStateBelow     ! derivative in the thermal conductivity w.r.t. water state in the layer above
  real(rkind)                      :: dThermalC_dNrgStateAbove     ! derivative in the thermal conductivity w.r.t. energy state in the layer above
  real(rkind)                      :: dThermalC_dNrgStateBelow     ! derivative in the thermal conductivity w.r.t. energy state in the layer above
  ! ------------------------------------------------------------------------------------------------------------------------------------------------------
  ! make association of local variables with information in the data structures
  associate(&
    ixDerivMethod           => model_decisions(iLookDECISIONS%fDerivMeth)%iDecision,      & ! intent(in): method used to calculate flux derivatives
    ix_bcUpprTdyn           => model_decisions(iLookDECISIONS%bcUpprTdyn)%iDecision,      & ! intent(in): method used to calculate the upper boundary condition for thermodynamics
    ix_bcLowrTdyn           => model_decisions(iLookDECISIONS%bcLowrTdyn)%iDecision,      & ! intent(in): method used to calculate the lower boundary condition for thermodynamics
    ixRichards                 => model_decisions(iLookDECISIONS%f_Richards)%iDecision,   & ! intent(in): index of the form of Richards' equation
    ixThCondSnow            => model_decisions(iLookDECISIONS%thCondSnow)%iDecision,      & ! intent(in): choice of method for thermal conductivity of snow
    ixThCondSoil            => model_decisions(iLookDECISIONS%thCondSoil)%iDecision,      & ! intent(in): choice of method for thermal conductivity of soil
    ! input: model coordinates
    nSnow                   => indx_data%var(iLookINDEX%nSnow)%dat(1),               & ! intent(in): number of snow layers
    nLayers                 => indx_data%var(iLookINDEX%nLayers)%dat(1),             & ! intent(in): total number of layers
    layerType               => indx_data%var(iLookINDEX%layerType)%dat,              & ! intent(in): layer type (iname_soil or iname_snow)
    ixLayerState            => indx_data%var(iLookINDEX%ixLayerState)%dat,           & ! intent(in): list of indices for all model layers
    ixSnowSoilNrg           => indx_data%var(iLookINDEX%ixSnowSoilNrg)%dat,          & ! intent(in): index in the state subset for energy state variables in the snow+soil domain
    ! input: thermal properties
    mLayerDepth             => prog_data%var(iLookPROG%mLayerDepth)%dat,             & ! intent(in): depth of each layer (m)
    mLayerHeight            => prog_data%var(iLookPROG%mLayerHeight)%dat,            & ! intent(in): height at the mid-point of each layer (m)
    upperBoundTemp          => mpar_data%var(iLookPARAM%upperBoundTemp)%dat(1),      & ! intent(in): temperature of the upper boundary (K)
    lowerBoundTemp          => mpar_data%var(iLookPARAM%lowerBoundTemp)%dat(1),      & ! intent(in): temperature of the lower boundary (K)
    iLayerHeight            => prog_data%var(iLookPROG%iLayerHeight)%dat,            & ! intent(in): height at the interface of each layer (m)
    fixedThermalCond_snow   => mpar_data%var(iLookPARAM%fixedThermalCond_snow)%dat(1),    & ! intent(in): temporally constant thermal conductivity of snow (W m-1 K-1)
    iLayerThermalC          => diag_data%var(iLookDIAG%iLayerThermalC)%dat,          & ! intent(inout): thermal conductivity at the interface of each layer (W m-1 K-1)
    ! input: depth varying soil parameters
    iden_soil               => mpar_data%var(iLookPARAM%soil_dens_intr)%dat,              & ! intent(in): intrinsic density of soil (kg m-3)
    thCond_soil             => mpar_data%var(iLookPARAM%thCond_soil)%dat,                 & ! intent(in): thermal conductivity of soil (W m-1 K-1)
    theta_sat               => mpar_data%var(iLookPARAM%theta_sat)%dat,                   & ! intent(in): soil porosity (-)
    frac_sand               => mpar_data%var(iLookPARAM%frac_sand)%dat,                   & ! intent(in): fraction of sand (-)
    frac_clay               => mpar_data%var(iLookPARAM%frac_clay)%dat,                   & ! intent(in): fraction of clay (-)
    vGn_m                   => diag_data%var(iLookDIAG%scalarVGn_m)%dat,                  & ! intent(in):  [dp(:)] van Genutchen "m" parameter (-)
    vGn_n                   => mpar_data%var(iLookPARAM%vGn_n)%dat,                       & ! intent(in):  [dp(:)] van Genutchen "n" parameter (-)
    vGn_alpha               => mpar_data%var(iLookPARAM%vGn_alpha)%dat,                   & ! intent(in):  [dp(:)] van Genutchen "alpha" parameter (m-1)
    theta_res               => mpar_data%var(iLookPARAM%theta_res)%dat,                   & ! intent(in):  [dp(:)] soil residual volumetric water content (-)
    ! input: snow parameters
    snowfrz_scale           => mpar_data%var(iLookPARAM%snowfrz_scale)%dat(1),            & ! intent(in):  [dp] scaling parameter for the snow freezing curve (K-1)
    ! output: diagnostic fluxes
    iLayerConductiveFlux => flux_data%var(iLookFLUX%iLayerConductiveFlux)%dat,    & ! intent(out): conductive energy flux at layer interfaces at end of time step (W m-2)
    iLayerAdvectiveFlux  => flux_data%var(iLookFLUX%iLayerAdvectiveFlux)%dat      & ! intent(out): advective energy flux at layer interfaces at end of time step (W m-2)
    )  ! association of local variables with information in the data structures
    ! ------------------------------------------------------------------------------------------------------------------------------------------------------
    ! initialize error control
    err=0; message='ssdNrgFlux/'

    ! set conductive and advective fluxes to missing in the upper boundary
    ! NOTE: advective flux at the upper boundary is included in the ground heat flux
    iLayerConductiveFlux(0) = realMissing
    iLayerAdvectiveFlux(0)  = realMissing

    ! check the need to compute numerical derivatives
    if(ixDerivMethod==numerical)then
      nFlux=5  ! compute the derivatives and cross derivates using one-sided finite differences
    else
      nFlux=0  ! compute analytical derivatives
    end if

    ! get the indices for the snow+soil layers
    if(scalarSolution)then
      ixLayerDesired = pack(ixLayerState, ixSnowSoilNrg/=integerMissing)
      ixTop = ixLayerDesired(1)
      ixBot = ixLayerDesired(1)
    else
      ixTop = 1
      ixBot = nLayers
    endif

    ! -------------------------------------------------------------------------------------------------------------------
    ! ***** compute the derivative in fluxes at layer interfaces w.r.t state in the layer above and the layer below *****
    ! -------------------------------------------------------------------------------------------------------------------

    ! initialize un-used elements
    ! ***** the upper boundary
    dFlux_dTempAbove(0) = 0._rkind ! this will be in canopy
    dFlux_dWatAbove(0) = 0._rkind ! this will be in canopy

    ! ***** the lower boundary
    dFlux_dTempBelow(nLayers) = -huge(lowerBoundTemp)  ! don't expect this to be used, so deliberately set to a ridiculous value to cause problems
    dFlux_dWatBelow(nLayers) = -huge(lowerBoundTemp)  ! don't expect this to be used, so deliberately set to a ridiculous value to cause problems

    ! loop through INTERFACES...
    do iLayer0=ixTop-1,ixBot
      if(ixTop-1/=0 .and. iLayer==ixTop-1)then
       iLayer = 0 !need to always do layer 0
      else
       iLayer = iLayer0
      endif

      ! either one or multiple flux calls, depending on if using analytical or numerical derivatives
      do itry=nFlux,0,-1  ! (work backwards to ensure all computed fluxes come from the un-perturbed case)

          ! =====
          ! determine layer to perturb
          ! ==========================
        select case(itry)
          ! skip undesired perturbations
          case(perturbState); cycle       ! perturbing the layers above and below the flux at the interface
          ! identify the index for the perturbation
          case(unperturbed);       ixPerturb = 0
          case(perturbStateTempAbove)
            if(iLayer==0) cycle ! cannot perturb state above (does not exist) -- so keep cycling
            ixPerturb = 1
          case(perturbStateTempBelow)
            if(iLayer==nLayers) cycle  ! cannot perturb state below (does not exist) -- so keep cycling
            ixPerturb = 2
          case(perturbStateWatAbove)
            if(iLayer==0) cycle ! cannot perturb state above (does not exist) -- so keep cycling
            ixPerturb = 3
          case(perturbStateWatBelow)
            if(iLayer==nLayers) cycle  ! cannot perturb state below (does not exist) -- so keep cycling
            ixPerturb = 4
          case default; err=10; message=trim(message)//"unknown perturbation"; return
        end select ! (identifying layer to of perturbation)
        ! determine the index in the original vector
        ixOriginal = iLayer + (ixPerturb-1)

        ! =====
        ! set indices and parameters needed for layer perturbation
        ! ========================================================
        mLayer_ind(1) = iLayer
        mLayer_ind(2) = iLayer+1
        if (iLayer==0 ) mLayer_ind(1) = 1
        if (iLayer==nLayers ) mLayer_ind(2) = nLayers
        ! indices of interface are different at top layer since interface 0 exists
        iLayer_ind = mLayer_ind
        if (iLayer==0 ) iLayer_ind(1) = 0

        ! =====
        ! get input state variables...
        ! ============================
        ! start with the un-perturbed case
        vectorVolFracLiqTrial(1:2) = mLayerVolFracLiqTrial(mLayer_ind)
        ! need to protect against negative indexes
        if ( mLayer_ind(1) > nSnow .and. mLayer_ind(2) > nSnow)then
         vectorMatricHeadTrial(1:2) = mLayerMatricHeadTrial(mLayer_ind-nSnow)
        else if (mLayer_ind(2) > nSnow)then !mLayer_ind(1) == nSnow
         vectorMatricHeadTrial(1) = realMissing
         vectorMatricHeadTrial(2) = mLayerMatricHeadTrial(mLayer_ind(2)-nSnow)
        else !snow layer
         vectorMatricHeadTrial(1:2) = realMissing
        end if
        vectorTempTrial(1:2) = mLayerTempTrial(mLayer_ind)
        vectorVolFracIceTrial(1:2) = mLayerVolFracIceTrial(mLayer_ind)
        ! make appropriate perturbations,
        if(ixPerturb > 2)then
          vectorMatricHeadTrial(ixPerturb-2) = vectorMatricHeadTrial(ixPerturb-2) + dx
          vectorVolFracLiqTrial(ixPerturb-2) = vectorVolFracLiqTrial(ixPerturb-2) + dx
        else if(ixPerturb > 0)then
          vectorTempTrial(ixPerturb) = vectorTempTrial(ixPerturb) + dx
        endif

        ! *****
        ! * compute the volumetric fraction of liquid, ice, and air in each layer in response to perturbation ...
        ! *******************************************************************************************************

        do i = 1,2 !(layer above and below)
          select case(layerType(mLayer_ind(i))) !(snow or soil)

            case(iname_soil)
              j = mLayer_ind(i)-nSnow !soil layer
              if(ixPerturb > 0)then ! only recompute these if perturbed
                select case(ixRichards)  ! (form of Richards' equation)
                  case(moisture) !
                    vectorMatricHeadTrial(i) = matricHead(vectorVolFracLiqTrial(i),vGn_alpha(j),theta_res(j),theta_sat(j),vGn_n(j),vGn_m(j))
                    Tcrit = crit_soilT(vectorMatricHeadTrial(i))
                    !if change temp and below critical, it changes the state variable, seems like a problem FIX
                    if(vectorTempTrial(i) < Tcrit) then !if do not perturb temperature, this should not change
                      matricFHead = (LH_fus/gravity)*(vectorTempTrial(i) - Tfreeze)/Tfreeze
                      vectorVolFracLiqTrial(i) = volFracLiq(matricFHead,vGn_alpha(j),theta_res(j),theta_sat(j),vGn_n(j),vGn_m(j))
                    endif
                  case(mixdform)
                    Tcrit = crit_soilT(vectorMatricHeadTrial(i))
                    if(vectorTempTrial(i) < Tcrit) then !if do not perturb temperature, this should not change, but matricHeadTrial will have changed
                      matricFHead = (LH_fus/gravity)*(vectorTempTrial(i) - Tfreeze)/Tfreeze
                      vectorVolFracLiqTrial(i) = volFracLiq(matricFHead,vGn_alpha(j),theta_res(j),theta_sat(j),vGn_n(j),vGn_m(j))
                    else
                      vectorVolFracLiqTrial(i) = volFracLiq(vectorMatricHeadTrial(i),vGn_alpha(j),theta_res(j),theta_sat(j),vGn_n(j),vGn_m(j))
                    endif
                end select ! (form of Richards' equation)
                vectorVolFracIceTrial(i) = volFracLiq(vectorMatricHeadTrial(i),vGn_alpha(j),theta_res(j),theta_sat(j),vGn_n(j),vGn_m(j)) - vectorVolFracLiqTrial(i)
              endif ! (recompute if perturbed)
              ! derivatives, these need to be computed because they are computed only in soilLiqFlx which is called after this
              vectordPsi_dTheta(i) = dPsi_dTheta(vectorVolFracLiqTrial(i),vGn_alpha(j),theta_res(j),theta_sat(j),vGn_n(j),vGn_m(j))
              vectordTheta_dPsi(i) = dTheta_dPsi(vectorMatricHeadTrial(i),vGn_alpha(j),theta_res(j),theta_sat(j),vGn_n(j),vGn_m(j))
              vectorFracLiqSnow(i) = realMissing
              ! soil parameters
              vectortheta_sat(i) = theta_sat(j)
              vectoriden_soil(i) = iden_soil(j)
              vectorthCond_soil(i) = thCond_soil(j)
              vectorfrac_sand(i) = frac_sand(j)
              vectorfrac_clay(i) = frac_clay(j)

            case(iname_snow)
              fLiq = fracliquid(vectorTempTrial(i),snowfrz_scale) ! fraction of liquid water
              if(ixPerturb > 0) vectorVolFracIceTrial(i) = ( vectorVolFracLiqTrial(i) / fLiq - vectorVolFracLiqTrial(i) )*(iden_water/iden_ice) ! use perturbed nodeVolTotWatTrial
              ! derivatives
              vectordPsi_dTheta(i) = realMissing
              vectordTheta_dPsi(i) = realMissing
              vectorFracLiqSnow(i) = mLayerFracLiqSnow(mLayer_ind(i))
              ! soil parameters do not exist
              vectortheta_sat(i) = realMissing
              vectoriden_soil(i) = realMissing
              vectorthCond_soil(i) = realMissing
              vectorfrac_sand(i) = realMissing
              vectorfrac_clay(i) = realMissing

            case default; err=20; message=trim(message)//'unable to identify type of layer (snow or soil) to compute volumetric fraction of air'; return

          end select !(snow or soil)
        enddo !(layer above and below)

        ! =====
        ! get thermal conductivity at layer interface and its derivative w.r.t. the state above and the state below...
        ! ============================================================================================================
        call iLayerThermalConduct(&
                            ! input: model control
                            deriv_desired,             & ! intent(in): flag indicating if derivatives are desired
                            ixRichards,                         & ! intent(in): index defining the form of Richards' equation (moisture or mixdform)
                            ixThCondSnow,                        & ! intent(in): choice of method for thermal conductivity of snow
                            ixThCondSoil,                        & ! intent(in): choice of method for thermal conductivity of soil
                            ! input: coordinate variables
                            nLayers,                             & ! intent(in): number of layers
                            iLayer,                              & ! intent(in): layer index for output
                            layerType(mLayer_ind),               & ! intent(in): layer type (iname_soil or iname_snow)
                            ! input: state variables (adjacent layers)
                            vectorMatricHeadTrial,               & ! intent(in): matric head at the nodes (m)
                            vectorVolFracLiqTrial,               & ! intent(in): volumetric liquid water at the nodes (m)
                            vectorVolFracIceTrial,               & ! intent(in): volumetric ice at the nodes (m)
                            vectorTempTrial,                     & ! intent(in): temperature at the nodes (m)
                            ! input: pre-computed derivatives
                            mLayerdTheta_dTk(mLayer_ind),        & ! intent(in): derivative in volumetric liquid water content w.r.t. temperature (K-1)
                            vectorFracLiqSnow,                   & ! intent(in): fraction of liquid water (-)
                            vectordTheta_dPsi,                   & ! intent(in): derivative in the soil water characteristic w.r.t. psi (m-1)
                            vectordPsi_dTheta,                   & ! intent(in): derivative in the soil water characteristic w.r.t. theta (m)
                            ! input: model coordinate variables (adjacent layers)
                            mLayerHeight(mLayer_ind),            & ! intent(in): height at the mid-point of the node (m)
                            iLayerHeight(iLayer_ind),            & ! intent(in): height at the interface of the nodes (m)
                            ! input: soil parameters
                            vectortheta_sat,                     & ! intent(in): soil porosity (-)
                            vectoriden_soil,                     & ! intent(in): intrinsic density of soil (kg m-3)
                            vectorthCond_soil,                   & ! intent(in): thermal conductivity of soil (W m-1 K-1)
                            vectorfrac_sand,                     & ! intent(in): fraction of sand (-)
                            vectorfrac_clay,                     & ! intent(in): fraction of clay (-)
                            ! input: snow parameters
                            fixedThermalCond_snow,               & ! intent(in): temporally constant thermal conductivity of snow (W m-1 K-1)
                            ! output: conductivity at the layer interface (scalars)
                            iLayerThermalC(iLayer),       & ! intent(inout): thermal conductivity at the interface of each layer (W m-1 K-1)
                            ! output: derivatives in thermal conductivity w.r.t. state variables -- matric head or volumetric lquid water -- in the layer above and layer below
                            dThermalC_dHydStateAbove,    & ! intent(out): derivative in the thermal conductivity w.r.t. water state in the layer above
                            dThermalC_dHydStateBelow,    & ! intent(out): derivative in the thermal conductivity w.r.t. water state in the layer above
                            ! output: derivatives in thermal conductivity w.r.t. energy state variables -- now just temperature -- in the layer above and layer below (W m-1 K-2)
                            dThermalC_dNrgStateAbove,    & ! intent(out): derivative in the thermal conductivity w.r.t. energy state in the layer above
                            dThermalC_dNrgStateBelow,    & ! intent(out): derivative in the thermal conductivity w.r.t. energy state in the layer above
                            ! output: error control
                            err,cmessage)               ! intent(out): error control

        if(err/=0)then; message=trim(message)//trim(cmessage); return; end if

        ! compute total vertical flux, to compute derivatives
        if(deriv_desired .and. ixDerivMethod==numerical)then
          select case(itry)
            case(unperturbed);           scalarThermCFlux            = iLayerThermalC(iLayer)
            case(perturbStateTempAbove); scalarThermCFlux_dTempAbove = iLayerThermalC(iLayer)
            case(perturbStateTempBelow); scalarThermCFlux_dTempBelow = iLayerThermalC(iLayer)
            case(perturbStateWatAbove);  scalarThermCFlux_dWatAbove  = iLayerThermalC(iLayer)
            case(perturbStateWatBelow);  scalarThermCFlux_dWatBelow  = iLayerThermalC(iLayer)
            case default; err=10; message=trim(message)//"unknown perturbation"; return
          end select
        end if

      end do  ! (looping through different flux calculations -- one or multiple calls depending if desire for numerical or analytical derivatives)


      ! ***** the upper boundary
      if(iLayer==0)then  ! (upper boundary)

        ! identify the upper boundary condition
        select case(ix_bcUpprTdyn)

          ! * prescribed temperature at the upper boundary
          case(prescribedTemp)
            dz = mLayerHeight(iLayer+1)*0.5_rkind
            if(ixDerivMethod==analytical)then    ! ** analytical derivatives
              dFlux_dWatBelow(iLayer)  = -dThermalC_dHydStateBelow * ( mLayerTempTrial(iLayer+1) - upperBoundTemp )/dz
              dFlux_dTempBelow(iLayer) = -dThermalC_dNrgStateBelow * ( mLayerTempTrial(iLayer+1) - upperBoundTemp )/dz - iLayerThermalC(iLayer)/dz
            else                              ! ** numerical derivatives
              flux0 = -scalarThermCFlux          *( mLayerTempTrial(iLayer+1) - upperBoundTemp ) / dz
              flux2 = -scalarThermCFlux_dWatBelow*( mLayerTempTrial(iLayer+1) - upperBoundTemp ) / dz
              dFlux_dWatBelow(iLayer) = (flux2 - flux0)/dx
              flux0 = -scalarThermCFlux           *( mLayerTempTrial(iLayer+1)     - upperBoundTemp ) / dz
              flux2 = -scalarThermCFlux_dTempBelow*((mLayerTempTrial(iLayer+1)+dx) - upperBoundTemp ) / dz
              dFlux_dTempBelow(iLayer) = (flux2 - flux0)/dx
            end if

          ! * zero flux at the upper boundary
          case(zeroFlux)
            dFlux_dWatBelow(iLayer) = 0._rkind
            dFlux_dTempBelow(iLayer) = 0._rkind

          ! * compute flux inside vegetation energy flux routine, use here
          case(energyFlux)
            dFlux_dWatBelow(iLayer) = 0._rkind !dGroundNetFlux_dGroundWat, does not exist in vegNrgFlux
            dFlux_dTempBelow(iLayer) = dGroundNetFlux_dGroundTemp

          case default; err=20; message=trim(message)//'unable to identify upper boundary condition for thermodynamics'; return

        end select  ! (identifying the upper boundary condition for thermodynamics)
        !dGroundNetFlux_dGroundWat  = dFlux_dWatBelow(iLayer) ! this is true, but since not used in vegNrgFlux do not define
        dGroundNetFlux_dGroundTemp = dFlux_dTempBelow(iLayer) ! need this in vegNrgFlux

      ! ***** the lower boundary
      else if(iLayer==nLayers)then  ! (lower boundary)

        ! identify the lower boundary condition
        select case(ix_bcLowrTdyn)

          ! * prescribed temperature at the lower boundary
          case(prescribedTemp)
          dz = mLayerDepth(iLayer)*0.5_rkind
          if(ixDerivMethod==analytical)then    ! ** analytical derivatives
            dFlux_dWatAbove(iLayer)  = -dThermalC_dHydStateAbove * ( lowerBoundTemp - mLayerTempTrial(iLayer) )/dz
            dFlux_dTempAbove(iLayer) = -dThermalC_dNrgStateAbove * ( lowerBoundTemp - mLayerTempTrial(iLayer) )/dz + iLayerThermalC(iLayer)/dz
          else                              ! ** numerical derivatives
            flux0 = -scalarThermCFlux           * ( lowerBoundTemp - mLayerTempTrial(iLayer) )/dz
            flux1 = -scalarThermCFlux_dWatAbove * ( lowerBoundTemp - mLayerTempTrial(iLayer) )/dz
            dFlux_dWatAbove(iLayer) = (flux1 - flux0)/dx
            flux0 = -scalarThermCFlux            * ( lowerBoundTemp -  mLayerTempTrial(iLayer)     )/dz
            flux1 = -scalarThermCFlux_dTempAbove * ( lowerBoundTemp - (mLayerTempTrial(iLayer)+dx) )/dz
            dFlux_dTempAbove(iLayer) = (flux1 - flux0)/dx
          end if

          ! * zero flux at the lower boundary
          case(zeroFlux)
          dFlux_dWatAbove(iLayer) = 0._rkind
          dFlux_dTempAbove(iLayer) = 0._rkind

          case default; err=20; message=trim(message)//'unable to identify lower boundary condition for thermodynamics'; return

        end select  ! (identifying the lower boundary condition for thermodynamics)

        ! ***** internal layers

      else
        dz = (mLayerHeight(iLayer+1) - mLayerHeight(iLayer))
        if(ixDerivMethod==analytical)then    ! ** analytical derivatives
          dFlux_dWatAbove(iLayer)  = -dThermalC_dHydStateAbove * ( mLayerTempTrial(iLayer+1) - mLayerTempTrial(iLayer) )/dz
          dFlux_dWatBelow(iLayer)  = -dThermalC_dHydStateBelow * ( mLayerTempTrial(iLayer+1) - mLayerTempTrial(iLayer) )/dz
          dFlux_dTempAbove(iLayer) = -dThermalC_dNrgStateAbove * ( mLayerTempTrial(iLayer+1) - mLayerTempTrial(iLayer) )/dz + iLayerThermalC(iLayer)/dz
          dFlux_dTempBelow(iLayer) = -dThermalC_dNrgStateBelow * ( mLayerTempTrial(iLayer+1) - mLayerTempTrial(iLayer) )/dz - iLayerThermalC(iLayer)/dz
        else                              ! ** numerical derivatives
          flux0 = -scalarThermCFlux          *( mLayerTempTrial(iLayer+1) - mLayerTempTrial(iLayer) ) / dz
          flux1 = -scalarThermCFlux_dWatAbove*( mLayerTempTrial(iLayer+1) - mLayerTempTrial(iLayer) ) / dz
          flux2 = -scalarThermCFlux_dWatBelow*( mLayerTempTrial(iLayer+1) - mLayerTempTrial(iLayer) ) / dz
          dFlux_dWatAbove(iLayer) = (flux1 - flux0)/dx
          dFlux_dWatBelow(iLayer) = (flux2 - flux0)/dx
          flux0 = -scalarThermCFlux           *( mLayerTempTrial(iLayer+1)     -  mLayerTempTrial(iLayer)    ) / dz
          flux1 = -scalarThermCFlux_dTempAbove*( mLayerTempTrial(iLayer+1)     - (mLayerTempTrial(iLayer)+dx)) / dz
          flux2 = -scalarThermCFlux_dTempBelow*((mLayerTempTrial(iLayer+1)+dx) -  mLayerTempTrial(iLayer)    ) / dz
          dFlux_dTempAbove(iLayer) = (flux1 - flux0)/dx
          dFlux_dTempBelow(iLayer) = (flux2 - flux0)/dx
        end if

      end if  ! type of layer (upper, internal, or lower)

    end do  ! (looping through layers)

    ! -------------------------------------------------------------------------------------------------------------------------
    ! ***** compute the conductive fluxes at layer interfaces *****
    ! Compute flux after the derivatives, because need iLayerThermal as calculated above
    ! -------------------------------------------------------------------------------------------------------------------------
    do iLayer0=ixTop-1,ixBot
      if(ixTop-1/=0 .and. iLayer==ixTop-1)then
       iLayer = 0 !need to always do layer 0
      else
       iLayer = iLayer0
      endif

      if(iLayer==0)then  ! (upper boundary fluxes -- positive downwards)
        iLayerConductiveFlux(iLayer) = groundNetFlux !from vegNrgFlux module

      else if(iLayer==nLayers)then ! (lower boundary fluxes -- positive downwards)
      ! flux depends on the type of lower boundary condition
        select case(ix_bcLowrTdyn) ! (identify the lower boundary condition for thermodynamics
          case(prescribedTemp); iLayerConductiveFlux(iLayer) = -iLayerThermalC(iLayer)*(lowerBoundTemp - mLayerTempTrial(iLayer))/(mLayerDepth(iLayer)*0.5_rkind)
          case(zeroFlux);       iLayerConductiveFlux(iLayer) = 0._rkind
        end select  ! (identifying the lower boundary condition for thermodynamics)

      else ! (domain boundary fluxes -- positive downwards)
        iLayerConductiveFlux(iLayer)  = -iLayerThermalC(iLayer)*(mLayerTempTrial(iLayer+1) - mLayerTempTrial(iLayer)) / &
                                        (mLayerHeight(iLayer+1) - mLayerHeight(iLayer))

        !write(*,'(a,i4,1x,2(f9.3,1x))') 'iLayer, iLayerConductiveFlux(iLayer), iLayerThermalC(iLayer) = ', iLayer, iLayerConductiveFlux(iLayer), iLayerThermalC(iLayer)
      end if ! (the type of layer)
    end do  ! looping through layers

    ! -------------------------------------------------------------------------------------------------------------------------
    ! ***** compute the advective fluxes at layer interfaces *****
    ! -------------------------------------------------------------------------------------------------------------------------
    do iLayer0=ixTop-1,ixBot
      if(ixTop-1/=0 .and. iLayer==ixTop-1)then
       iLayer = 0 !need to always do layer 0
      else
       iLayer = iLayer0
      endif

      if (iLayer==0) then
        iLayerAdvectiveFlux(iLayer) = realMissing !advective flux at the upper boundary is included in the ground heat flux
      else ! get the liquid flux at layer interfaces
        select case(layerType(iLayer))
          case(iname_snow); qFlux = iLayerLiqFluxSnow(iLayer)
          case(iname_soil); qFlux = iLayerLiqFluxSoil(iLayer-nSnow)
          case default; err=20; message=trim(message)//'unable to identify layer type'; return
        end select
        ! compute fluxes at the lower boundary -- positive downwards
        if(iLayer==nLayers)then
          iLayerAdvectiveFlux(iLayer) = -Cp_water*iden_water*qFlux*(lowerBoundTemp - mLayerTempTrial(iLayer))
        ! compute fluxes within the domain -- positive downwards
        else
          iLayerAdvectiveFlux(iLayer) = -Cp_water*iden_water*qFlux*(mLayerTempTrial(iLayer+1) - mLayerTempTrial(iLayer))
        end if
      end if ! (all layers except surface)
    end do  ! looping through layers

    ! -------------------------------------------------------------------------------------------------------------------------
    ! ***** compute the total fluxes at layer interfaces *****
    ! -------------------------------------------------------------------------------------------------------------------------
    ! NOTE: ignore advective fluxes for now
    iLayerNrgFlux(0)           = iLayerConductiveFlux(0)
    iLayerNrgFlux(ixTop:ixBot) = iLayerConductiveFlux(ixTop:ixBot)

  ! end association of local variables with information in the data structures
  end associate

end subroutine ssdNrgFlux


! *************************************************************************************************************************************
! private subroutine iLayerThermalConduct: compute diagnostic energy variables (thermal conductivity and heat capacity) and derivatives
! *************************************************************************************************************************************
subroutine iLayerThermalConduct(&
                      ! input: model control
                      deriv_desired,             & ! intent(in): flag indicating if derivatives are desired
                      ixRichards,                & ! intent(in):    choice of option for Richards' equation
                      ixThCondSnow,              & ! intent(in): choice of method for thermal conductivity of snow
                      ixThCondSoil,              & ! intent(in): choice of method for thermal conductivity of soil
                      ! input: coordinate variables
                      nLayers,                   & ! intent(in): number of layers
                      ixLayerDesired,            & ! intent(in): layer index for output
                      layerType,                 & ! intent(in): layer type (iname_soil or iname_snow)
                      ! input: state variables (adjacent layers)
                      nodeMatricHead,            & ! intent(in): matric head at the nodes (m)
                      nodeVolFracLiq,            & ! intent(in): volumetric liquid water content at the nodes (m)
                      nodeVolFracIce,            & ! intent(in): volumetric ice at the nodes (m)
                      nodeTemp,                  & ! intent(in): temperature at the nodes (m)
                      ! input: pre-computed derivatives
                      dTheta_dTk,                & ! intent(in): derivative in volumetric liquid water content w.r.t. temperature (K-1)
                      fracLiqSnow,               & ! intent(in)  : fraction of liquid water (-)
                      dTheta_dPsi,               & ! intent(in): derivative in the soil water characteristic w.r.t. psi (m-1)
                      dPsi_dTheta,               & ! intent(in): derivative in the soil water characteristic w.r.t. theta (m)
                      ! input: model coordinate variables (adjacent layers)
                      nodeHeight,                & ! intent(in): height at the mid-point of the node (m)
                      node_iHeight,              & ! intent(in): height at the interface of the nodes (m)
                      ! input: soil parameters at nodes
                      theta_sat,                 & ! intent(in): soil porosity (-)
                      iden_soil,                 & !intrinsic density of soil (kg m-3)
                      thCond_soil,               & ! thermal conductivity of soil (W m-1 K-1)
                      frac_sand,                 & ! intent(in): fraction of sand (-)
                      frac_clay,                 & ! fraction of clay (-)
                      ! input: snow parameters
                      fixedThermalCond_snow,     & ! intent(in): temporally constant thermal conductivity of snow (W m-1 K-1)
                      ! output: conductivity at the layer interface (scalars)
                      iLayerThermalC,            & ! intent(inout) thermal conductivity at the interface of each layer (W m-1 K-1)
                      ! output: derivatives in thermal conductivity w.r.t. state variables -- matric head or volumetric lquid water -- in the layer above and layer below
                      dThermalC_dHydStateAbove,  & ! intent(out): derivative in the thermal conductivity w.r.t. water state in the layer above
                      dThermalC_dHydStateBelow,  & ! intent(out): derivative in the thermal conductivity w.r.t. water state in the layer above
                    ! output: derivatives in thermal conductivity w.r.t. energy state variables -- now just temperature -- in the layer above and layer below (W m-1 K-2)
                      dThermalC_dNrgStateAbove,  & ! intent(out): derivative in the thermal conductivity w.r.t. energy state in the layer above
                      dThermalC_dNrgStateBelow,  & ! intent(out): derivative in the thermal conductivity w.r.t. energy state in the layer above
                      ! output: error control
                      err,message)                 ! intent(out): error control

  USE snow_utils_module,only:tcond_snow     ! compute thermal conductivity of snow
  USE soil_utils_module,only:crit_soilT     ! compute critical temperature below which ice exists
  ! constants
  USE multiconst, only: gravity, &                          ! gravitational acceleration (m s-1)
                        Tfreeze, &                          ! freezing point of water (K)
                        iden_water,iden_ice,&               ! intrinsic density of water and ice (kg m-3)
                        LH_fus                              ! latent heat of fusion (J kg-1)

  implicit none
  ! --------------------------------------------------------------------------------------------------------------------------------------
  ! input: model control
  logical(lgt),intent(in)          :: deriv_desired             ! flag to indicate if derivatives are desired
  integer(i4b),intent(in)          :: ixRichards                ! choice of option for Richards' equation
  integer(i4b),intent(in)          :: ixThCondSnow              ! choice of method for thermal conductivity of snow
  integer(i4b),intent(in)          :: ixThCondSoil              ! choice of method for thermal conductivity of soil
  ! input: coordinate variables
  integer(i4b),intent(in)          :: nLayers                   ! number of layers
  integer(i4b),intent(in)          :: ixLayerDesired            ! layer index for output
  integer(i4b),intent(in)          :: layerType(:)              ! layer type (iname_soil or iname_snow)
  ! input: state variables
  real(rkind),intent(in)           :: nodeMatricHead(:)         ! trial vector of total water matric potential (m)
  real(rkind),intent(in)           :: nodeVolFracLiq(:)         ! trial vector of volumetric liquid water content, recomputed with perturbed water state(-)
  real(rkind),intent(in)           :: nodeVolFracIce(:)         ! trial vector of ice content, recomputed with perturbed water state(-)
  real(rkind),intent(in)           :: nodeTemp(:)               ! trial vector of temperature (K)
  ! input: pre-computed derivatives
  real(rkind),intent(in)           :: dTheta_dTk(:)             ! derivative in volumetric liquid water content w.r.t. temperature (K-1)
  real(rkind),intent(in)           :: fracLiqSnow(:)            ! fraction of liquid water (-)
  real(rkind),intent(in)           :: dTheta_dPsi(:)            ! derivative in the soil water characteristic w.r.t. psi (m-1)
  real(rkind),intent(in)           :: dPsi_dTheta(:)            ! derivative in the soil water characteristic w.r.t. theta (m)
  ! input: model coordinate variables
  real(rkind),intent(in)           :: nodeHeight(:)             ! height at the mid-point of the lower node (m)
  real(rkind),intent(in)           :: node_iHeight(:)           ! height at the interface of each node (m)
  ! input: soil parameters
  real(rkind),intent(in)           :: theta_sat(:)              ! soil porosity (-)
  real(rkind),intent(in)           :: iden_soil(:)              ! intrinsic density of soil (kg m-3)
  real(rkind),intent(in)           :: thCond_soil(:)            ! thermal conductivity of soil (W m-1 K-1)
  real(rkind),intent(in)           :: frac_sand(:)              ! intent(in): fraction of sand (-)
  real(rkind),intent(in)           :: frac_clay(:)              ! fraction of clay (-)
  ! input: snow parameters
  real(rkind),intent(in)           :: fixedThermalCond_snow     ! intent(in): temporally constant thermal conductivity of snow (W m-1 K-1)
  ! output: thermal conductivity at layer interfaces
  real(rkind),intent(inout)        :: iLayerThermalC            ! thermal conductivity at the interface of each layer (W m-1 K-1)
  ! output: thermal conductivity derivatives at all layer interfaces
  real(rkind),intent(out)          :: dThermalC_dHydStateAbove  ! derivatives in the thermal conductivity w.r.t. matric head or volumetric liquid water in the layer above
  real(rkind),intent(out)          :: dThermalC_dHydStateBelow  ! derivatives in the thermal conductivity w.r.t. matric head or volumetric liquid water in the layer below
  real(rkind),intent(out)          :: dThermalC_dNrgStateAbove  ! derivatives in the thermal conductivity w.r.t. temperature in the layer above (W m-1 K-2)
  real(rkind),intent(out)          :: dThermalC_dNrgStateBelow  ! derivatives in the thermal conductivity w.r.t. temperature in the layer below (W m-1 K-2)
  ! output: error control
  integer(i4b),intent(out)         :: err                       ! error code
  character(*),intent(out)         :: message                   ! error message
  ! --------------------------------------------------------------------------------------------------------------------------------
  ! local variables (named variables to provide index of 2-element vectors)
  integer(i4b),parameter           :: ixUpper=1                 ! index of upper node in the 2-element vectors
  integer(i4b),parameter           :: ixLower=2                 ! index of lower node in the 2-element vectors
  character(LEN=256)               :: cmessage                  ! error message of downwind routine
  integer(i4b)                     :: iLayer                    ! index of model layer
  real(rkind)                      :: TCn                       ! thermal conductivity below the layer interface (W m-1 K-1)
  real(rkind)                      :: TCp                       ! thermal conductivity above the layer interface (W m-1 K-1)
  real(rkind)                      :: zdn                       ! height difference between interface and lower value (m)
  real(rkind)                      :: zdp                       ! height difference between interface and upper value (m)
  real(rkind)                      :: bulkden_soil              ! bulk density of soil (kg m-3)
  real(rkind)                      :: lambda_drysoil            ! thermal conductivity of dry soil (W m-1)
  real(rkind)                      :: lambda_wetsoil            ! thermal conductivity of wet soil (W m-1)
  real(rkind)                      :: lambda_wet                ! thermal conductivity of the wet material
  real(rkind)                      :: relativeSat               ! relative saturation (-)
  real(rkind)                      :: kerstenNum                ! the Kersten number (-), defining weight applied to conductivity of the wet medium
  real(rkind)                      :: den                       ! denominator in the thermal conductivity calculations
  real(rkind)                      :: dThermalC_dWat(2)   ! derivative in thermal conductivity w.r.t. matric head or volumetric liquid water
  real(rkind)                      :: dThermalC_dNrg(2)   ! derivative in thermal conductivity w.r.t. temperature
  real(rkind)                      :: Tcrit                     ! temperature where all water is unfrozen (K)
  real(rkind)                      :: dlambda_wet_dWat          ! derivative in thermal conductivity of wet material w.r.t.soil water state variable
  real(rkind)                      :: dlambda_wet_dTk           ! derivative in thermal conductivity of wet material w.r.t. temperature
  real(rkind)                      :: dkerstenNum_dWat          ! derivative in Kersten number w.r.t. soil water state variable
  real(rkind)                      :: mLayerThermalC(2)         ! thermal conductivity of each layer (W m-1 K-1)
  real(rkind)                      :: dVolFracLiq_dWat    ! derivative in vol fraction of liquid w.r.t. water state variable
  real(rkind)                      :: dVolFracIce_dWat    ! derivative in vol fraction of ice w.r.t. water state variable
  real(rkind)                      :: dVolFracLiq_dTk     ! derivative in vol fraction of liquid w.r.t. temperature
  real(rkind)                      :: dVolFracIce_dTk     ! derivative in vol fraction of ice w.r.t. temperature
  ! local variables to reproduce the thermal conductivity of Hansson et al. VZJ 2005
  real(rkind),parameter            :: c1=0.55_rkind             ! optimized parameter from Hansson et al. VZJ 2005 (W m-1 K-1)
  real(rkind),parameter            :: c2=0.8_rkind              ! optimized parameter from Hansson et al. VZJ 2005 (W m-1 K-1)
  real(rkind),parameter            :: c3=3.07_rkind             ! optimized parameter from Hansson et al. VZJ 2005 (-)
  real(rkind),parameter            :: c4=0.13_rkind             ! optimized parameter from Hansson et al. VZJ 2005 (W m-1 K-1)
  real(rkind),parameter            :: c5=4._rkind               ! optimized parameter from Hansson et al. VZJ 2005 (-)
  real(rkind),parameter            :: f1=13.05_rkind            ! optimized parameter from Hansson et al. VZJ 2005 (-)
  real(rkind),parameter            :: f2=1.06_rkind             ! optimized parameter from Hansson et al. VZJ 2005 (-)
  real(rkind)                      :: fArg,xArg                 ! temporary variables (see Hansson et al. VZJ 2005 for details)
  real(rkind)                      :: dxArg_dWat,dxArg_dTk      ! derivates of the temporary variables with respect to soil water state variable and temperature

  ! --------------------------------------------------------------------------------------------------------------------------------
  ! initialize error control
  err=0; message="iLayerThermalConduct/"

  ! loop through layers
  do iLayer=ixUpper,ixLower

    ! compute the thermal conductivity of dry and wet soils (W m-1)
    ! NOTE: this is actually constant over the simulation, and included here for clarity
    if(ixThCondSoil==funcSoilWet .and. layerType(iLayer)==iname_soil)then
    bulkden_soil   = iden_soil(iLayer)*( 1._rkind - theta_sat(iLayer) )
    lambda_drysoil = (0.135_rkind*bulkden_soil + 64.7_rkind) / (iden_soil(iLayer) - 0.947_rkind*bulkden_soil)
    lambda_wetsoil = (8.80_rkind*frac_sand(iLayer) + 2.92_rkind*frac_clay(iLayer)) / (frac_sand(iLayer) + frac_clay(iLayer))
    end if

    ! *****
    ! * compute the thermal conductivity of snow and soil and derivates at the mid-point of each layer...
    ! ***************************************************************************************************
    dThermalC_dWat(iLayer) = 0._rkind
    dThermalC_dNrg(iLayer) = 0._rkind

    select case(layerType(iLayer))

      ! ***** soil
      case(iname_soil)

        ! (process derivatives)
        dVolFracLiq_dWat = 0._rkind
        dVolFracIce_dWat = 0._rkind
        dVolFracLiq_dTk  = 0._rkind
        dVolFracIce_dTk  = 0._rkind
        if(deriv_desired)then
          select case(ixRichards)  ! (form of Richards' equation)
            case(moisture)
              dVolFracLiq_dWat = 1._rkind
              dVolFracIce_dWat = dPsi_dTheta(iLayer) - 1._rkind
            case(mixdform)
              Tcrit = crit_soilT(nodeMatricHead(iLayer) )
              if(nodeTemp(iLayer) < Tcrit) then
                dVolFracLiq_dWat = 0._rkind
                dVolFracIce_dWat = dTheta_dPsi(iLayer)
              else
                dVolFracLiq_dWat = dTheta_dPsi(iLayer)
                dVolFracIce_dWat = 0._rkind
              endif
          end select
          dVolFracLiq_dTk = dTheta_dTk(iLayer) !already zeroed out if not below critical temperature
          dVolFracIce_dTk = -dVolFracLiq_dTk !often can and will simplify one of these terms out
        endif

        ! select option for thermal conductivity of soil
        select case(ixThCondSoil)

          ! ** function of soil wetness
          case(funcSoilWet)

            ! compute the thermal conductivity of the wet material (W m-1)
            lambda_wet  = lambda_wetsoil**( 1._rkind - theta_sat(iLayer) ) * lambda_water**theta_sat(iLayer) * lambda_ice**(theta_sat(iLayer) - nodeVolFracLiq(iLayer))
            dlambda_wet_dWat = -lambda_wet * log(lambda_ice) * dVolFracLiq_dWat
            dlambda_wet_dTk  = -lambda_wet * log(lambda_ice) * dVolFracLiq_dTk

            relativeSat = (nodeVolFracIce(iLayer) + nodeVolFracLiq(iLayer))/theta_sat(iLayer)  ! relative saturation
            ! drelativeSat_dWat = dPsi0_dWat/theta_sat(iLayer), and drelativeSat_dTk = 0 (so dkerstenNum_dTk = 0)
            ! compute the Kersten number (-)
            if(relativeSat > 0.1_rkind)then ! log10(0.1) = -1
              kerstenNum = log10(relativeSat) + 1._rkind
              dkerstenNum_dWat = (dVolFracIce_dWat + dVolFracLiq_dWat) / ( theta_sat(iLayer) * relativeSat * log(10._rkind) )
            else
              kerstenNum = 0._rkind  ! dry thermal conductivity
              dkerstenNum_dWat = 0._rkind
            endif
            ! ...and, compute the thermal conductivity
            mLayerThermalC(iLayer) = kerstenNum*lambda_wet + (1._rkind - kerstenNum)*lambda_drysoil

            ! compute derivatives
            dThermalC_dWat(iLayer) = dkerstenNum_dWat * ( lambda_wet - lambda_drysoil ) + kerstenNum*dlambda_wet_dWat
            dThermalC_dNrg(iLayer) = kerstenNum*dlambda_wet_dTk

          ! ** mixture of constituents
          case(mixConstit)
            mLayerThermalC(iLayer) = thCond_soil(iLayer) * ( 1._rkind - theta_sat(iLayer) ) + & ! soil component
                                    lambda_ice         * nodeVolFracIce(iLayer)     + & ! ice component
                                    lambda_water       * nodeVolFracLiq(iLayer)     + & ! liquid water component
                                    lambda_air         * ( theta_sat(iLayer) - (nodeVolFracIce(iLayer) + nodeVolFracLiq(iLayer)) ) ! air component
            ! compute derivatives
            dThermalC_dWat(iLayer) = lambda_ice*dVolFracIce_dWat + lambda_water*dVolFracLiq_dWat + lambda_air*(-dVolFracIce_dWat - dVolFracLiq_dWat)
            dThermalC_dNrg(iLayer) = (lambda_ice - lambda_water) * dVolFracIce_dTk

          ! ** test case for the mizoguchi lab experiment, Hansson et al. VZJ 2004
          case(hanssonVZJ)
            fArg  = 1._rkind + f1*nodeVolFracIce(iLayer)**f2
            xArg  = nodeVolFracLiq(iLayer) + fArg*nodeVolFracIce(iLayer)
            dxArg_dWat = dVolFracLiq_dWat + dVolFracIce_dWat * (1._rkind + f1*(f2+1)*nodeVolFracIce(iLayer)**f2)
            dxArg_dTk  = dVolFracIce_dTk * f1*(f2+1)*nodeVolFracIce(iLayer)**f2
            ! ...and, compute the thermal conductivity
            mLayerThermalC(iLayer) = c1 + c2*xArg + (c1 - c4)*exp(-(c3*xArg)**c5)

            ! compute derivatives
            dThermalC_dWat(iLayer) = ( c2 - c5*c3*(c3*xArg)**(c5-1)*(c1 - c4)*exp(-(c3*xArg)**c5) ) * dxArg_dWat
            dThermalC_dNrg(iLayer) = ( c2 - c5*c3*(c3*xArg)**(c5-1)*(c1 - c4)*exp(-(c3*xArg)**c5) ) * dxArg_dTk

          ! ** check
          case default; err=20; message=trim(message)//'unable to identify option for thermal conductivity of soil'; return

        end select  ! option for the thermal conductivity of soil

      ! ***** snow
      case(iname_snow)
        dVolFracIce_dWat = ( 1._rkind - fracLiqSnow(iLayer) )*(iden_water/iden_ice)
        dVolFracIce_dTk = -dTheta_dTk(iLayer)*(iden_water/iden_ice)

        ! temporally constant thermal conductivity
        if(ixThCondSnow==Smirnova2000)then
          mLayerThermalC(iLayer) = fixedThermalCond_snow
          dThermalC_dWat(iLayer) = 0._rkind
          dThermalC_dNrg(iLayer) = 0._rkind
          ! thermal conductivity as a function of snow density
        else
          call tcond_snow(nodeVolFracIce(iLayer)*iden_ice,  & ! input: snow density (kg m-3)
                          mLayerThermalC(iLayer),             & ! output: thermal conductivity (W m-1 K-1)
                          err,cmessage)                         ! output: error control
          if(err/=0)then; message=trim(message)//trim(cmessage); return; end if

          select case(ixThCondSnow)
            case(Yen1965)
              dThermalC_dWat(iLayer) = 2._rkind * 3.217d-6 * nodeVolFracIce(iLayer) * iden_ice**2._rkind * dVolFracIce_dWat
              dThermalC_dNrg(iLayer) = 2._rkind * 3.217d-6 * nodeVolFracIce(iLayer) * iden_ice**2._rkind * dVolFracIce_dTk
            case(Mellor1977)
              dThermalC_dWat(iLayer) = 2._rkind * 2.576d-6 * nodeVolFracIce(iLayer) * iden_ice**2._rkind * dVolFracIce_dWat
              dThermalC_dNrg(iLayer) = 2._rkind * 2.576d-6 * nodeVolFracIce(iLayer) * iden_ice**2._rkind * dVolFracIce_dTk
            case(Jordan1991)
              dThermalC_dWat(iLayer) = ( 7.75d-5 + 2._rkind * 1.105d-6 * nodeVolFracIce(iLayer) * iden_ice ) * (lambda_ice-lambda_air) * iden_ice * dVolFracIce_dWat
              dThermalC_dNrg(iLayer) = ( 7.75d-5 + 2._rkind * 1.105d-6 * nodeVolFracIce(iLayer) * iden_ice ) * (lambda_ice-lambda_air) * iden_ice * dVolFracIce_dTk
          end select  ! option for the thermal conductivity of snow
        end if

      ! * error check
      case default; err=20; message=trim(message)//'unable to identify type of layer (snow or soil) to compute thermal conductivity'; return

    end select

  end do  ! looping through layers

  ! *****
  ! * compute the thermal conductivity of snow at the interface of each layer...
  ! ****************************************************************************
  if (ixLayerDesired==0) then
    ! special case of hansson
    if(ixThCondSoil==hanssonVZJ)then
      iLayerThermalC = 28._rkind*(0.5_rkind*(node_iHeight(ixLower) - node_iHeight(ixUpper))) ! these are indices 1,0 since was passed with 0:1
      dThermalC_dHydStateBelow = 0._rkind
      dThermalC_dNrgStateBelow = 0._rkind
    else
      iLayerThermalC = mLayerThermalC(ixUpper) ! index was passed with 1:1
      dThermalC_dHydStateBelow = dThermalC_dWat(ixUpper)
      dThermalC_dNrgStateBelow = dThermalC_dNrg(ixUpper)
    end if
    dThermalC_dHydStateAbove = realMissing
    dThermalC_dNrgStateAbove = realMissing
  else if (ixLayerDesired==nLayers ) then
    ! assume the thermal conductivity at the domain boundaries is equal to the thermal conductivity of the layer
    iLayerThermalC = mLayerThermalC(ixLower) ! index was passed with iLayers:iLayers
    dThermalC_dHydStateAbove = dThermalC_dWat(ixLower)
    dThermalC_dNrgStateAbove = dThermalC_dNrg(ixLower)
    dThermalC_dHydStateBelow = realMissing
    dThermalC_dNrgStateBelow = realMissing
  else
    ! get temporary variables
    TCn = mLayerThermalC(ixUpper)    ! thermal conductivity below the layer interface (W m-1 K-1)
    TCp = mLayerThermalC(ixLower)  ! thermal conductivity above the layer interface (W m-1 K-1)
    zdn = node_iHeight(ixUpper)   - nodeHeight(ixUpper) ! height difference between interface and lower value (m)
    zdp = nodeHeight(ixLower) - node_iHeight(ixUpper) ! height difference between interface and upper value (m)
    den = TCn*zdp + TCp*zdn  ! denominator
    ! compute thermal conductivity
    if(TCn+TCp > epsilon(TCn))then
      iLayerThermalC = (TCn*TCp*(zdn + zdp)) / den
      dThermalC_dHydStateBelow = ( TCn*(zdn + zdp) - iLayerThermalC*zdn ) / den * dThermalC_dWat(ixLower)
      dThermalC_dHydStateAbove = ( TCp*(zdn + zdp) - iLayerThermalC*zdp ) / den * dThermalC_dWat(ixUpper)
      dThermalC_dNrgStateBelow = ( TCn*(zdn + zdp) - iLayerThermalC*zdn ) / den * dThermalC_dNrg(ixLower)
      dThermalC_dNrgStateAbove = ( TCp*(zdn + zdp) - iLayerThermalC*zdp ) / den * dThermalC_dNrg(ixUpper)
    else
      iLayerThermalC = (TCn*zdn +  TCp*zdp) / (zdn + zdp)
      dThermalC_dHydStateBelow = zdp / (zdn + zdp) * dThermalC_dWat(ixLower)
      dThermalC_dHydStateAbove = zdn / (zdn + zdp) * dThermalC_dWat(ixUpper)
      dThermalC_dNrgStateBelow = zdp / (zdn + zdp) * dThermalC_dNrg(ixLower)
      dThermalC_dNrgStateAbove = zdn / (zdn + zdp) * dThermalC_dNrg(ixUpper)
    end if
  endif

end subroutine iLayerThermalConduct



end module ssdNrgFlux_module


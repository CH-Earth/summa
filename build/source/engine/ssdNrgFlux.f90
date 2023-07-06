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
                    iden_water,  &  ! intrinsic density of water    (kg m-3)
                    ! specific heat
                    Cp_water        ! specific heat of liquid water (J kg-1 K-1)

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
 zeroFlux                         ! zero flux

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
                      ! input: fluxes and derivatives at the upper boundary
                      groundNetFlux,                      & ! intent(in):    total flux at the ground surface (W m-2)
                      dGroundNetFlux_dGroundTemp,         & ! intent(in):    derivative in total ground surface flux w.r.t. ground temperature (W m-2 K-1)
                      ! input: liquid water fluxes
                      iLayerLiqFluxSnow,                  & ! intent(in):    liquid flux at the interface of each snow layer (m s-1)
                      iLayerLiqFluxSoil,                  & ! intent(in):    liquid flux at the interface of each soil layer (m s-1)
                      ! input: trial value of model state variables
                      mLayerTempTrial,                    & ! intent(in):    trial temperature at the current iteration (K)
                      ! input: derivatives
                      dThermalC_dWatAbove,                & ! intent(in):    derivative in the thermal conductivity w.r.t. water state in the layer above
                      dThermalC_dWatBelow,                & ! intent(in):    derivative in the thermal conductivity w.r.t. water state in the layer above
                      dThermalC_dTempAbove,               & ! intent(in):    derivative in the thermal conductivity w.r.t. energy state in the layer above
                      dThermalC_dTempBelow,               & ! intent(in):    derivative in the thermal conductivity w.r.t. energy state in the layer above
                      ! input-output: data structures
                      mpar_data,                          & ! intent(in):    model parameters
                      indx_data,                          & ! intent(in):    model indices
                      prog_data,                          & ! intent(in):    model prognostic variables for a local HRU
                      diag_data,                          & ! intent(in):    model diagnostic variables for a local HRU
                      flux_data,                          & ! intent(inout): model fluxes for a local HRU
                      ! output: fluxes and derivatives at all layer interfaces
                      iLayerNrgFlux,                      & ! intent(out):   energy flux at the layer interfaces (W m-2)
                      dFlux_dTempAbove,                   & ! intent(out):   derivatives in the flux w.r.t. temperature in the layer above (W m-2 K-1)
                      dFlux_dTempBelow,                   & ! intent(out):   derivatives in the flux w.r.t. temperature in the layer below (W m-2 K-1)
                      dFlux_dWatAbove,                    & ! intent(out):   derivatives in the flux w.r.t. water state in the layer above (W m-2 K-1)
                      dFlux_dWatBelow,                    & ! intent(out):   derivatives in the flux w.r.t. water state in the layer below (W m-2 K-1)
                      ! output: error control
                      err,message)                          ! intent(out):   error control

  ! -------------------------------------------------------------------------------------------------------------------------------------------------
  implicit none
  ! input: model control
  logical(lgt),intent(in)             :: scalarSolution             ! flag to denote if implementing the scalar solution
  ! input: fluxes and derivatives at the upper boundary
  real(rkind),intent(in)              :: groundNetFlux              ! net energy flux for the ground surface (W m-2)
  real(rkind),intent(inout)           :: dGroundNetFlux_dGroundTemp ! derivative in net ground flux w.r.t. ground temperature (W m-2 K-1)
  ! input: liquid water fluxes
  real(rkind),intent(in)              :: iLayerLiqFluxSnow(0:)      ! liquid flux at the interface of each snow layer (m s-1)
  real(rkind),intent(in)              :: iLayerLiqFluxSoil(0:)      ! liquid flux at the interface of each soil layer (m s-1)
  ! input: trial model state variables
  real(rkind),intent(in)              :: mLayerTempTrial(:)         ! temperature in each layer at the current iteration (m)
  ! input: derivatives
  real(rkind),intent(in)              :: dThermalC_dWatAbove(0:)    ! derivative in the thermal conductivity w.r.t. water state in the layer above
  real(rkind),intent(in)              :: dThermalC_dWatBelow(0:)    ! derivative in the thermal conductivity w.r.t. water state in the layer above
  real(rkind),intent(in)              :: dThermalC_dTempAbove(0:)   ! derivative in the thermal conductivity w.r.t. energy state in the layer above
  real(rkind),intent(in)              :: dThermalC_dTempBelow(0:)   ! derivative in the thermal conductivity w.r.t. energy state in the layer above
  ! input-output: data structures
  type(var_dlength),intent(in)        :: mpar_data                  ! model parameters
  type(var_ilength),intent(in)        :: indx_data                  ! state vector geometry
  type(var_dlength),intent(in)        :: prog_data                  ! prognostic variables for a local HRU
  type(var_dlength),intent(in)        :: diag_data                  ! diagnostic variables for a local HRU
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
  character(LEN=256)                  :: cmessage                   ! error message of downwind routine
  integer(i4b)                        :: iLayer                     ! index of model layers
  integer(i4b)                        :: ixLayerDesired(1)          ! layer desired (scalar solution)
  integer(i4b)                        :: ixTop                      ! top layer in subroutine call
  integer(i4b)                        :: ixBot                      ! bottom layer in subroutine call
  real(rkind)                         :: qFlux                      ! liquid flux at layer interfaces (m s-1)
  real(rkind)                         :: dz                         ! height difference (m)
  ! ------------------------------------------------------------------------------------------------------------------------------------------------------
  ! make association of local variables with information in the data structures
  associate(&
    ix_bcUpprTdyn           => model_decisions(iLookDECISIONS%bcUpprTdyn)%iDecision, & ! intent(in): method used to calculate the upper boundary condition for thermodynamics
    ix_bcLowrTdyn           => model_decisions(iLookDECISIONS%bcLowrTdyn)%iDecision, & ! intent(in): method used to calculate the lower boundary condition for thermodynamics
    ! input: coordinate variables
    nSnow                   => indx_data%var(iLookINDEX%nSnow)%dat(1),               & ! intent(in): number of snow layers
    nLayers                 => indx_data%var(iLookINDEX%nLayers)%dat(1),             & ! intent(in): total number of layers
    layerType               => indx_data%var(iLookINDEX%layerType)%dat,              & ! intent(in): layer type (iname_soil or iname_snow)
    ixLayerState            => indx_data%var(iLookINDEX%ixLayerState)%dat,           & ! intent(in): list of indices for all model layers
    ixSnowSoilNrg           => indx_data%var(iLookINDEX%ixSnowSoilNrg)%dat,          & ! intent(in): index in the state subset for energy state variables in the snow+soil domain
    mLayerDepth             => prog_data%var(iLookPROG%mLayerDepth)%dat,             & ! intent(in): depth of each layer (m)
    mLayerHeight            => prog_data%var(iLookPROG%mLayerHeight)%dat,            & ! intent(in): height at the mid-point of each layer (m)
    ! input: thermal properties
    upperBoundTemp          => mpar_data%var(iLookPARAM%upperBoundTemp)%dat(1),      & ! intent(in): temperature of the upper boundary (K)
    lowerBoundTemp          => mpar_data%var(iLookPARAM%lowerBoundTemp)%dat(1),      & ! intent(in): temperature of the lower boundary (K)
    iLayerThermalC          => diag_data%var(iLookDIAG%iLayerThermalC)%dat,          & ! intent(in): thermal conductivity at the interface of each layer (W m-1 K-1)
     ! output: diagnostic fluxes
    iLayerConductiveFlux => flux_data%var(iLookFLUX%iLayerConductiveFlux)%dat,       & ! intent(out): conductive energy flux at layer interfaces at end of time step (W m-2)
    iLayerAdvectiveFlux  => flux_data%var(iLookFLUX%iLayerAdvectiveFlux)%dat         & ! intent(out): advective energy flux at layer interfaces at end of time step (W m-2)
    )  ! association of local variables with information in the data structures
    ! ------------------------------------------------------------------------------------------------------------------------------------------------------
    ! initialize error control
    err=0; message='ssdNrgFlux/'

    ! set conductive and advective fluxes to missing in the upper boundary
    iLayerConductiveFlux(0) = realMissing
    iLayerAdvectiveFlux(0)  = realMissing  !included in the ground heat flux

    ! get the indices for the snow+soil layers
    if(scalarSolution)then
      ixLayerDesired = pack(ixLayerState, ixSnowSoilNrg/=integerMissing)
      ixTop = ixLayerDesired(1)
      ixBot = ixLayerDesired(1)
    else
      ixTop = 1
      ixBot = nLayers
    endif

    ! -------------------------------------------------------------------------------------------------------------------------
    ! ***** compute the conductive fluxes at layer interfaces *****
    ! -------------------------------------------------------------------------------------------------------------------------
    do iLayer=ixTop,ixBot
      if(iLayer==nLayers)then ! (lower boundary fluxes -- positive downwards)
      ! flux depends on the type of lower boundary condition
        select case(ix_bcLowrTdyn) ! (identify the lower boundary condition for thermodynamics
          case(prescribedTemp); iLayerConductiveFlux(iLayer) = -iLayerThermalC(iLayer)*(lowerBoundTemp - mLayerTempTrial(iLayer))/(mLayerDepth(iLayer)*0.5_rkind)
          case(zeroFlux);       iLayerConductiveFlux(iLayer) = 0._rkind
        end select  ! (identifying the lower boundary condition for thermodynamics)
      else ! (domain boundary fluxes -- positive downwards)
        iLayerConductiveFlux(iLayer)  = -iLayerThermalC(iLayer)*(mLayerTempTrial(iLayer+1) - mLayerTempTrial(iLayer)) / &
                                        (mLayerHeight(iLayer+1) - mLayerHeight(iLayer))
      end if ! (the type of layer)
    end do  ! looping through layers

    ! -------------------------------------------------------------------------------------------------------------------------
    ! ***** compute the advective fluxes at layer interfaces *****
    ! -------------------------------------------------------------------------------------------------------------------------
    do iLayer=ixTop,ixBot
      ! get the liquid flux at layer interfaces
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
    end do  ! looping through layers

    ! -------------------------------------------------------------------------------------------------------------------------
    ! ***** compute the total fluxes at layer interfaces *****
    ! -------------------------------------------------------------------------------------------------------------------------
    ! NOTE: ignore advective fluxes for now
    iLayerNrgFlux(0)           = groundNetFlux !from vegNrgFlux module
    iLayerNrgFlux(ixTop:ixBot) = iLayerConductiveFlux(ixTop:ixBot)

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

    ! ***** the upper boundary, always do
  	select case(ix_bcUpprTdyn)

  	  ! * prescribed temperature at the upper boundary
  	  case(prescribedTemp)
  		dz = mLayerHeight(1)*0.5_rkind
  		dFlux_dWatBelow(0)  = -dThermalC_dWatBelow(0) * ( mLayerTempTrial(1) - upperBoundTemp )/dz
   		dFlux_dTempBelow(0) = -dThermalC_dTempBelow(0) * ( mLayerTempTrial(1) - upperBoundTemp )/dz - iLayerThermalC(0)/dz

  	  ! * zero flux at the upper boundary
  	  case(zeroFlux)
  		dFlux_dWatBelow(0) = 0._rkind
  		dFlux_dTempBelow(0) = 0._rkind

  	  ! * compute flux inside vegetation energy flux routine, use here
  	  case(energyFlux)
  		dFlux_dWatBelow(0) = 0._rkind !dGroundNetFlux_dGroundWat, does not exist in vegNrgFlux
  		dFlux_dTempBelow(0) = dGroundNetFlux_dGroundTemp

  	  case default; err=20; message=trim(message)//'unable to identify upper boundary condition for thermodynamics'; return

  	end select  ! (identifying the upper boundary condition for thermodynamics)
  	!dGroundNetFlux_dGroundWat  = dFlux_dWatBelow(0) ! this is true, but since not used in vegNrgFlux do not define
  	dGroundNetFlux_dGroundTemp = dFlux_dTempBelow(0) ! need this in vegNrgFlux

    ! loop through INTERFACES...
    do iLayer=ixTop,ixBot
      ! ***** the lower boundary
      if(iLayer==nLayers)then  ! (lower boundary)

        ! identify the lower boundary condition
        select case(ix_bcLowrTdyn) !prescribed temperature at the lower boundary
          case(prescribedTemp)
            dz = mLayerDepth(iLayer)*0.5_rkind
            dFlux_dWatAbove(iLayer)  = -dThermalC_dWatAbove(iLayer) * ( lowerBoundTemp - mLayerTempTrial(iLayer) )/dz
            dFlux_dTempAbove(iLayer) = -dThermalC_dTempAbove(iLayer) * ( lowerBoundTemp - mLayerTempTrial(iLayer) )/dz + iLayerThermalC(iLayer)/dz
          case(zeroFlux)  ! zero flux at the lower boundary
            dFlux_dWatAbove(iLayer) = 0._rkind
            dFlux_dTempAbove(iLayer) = 0._rkind
          case default; err=20; message=trim(message)//'unable to identify lower boundary condition for thermodynamics'; return
        end select  ! (identifying the lower boundary condition for thermodynamics)

      ! ***** internal layers
      else
        dz = (mLayerHeight(iLayer+1) - mLayerHeight(iLayer))
        dFlux_dWatAbove(iLayer)  = -dThermalC_dWatAbove(iLayer) * ( mLayerTempTrial(iLayer+1) - mLayerTempTrial(iLayer) )/dz
        dFlux_dWatBelow(iLayer)  = -dThermalC_dWatBelow(iLayer) * ( mLayerTempTrial(iLayer+1) - mLayerTempTrial(iLayer) )/dz
        dFlux_dTempAbove(iLayer) = -dThermalC_dTempAbove(iLayer) * ( mLayerTempTrial(iLayer+1) - mLayerTempTrial(iLayer) )/dz + iLayerThermalC(iLayer)/dz
        dFlux_dTempBelow(iLayer) = -dThermalC_dTempBelow(iLayer) * ( mLayerTempTrial(iLayer+1) - mLayerTempTrial(iLayer) )/dz - iLayerThermalC(iLayer)/dz
      end if  ! type of layer (upper, internal, or lower)

    end do  ! (looping through layers)

  ! end association of local variables with information in the data structures
  end associate

end subroutine ssdNrgFlux


end module ssdNrgFlux_module


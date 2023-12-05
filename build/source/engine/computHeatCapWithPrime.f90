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

module computHeatCapWithPrime_module

! data types
USE nrtype

! derived types to define the data structures
USE data_types,only:&
                    var_d,            & ! data vector (rkind)
                    var_ilength,      & ! data vector with variable length dimension (i4b)
                    var_dlength         ! data vector with variable length dimension (rkind)

! named variables defining elements in the data structures
USE var_lookup,only:iLookPARAM,iLookDIAG,iLookINDEX  ! named variables for structure elements

! physical constants
USE multiconst,only:&
                    Tfreeze,     & ! freezing point of water (K)
                    iden_air,    & ! intrinsic density of air      (kg m-3)
                    iden_ice,    & ! intrinsic density of ice      (kg m-3)
                    iden_water,  & ! intrinsic density of water    (kg m-3)
                    ! specific heat
                    Cp_air,      & ! specific heat of air          (J kg-1 K-1)
                    Cp_ice,      & ! specific heat of ice          (J kg-1 K-1)
                    Cp_soil,     & ! specific heat of soil         (J kg-1 K-1)
                    Cp_water       ! specific heat of liquid water (J kg-1 K-1)
! named variables to describe the state variable type
USE globalData,only:iname_nrgCanair  ! named variable defining the energy of the canopy air space
USE globalData,only:iname_nrgCanopy  ! named variable defining the energy of the vegetation canopy
USE globalData,only:iname_watCanopy  ! named variable defining the mass of total water on the vegetation canopy
USE globalData,only:iname_liqCanopy  ! named variable defining the mass of liquid water on the vegetation canopy
USE globalData,only:iname_nrgLayer   ! named variable defining the energy state variable for snow+soil layers
USE globalData,only:iname_watLayer   ! named variable defining the total water state variable for snow+soil layers
USE globalData,only:iname_liqLayer   ! named variable defining the liquid  water state variable for snow+soil layers
USE globalData,only:iname_matLayer   ! named variable defining the matric head state variable for soil layers
USE globalData,only:iname_lmpLayer   ! named variable defining the liquid matric potential state variable for soil layers
USE globalData,only:iname_watAquifer ! named variable defining the water storage in the aquifer

! missing values
USE globalData,only:integerMissing ! missing integer
USE globalData,only:realMissing    ! missing real

! named variables that define the layer type
USE globalData,only:iname_snow     ! snow
USE globalData,only:iname_soil     ! soil


! privacy
implicit none
private
public::computHeatCapWithPrime

contains


! **********************************************************************************************************
! public subroutine computHeatCapWithPrime: compute diagnostic energy variables (heat capacity)
! **********************************************************************************************************
subroutine computHeatCapWithPrime(&
                     ! input: control variables
                     nLayers,                     & ! intent(in):    number of layers (soil+snow)
                     computeVegFlux,              & ! intent(in):    flag to denote if computing the vegetation flux
                     canopyDepth,                 & ! intent(in):    canopy depth (m)
                     ! input output data structures
                     mpar_data,                   & ! intent(in):    model parameters
                     indx_data,                   & ! intent(in):    model layer indices
                     diag_data,                   & ! intent(inout): model diagnostic variables for a local HRU
                     ! input: state variables
                     scalarCanopyIce,             & ! intent(in):    trial value for mass of ice on the vegetation canopy (kg m-2)
                     scalarCanopyLiquid,          & ! intent(in):    trial value for the liquid water on the vegetation canopy (kg m-2)
                     scalarCanopyTempTrial,       & ! intent(in):    trial value of canopy temperature
                     scalarCanopyTempPrime,       & ! intent(in):    prime value of canopy temperature (K)
                     scalarCanopyEnthalpyPrime,   & ! intent(in):    prime enthalpy of the vegetation canopy (J m-3)
                     mLayerVolFracIce,            & ! intent(in):    volumetric fraction of ice at the start of the sub-step (-)
                     mLayerVolFracLiq,            & ! intent(in):    volumetric fraction of liquid water at the start of the sub-step (-)
                     mLayerTempTrial,             & ! intent(in):    trial temperature
                     mLayerTempPrime,             & ! intent(in):    prime temperature
                     mLayerEnthalpyPrime,         & ! intent(in):    prime enthalpy for snow and soil
                     mLayerMatricHeadTrial,       & ! intent(in):    trial total water matric potential (m)
                     ! input: pre-computed derivatives
                     dTheta_dTkCanopy,            & ! intent(in):    derivative in canopy volumetric liquid water content w.r.t. temperature (K-1)
                     scalarFracLiqVeg,            & ! intent(in):    fraction of canopy liquid water (-)
                     mLayerdTheta_dTk,            & ! intent(in):    derivative of volumetric liquid water content w.r.t. temperature (K-1)
                     mLayerFracLiqSnow,           & ! intent(in):    fraction of liquid water (-)
                     dVolTot_dPsi0,               & ! intent(in):    derivative in total water content w.r.t. total water matric potential (m-1)
                     ! output
                     heatCapVeg,                  & ! intent(out):   heat capacity for canopy
                     mLayerHeatCap,               & ! intent(out):   heat capacity for snow and soil
                     dVolHtCapBulk_dPsi0,         & ! intent(out):   derivative in bulk heat capacity w.r.t. matric potential
                     dVolHtCapBulk_dTheta,        & ! intent(out):   derivative in bulk heat capacity w.r.t. volumetric water content
                     dVolHtCapBulk_dCanWat,       & ! intent(out):   derivative in bulk heat capacity w.r.t. volumetric water content
                     dVolHtCapBulk_dTk,           & ! intent(out):   derivative in bulk heat capacity w.r.t. temperature
                     dVolHtCapBulk_dTkCanopy,     & ! intent(out):   derivative in bulk heat capacity w.r.t. temperature         
                     ! output: error control
                     err,message)                   ! intent(out): error control
  ! --------------------------------------------------------------------------------------------------------------------------------------
  ! provide access to external subroutines
  USE soil_utils_module,only:crit_soilT     ! compute critical temperature below which ice exists
  ! --------------------------------------------------------------------------------------------------------------------------------
  ! input: control variables
  integer(i4b),intent(in)         :: nLayers                   ! number of layers (soil+snow)
  logical(lgt),intent(in)         :: computeVegFlux            ! logical flag to denote if computing the vegetation flux
  real(rkind),intent(in)          :: canopyDepth               ! depth of the vegetation canopy (m)
  ! input/output: data structures
  type(var_dlength),intent(in)    :: mpar_data                 ! model parameters
  type(var_ilength),intent(in)    :: indx_data                 ! model layer indices
  type(var_dlength),intent(inout) :: diag_data                 ! diagnostic variables for a local HRU
  ! input: state variables
  real(rkind),intent(in)          :: scalarCanopyIce           ! trial value of canopy ice content (kg m-2)
  real(rkind),intent(in)          :: scalarCanopyLiquid        ! trial value for the liquid water on the vegetation canopy (kg m-2)
  real(rkind),intent(in)          :: scalarCanopyTempTrial     ! trial value of canopy temperature
  real(rkind),intent(in)          :: scalarCanopyTempPrime     ! prime value of canopy temperature
  real(rkind),intent(in)          :: scalarCanopyEnthalpyPrime ! prime enthalpy of the vegetation canopy (J m-3)
  real(rkind),intent(in)          :: mLayerVolFracLiq(:)       ! trial vector of volumetric liquid water content (-)
  real(rkind),intent(in)          :: mLayerVolFracIce(:)       ! trial vector of volumetric ice water content (-)
  real(rkind),intent(in)          :: mLayerTempTrial(:)        ! trial temperature
  real(rkind),intent(in)          :: mLayerTempPrime(:)        ! prime temperature
  real(rkind),intent(in)          :: mLayerEnthalpyPrime(:)    ! prime enthalpy for snow and soil
  real(rkind),intent(in)          :: mLayerMatricHeadTrial(:)  ! vector of total water matric potential (m)
  ! input: pre-computed derivatives
  real(rkind),intent(in)          :: dTheta_dTkCanopy          ! derivative in canopy volumetric liquid water content w.r.t. temperature (K-1)
  real(rkind),intent(in)          :: scalarFracLiqVeg          ! fraction of canopy liquid water (-)
  real(rkind),intent(in)          :: mLayerdTheta_dTk(:)       ! derivative of volumetric liquid water content w.r.t. temperature (K-1)
  real(rkind),intent(in)          :: mLayerFracLiqSnow(:)      ! fraction of liquid water (-)
  real(rkind),intent(in)          :: dVolTot_dPsi0(:)          ! derivative in total water content w.r.t. total water matric potential (m-1)
  ! output:
  real(qp),intent(out)            :: heatCapVeg                ! heat capacity for canopy
  real(qp),intent(out)            :: mLayerHeatCap(:)          ! heat capacity for snow and soil
  real(rkind),intent(out)         :: dVolHtCapBulk_dPsi0(:)    ! derivative in bulk heat capacity w.r.t. matric potential
  real(rkind),intent(out)         :: dVolHtCapBulk_dTheta(:)   ! derivative in bulk heat capacity w.r.t. volumetric water content
  real(rkind),intent(out)         :: dVolHtCapBulk_dCanWat     ! derivative in bulk heat capacity w.r.t. volumetric water content
  real(rkind),intent(out)         :: dVolHtCapBulk_dTk(:)      ! derivative in bulk heat capacity w.r.t. temperature
  real(rkind),intent(out)         :: dVolHtCapBulk_dTkCanopy   ! derivative in bulk heat capacity w.r.t. temperature
  ! output: error control
  integer(i4b),intent(out)        :: err                       ! error code
  character(*),intent(out)        :: message                   ! error message
  ! --------------------------------------------------------------------------------------------------------------------------------
  ! local variables
  integer(i4b)                    :: iLayer                    ! index of model layer
  integer(i4b)                    :: iSoil                     ! index of soil layer
  real(rkind)                     :: delT                      ! temperature change
  real(rkind)                     :: delEnt                    ! enthalpy change
  real(rkind)                     :: fLiq                      ! fraction of liquid water
  real(rkind)                     :: Tcrit                     ! temperature where all water is unfrozen (K)
  ! --------------------------------------------------------------------------------------------------------------------------------
  ! associate variables in data structure
  associate(&
    ! input: coordinate variables
    nSnow                   => indx_data%var(iLookINDEX%nSnow)%dat(1),             & ! intent(in): number of snow layers
    layerType               => indx_data%var(iLookINDEX%layerType)%dat,            & ! intent(in): layer type (iname_soil or iname_snow)
    ! input: heat capacity and thermal conductivity
    specificHeatVeg         => mpar_data%var(iLookPARAM%specificHeatVeg)%dat(1),   & ! intent(in): specific heat of vegetation (J kg-1 K-1)
    maxMassVegetation       => mpar_data%var(iLookPARAM%maxMassVegetation)%dat(1), & ! intent(in): maximum mass of vegetation (kg m-2)
    ! input: depth varying soil parameters
    iden_soil               => mpar_data%var(iLookPARAM%soil_dens_intr)%dat,       & ! intent(in): intrinsic density of soil (kg m-3)
    theta_sat               => mpar_data%var(iLookPARAM%theta_sat)%dat             & ! intent(in): soil porosity (-)
    )  ! end associate statemen
    ! initialize error control
    err=0; message="computHeatCapWithPrime/"

    ! initialize the soil layer
    iSoil=integerMissing

    ! compute the bulk volumetric heat capacity of vegetation (J m-3 K-1)
    if(computeVegFlux)then
      delT = scalarCanopyTempPrime
      if(abs(delT) <= 1e-2_rkind)then
        heatCapVeg = specificHeatVeg*maxMassVegetation/canopyDepth + & ! vegetation component
                     Cp_water*scalarCanopyLiquid/canopyDepth       + & ! liquid water component
                     Cp_ice*scalarCanopyIce/canopyDepth                ! ice component
      else
        delEnt = scalarCanopyEnthalpyPrime
        heatCapVeg = delEnt / delT
      endif

      ! derivatives for Jacobian are from analytical solution
      fLiq = scalarFracLiqVeg
      dVolHtCapBulk_dCanWat = ( -Cp_ice*( fLiq-1._rkind ) + Cp_water*fLiq )/canopyDepth !this is iden_water/(iden_water*canopyDepth)
      if(scalarCanopyTempTrial < Tfreeze)then
        dVolHtCapBulk_dTkCanopy = iden_water * (-Cp_ice + Cp_water) * dTheta_dTkCanopy ! no derivative in air
      else
        dVolHtCapBulk_dTkCanopy = 0._rkind
      endif
    endif

    ! loop through layers
    do iLayer=1,nLayers
      delT = mLayerTempPrime(iLayer)
      if(abs(delT) <= 1e-2_rkind)then
        ! get the soil layer
        if(iLayer>nSnow) iSoil = iLayer-nSnow
          select case(layerType(iLayer))
            ! * soil
            case(iname_soil)
              mLayerHeatCap(iLayer) = iden_soil(iSoil) * Cp_soil  * ( 1._rkind - theta_sat(iSoil) ) + & ! soil component
                                      iden_ice         * Cp_ice   * mLayerVolFracIce(iLayer)        + & ! ice component
                                      iden_water       * Cp_water * mLayerVolFracLiq(iLayer)        + & ! liquid water component
                                      iden_air         * Cp_air   * ( theta_sat(iSoil) - (mLayerVolFracIce(iLayer) + mLayerVolFracLiq(iLayer)) )! air component
            ! * snow
            case(iname_snow)
              mLayerHeatCap(iLayer) = iden_ice         * Cp_ice   * mLayerVolFracIce(iLayer)     + & ! ice component
                                      iden_water       * Cp_water * mLayerVolFracLiq(iLayer)     + & ! liquid water component
                                      iden_air         * Cp_air   * ( 1._rkind - (mLayerVolFracIce(iLayer) + mLayerVolFracLiq(iLayer)) )   ! air component

           case default; err=20; message=trim(message)//'unable to identify type of layer (snow or soil) to compute olumetric heat capacity'; return
          end select
      else
        delEnt = mLayerEnthalpyPrime(iLayer)
        mLayerHeatCap(iLayer) = delEnt / delT
      endif

      ! derivatives for Jacobian are from analytical solution
      select case(layerType(iLayer))
        ! * soil
        case(iname_soil)
          dVolHtCapBulk_dTheta(iLayer) = realMissing ! do not use
          Tcrit = crit_soilT( mLayerMatricHeadTrial(iSoil) )
          if( mLayerTempTrial(iLayer) < Tcrit)then
            dVolHtCapBulk_dPsi0(iSoil) = (iden_ice * Cp_ice   - iden_air * Cp_air) * dVolTot_dPsi0(iSoil)
            dVolHtCapBulk_dTk(iLayer) = (-iden_ice * Cp_ice + iden_water * Cp_water) * mLayerdTheta_dTk(iLayer)
          else
            dVolHtCapBulk_dPsi0(iSoil) = (iden_water*Cp_water - iden_air * Cp_air) * dVolTot_dPsi0(iSoil)
            dVolHtCapBulk_dTk(iLayer) = 0._rkind
          endif
        ! * snow
        case(iname_snow)
          fLiq = mLayerFracLiqSnow(iLayer)
          dVolHtCapBulk_dTheta(iLayer) = iden_water * ( -Cp_ice*( fLiq-1._rkind ) + Cp_water*fLiq ) + iden_air * ( ( fLiq-1._rkind )*iden_water/iden_ice - fLiq ) * Cp_air
          if( mLayerTempTrial(iLayer) < Tfreeze)then
            dVolHtCapBulk_dTk(iLayer) = ( iden_water * (-Cp_ice + Cp_water) + iden_air * (iden_water/iden_ice - 1._rkind) * Cp_air ) * mLayerdTheta_dTk(iLayer)
          else
            dVolHtCapBulk_dTk(iLayer) = 0._rkind
          endif
        end select

    end do  ! looping through layers

  end associate

end subroutine computHeatCapWithPrime


end module computHeatCapWithPrime_module
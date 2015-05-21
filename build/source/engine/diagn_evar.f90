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

module diagn_evar_module
! data types
USE nrtype

! physical constants
USE multiconst,only:&
                    iden_air,    & ! intrinsic density of air      (kg m-3)
                    iden_ice,    & ! intrinsic density of ice      (kg m-3)
                    iden_water,  & ! intrinsic density of water    (kg m-3)
                    ! specific heat
                    Cp_air,      & ! specific heat of air          (J kg-1 K-1)
                    Cp_ice,      & ! specific heat of ice          (J kg-1 K-1)
                    Cp_soil,     & ! specific heat of soil         (J kg-1 K-1)
                    Cp_water,    & ! specific heat of liquid water (J kg-1 K-1)
                    ! thermal conductivity
                    lambda_air,  & ! thermal conductivity of air   (J s-1 m-1)
                    lambda_ice,  & ! thermal conductivity of ice   (J s-1 m-1)
                    lambda_water   ! thermal conductivity of water (J s-1 m-1)

! access the number of snow and soil layers
USE data_struc,only:&
                    nSnow,        & ! number of snow layers
                    nSoil,        & ! number of soil layers
                    nLayers         ! total number of layers

! named variables that define the layer type
USE data_struc,only:ix_soil        ! soil
USE data_struc,only:ix_snow        ! snow
implicit none
private
public::diagn_evar
! algorithmic parameters
real(dp),parameter     :: valueMissing=-9999._dp  ! missing value, used when diagnostic or state variables are undefined
real(dp),parameter     :: verySmall=1.e-6_dp   ! used as an additive constant to check if substantial difference among real numbers
real(dp),parameter     :: mpe=1.e-6_dp         ! prevents overflow error if division by zero
real(dp),parameter     :: dx=1.e-6_dp          ! finite difference increment
contains


 ! **********************************************************************************************************
 ! public subroutine diagn_evar: compute diagnostic energy variables (thermal conductivity and heat capacity)
 ! **********************************************************************************************************
 subroutine diagn_evar(&
                       ! input: control variables
                       computeVegFlux,          & ! intent(in): flag to denote if computing the vegetation flux
                       canopyDepth,             & ! intent(in): canopy depth (m)
                       ! input/output: data structures
                       mpar_data,               & ! intent(in):    model parameters
                       indx_data,               & ! intent(in):    model layer indices
                       mvar_data,               & ! intent(inout): model variables for a local HRU
                       ! output: error control
                       err,message)               ! intent(out): error control
 ! --------------------------------------------------------------------------------------------------------------------------------------
 ! provide access to the derived types to define the data structures
 USE data_struc,only:&
                     var_d,            & ! data vector (dp)
                     var_ilength,      & ! data vector with variable length dimension (i4b)
                     var_dlength         ! data vector with variable length dimension (dp)
 ! provide access to named variables defining elements in the data structures
 USE var_lookup,only:iLookTIME,iLookTYPE,iLookATTR,iLookFORCE,iLookPARAM,iLookMVAR,iLookBVAR,iLookINDEX  ! named variables for structure elements
 USE var_lookup,only:iLookDECISIONS               ! named variables for elements of the decision structure
 ! provide access to named variables for thermal conductivity of soil
 USE data_struc,only:model_decisions        ! model decision structure
 USE mDecisions_module,only: funcSoilWet, & ! function of soil wetness
                             mixConstit,  & ! mixture of constituents
                             hanssonVZJ     ! test case for the mizoguchi lab experiment, Hansson et al. VZJ 2004
 ! provide access to external subroutines
 USE snow_utils_module,only:tcond_snow            ! compute thermal conductivity of snow
 ! --------------------------------------------------------------------------------------------------------------------------------------
 ! input: model control
 logical(lgt),intent(in)         :: computeVegFlux         ! logical flag to denote if computing the vegetation flux
 real(dp),intent(in)             :: canopyDepth            ! depth of the vegetation canopy (m)
 ! input/output: data structures
 type(var_d),intent(in)          :: mpar_data              ! model parameters
 type(var_ilength),intent(inout) :: indx_data              ! model layer indices
 type(var_dlength),intent(inout) :: mvar_data              ! model variables for a local HRU
 ! output: error control
 integer(i4b),intent(out)        :: err                    ! error code
 character(*),intent(out)        :: message                ! error message
 ! --------------------------------------------------------------------------------------------------------------------------------
 ! local variables
 character(LEN=256)                :: cmessage               ! error message of downwind routine
 integer(i4b)                      :: iLayer                 ! index of model layer
 real(dp)                          :: TCn                    ! thermal conductivity below the layer interface (W m-1 K-1)
 real(dp)                          :: TCp                    ! thermal conductivity above the layer interface (W m-1 K-1)
 real(dp)                          :: zdn                    ! height difference between interface and lower value (m)
 real(dp)                          :: zdp                    ! height difference between interface and upper value (m)
 real(dp)                          :: bulkden_soil           ! bulk density of soil (kg m-3)
 real(dp)                          :: lambda_drysoil         ! thermal conductivity of dry soil (W m-1)
 real(dp)                          :: lambda_wetsoil         ! thermal conductivity of wet soil (W m-1)
 real(dp)                          :: lambda_wet             ! thermal conductivity of the wet material
 real(dp)                          :: kerstenNum             ! the Kersten number (-), defining weight applied to conductivity of the wet medium
 ! local variables to reproduce the thermal conductivity of Hansson et al. VZJ 2005
 real(dp),parameter                :: c1=0.55_dp             ! optimized parameter from Hansson et al. VZJ 2005 (W m-1 K-1)
 real(dp),parameter                :: c2=0.8_dp              ! optimized parameter from Hansson et al. VZJ 2005 (W m-1 K-1)
 real(dp),parameter                :: c3=3.07_dp             ! optimized parameter from Hansson et al. VZJ 2005 (-)
 real(dp),parameter                :: c4=0.13_dp             ! optimized parameter from Hansson et al. VZJ 2005 (W m-1 K-1)
 real(dp),parameter                :: c5=4._dp               ! optimized parameter from Hansson et al. VZJ 2005 (-)
 real(dp),parameter                :: f1=13.05_dp            ! optimized parameter from Hansson et al. VZJ 2005 (-)
 real(dp),parameter                :: f2=1.06_dp             ! optimized parameter from Hansson et al. VZJ 2005 (-)
 real(dp)                          :: fArg,xArg              ! temporary variables (see Hansson et al. VZJ 2005 for details)
 ! --------------------------------------------------------------------------------------------------------------------------------
 ! initialize error control
 err=0; message="diagn_evar/"
 ! associate variables in data structure
 associate(&
 ! input: model decisions
 ixThCondSoil            => model_decisions(iLookDECISIONS%thCondSoil)%iDecision,      & ! intent(in): choice of method for thermal conductivity of soil
 ! input: state variables
 scalarCanopyIce         => mvar_data%var(iLookMVAR%scalarCanopyIce)%dat(1),           & ! intent(in): canopy ice content (kg m-2)
 scalarCanopyLiquid      => mvar_data%var(iLookMVAR%scalarCanopyLiq)%dat(1),           & ! intent(in): canopy liquid water content (kg m-2)
 mLayerVolFracIce        => mvar_data%var(iLookMVAR%mLayerVolFracIce)%dat,             & ! intent(in): volumetric fraction of ice at the start of the sub-step (-)
 mLayerVolFracLiq        => mvar_data%var(iLookMVAR%mLayerVolFracLiq)%dat,             & ! intent(in): volumetric fraction of liquid water at the start of the sub-step (-)
 ! input: coordinate variables
 layerType               => indx_data%var(iLookINDEX%layerType)%dat,                   & ! intent(in): layer type (ix_soil or ix_snow)
 mLayerHeight            => mvar_data%var(iLookMVAR%mLayerHeight)%dat,                 & ! intent(in): height at the mid-point of each layer (m)
 iLayerHeight            => mvar_data%var(iLookMVAR%iLayerHeight)%dat,                 & ! intent(in): height at the interface of each layer (m)
 ! input: heat capacity and thermal conductivity
 specificHeatVeg         => mpar_data%var(iLookPARAM%specificHeatVeg),                 & ! intent(in): specific heat of vegetation (J kg-1 K-1)
 maxMassVegetation       => mpar_data%var(iLookPARAM%maxMassVegetation),               & ! intent(in): maximum mass of vegetation (kg m-2)
 iden_soil               => mpar_data%var(iLookPARAM%soil_dens_intr),                  & ! intent(in): intrinsic density of soil (kg m-3)
 thCond_soil             => mpar_data%var(iLookPARAM%thCond_soil),                     & ! intent(in): thermal conductivity of soil (W m-1 K-1)
 theta_sat               => mpar_data%var(iLookPARAM%theta_sat),                       & ! intent(in): soil porosity (-)
 frac_sand               => mpar_data%var(iLookPARAM%frac_sand),                       & ! intent(in): fraction of sand (-)
 frac_silt               => mpar_data%var(iLookPARAM%frac_silt),                       & ! intent(in): fraction of silt (-)
 frac_clay               => mpar_data%var(iLookPARAM%frac_clay),                       & ! intent(in): fraction of clay (-)
 ! output: diagnostic variables
 scalarBulkVolHeatCapVeg => mvar_data%var(iLookMVAR%scalarBulkVolHeatCapVeg)%dat(1),   & ! intent(out): volumetric heat capacity of the vegetation (J m-3 K-1)
 mLayerVolHtCapBulk      => mvar_data%var(iLookMVAR%mLayerVolHtCapBulk)%dat,           & ! intent(out): volumetric heat capacity in each layer (J m-3 K-1)
 mLayerThermalC          => mvar_data%var(iLookMVAR%mLayerThermalC)%dat,               & ! intent(out): thermal conductivity at the mid-point of each layer (W m-1 K-1)
 iLayerThermalC          => mvar_data%var(iLookMVAR%iLayerThermalC)%dat,               & ! intent(out): thermal conductivity at the interface of each layer (W m-1 K-1)
 mLayerVolFracAir        => mvar_data%var(iLookMVAR%mLayerVolFracAir)%dat              & ! intent(out): volumetric fraction of air in each layer (-)
 )  ! end associate statement
 ! --------------------------------------------------------------------------------------------------------------------------------

 ! compute the bulk volumetric heat capacity of vegetation (J m-3 K-1)
 if(computeVegFlux)then
  scalarBulkVolHeatCapVeg = specificHeatVeg*maxMassVegetation/canopyDepth + & ! vegetation component
                            Cp_water*scalarCanopyLiquid/canopyDepth       + & ! liquid water component
                            Cp_ice*scalarCanopyIce/canopyDepth                ! ice component
 else
  scalarBulkVolHeatCapVeg = valueMissing
 endif
 !print*, 'diagn_evar: scalarBulkVolHeatCapVeg = ', scalarBulkVolHeatCapVeg

 ! compute the thermal conductivity of dry and wet soils (W m-1)
 ! NOTE: this is actually constant over the simulation, and included here for clarity
 if(ixThCondSoil == funcSoilWet)then
  bulkden_soil   = iden_soil*(1._dp - theta_sat)
  lambda_drysoil = (0.135_dp*bulkden_soil + 64.7_dp) / (iden_soil - 0.947_dp*bulkden_soil)
  lambda_wetsoil = (8.80_dp*frac_sand + 2.92_dp*frac_clay) / (frac_sand + frac_clay)
 endif

 ! loop through layers
 do iLayer=1,nLayers

  ! *****
  ! * compute the volumetric fraction of air in each layer...
  ! *********************************************************
  select case(layerType(iLayer))
   case(ix_soil); mLayerVolFracAir(iLayer) = theta_sat - (mLayerVolFracIce(iLayer) + mLayerVolFracLiq(iLayer))
   case(ix_snow); mLayerVolFracAir(iLayer) = 1._dp - (mLayerVolFracIce(iLayer) + mLayerVolFracLiq(iLayer))
   case default; err=20; message=trim(message)//'unable to identify type of layer (snow or soil) to compute volumetric fraction of air'; return
  end select

  ! *****
  ! * compute the volumetric heat capacity of each layer (J m-3 K-1)...
  ! *******************************************************************
  select case(layerType(iLayer))
   ! * soil
   case(ix_soil)
    mLayerVolHtCapBulk(iLayer) = iden_soil  * Cp_soil  * (1._dp - theta_sat)      + & ! soil component
                                 iden_ice   * Cp_Ice   * mLayerVolFracIce(iLayer) + & ! ice component
                                 iden_water * Cp_water * mLayerVolFracLiq(iLayer) + & ! liquid water component
                                 iden_air   * Cp_air   * mLayerVolFracAir(iLayer)     ! air component
   ! * snow
   case(ix_snow)
    mLayerVolHtCapBulk(iLayer) = iden_ice   * Cp_ice   * mLayerVolFracIce(iLayer) + & ! ice component
                                 iden_water * Cp_water * mLayerVolFracLiq(iLayer) + & ! liquid water component
                                 iden_air   * Cp_air   * mLayerVolFracAir(iLayer)     ! air component
   case default; err=20; message=trim(message)//'unable to identify type of layer (snow or soil) to compute olumetric heat capacity'; return
  end select

  ! *****
  ! * compute the thermal conductivity of snow and soil at the mid-point of each layer...
  ! *************************************************************************************
  select case(layerType(iLayer))

   ! ***** soil
   case(ix_soil)

    ! select option for thermal conductivity of soil
    select case(ixThCondSoil)

     ! ** function of soil wetness
     case(funcSoilWet)

      ! compute the thermal conductivity of the wet material (W m-1)
      lambda_wet = lambda_wetsoil**(1._dp - theta_sat) * lambda_water**theta_sat * lambda_ice**(theta_sat - mLayerVolFracLiq(iLayer))
      ! compute the Kersten number (-)
      kerstenNum = log10( (mLayerVolFracIce(iLayer) + mLayerVolFracLiq(iLayer))/theta_sat ) + 1._dp
      ! ...and, compute the thermal conductivity
      mLayerThermalC(iLayer) = kerstenNum*lambda_wet + (1._dp - kerstenNum)*lambda_drysoil

     ! ** mixture of constituents
     case(mixConstit)
      mLayerThermalC(iLayer) = thCond_soil * (1._dp - theta_sat)      + & ! soil component
                               lambda_ice  * mLayerVolFracIce(iLayer) + & ! ice component
                               lambda_water* mLayerVolFracLiq(iLayer) + & ! liquid water component
                               lambda_air  * mLayerVolFracAir(iLayer)     ! air component

     ! ** test case for the mizoguchi lab experiment, Hansson et al. VZJ 2004
     case(hanssonVZJ)
      fArg  = 1._dp + f1*mLayerVolFracIce(iLayer)**f2
      xArg  = mLayerVolFracLiq(iLayer) + fArg*mLayerVolFracIce(iLayer)
      mLayerThermalC(iLayer) = c1 + c2*xArg + (c1 - c4)*exp(-(c3*xArg)**c5)

     ! ** check
     case default; err=20; message=trim(message)//'unable to identify option for thermal conductivity of soil'; return

    end select  ! option for the thermal conductivity of soil
    
   ! ***** snow
   case(ix_snow)
    call tcond_snow(mLayerVolFracIce(iLayer)*iden_ice,mLayerThermalC(iLayer),err,cmessage)
    if(err/=0)then; message=trim(message)//trim(cmessage); return; endif
   ! * error check
   case default; err=20; message=trim(message)//'unable to identify type of layer (snow or soil) to compute thermal conductivity'; return
  end select
  !print*, 'iLayer, mLayerThermalC(iLayer) = ', iLayer, mLayerThermalC(iLayer)

 end do  ! looping through layers
 !pause

 ! *****
 ! * compute the thermal conductivity of snow at the interface of each layer...
 ! ****************************************************************************
 do iLayer=1,nLayers-1  ! (loop through layers)
  TCn = mLayerThermalC(iLayer)    ! thermal conductivity below the layer interface (W m-1 K-1)
  TCp = mLayerThermalC(iLayer+1)  ! thermal conductivity above the layer interface (W m-1 K-1)
  zdn = iLayerHeight(iLayer)   - mLayerHeight(iLayer) ! height difference between interface and lower value (m)
  zdp = mLayerHeight(iLayer+1) - iLayerHeight(iLayer) ! height difference between interface and upper value (m)
  iLayerThermalC(iLayer) = (TCn*TCp*(zdn + zdp)) / (TCn*zdp + TCp*zdn)
  !write(*,'(a,1x,i4,1x,10(f9.3,1x))') 'iLayer, TCn, TCp, zdn, zdp, iLayerThermalC(iLayer) = ', iLayer, TCn, TCp, zdn, zdp, iLayerThermalC(iLayer)
 end do

 ! special case of hansson
 if(ixThCondSoil==hanssonVZJ)then
  iLayerThermalC(0) = 28._dp*(0.5_dp*(iLayerHeight(1) - iLayerHeight(0)))
  print*, 'iLayerThermalC(0) = ', iLayerThermalC(0)
 else
  iLayerThermalC(0) = mLayerThermalC(1)
 endif

 ! assume the thermal conductivity at the domain boundaries is equal to the thermal conductivity of the layer
 iLayerThermalC(nLayers) = mLayerThermalC(nLayers)

 ! end association to variables in the data structure
 end associate

 end subroutine diagn_evar


end module diagn_evar_module

module diagn_evar_module
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

 ! ************************************************************************************************
 ! new subroutine: compute diagnostic energy variables (thermal conductivity and heat capacity) 
 ! ************************************************************************************************
 subroutine diagn_evar(&
                       ! input: control variables
                       computeVegFlux,          & ! intent(in): flag to denote if computing the vegetation flux
                       ! input: state variables
                       scalarCanopyIce,         & ! intent(in): mass of ice on the vegetation canopy (kg m-2)
                       scalarCanopyLiquid,      & ! intent(in): mass of liquid water on the vegetation canopy (kg m-2)
                       mLayerVolFracIce,        & ! intent(in): volumetric fraction of ice in each layer (-)
                       mLayerVolFracLiq,        & ! intent(in): volumetric fraction of liquid water in each layer (-)
                       ! input: coordinate variables
                       layerType,               & ! intent(in): type of the layer (snow or soil)
                       mLayerHeight,            & ! intent(in): height of the layer mid-point (top of soil = 0)
                       iLayerHeight,            & ! intent(in): height of the layer interface (top of soil = 0)
                       ! input: model parameters
                       heightCanopyTop,         & ! intent(in): height at the top of the veg canopy (m)
                       heightCanopyBottom,      & ! intent(in): height at the bottom of the veg canopy (m)
                       specificHeatVeg,         & ! intent(in): specific heat of vegetation (J kg-1 K-1)
                       maxMassVegetation,       & ! intent(in): maximum mass of vegetation (full foliage) (kg m-2)
                       iden_soil,               & ! intent(in): intrinsic density of soil (kg m-3)
                       thCond_soil,             & ! intent(in): thermal conductivity of soil (W m-1 K-1)
                       theta_sat,               & ! intent(in): soil porosity (-)
                       frac_sand,               & ! intent(in): fraction of sand (-)
                       frac_silt,               & ! intent(in): fraction of silt (-)
                       frac_clay,               & ! intent(in): fraction of clay (-)
                       ! output: diagnostic variables
                       scalarBulkVolHeatCapVeg, & ! intent(out): bulk volumetric heat capacity of vegetation (J m-3 K-1)
                       mLayerVolHtCapBulk,      & ! intent(out): volumetric heat capacity in each layer (J m-3 K-1)
                       mLayerThermalC,          & ! intent(out): thermal conductivity at the mid-point of each layer (W m-1 K-1)
                       iLayerThermalC,          & ! intent(out): thermal conductivity at the interface of each layer (W m-1 K-1)
                       mLayerVolFracAir,        & ! intent(out): volumetric fraction of air in each layer (-)
                       ! output: error control
                       err,message)               ! intent(out): error control
 ! --------------------------------------------------------------------------------------------------------------------------------
 USE snow_utils_module,only:tcond_snow            ! compute thermal conductivity of snow
 implicit none
 ! --------------------------------------------------------------------------------------------------------------------------------
 ! input: control variables
 logical(lgt),intent(in)       :: computeVegFlux          ! flag to indicate if we are computing fluxes over vegetation (.false. means veg is buried with snow)
 ! input: state variables
 real(dp),intent(in)           :: scalarCanopyIce         ! mass of ice on the vegetation canopy (kg m-2)
 real(dp),intent(in)           :: scalarCanopyLiquid      ! mass of liquid water on the vegetation canopy (kg m-2)
 real(dp),intent(in)           :: mLayerVolFracIce(:)     ! volumetric fraction of ice in each layer (-)
 real(dp),intent(in)           :: mLayerVolFracLiq(:)     ! volumetric fraction of liquid water in each layer (-)
 ! input: coordinate variables
 integer(i4b),intent(in)       :: layerType(:)            ! type of the layer (snow or soil)
 real(dp),intent(in)           :: mLayerHeight(:)         ! height of the layer mid-point (top of soil = 0)
 real(dp),intent(in)           :: iLayerHeight(0:)        ! height of the layer interface (top of soil = 0)
 ! input: model parameters
 real(dp),intent(in)           :: heightCanopyTop         ! height at the top of the veg canopy
 real(dp),intent(in)           :: heightCanopyBottom      ! height at the bottom of the veg canopy
 real(dp),intent(in)           :: specificHeatVeg         ! specific heat of vegetation (J kg-1 K-1)
 real(dp),intent(in)           :: maxMassVegetation       ! maximum mass of vegetation (full foliage) (kg m-2)
 real(dp),intent(in)           :: iden_soil               ! intrinsic density of soil (kg m-3)
 real(dp),intent(in)           :: thCond_soil             ! thermal conductivity of soil (W m-1 K-1)
 real(dp),intent(in)           :: theta_sat               ! soil porosity (-)
 real(dp),intent(in)           :: frac_sand               ! fraction of sand (-)
 real(dp),intent(in)           :: frac_silt               ! fraction of silt (-)
 real(dp),intent(in)           :: frac_clay               ! fraction of clay (-)
 ! --------------------------------------------------------------------------------------------------------------------------------
 ! output: diagnostic variables
 real(dp),intent(out)          :: scalarBulkVolHeatCapVeg ! bulk volumetric heat capacity of vegetation (J m-3 K-1)
 real(dp),intent(out)          :: mLayerVolHtCapBulk(:)   ! volumetric heat capacity in each layer (J m-3 K-1)
 real(dp),intent(out),target   :: mLayerThermalC(:)       ! thermal conductivity at the mid-point of each layer (W m-1 K-1)
 real(dp),intent(out)          :: iLayerThermalC(0:)      ! thermal conductivity at the interface of each layer (W m-1 K-1)
 real(dp),intent(out)          :: mLayerVolFracAir(:)     ! volumetric fraction of air in each layer (-)
 ! output: error control
 integer(i4b),intent(out)      :: err                     ! error code
 character(*),intent(out)      :: message                 ! error message
 ! --------------------------------------------------------------------------------------------------------------------------------
 ! --------------------------------------------------------------------------------------------------------------------------------
 ! local variables
 character(LEN=256)            :: cmessage               ! error message of downwind routine
 integer(i4b)                  :: nLayers                ! number of model layers
 integer(i4b)                  :: iLayer                 ! index of model layer
 real(dp),pointer              :: TCn                    ! thermal conductivity below the layer interface (W m-1 K-1)
 real(dp),pointer              :: TCp                    ! thermal conductivity above the layer interface (W m-1 K-1)
 real(dp)                      :: zdn                    ! height difference between interface and lower value (m)
 real(dp)                      :: zdp                    ! height difference between interface and upper value (m)
 real(dp)                      :: lambda_wet             ! thermal conductivity of the wet material
 !real(dp)                      :: kerstenNum             ! the Kersten number (-), defining weight applied to conductivity of the wet medium
 real(dp)                      :: bulkden_soil           ! bulk density of soil (kg m-3)
 real(dp)                      :: lambda_drysoil         ! thermal conductivity of dry soil (W m-1)
 real(dp)                      :: lambda_wetsoil         ! thermal conductivity of wet soil (W m-1)
 real(dp)                      :: canopyDepth            ! depth of the vegetation canopy (m)

 ! --------------------------------------------------------------------------------------------------------------------------------
 ! initialize error control
 err=0; message="diagn_evar/"

 ! compute the number of layers
 nLayers = size(mLayerHeight)

 ! define the canopy depth (m)
 canopyDepth = heightCanopyTop - heightCanopyBottom
 if(heightCanopyBottom > heightCanopyTop)then
  err=20; message=trim(message)//'height of the bottom of the canopy > top of the canopy'; return
 endif

 ! compute the bulk volumetric heat capacity of vegetation (J m-3 K-1)
 if(computeVegFlux)then
  scalarBulkVolHeatCapVeg = specificHeatVeg*maxMassVegetation/canopyDepth + & ! vegetation component
                            Cp_water*scalarCanopyLiquid/canopyDepth       + & ! liquid water component
                            Cp_ice*scalarCanopyIce/canopyDepth                ! ice component
 else
  scalarBulkVolHeatCapVeg = valueMissing
 endif 

 ! compute the thermal conductivity of dry and wet soils (W m-1)
 ! NOTE: this is actually constant over the simulation, and included here for clarity
 bulkden_soil   = iden_soil*(1._dp - theta_sat)
 lambda_drysoil = (0.135_dp*bulkden_soil + 64.7_dp) / (iden_soil - 0.947_dp*bulkden_soil)
 lambda_wetsoil = (8.80_dp*frac_sand + 2.92_dp*frac_clay) / (frac_sand + frac_clay)

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
   ! * soil
   case(ix_soil)
    mLayerThermalC(iLayer) = thCond_soil * (1._dp - theta_sat)      + & ! soil component
                             lambda_ice  * mLayerVolFracIce(iLayer) + & ! ice component
                             lambda_water* mLayerVolFracLiq(iLayer) + & ! liquid water component
                             lambda_air  * mLayerVolFracAir(iLayer)     ! air component
    ! compute the thermal conductivity of the wet material (W m-1)
    !lambda_wet = lambda_wetsoil**(1._dp - theta_sat) * lambda_water**theta_sat * lambda_ice**(theta_sat - mLayerVolFracLiq(iLayer))
    ! compute the Kersten number (-)
    !kerstenNum = log10( (mLayerVolFracIce(iLayer) + mLayerVolFracLiq(iLayer))/theta_sat ) + 1._dp
    ! ...and, compute the thermal conductivity
    !mLayerThermalC(iLayer) = kerstenNum*lambda_wet + (1._dp - kerstenNum)*lambda_drysoil
   ! * snow
   case(ix_snow)
    call tcond_snow(mLayerVolFracIce(iLayer)*iden_ice,mLayerThermalC(iLayer),err,cmessage)
    if(err/=0)then; message=trim(message)//trim(cmessage); return; endif
   ! * error check
   case default; err=20; message=trim(message)//'unable to identify type of layer (snow or soil) to compute thermal conductivity'; return
  end select

 end do  ! looping through layers

 ! *****
 ! * compute the thermal conductivity of snow at the interface of each layer...
 ! ****************************************************************************
 do iLayer=1,nLayers-1  ! (loop through layers)
  TCn => mLayerThermalC(iLayer)    ! thermal conductivity below the layer interface (W m-1 K-1)
  TCp => mLayerThermalC(iLayer+1)  ! thermal conductivity above the layer interface (W m-1 K-1)
  zdn =  iLayerHeight(iLayer)   - mLayerHeight(iLayer) ! height difference between interface and lower value (m)
  zdp =  mLayerHeight(iLayer+1) - iLayerHeight(iLayer) ! height difference between interface and upper value (m)
  iLayerThermalC(iLayer) = (TCn*TCp*(zdn + zdp)) / (TCn*zdp + TCp*zdn)
  !write(*,'(a,1x,i4,1x,10(f9.3,1x))') 'iLayer, TCn, TCp, zdn, zdp, iLayerThermalC(iLayer) = ', iLayer, TCn, TCp, zdn, zdp, iLayerThermalC(iLayer)
 end do
 ! assume the thermal conductivity at the domain boundaries is equal to the thermal conductivity of the layer
 iLayerThermalC(0)       = mLayerThermalC(1)
 iLayerThermalC(nLayers) = mLayerThermalC(nLayers)

 end subroutine diagn_evar

end module diagn_evar_module

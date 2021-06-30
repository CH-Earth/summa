
module computEnthalpy_module

! data types
USE nrtype

! derived types to define the data structures
USE data_types,only:&
                    var_ilength,  & ! data vector with variable length dimension (i4b)
                    var_dlength     ! data vector with variable length dimension (rkind)

! named variables
USE var_lookup,only:iLookPROG       ! named variables for structure elements
USE var_lookup,only:iLookDIAG       ! named variables for structure elements
USE var_lookup,only:iLookFLUX       ! named variables for structure elements
USE var_lookup,only:iLookINDEX      ! named variables for structure elements
USE var_lookup,only:iLookPARAM      ! named variables for structure elements

! access the global print flag
USE globalData,only:globalPrintFlag

! access missing values
USE globalData,only:integerMissing  ! missing integer
USE globalData,only:realMissing     ! missing real number
USE globalData,only:iname_snow      ! named variables for snow
USE globalData,only:iname_soil      ! named variables for soil

! constants
USE multiconst,only:&
                    Tfreeze,      & ! freezing temperature                 (K)
                    LH_fus,       & ! latent heat of fusion                (J kg-1)
                    LH_vap,       & ! latent heat of vaporization          (J kg-1)
                    iden_ice,     & ! intrinsic density of ice             (kg m-3)
                    iden_water,   & ! intrinsic density of liquid water    (kg m-3)
                    iden_air,     &
                    ! specific heat
                    Cp_air,      & ! specific heat of air          (J kg-1 K-1)
                    Cp_ice,      & ! specific heat of ice          (J kg-1 K-1)
                    Cp_soil,     & ! specific heat of soil         (J kg-1 K-1)
                    Cp_water       ! specific heat of liquid water (J kg-1 K-1)
! privacy
implicit none
private
public::computEnthalpy
public::computEnthalpyPrime
contains
 
  ! **********************************************************************************************************
 ! public subroutine computEnthalpy
 ! **********************************************************************************************************
 subroutine computEnthalpy(&
                        ! input
                        indx_data,                    &
                        nLayers,                      &
                        mLayerTemp,                   &
                        mLayerVolFracIce,             &
                        mLayerHeatCap,                  & 
                        ! output
                        mLayerEnthalpy                &
                        )                
 ! --------------------------------------------------------------------------------------------------------------------------------
 implicit none
 ! input: model control
 type(var_ilength),intent(in)    :: indx_data                 ! indices defining model states and layers
 integer(i4b),intent(in)         :: nLayers                     ! number of snow layers
 real(rkind),intent(in)             :: mLayerTemp(:)             ! temperature of each snow/soil layer (K)
 real(rkind),intent(in)             :: mLayerVolFracIce(:)       ! volumetric fraction of ice (-)
 real(rkind),intent(in)             :: mLayerHeatCap(:)
 real(rkind),intent(out)            :: mLayerEnthalpy(:)
 
 ! local variables
 integer(i4b)                    :: iLayer

 ! --------------------------------------------------------------------------------------------------------------------------------
 
 associate(&
   layerType               => indx_data%var(iLookINDEX%layerType)%dat               ,& ! intent(in): [i4b(:)] named variables defining the type of layer
   ixSnowSoilNrg           => indx_data%var(iLookINDEX%ixSnowSoilNrg)%dat            & ! intent(in): [i4b(:)] indices for energy states
 )
  ! (loop through non-missing energy state variables in the snow+soil domain)
  do concurrent (iLayer=1:nLayers,ixSnowSoilNrg(iLayer)/=integerMissing)   
   select case( layerType(iLayer) )
    case(iname_snow)
         mLayerEnthalpy(iLayer) = mLayerHeatCap(iLayer)*mLayerTemp(iLayer) - LH_fus*iden_ice * mLayerVolFracIce(iLayer)
    case(iname_soil)
         mLayerEnthalpy(iLayer) = mLayerHeatCap(iLayer)*mLayerTemp(iLayer) - LH_fus*iden_water * mLayerVolFracIce(iLayer)  
   end select
  end do  ! looping through non-missing energy state variables in the snow+soil domain
  
 end associate

 end subroutine computEnthalpy
 
   ! **********************************************************************************************************
 ! public subroutine computEnthalpyPrime
 ! **********************************************************************************************************
 subroutine computEnthalpyPrime(&
                        ! input
                        computeVegFlux,				  &
                        indx_data,                    &
                        nLayers,                      &
                        canopyDepth,               	  & ! intent(in): canopy depth (m)
                        scalarCanopyTempPrime,        & ! intent(in):    Prime value for the temperature of the vegetation canopy (K)
                        scalarCanopyIcePrime,         & ! intent(in):    Prime value for the ice on the vegetation canopy (kg m-2)
                        mLayerTempPrime,              &
                        mLayerVolFracIcePrime,        &
                        heatCapVeg,					  &
                        mLayerHeatCap,                & 
                        ! output
                        scalarCanopyEnthalpyPrime,	  &
                        mLayerEnthalpyPrime           &
                        )                
 ! --------------------------------------------------------------------------------------------------------------------------------
 implicit none
 ! input: model control
 logical(lgt),intent(in)         :: computeVegFlux         		! logical flag to denote if computing the vegetation flux
 type(var_ilength),intent(in)    :: indx_data                 	! indices defining model states and layers
 integer(i4b),intent(in)         :: nLayers                     ! number of snow layers
 real(rkind),intent(in)			 :: canopyDepth					! canopy depth (m)
 real(rkind),intent(in)			 :: scalarCanopyTempPrime       ! Prime value for the temperature of the vegetation canopy (K)
 real(rkind),intent(in)			 :: scalarCanopyIcePrime		! Prime value for the ice on the vegetation canopy (kg m-2)
 real(rkind),intent(in)			 :: heatCapVeg
 real(rkind),intent(in)             :: mLayerTempPrime(:)          ! temperature of each snow/soil layer (K)
 real(rkind),intent(in)             :: mLayerVolFracIcePrime(:)    ! volumetric fraction of ice (-)
 real(rkind),intent(in)             :: mLayerHeatCap(:)
 real(rkind),intent(out)			 :: scalarCanopyEnthalpyPrime
 real(rkind),intent(out)            :: mLayerEnthalpyPrime(:)
 
 ! local variables
 integer(i4b)                    :: iLayer

 ! --------------------------------------------------------------------------------------------------------------------------------
 
 associate(&
   layerType               => indx_data%var(iLookINDEX%layerType)%dat               ,& ! intent(in): [i4b(:)] named variables defining the type of layer
   ixSnowSoilNrg           => indx_data%var(iLookINDEX%ixSnowSoilNrg)%dat            & ! intent(in): [i4b(:)] indices for energy states
 )
 
 if(computeVegFlux)then
	scalarCanopyEnthalpyPrime = heatCapVeg * scalarCanopyTempPrime - LH_fus*scalarCanopyIcePrime/canopyDepth
 end if
  ! (loop through non-missing energy state variables in the snow+soil domain)
  do concurrent (iLayer=1:nLayers,ixSnowSoilNrg(iLayer)/=integerMissing)   
   select case( layerType(iLayer) )
    case(iname_snow)
         mLayerEnthalpyPrime(iLayer) = mLayerHeatCap(iLayer)*mLayerTempPrime(iLayer) - LH_fus*iden_ice * mLayerVolFracIcePrime(iLayer)
    case(iname_soil)
         mLayerEnthalpyPrime(iLayer) = mLayerHeatCap(iLayer)*mLayerTempPrime(iLayer) - LH_fus*iden_water * mLayerVolFracIcePrime(iLayer)  
   end select
  end do  ! looping through non-missing energy state variables in the snow+soil domain
  
 end associate

 end subroutine computEnthalpyPrime

end module computEnthalpy_module

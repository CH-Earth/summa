module pOverwrite_module
USE nrtype
implicit none
private
public::pOverwrite
contains

 ! ************************************************************************************************
 ! (1) new subroutine: use Noah tables to overwrite default model parameters
 ! ************************************************************************************************
 subroutine pOverwrite(err,message)
 ! FUSE data structures
 USE data_struc,only:parFallback            ! data structures for default values and constraints for model parameters
 USE data_struc,only:type_data              ! data structures for local classification of veg, soil, etc.
 USE var_lookup,only:iLookTYPE,iLookPARAM   ! named variables for elements of the data structures
 ! Noah table dimensions
 USE module_sf_noahlsm, only: LUCATS        ! dimension of the vegetation tables (number of land use catagories)
 USE module_sf_noahlsm, only: NSLTYPE       ! dimension of the soil tables
 ! Noah vegetation tables
 USE NOAHMP_VEG_PARAMETERS, only: Z0MVT     ! Noah-MP: momentum roughness length (m)
 USE NOAHMP_VEG_PARAMETERS, only: HVT       ! Noah-MP: height at top of canopy (m)
 USE NOAHMP_VEG_PARAMETERS, only: DLEAF     ! Noah-MP: characteristic leaf dimension (m)
 ! Noah soil tables
 USE module_sf_noahlsm, only: theta_res, theta_sat, vGn_alpha, vGn_n, k_soil  ! van Genutchen soil parameters 
 USE module_sf_noahlsm, only: REFSMC        ! Noah-MP: reference volumetric soil moisture content (-)
 USE module_sf_noahlsm, only: WLTSMC        ! Noah-MP: volumetric soil moisture content when plants are wilting (-)
 implicit none
 ! define output
 integer(i4b),intent(out)             :: err         ! error code
 character(*),intent(out)             :: message     ! error message
 ! define general variables
 integer(i4b)                         :: ixVeg       ! veg type index
 integer(i4b)                         :: ixSoil      ! soil type index
 ! Start procedure here
 err=0; message="pOverwrite/"

 ! define vegetation class
 ixVeg = type_data%var(iLookTYPE%vegTypeIndex)
 if(ixVeg < 1)then; err=20; message=trim(message)//'index for vegetation type must be > 0'; return; endif
 if(ixVeg > LUCATS)then; err=20; message=trim(message)//'index for vegetation type is greater than dimension of vegetation table'; return; endif

 ! define soil class
 ixSoil = type_data%var(iLookTYPE%soilTypeIndex)
 if(ixSoil < 1)then; err=20; message=trim(message)//'index for soil type must be > 0'; return; endif
 if(ixSoil > NSLTYPE)then; err=20; message=trim(message)//'index for soil type is greater than dimension of soil table'; return; endif

 ! include parameters from the vegetation tables
 parFallback(iLookPARAM%canopyHeight)%default_val        = HVT(ixVeg)          ! Noah-MP: height at top of canopy (m)
 parFallback(iLookPARAM%z0Canopy)%default_val            = Z0MVT(ixVeg)        ! Noah-MP: momentum roughness length (m)
 parFallback(iLookPARAM%leafDimension)%default_val       = DLEAF(ixVeg)        ! Noah-MP: characteristic leaf dimension (m)

 ! include parameters from the soil tables
 parFallback(iLookPARAM%k_soil)%default_val              = k_soil(ixSoil)      ! hydraulic conductivity (m s-1)
 parFallback(iLookPARAM%theta_res)%default_val           = theta_res(ixSoil)   ! residual volumetric liquid water content (-)
 parFallback(iLookPARAM%theta_sat)%default_val           = theta_sat(ixSoil)   ! soil porosity (-)
 parFallback(iLookPARAM%vGn_alpha)%default_val           = vGn_alpha(ixSoil)   ! van Genutchen "alpha" parameter (m-1)
 parFallback(iLookPARAM%vGn_n)%default_val               = vGn_n(ixSoil)       ! van Genutchen "n" parameter (-)
 parFallback(iLookPARAM%critSoilTranspire)%default_val   = REFSMC(ixSoil)      ! Noah-MP: reference volumetric soil moisture content (-)
 parFallback(iLookPARAM%critSoilWilting)%default_val     = WLTSMC(ixSoil)      ! Noah-MP: volumetric soil moisture content when plants are wilting (-)

 end subroutine pOverwrite

end module pOverwrite_module

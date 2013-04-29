module pOverwrite_module
USE nrtype
implicit none
private
public::pOverwrite
contains

 ! ************************************************************************************************
 ! (1) new subroutine: use Noah tables to overwrite default model parameters
 ! ************************************************************************************************
 subroutine pOverwrite(ixVeg,ixSoil,err,message)
 ! FUSE data structures
 USE data_struc,only:localParFallback       ! data structures for default values and constraints for model parameters
 USE var_lookup,only:iLookPARAM             ! named variables for elements of the data structures
 ! Noah table dimensions
 USE module_sf_noahlsm, only: LUCATS        ! dimension of the vegetation tables (number of land use catagories)
 USE module_sf_noahlsm, only: NSLTYPE       ! dimension of the soil tables
 ! Noah vegetation tables
 USE NOAHMP_VEG_PARAMETERS, only: Z0MVT     ! Noah-MP: momentum roughness length (m)
 USE NOAHMP_VEG_PARAMETERS, only: HVT       ! Noah-MP: height at top of canopy (m)
 USE NOAHMP_VEG_PARAMETERS, only: HVB       ! Noah-MP: height at bottom of canopy (m)
 USE NOAHMP_VEG_PARAMETERS, only: DLEAF     ! Noah-MP: characteristic leaf dimension (m)
 ! Noah soil tables
 USE module_sf_noahlsm, only: theta_res, theta_sat, vGn_alpha, vGn_n, k_soil  ! van Genutchen soil parameters 
 USE module_sf_noahlsm, only: REFSMC        ! Noah-MP: reference volumetric soil moisture content (-)
 USE module_sf_noahlsm, only: WLTSMC        ! Noah-MP: volumetric soil moisture content when plants are wilting (-)
 implicit none
 ! define input
 integer(i4b),intent(in)              :: ixVeg       ! vegetation category
 integer(i4b),intent(in)              :: ixSoil      ! soil category
 ! define output
 integer(i4b),intent(out)             :: err         ! error code
 character(*),intent(out)             :: message     ! error message
 ! Start procedure here
 err=0; message="pOverwrite/"

 ! define vegetation class
 if(ixVeg < 1)then; err=20; message=trim(message)//'index for vegetation type must be > 0'; return; endif
 if(ixVeg > LUCATS)then
  write(message,'(2(a,i0),a)')trim(message)//'index for vegetation type is greater than dimension of vegetation table [ixVeg = ', ixVeg, &
                            '; LUCATS = ', LUCATS, ']'
  err=20; return
 endif

 ! define soil class
 if(ixSoil < 1)then; err=20; message=trim(message)//'index for soil type must be > 0'; return; endif
 if(ixSoil > NSLTYPE)then
  write(message,'(2(a,i0),a)')trim(message)//'index for soil type is greater than dimension of soil table [ixSoil = ', ixSoil, &
                            '; NSLTYPE = ', NSLTYPE, ']'
  err=20; return
 endif

 ! include parameters from the vegetation tables
 localParFallback(iLookPARAM%heightCanopyTop)%default_val     = HVT(ixVeg)          ! Noah-MP: height at top of canopy (m)
 localParFallback(iLookPARAM%heightCanopyBottom)%default_val  = HVB(ixVeg)          ! Noah-MP: height at bottom of canopy (m)
 localParFallback(iLookPARAM%z0Canopy)%default_val            = Z0MVT(ixVeg)        ! Noah-MP: momentum roughness length (m)
 localParFallback(iLookPARAM%leafDimension)%default_val       = DLEAF(ixVeg)        ! Noah-MP: characteristic leaf dimension (m)

 ! include parameters from the soil tables
 localParFallback(iLookPARAM%k_soil)%default_val              = k_soil(ixSoil)      ! hydraulic conductivity (m s-1)
 localParFallback(iLookPARAM%theta_res)%default_val           = theta_res(ixSoil)   ! residual volumetric liquid water content (-)
 localParFallback(iLookPARAM%theta_sat)%default_val           = theta_sat(ixSoil)   ! soil porosity (-)
 localParFallback(iLookPARAM%vGn_alpha)%default_val           = vGn_alpha(ixSoil)   ! van Genutchen "alpha" parameter (m-1)
 localParFallback(iLookPARAM%vGn_n)%default_val               = vGn_n(ixSoil)       ! van Genutchen "n" parameter (-)
 localParFallback(iLookPARAM%critSoilTranspire)%default_val   = REFSMC(ixSoil)      ! Noah-MP: reference volumetric soil moisture content (-)
 localParFallback(iLookPARAM%critSoilWilting)%default_val     = WLTSMC(ixSoil)      ! Noah-MP: volumetric soil moisture content when plants are wilting (-)

 end subroutine pOverwrite

end module pOverwrite_module

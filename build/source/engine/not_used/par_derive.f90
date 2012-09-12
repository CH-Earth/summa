module par_derive_module
USE nrtype
implicit none
contains

 ! **********************************************************************************************************
 ! new subroutine: compute derived model parameters
 ! **********************************************************************************************************
 subroutine par_derive(err,message)
 ! used to compute derived model parameters
 USE multiconst, only: vkc                                       ! von Karman's constant
 USE var_lookup, only: iLookPARAM                                ! named variables that refer to elements of data structure
 USE data_struc, only: mpar_data                                 ! model parameters
 implicit none
 ! declare dummy variables
 integer(i4b),intent(out) :: err            ! error code
 character(*),intent(out) :: message        ! error message
 ! declare local variables
 real(dp),pointer         :: zon            ! roughness length (m)
 real(dp),pointer         :: mheight        ! measurement height (m)
 real(dp),pointer         :: bparam         ! parameter in Louis (1979) stability function
 real(dp),pointer         :: c_star         ! parameter in Louis (1979) stability function
 real(dp),pointer         :: ExNeut         ! exchange coefficient in neutral conditions
 real(dp),pointer         :: bprime         ! used in Louis (1979) stability function
 real(dp),pointer         :: cparam         ! used in Louis (1979) stability function
 ! initialize error control
 err=0; message='f-par_derive/OK'
 ! assign local pointers to the values in the data structures
 zon    =>mpar_data%var(iLookPARAM%zon)     ! roughness length (m)
 mheight=>mpar_data%var(iLookPARAM%mheight) ! measurement height (m)
 bparam =>mpar_data%var(iLookPARAM%bparam)  ! parameter in Louis (1979) stability function
 c_star =>mpar_data%var(iLookPARAM%c_star)  ! parameter in Louis (1979) stability function
 ExNeut =>mpar_data%var(iLookPARAM%ExNeut)  ! exchange coefficient in neutral conditions
 bprime =>mpar_data%var(iLookPARAM%bprime)  ! used in Louis (1979) stability function
 cparam =>mpar_data%var(iLookPARAM%cparam)  ! used in Louis (1979) stability function
 ! compute derived parameters for turbulent heat transfer
 ExNeut = (vkc**2._dp) / (log((mheight/zon))**2._dp)           ! exchange coefficient in neutral conditions
 bprime = bparam / 2._dp                                       ! used in Louis (1979) stability function
 cparam = c_star * ExNeut * bparam * ((mheight/zon)**0.5_dp)   ! used in Louis (1979) stability function
 ! compute volumetric heat capacity
 iden_soil    = soil_dens_bulk / (1._dp - theta_sat)
 volHtCap_dry = iden_soil  * Cp_soil
 volHtCap_air = iden_air   * Cp_air
 volHtCap_ice = iden_ice   * Cp_Ice
 volHtCap_wat = iden_water * Cp_water

 ! compute the volumetric latent heat of fusion (J m-3 K-1)
 volLatHt_fus = iden_ice   * LH_fus


 end subroutine par_derive

end module par_derive_module

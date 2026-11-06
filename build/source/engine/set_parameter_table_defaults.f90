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

module parameter_table_defaults_module
USE nr_type
implicit none
private
public::set_parameter_table_defaults
contains


 ! ************************************************************************************************
 ! public subroutine set_parameter_table_defaults: use Noah tables to overwrite default model parameters
 ! ************************************************************************************************
 subroutine set_parameter_table_defaults(ixVeg,ixSoil,soilCatTbl, defaultParam,err,message)
 ! SUMMA dictionary
 USE var_lookup,only:iLookPARAM             ! named variables for elements of the data structures
 ! Noah table dimensions
 USE module_sf_noahlsm, only: LUCATS        ! dimension of the vegetation tables (number of land use catagories)
 USE module_sf_noahlsm, only: NSLTYPE       ! dimension of the soil tables
 ! Noah vegetation tables
 USE NOAHMP_VEG_PARAMETERS, only: Z0MVT     ! Noah-MP: momentum roughness length (m)
 USE NOAHMP_VEG_PARAMETERS, only: HVT       ! Noah-MP: height at top of canopy (m)
 USE NOAHMP_VEG_PARAMETERS, only: HVB       ! Noah-MP: height at bottom of canopy (m)
 USE NOAHMP_VEG_PARAMETERS, only: DLEAF     ! Noah-MP: characteristic leaf dimension (m)
 USE NOAHMP_VEG_PARAMETERS, only: VCMX25    ! Noah-MP: maximum Rubisco carboxylation rate (umol m-2 s-1)
 USE NOAHMP_VEG_PARAMETERS, only: MP        ! Noah-MP: slope of conductance-photosynthesis relationship (-)
 ! Noah soil tables
 USE module_sf_noahlsm, only: theta_res, theta_sat, vGn_alpha, vGn_n, k_soil  ! van Genutchen soil parameters, only ROSETTA
 USE module_sf_noahlsm, only: REFSMC        ! Noah-MP: reference volumetric soil moisture content (-)
 USE module_sf_noahlsm, only: WLTSMC        ! Noah-MP: volumetric soil moisture content when plants are wilting (-)
 USE module_sf_noahlsm, only: SATDK         ! Noah-MP: saturated hydraulic conductivity (m s-1)
 USE module_sf_noahlsm, only: DRYSMC        ! Noah-MP: dry soil moisture content (-)
 USE module_sf_noahlsm, only: MAXSMC        ! Noah-MP: soil porosity (-)
 USE module_sf_noahlsm, only: SATPSI        ! Noah-MP: saturated soil matric potential (m)
 USE module_sf_noahlsm, only: BB            ! Noah-MP: Clapp-Hornberger exponent (-)
 implicit none
 ! define input
 integer(i4b),intent(in)              :: ixVeg           ! vegetation category
 integer(i4b),intent(in)              :: ixSoil          ! soil category
 character(*), intent(in)             :: soilCatTbl      ! soil category table name
 ! define output
 real(rkind),intent(inout)            :: defaultParam(:) ! default model parameters
 integer(i4b),intent(out)             :: err             ! error code
 character(*),intent(out)             :: message         ! error message
 ! local variables
 real(rkind)                          :: p               ! temporary variable for van Genuchten parameter calculation
 ! Start procedure here
 err=0; message="set_parameter_table_defaults/"

 ! define vegetation class
 if(ixVeg < 1)then; err=20; message=trim(message)//'index for vegetation type must be > 0'; return; end if
 if(ixVeg > LUCATS)then
  write(message,'(2(a,i0),a)')trim(message)//'index for vegetation type is greater than dimension of vegetation table [ixVeg = ', ixVeg, &
                            '; LUCATS = ', LUCATS, ']'
  err=20; return
 end if

 ! define soil class
 if(ixSoil < 1)then; err=20; message=trim(message)//'index for soil type must be > 0'; return; end if
 if(ixSoil > NSLTYPE)then
  write(message,'(2(a,i0),a)')trim(message)//'index for soil type is greater than dimension of soil table [ixSoil = ', ixSoil, &
                            '; NSLTYPE = ', NSLTYPE, ']'
  err=20; return
 end if

 ! include parameters from the vegetation tables
 defaultParam(iLookPARAM%heightCanopyTop)     = HVT(ixVeg)          ! Noah-MP: height at top of canopy (m)
 defaultParam(iLookPARAM%heightCanopyBottom)  = HVB(ixVeg)          ! Noah-MP: height at bottom of canopy (m)
 defaultParam(iLookPARAM%z0Canopy)            = Z0MVT(ixVeg)        ! Noah-MP: momentum roughness length (m)
 defaultParam(iLookPARAM%leafDimension)       = DLEAF(ixVeg)        ! Noah-MP: characteristic leaf dimension (m)
 defaultParam(iLookPARAM%vcmax25_canopyTop)   = VCMX25(ixVeg)       ! Noah-MP: maximum Rubisco carboxylation rate (umol m-2 s-1)
 defaultParam(iLookPARAM%cond2photo_slope)    = MP(ixVeg)           ! Noah-MP: slope of conductance-photosynthesis relationship (-)

 ! include parameters from the soil tables
 select case(trim(soilCatTbl))
  case('ROSETTA')
   ! use van Genutchen parameters only available for ROSETTA soil categories
   defaultParam(iLookPARAM%k_soil)          = k_soil(ixSoil)          ! hydraulic conductivity (m s-1)
   defaultParam(iLookPARAM%theta_res)       = theta_res(ixSoil)       ! residual volumetric liquid water content (-)
   defaultParam(iLookPARAM%theta_sat)       = theta_sat(ixSoil)       ! soil porosity (-)
   defaultParam(iLookPARAM%vGn_alpha)       = vGn_alpha(ixSoil)       ! van Genutchen "alpha" parameter (m-1)
   defaultParam(iLookPARAM%vGn_n)           = vGn_n(ixSoil)           ! van Genutchen "n" parameter (-)
  case('STAS','STAS-RUC')
   defaultParam(iLookPARAM%k_soil)          = SATDK(ixSoil)           ! Noah-MP: saturated hydraulic conductivity (m s-1)
   defaultParam(iLookPARAM%theta_res)       = DRYSMC(ixSoil)          ! Noah-MP: minumum soil moisture content used as residual volumetric liquid water content (-)
   defaultParam(iLookPARAM%theta_sat)       = MAXSMC(ixSoil)          ! Noah-MP: soil porosity (-)
   ! following Morel-Seytouxet al., 1996, use empirical relations to estimate van Genuchten parameters from Clapp-Hornberger parameters
   defaultParam(iLookPARAM%vGn_n)          = (BB(ixSoil)+1._rkind)/BB(ixSoil)        ! Noah-MP: van Genuchten n parameter (-)
   p = 2.0_rkind*BB(ixSoil)+3.0_rkind
   defaultParam(iLookPARAM%vGn_alpha)      = -1._rkind/SATPSI(ixSoil) * (p + 3._rkind)/(2._rkind*p* (p - 1._rkind)) * &
                                             (147.8_rkind + 8.1_rkind*p + 0.092_rkind * p**2_i4b) / (55.6_rkind + 7.4_rkind*p + p**2_i4b) ! Noah-MP: van Genuchten alpha parameter (m-1)
 end select
 defaultParam(iLookPARAM%critSoilTranspire) = REFSMC(ixSoil)          ! Noah-MP: reference volumetric soil moisture content (-)
 defaultParam(iLookPARAM%critSoilWilting)   = WLTSMC(ixSoil)          ! Noah-MP: volumetric soil moisture content when plants are wilting (-)

 end subroutine set_parameter_table_defaults


end module parameter_table_defaults_module

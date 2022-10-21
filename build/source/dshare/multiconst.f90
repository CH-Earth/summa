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

MODULE multiconst
 USE nrtype
 ! define physical constants
 real(rkind), PARAMETER           :: ave_slp        =  101325.0_rkind ! mean sea level pressure              (Pa)
 real(rkind), PARAMETER           :: vkc            =  0.4_rkind      ! von Karman constant                  (-)
 real(rkind), PARAMETER           :: satvpfrz       =  610.8_rkind    ! sat vapour pressure at 273.16K       (Pa)
 real(rkind), PARAMETER           :: w_ratio        =  0.622_rkind    ! molecular ratio water to dry air     (-)
 real(rkind), PARAMETER           :: R_da           =  287.053_rkind  ! gas constant for dry air             (Pa K-1 m3 kg-1; J kg-1 K-1)
 real(rkind), PARAMETER           :: R_wv           = 461.285_rkind   ! gas constant for water vapor         (Pa K-1 m3 kg-1; J kg-1 K-1)
 real(rkind), PARAMETER           :: Rgas           = 8.314_rkind     ! universal gas constant               (J mol-1 K-1)
 real(rkind), PARAMETER           :: gravity        = 9.80616_rkind   ! acceleration of gravity              (m s-2)
 real(rkind), PARAMETER           :: Cp_air         = 1005._rkind     ! specific heat of air                 (J kg-1 K-1)
 real(rkind), PARAMETER           :: Cp_ice         = 2114._rkind     ! specific heat of ice                 (J kg-1 K-1)
 real(rkind), PARAMETER           :: Cp_soil        = 850._rkind      ! specific heat of soil                (J kg-1 K-1)
 real(rkind), PARAMETER           :: Cp_water       = 4181._rkind     ! specific heat of liquid water        (J kg-1 K-1)
 real(rkind), PARAMETER           :: Tfreeze        = 273.16_rkind    ! temperature at freezing              (K)
 real(rkind), PARAMETER           :: TriplPt        = 273.16_rkind    ! triple point of water                (K)
 real(rkind), PARAMETER           :: LH_fus         = 333700.0_rkind  ! latent heat of fusion                (J kg-1)
 real(rkind), PARAMETER           :: LH_vap         = 2501000.0_rkind ! latent heat of vaporization          (J kg-1)
 real(rkind), PARAMETER           :: LH_sub         = 2834700.0_rkind ! latent heat of sublimation           (J kg-1)
 real(rkind), PARAMETER           :: sb             = 5.6705d-8    ! Stefan Boltzman constant             (W m-2 K-4)
 real(rkind), PARAMETER           :: em_sno         = 0.99_rkind      ! emissivity of snow                   (-)
 real(rkind), PARAMETER           :: lambda_air     = 0.026_rkind     ! thermal conductivity of air          (W m-1 K-1)
 real(rkind), PARAMETER           :: lambda_ice     = 2.50_rkind      ! thermal conductivity of ice          (W m-1 K-1)
 real(rkind), PARAMETER           :: lambda_water   = 0.60_rkind      ! thermal conductivity of liquid water (W m-1 K-1)
 real(rkind), PARAMETER           :: iden_air       = 1.293_rkind     ! intrinsic density of air             (kg m-3)
 real(rkind), PARAMETER           :: iden_ice       = 917.0_rkind     ! intrinsic density of ice             (kg m-3)
 real(rkind), PARAMETER           :: iden_water     = 1000.0_rkind    ! intrinsic density of liquid water    (kg m-3)
 real(rkind), PARAMETER           :: secprday       = 86400._rkind    ! number of seconds in a day
 real(rkind), PARAMETER           :: secprhour      = 3600._rkind     ! number of seconds in an hour
 real(rkind), PARAMETER           :: secprmin       = 60._rkind       ! number of seconds in a minute
 real(rkind), PARAMETER           :: minprhour      = 60._rkind       ! number of minutes in an hour
END MODULE multiconst

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
 real(rk), PARAMETER           :: ave_slp        =  101325.0_rk ! mean sea level pressure              (Pa)
 real(rk), PARAMETER           :: vkc            =  0.4_rk      ! von Karman constant                  (-)
 real(rk), PARAMETER           :: satvpfrz       =  610.8_rk    ! sat vapour pressure at 273.16K       (Pa)
 real(rk), PARAMETER           :: w_ratio        =  0.622_rk    ! molecular ratio water to dry air     (-)
 real(rk), PARAMETER           :: R_da           =  287.053_rk  ! gas constant for dry air             (Pa K-1 m3 kg-1; J kg-1 K-1)
 real(rk), PARAMETER           :: R_wv           = 461.285_rk   ! gas constant for water vapor         (Pa K-1 m3 kg-1; J kg-1 K-1)
 real(rk), PARAMETER           :: Rgas           = 8.314_rk     ! universal gas constant               (J mol-1 K-1)
 real(rk), PARAMETER           :: gravity        = 9.80616_rk   ! acceleration of gravity              (m s-2)
 real(rk), PARAMETER           :: Cp_air         = 1005._rk     ! specific heat of air                 (J kg-1 K-1)
 real(rk), PARAMETER           :: Cp_ice         = 2114._rk     ! specific heat of ice                 (J kg-1 K-1)
 real(rk), PARAMETER           :: Cp_soil        = 850._rk      ! specific heat of soil                (J kg-1 K-1)
 real(rk), PARAMETER           :: Cp_water       = 4181._rk     ! specific heat of liquid water        (J kg-1 K-1)
 real(rk), PARAMETER           :: Tfreeze        = 273.16_rk    ! temperature at freezing              (K)
 real(rk), PARAMETER           :: TriplPt        = 273.16_rk    ! triple point of water                (K)
 real(rk), PARAMETER           :: LH_fus         = 333700.0_rk  ! latent heat of fusion                (J kg-1)
 real(rk), PARAMETER           :: LH_vap         = 2501000.0_rk ! latent heat of vaporization          (J kg-1)
 real(rk), PARAMETER           :: LH_sub         = 2834700.0_rk ! latent heat of sublimation           (J kg-1)
 real(rk), PARAMETER           :: sb             = 5.6705d-8    ! Stefan Boltzman constant             (W m-2 K-4)
 real(rk), PARAMETER           :: em_sno         = 0.99_rk      ! emissivity of snow                   (-)
 real(rk), PARAMETER           :: lambda_air     = 0.026_rk     ! thermal conductivity of air          (W m-1 K-1)
 real(rk), PARAMETER           :: lambda_ice     = 2.50_rk      ! thermal conductivity of ice          (W m-1 K-1)
 real(rk), PARAMETER           :: lambda_water   = 0.60_rk      ! thermal conductivity of liquid water (W m-1 K-1)
 real(rk), PARAMETER           :: iden_air       = 1.293_rk     ! intrinsic density of air             (kg m-3)
 real(rk), PARAMETER           :: iden_ice       = 917.0_rk     ! intrinsic density of ice             (kg m-3)
 real(rk), PARAMETER           :: iden_water     = 1000.0_rk    ! intrinsic density of liquid water    (kg m-3)
 real(rk), PARAMETER           :: secprday       = 86400._rk    ! number of seconds in a day
 real(rk), PARAMETER           :: secprhour      = 3600._rk     ! number of seconds in an hour
 real(rk), PARAMETER           :: secprmin       = 60._rk       ! number of seconds in a minute
 real(rk), PARAMETER           :: minprhour      = 60._rk       ! number of minutes in an hour
END MODULE multiconst

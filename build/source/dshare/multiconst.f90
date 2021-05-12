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
 real(summa_prec), PARAMETER           :: ave_slp        =  101325.0_summa_prec ! mean sea level pressure              (Pa)
 real(summa_prec), PARAMETER           :: vkc            =  0.4_summa_prec      ! von Karman constant                  (-)
 real(summa_prec), PARAMETER           :: satvpfrz       =  610.8_summa_prec    ! sat vapour pressure at 273.16K       (Pa)
 real(summa_prec), PARAMETER           :: w_ratio        =  0.622_summa_prec    ! molecular ratio water to dry air     (-)
 real(summa_prec), PARAMETER           :: R_da           =  287.053_summa_prec  ! gas constant for dry air             (Pa K-1 m3 kg-1; J kg-1 K-1)
 real(summa_prec), PARAMETER           :: R_wv           = 461.285_summa_prec   ! gas constant for water vapor         (Pa K-1 m3 kg-1; J kg-1 K-1)
 real(summa_prec), PARAMETER           :: Rgas           = 8.314_summa_prec     ! universal gas constant               (J mol-1 K-1)
 real(summa_prec), PARAMETER           :: gravity        = 9.80616_summa_prec   ! acceleration of gravity              (m s-2)
 real(summa_prec), PARAMETER           :: Cp_air         = 1005._summa_prec     ! specific heat of air                 (J kg-1 K-1)
 real(summa_prec), PARAMETER           :: Cp_ice         = 2114._summa_prec     ! specific heat of ice                 (J kg-1 K-1)
 real(summa_prec), PARAMETER           :: Cp_soil        = 850._summa_prec      ! specific heat of soil                (J kg-1 K-1)
 real(summa_prec), PARAMETER           :: Cp_water       = 4181._summa_prec     ! specific heat of liquid water        (J kg-1 K-1)
 real(summa_prec), PARAMETER           :: Tfreeze        = 273.16_summa_prec    ! temperature at freezing              (K)
 real(summa_prec), PARAMETER           :: TriplPt        = 273.16_summa_prec    ! triple point of water                (K)
 real(summa_prec), PARAMETER           :: LH_fus         = 333700.0_summa_prec  ! latent heat of fusion                (J kg-1)
 real(summa_prec), PARAMETER           :: LH_vap         = 2501000.0_summa_prec ! latent heat of vaporization          (J kg-1)
 real(summa_prec), PARAMETER           :: LH_sub         = 2834700.0_summa_prec ! latent heat of sublimation           (J kg-1)
 real(summa_prec), PARAMETER           :: sb             = 5.6705d-8    ! Stefan Boltzman constant             (W m-2 K-4)
 real(summa_prec), PARAMETER           :: em_sno         = 0.99_summa_prec      ! emissivity of snow                   (-)
 real(summa_prec), PARAMETER           :: lambda_air     = 0.026_summa_prec     ! thermal conductivity of air          (W m-1 K-1)
 real(summa_prec), PARAMETER           :: lambda_ice     = 2.50_summa_prec      ! thermal conductivity of ice          (W m-1 K-1)
 real(summa_prec), PARAMETER           :: lambda_water   = 0.60_summa_prec      ! thermal conductivity of liquid water (W m-1 K-1)
 real(summa_prec), PARAMETER           :: iden_air       = 1.293_summa_prec     ! intrinsic density of air             (kg m-3)
 real(summa_prec), PARAMETER           :: iden_ice       = 917.0_summa_prec     ! intrinsic density of ice             (kg m-3)
 real(summa_prec), PARAMETER           :: iden_water     = 1000.0_summa_prec    ! intrinsic density of liquid water    (kg m-3)
 real(summa_prec), PARAMETER           :: secprday       = 86400._summa_prec    ! number of seconds in a day
 real(summa_prec), PARAMETER           :: secprhour      = 3600._summa_prec     ! number of seconds in an hour
 real(summa_prec), PARAMETER           :: secprmin       = 60._summa_prec       ! number of seconds in a minute
 real(summa_prec), PARAMETER           :: minprhour      = 60._summa_prec       ! number of minutes in an hour
END MODULE multiconst

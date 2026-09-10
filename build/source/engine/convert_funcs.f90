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

module convert_funcs_module
USE nr_type                                ! variable types
USE multiconst                             ! fixed parameters (lh vapzn, etc.)
implicit none
private
public::RELHM2SPHM,SPHM2RELHM,WETBULBTMP,satVapPress,vapPress,getLatentHeatValue
public::MSLP2AIRP,AIRP2MSLP
contains

! ----------------------------------------------------------------------
! series of functions to convert one thing to another
! (partially courtesy of Drew Slater)
! ----------------------------------------------------------------------

! ***************************************************************************************************************
! public function getLatentHeatValue: get appropriate latent heat of sublimation/vaporization for a given surface
! ***************************************************************************************************************
function getLatentHeatValue(T)
implicit none
real(rkind),intent(in)   :: T                    ! temperature (K)
real(rkind)              :: getLatentHeatValue   ! latent heat of sublimation/vaporization (J kg-1)
!---------------------------------------------------------------------------------------------------
if(T > Tfreeze)then
 getLatentHeatValue = LH_vap     ! latent heat of vaporization          (J kg-1)
else
 getLatentHeatValue = LH_sub     ! latent heat of sublimation           (J kg-1)
end if
end function getLatentHeatValue

! ***************************************************************************************************************
! public function vapPress: convert specific humidity (g g-1) to vapor pressure (Pa)
! ***************************************************************************************************************
function vapPress(q,p)
implicit none
real(rkind),intent(in)   :: q        ! specific humidity (g g-1)
real(rkind),intent(in)   :: p        ! pressure (Pa)
real(rkind)              :: vapPress ! vapor pressure (Pa)
real(rkind)              :: w        ! mixing ratio
!---------------------------------------------------------------------------------------------------
w = q / (1._rkind - q)             ! mixing ratio (-)
vapPress = (w/(w + w_ratio))*p     ! vapor pressure (Pa)
end function vapPress

! ***************************************************************************************************************
! public subroutine satVapPress: Uses Teten's formula to compute saturated vapor pressure (Pa)
! ***************************************************************************************************************
! NOTE: temperature units are degC !!!!
! ***************************************************************************************************************
subroutine satVapPress(TC, SVP, dSVP_dT)
implicit none
real(rkind), intent(in)            :: TC       ! temperature (C)
real(rkind), intent(out)           :: SVP      ! saturation vapor pressure (Pa)
real(rkind), intent(out)           :: dSVP_dT  ! d(SVP)/dT
real(rkind), parameter             :: X1 = 17.27_rkind
real(rkind), parameter             :: X2 = 237.30_rkind
!---------------------------------------------------------------------------------------------------
! Units note :              Pa = N m-2 = kg m-1 s-2
! SATVPFRZ=     610.8       ! Saturation water vapour pressure at 273.16K (Pa)
if(X2 + TC <= 0.0_rkind)then ! will fail if divide by 0, but will blow up if negative top and bottom of fraction
 SVP     = tiny(1.0_rkind)
 dSVP_dT = tiny(1.0_rkind)
else
 SVP     = SATVPFRZ * EXP( (X1*TC)/(X2 + TC) ) ! Saturated Vapour Press (Pa)
 dSVP_dT = SVP * (X1/(X2 + TC) - X1*TC/(X2 + TC)**2_i4b)
end if
end subroutine satVapPress

! ***************************************************************************************************************
! public function MSLP2AIRP: compute air pressure using mean sea level pressure and elevation
! ***************************************************************************************************************
! (after Shuttleworth, 1993)
!
! -- actually returns MSLP2AIRP in the same units as MSLP, because
!    ( (293.-0.0065*ELEV) / 293. )**5.256 is dimensionless
! ***************************************************************************************************************
function MSLP2AIRP(MSLP, ELEV)
implicit none
real(rkind),intent(in)        :: MSLP      ! base pressure (Pa)
real(rkind),intent(in)        :: ELEV      ! elevation difference from base (m)
real(rkind)                   :: MSLP2AIRP ! Air pressure (Pa)
!---------------------------------------------------------------------------------------------------
MSLP2AIRP = MSLP * ( (293.-0.0065*ELEV) / 293. )**5.256
end function MSLP2AIRP

! ***************************************************************************************************************
! public function AIRP2MSLP: compute mean sea level pressure using air pressure and elevation
! ***************************************************************************************************************
! (after Shuttleworth, 1993)
!
! -- actually returns AIRP2MSLP in the same units as AIRP, because
!    ( (293.-0.0065*ELEV) / 293. )**5.256 is dimensionless
! ***************************************************************************************************************
function AIRP2MSLP(AIRP, ELEV)
implicit none
real(rkind),intent(in)        :: AIRP      ! air pressure (Pa)
real(rkind),intent(in)        :: ELEV      ! elevation difference from base (m)
real(rkind)                   :: AIRP2MSLP ! base pressure (Pa)
!---------------------------------------------------------------------------------------------------
AIRP2MSLP = AIRP / ( (293.-0.0065*ELEV) / 293. )**5.256
end function AIRP2MSLP

! ***************************************************************************************************************
! private function RLHUM2DEWPT: compute dewpoint temperature from relative humidity
! ***************************************************************************************************************
! ---- This is done with respect to water ONLY ----
!
! All units are SI standard - i.e. Kelvin and pascals
! Based on Tetens' formula (1930)
! ***************************************************************************************************************
function RLHUM2DEWPT(T, RLHUM)
implicit none
real(rkind),intent(in)        :: T           ! Temperature           (K)
real(rkind),intent(in)        :: RLHUM       ! Relative Humidity     (%)
real(rkind)                   :: RLHUM2DEWPT ! Dewpoint Temp   (K)
real(rkind)                   :: VPSAT       ! Sat. vapour pressure at T (Pa)
real(rkind)                   :: TDCEL       ! Dewpoint temp Celcius (C)
!---------------------------------------------------------------------------------------------------
! Units note :              Pa = N m-2 = kg m-1 s-2
! SATVPFRZ=     610.8       ! Saturation water vapour pressure at 273.16K (Pa)
! W_RATIO =       0.622     ! molecular weight ratio of water to dry air (-)
VPSAT = SATVPFRZ * EXP( (17.27*(T-TFREEZE)) / (237.30 + (T-TFREEZE)) ) ! sat vapor press at grid cell (Pa)
TDCEL = 237.30 * LOG( (VPSAT/SATVPFRZ)*(RLHUM/100.) ) / &              ! dewpoint temperature         (C)
        (17.27 - LOG( (VPSAT/SATVPFRZ)*(RLHUM/100.) ) )
RLHUM2DEWPT = TDCEL + TFREEZE
end function RLHUM2DEWPT

! ***************************************************************************************************************
! private function DEWPT2RLHUM: compute relative humidity from dewpoint temperature
! ***************************************************************************************************************
! ---- This is done with respect to water ONLY ----
!
! All units are SI standard - i.e. Kelvin and pascals
! Based on Tetens' formula (1930)
! ***************************************************************************************************************
function DEWPT2RLHUM(T, DEWPT)
implicit none
real(rkind),intent(in)        :: T           ! Temperature           (K)
real(rkind),intent(in)        :: DEWPT       ! Dewpoint temp         (K)
real(rkind)                   :: DEWPT2RLHUM ! Relative Humidity   (%)
real(rkind)                   :: VPSAT       ! Sat. vapour pressure at T (Pa)
real(rkind)                   :: TDCEL       ! Dewpt in celcius      (C)
!---------------------------------------------------------------------------------------------------
! Units note :              Pa = N m-2 = kg m-1 s-2
! SATVPFRZ=     610.8       ! Saturation water vapour pressure at 273.16K (Pa)
TDCEL = DEWPT-TFREEZE
VPSAT = SATVPFRZ * EXP( (17.27*(T-TFREEZE)) / (237.30 + (T-TFREEZE)) )      ! Sat vapor press (Pa)
DEWPT2RLHUM = 100. * (SATVPFRZ/VPSAT) * EXP((17.27*TDCEL)/(237.30+TDCEL))   ! Relative Humidity (%)
end function DEWPT2RLHUM

! ***************************************************************************************************************
! private function DEWPT2SPHM: compute specific humidity from dewpoint temperature
! ***************************************************************************************************************
! ---- This is done with respect to water ONLY ----
!
! All units are SI standard - i.e. Kelvin and pascals
! Based on Tetens' formula (1930)
! VPAIR is the current vapor pressure as it used dewpoint to compute staurated VP
! ***************************************************************************************************************
function DEWPT2SPHM(DEWPT, PRESS)
implicit none
real(rkind),intent(in)        :: DEWPT      ! Dewpoint temp         (K)
real(rkind),intent(in)        :: PRESS      ! Pressure              (Pa)
real(rkind)                   :: DEWPT2SPHM ! Specific Humidity    (g/g)
real(rkind)                   :: VPAIR      ! vapour pressure at T  (Pa)
real(rkind)                   :: TDCEL      ! Dewpt in celcius      (C)
!---------------------------------------------------------------------------------------------------
! Units note :              Pa = N m-2 = kg m-1 s-2
! SATVPFRZ=     610.8       ! Saturation water vapour pressure at 273.16K (Pa)
TDCEL = DEWPT-TFREEZE
VPAIR = SATVPFRZ * EXP( (17.27*TDCEL) / (237.30 + TDCEL) )        ! Vapour Press           (Pa)
DEWPT2SPHM = (VPAIR * W_RATIO)/(PRESS - (1.-W_RATIO)*VPAIR)       ! Specific humidity (g/g)
end function DEWPT2SPHM

! ***************************************************************************************************************
! private function DEWPT2VPAIR: compute vapor pressure of air from dewpoint temperature
! ***************************************************************************************************************
! ---- This is done with respect to water ONLY ----
!
! All units are SI standard - i.e. Kelvin and pascals
! Based on Tetens' formula (1930)
! VPAIR is the current vapor pressure as it used dewpoint to compute saturated VP
! ***************************************************************************************************************
function DEWPT2VPAIR(DEWPT)
implicit none
real(rkind),intent(in)        :: DEWPT       ! Dewpoint temp         (K)
real(rkind)                   :: TDCEL       ! Dewpt in celcius      (C)
real(rkind)                   :: DEWPT2VPAIR ! Vapour Press  (Pa)
!---------------------------------------------------------------------------------------------------
! Units note :              Pa = N m-2 = kg m-1 s-2
! SATVPFRZ=     610.8       ! Saturation water vapour pressure at 273.16K (Pa)
TDCEL = DEWPT-TFREEZE
DEWPT2VPAIR = SATVPFRZ * EXP( (17.27*TDCEL) / (237.30 + TDCEL) )   ! Vapour Press  (Pa)
end function DEWPT2VPAIR

! ***************************************************************************************************************
! public function SPHM2RELHM: compute relative humidity from specific humidity
! ***************************************************************************************************************
! ---- This is done with respect to water ONLY ----
!
! All units are SI standard - i.e. Kelvin and pascals
! Based on Tetens' formula (1930)
! ***************************************************************************************************************
function SPHM2RELHM(SPHM, PRESS, TAIR)
implicit none
real(rkind),intent(in)        :: SPHM       ! Specific Humidity (g/g)
real(rkind),intent(in)        :: PRESS      ! Pressure              (Pa)
real(rkind),intent(in)        :: TAIR       ! Air temp
real(rkind)                   :: SPHM2RELHM ! Dewpoint Temp (K)
real(rkind)                   :: VPSAT      ! vapour pressure at T  (Pa)
real(rkind)                   :: TDCEL      ! Dewpt in celcius      (C)
!---------------------------------------------------------------------------------------------------
! Units note :              Pa = N m-2 = kg m-1 s-2
! SATVPFRZ=     610.8       ! Saturation water vapour pressure at 273.16K (Pa)
TDCEL = TAIR-TFREEZE
VPSAT = SATVPFRZ * EXP( (17.27*TDCEL) / (237.30 + TDCEL) )       ! Vapour Press      (Pa)
SPHM2RELHM = (SPHM * PRESS)/(VPSAT * (W_RATIO + SPHM*(1.-W_RATIO)))
end function SPHM2RELHM

! ***************************************************************************************************************
! public function RELHM2SPHM: compute specific humidity from relative humidity
! ***************************************************************************************************************
! ---- This is done with respect to water ONLY ----
!
! All units are SI standard - i.e. Kelvin and pascals
! Based on Tetens' formula (1930)
! ***************************************************************************************************************
function RELHM2SPHM(RELHM, PRESS, TAIR)
implicit none
real(rkind),intent(in)        :: RELHM      ! Relative Humidity     (%)
real(rkind),intent(in)        :: PRESS      ! Pressure              (Pa)
real(rkind),intent(in)        :: TAIR       ! Air temp
real(rkind)                   :: RELHM2SPHM ! Specific Humidity (g/g)
real(rkind)                   :: PVP        ! Partial vapour pressure at T  (Pa)
real(rkind)                   :: TDCEL      ! Dewpt in celcius      (C)
!---------------------------------------------------------------------------------------------------
! Units note :              Pa = N m-2 = kg m-1 s-2
! SATVPFRZ=     610.8       ! Saturation water vapour pressure at 273.16K (Pa)
TDCEL = TAIR-TFREEZE
PVP = RELHM * SATVPFRZ * EXP( (17.27*TDCEL)/(237.30 + TDCEL) ) ! Partial Vapour Press (Pa)
RELHM2SPHM = (PVP * W_RATIO)/(PRESS - (1. - W_RATIO)*PVP)
end function RELHM2SPHM

! ***************************************************************************************************************
! public function WETBULBTMP: compute wet bulb temperature based on humidity and pressure
! ***************************************************************************************************************
function WETBULBTMP(TAIR, RELHM, PRESS)
implicit none
real(rkind),intent(in)        :: TAIR           ! Air temp              (K)
real(rkind),intent(in)        :: RELHM          ! Relative Humidity     (-)
real(rkind),intent(in)        :: PRESS          ! Pressure              (Pa)
real(rkind)                   :: WETBULBTMP     ! Wet bulb temperature (K)
real(rkind)                   :: Tcel           ! Temperature in celcius      (C)
real(rkind)                   :: PVP            ! Partial vapor pressure (Pa)
real(rkind)                   :: TWcel          ! Wet bulb temperature in celcius (C)
real(rkind),parameter         :: k=6.54E-4_DP   ! normalizing factor in wet bulb estimate (C-1)
real(rkind)                   :: Twet_trial0    ! trial value for wet bulb temperature (C)
real(rkind)                   :: Twet_trial1    ! trial value for wet bulb temperature (C)
real(rkind)                   :: f0,f1          ! function evaluations (C)
real(rkind)                   :: df_dT          ! derivative (-)
real(rkind)                   :: TWinc          ! wet bulb temperature increment (C)
INTEGER(I4B)                  :: iter           ! iterattion index
real(rkind),parameter         :: Xoff=1.E-5_DP  ! finite difference increment (C)
real(rkind),parameter         :: Xtol=1.E-8_DP  ! convergence tolerance (C)
INTEGER(I4B)                  :: maxiter=15     ! maximum number of iterations
!---------------------------------------------------------------------------------------------------
! convert temperature to Celcius
Tcel = TAIR-TFREEZE
! compute partial vapor pressure based on temperature (Pa)
PVP = RELHM * SATVPRESS(Tcel)
! define an initial trial value for wetbulb temperature
TWcel = Tcel - 5._rkind
! iterate until convergence
do iter=1,maxiter
 ! compute Twet estimates
 Twet_trial0 = Tcel - (SATVPRESS(TWcel)      - PVP)/(k*PRESS)
 Twet_trial1 = Tcel - (SATVPRESS(TWcel+Xoff) - PVP)/(k*PRESS)
 ! compute function evaluations
 f0 = Twet_trial0 - TWcel
 f1 = Twet_trial1 - (TWcel+Xoff)
 ! compute derivative and iteration increment
 df_dT = (f0 - f1)/Xoff
 TWinc = f0/df_dT
 ! compute new value of wet bulb temperature (C)
 TWcel = TWcel + TWinc
 ! check if achieved tolerance
 if(abs(f0) < Xtol) exit
 ! check convergence
 if(iter==maxiter)stop 'failed to converge in WETBULBTMP'
end do  ! (iterating)
! return value in K
WETBULBTMP = TWcel + TFREEZE
end function WETBULBTMP

! ***************************************************************************************************************
! private function SATVPRESS: compute saturated vapor pressure (Pa)
! ***************************************************************************************************************
! Units note :              Pa = N m-2 = kg m-1 s-2
! SATVPFRZ=     610.8       ! Saturation water vapour pressure at 273.16K (Pa)
! ***************************************************************************************************************
function SATVPRESS(TCEL)
implicit none
real(rkind),intent(in) :: TCEL      ! Temperature (C)
real(rkind)            :: SATVPRESS ! Saturated vapor pressure (Pa)
!---------------------------------------------------------------------------------------------------
SATVPRESS = SATVPFRZ * EXP( (17.27_rkind*TCEL)/(237.30_rkind + TCEL) ) ! Saturated Vapour Press (Pa)
end function SATVPRESS


end module convert_funcs_module

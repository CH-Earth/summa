! SUMMA - Structure for Unifying Multiple Modeling Alternatives
! Copyright (C) 2014-2015 NCAR/RAL
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

module conv_funcs_module
USE nrtype                                 ! variable types
USE multiconst                             ! fixed parameters (lh vapzn, etc.)
implicit none
private
public::RELHM2SPHM,SPHM2RELHM,WETBULBTMP,satVapPress,vapPress,getLatentHeatValue
contains

! ----------------------------------------------------------------------
! series of functions to convert one thing to another
! (partially courtesy of Drew Slater)
! ----------------------------------------------------------------------

function getLatentHeatValue(T)
! get appropriate latent heat of sublimation/vaporization for a given surface
implicit none
real(dp),intent(in)   :: T                    ! temperature (K)
real(dp)              :: getLatentHeatValue   ! latent heat of sublimation/vaporization (J kg-1)
if(T > Tfreeze)then
 getLatentHeatValue = LH_vap     ! latent heat of vaporization          (J kg-1)
else
 getLatentHeatValue = LH_sub     ! latent heat of sublimation           (J kg-1)
endif
end function getLatentHeatValue


! ----------------------------------------------------------------------
! ----------------------------------------------------------------------
function vapPress(q,p)
! convert specific humidity (g g-1) to vapor pressure (Pa)
implicit none
! input
real(dp),intent(in)   :: q        ! specific humidity (g g-1)
real(dp),intent(in)   :: p        ! pressure (Pa)
! output
real(dp)              :: vapPress ! vapor pressure (Pa)
! local
real(dp)              :: w        ! mixing ratio
!real(dp),parameter    :: w_ratio = 0.622_dp ! molecular weight ratio of water to dry air (-) 
w = q / (1._dp - q)                ! mixing ratio (-)
vapPress = (w/(w + w_ratio))*p     ! vapor pressure (Pa)
end function vapPress

! ----------------------------------------------------------------------
! ----------------------------------------------------------------------
subroutine satVapPress_noah(TC, SVP, dSVP_dT)
!---------------------------------------------------------------------------------------------------
! Modified from Noah-MP
! Use polynomials to calculate saturation vapor pressure and derivative with
! respect to temperature: over water when t > 0 c and over ice when t <= 0 c
! NOTE: temperature units are degC !!!!
! NOTE: analytical derivatives do not match analytical derivatives
!         --> there may be a problem in one of the constants....
IMPLICIT NONE
! input
real(dp), intent(in)            :: TC       ! temperature (C)
! output
real(dp), intent(out)           :: SVP      ! saturation vapor pressure (Pa)
real(dp), intent(out)           :: dSVP_dT  ! d(SVP)/dT
! local
REAL(DP) :: A0,A1,A2,A3,A4,A5,A6  !coefficients for esat over water
REAL(DP) :: B0,B1,B2,B3,B4,B5,B6  !coefficients for esat over ice
REAL(DP) :: C0,C1,C2,C3,C4,C5,C6  !coefficients for dsat over water
REAL(DP) :: D0,D1,D2,D3,D4,D5,D6  !coefficients for dsat over ice
! local (testing)
real(dp), parameter             :: X1 = 17.27_dp
real(dp), parameter             :: X2 = 237.30_dp
! parameters
PARAMETER (A0=6.107799961_dp    , A1=4.436518521E-01_dp,  &
           A2=1.428945805E-02_dp, A3=2.650648471E-04_dp,  &
           A4=3.031240396E-06_dp, A5=2.034080948E-08_dp,  &
           A6=6.136820929E-11_dp)
PARAMETER (B0=6.109177956_dp    , B1=5.034698970E-01_dp,  &
           B2=1.886013408E-02_dp, B3=4.176223716E-04_dp,  &
           B4=5.824720280E-06_dp, B5=4.838803174E-08_dp,  &
           B6=1.838826904E-10_dp)
PARAMETER (C0= 4.438099984E-01_dp, C1=2.857002636E-02_dp,  &
           C2= 7.938054040E-04_dp, C3=1.215215065E-05_dp,  &
           C4= 1.036561403E-07_dp, C5=3.532421810e-10_dp,  &
           C6=-7.090244804E-13_dp)
PARAMETER (D0=5.030305237E-01_dp, D1=3.773255020E-02_dp,  &
           D2=1.267995369E-03_dp, D3=2.477563108E-05_dp,  &
           D4=3.005693132E-07_dp, D5=2.158542548E-09_dp,  &
           D6=7.131097725E-12_dp)
! water
if(TC > 0._dp)then
 SVP     = 100._dp*(A0+TC*(A1+TC*(A2+TC*(A3+TC*(A4+TC*(A5+TC*A6))))))
 dSVP_dT = 100._dp*(C0+TC*(C1+TC*(C2+TC*(C3+TC*(C4+TC*(C5+TC*C6))))))
! ice
else
 SVP     = 100._dp*(B0+TC*(B1+TC*(B2+TC*(B3+TC*(B4+TC*(B5+TC*B6))))))
 dSVP_dT = 100._dp*(D0+TC*(D1+TC*(D2+TC*(D3+TC*(D4+TC*(D5+TC*D6))))))
endif
print*, 'noah: SVP, dSVP_dT = ', SVP, dSVP_dT
! test
SVP     = SATVPFRZ * EXP( (X1*TC)/(X2 + TC) ) ! Saturated Vapour Press (Pa)
dSVP_dT = SVP * (X1/(X2 + TC) - X1*TC/(X2 + TC)**2._dp)
print*, 'test: SVP, dSVP_dT = ', SVP, dSVP_dT
END SUBROUTINE satVapPress_noah

! ----------------------------------------------------------------------
! ----------------------------------------------------------------------
subroutine satVapPress(TC, SVP, dSVP_dT)
!---------------------------------------------------------------------------------------------------
! Uses Teten's formula to compute saturated vapor pressure (Pa)
! NOTE: temperature units are degC !!!!
IMPLICIT NONE
! input
real(dp), intent(in)            :: TC       ! temperature (C)
! output
real(dp), intent(out)           :: SVP      ! saturation vapor pressure (Pa)
real(dp), intent(out)           :: dSVP_dT  ! d(SVP)/dT
! local
real(dp), parameter             :: X1 = 17.27_dp
real(dp), parameter             :: X2 = 237.30_dp
! local (use to test derivative calculations)
real(dp),parameter              :: dx = 1.e-8_dp  ! finite difference increment
!---------------------------------------------------------------------------------------------------
! Units note :              Pa = N m-2 = kg m-1 s-2
! SATVPFRZ=     610.8       ! Saturation water vapour pressure at 273.16K (Pa)
SVP     = SATVPFRZ * EXP( (X1*TC)/(X2 + TC) ) ! Saturated Vapour Press (Pa)
dSVP_dT = SVP * (X1/(X2 + TC) - X1*TC/(X2 + TC)**2._dp)
!print*, 'dSVP_dT check... ', SVP, dSVP_dT, (SATVPRESS(TC+dx) - SVP)/dx
END SUBROUTINE satVapPress

! ----------------------------------------------------------------------
! ----------------------------------------------------------------------
FUNCTION MSLP2AIRP(MSLP, ELEV)
! compute air pressure using mean sea level pressure and elevation
! (after Shuttleworth, 1993)
!
! -- actually returns MSLP2AIRP in the same units as MSLP, because
!    ( (293.-0.0065*ELEV) / 293. )**5.256 is dimensionless
!
IMPLICIT NONE

REAL(DP), INTENT(IN)         :: MSLP      ! base pressure (Pa)
REAL(DP), INTENT(IN)         :: ELEV      ! elevation difference from base (m)

REAL(DP)                     :: MSLP2AIRP ! Air pressure (Pa)

MSLP2AIRP = MSLP * ( (293.-0.0065*ELEV) / 293. )**5.256

END FUNCTION MSLP2AIRP

! ----------------------------------------------------------------------
! ----------------------------------------------------------------------
FUNCTION RLHUM2DEWPT(T, RLHUM)
! Compute Dewpoint temperature from Relative Humidity
! ---- This is done with respect to water ONLY ----
!
! All units are SI standard - i.e. Kelvin and pascals
! Based on Tetens' formula (1930)
IMPLICIT NONE

REAL(DP), INTENT(IN)         :: T         ! Temperature           (K)
REAL(DP), INTENT(IN)         :: RLHUM     ! Relative Humidity     (%)


REAL(DP)                     :: RLHUM2DEWPT     ! Dewpoint Temp   (K)

REAL(DP)                     :: VPSAT     ! Sat. vapour pressure at T (Pa)
REAL(DP)                     :: TDCEL     ! Dewpoint temp Celcius (C)

! Units note :              Pa = N m-2 = kg m-1 s-2
! SATVPFRZ=     610.8       ! Saturation water vapour pressure at 273.16K (Pa)
! W_RATIO =       0.622     ! molecular weight ratio of water to dry air (-)

VPSAT = SATVPFRZ * EXP( (17.27*(T-TFREEZE)) / (237.30 + (T-TFREEZE)) ) ! sat vapor press at grid cell (Pa)
TDCEL = 237.30 * LOG( (VPSAT/SATVPFRZ)*(RLHUM/100.) ) / &              ! dewpoint temperature         (C)
        (17.27 - LOG( (VPSAT/SATVPFRZ)*(RLHUM/100.) ) ) 
RLHUM2DEWPT = TDCEL + TFREEZE

END FUNCTION RLHUM2DEWPT

! ----------------------------------------------------------------------
! ----------------------------------------------------------------------
FUNCTION DEWPT2RLHUM(T, DEWPT)
! Compute Relative humidity from dewpoint temperature
! ---- This is done with respect to water ONLY ----
!
! All units are SI standard - i.e. Kelvin and pascals
! Based on Tetens' formula (1930)
IMPLICIT NONE

REAL(DP), INTENT(IN)         :: T         ! Temperature           (K)
REAL(DP), INTENT(IN)         :: DEWPT     ! Dewpoint temp         (K)

REAL(DP)                     :: DEWPT2RLHUM ! Relative Humidity   (%)

REAL(DP)                     :: VPSAT     ! Sat. vapour pressure at T (Pa)
REAL(DP)                     :: TDCEL     ! Dewpt in celcius      (C)

! Units note :              Pa = N m-2 = kg m-1 s-2
! SATVPFRZ=     610.8       ! Saturation water vapour pressure at 273.16K (Pa)

TDCEL = DEWPT-TFREEZE
VPSAT = SATVPFRZ * EXP( (17.27*(T-TFREEZE)) / (237.30 + (T-TFREEZE)) )      ! Sat vapor press (Pa)
DEWPT2RLHUM = 100. * (SATVPFRZ/VPSAT) * EXP((17.27*TDCEL)/(237.30+TDCEL))   ! Relative Humidity (%)

END FUNCTION DEWPT2RLHUM


! ----------------------------------------------------------------------
! ----------------------------------------------------------------------
FUNCTION DEWPT2SPHM(DEWPT, PRESS)
! Compute specific humidity from dewpoint temp with respect to water
! ---- This is done with respect to water ONLY ----
!
! All units are SI standard - i.e. Kelvin and pascals
! Based on Tetens' formula (1930)
! VPAIR is the current vapor pressure as it used dewpoint to compute staurated VP
IMPLICIT NONE

REAL(DP), INTENT(IN)         :: DEWPT     ! Dewpoint temp         (K)
REAL(DP), INTENT(IN)         :: PRESS     ! Pressure              (Pa)

REAL(DP)                     :: DEWPT2SPHM ! Specific Humidity    (g/g)

REAL(DP)                     :: VPAIR     ! vapour pressure at T  (Pa)
REAL(DP)                     :: TDCEL     ! Dewpt in celcius      (C)

! Units note :              Pa = N m-2 = kg m-1 s-2
! SATVPFRZ=     610.8       ! Saturation water vapour pressure at 273.16K (Pa)

TDCEL = DEWPT-TFREEZE
VPAIR = SATVPFRZ * EXP( (17.27*TDCEL) / (237.30 + TDCEL) )        ! Vapour Press           (Pa)
DEWPT2SPHM = (VPAIR * W_RATIO)/(PRESS - (1.-W_RATIO)*VPAIR)       ! Specific humidity (g/g)

END FUNCTION DEWPT2SPHM

! ----------------------------------------------------------------------
! ----------------------------------------------------------------------
FUNCTION DEWPT2VPAIR(DEWPT)
! Compute vapor pressure of the air from dewpoint temp with respect to water
! ---- This is done with respect to water ONLY ----
!
! All units are SI standard - i.e. Kelvin and pascals
! Based on Tetens' formula (1930)
! VPAIR is the current vapor pressure as it used dewpoint to compute staurated VP
IMPLICIT NONE

REAL(DP), INTENT(IN)         :: DEWPT     ! Dewpoint temp         (K)
REAL(DP)                     :: TDCEL     ! Dewpt in celcius      (C)

REAL(DP)                     :: DEWPT2VPAIR ! Vapour Press  (Pa)

! Units note :              Pa = N m-2 = kg m-1 s-2
! SATVPFRZ=     610.8       ! Saturation water vapour pressure at 273.16K (Pa)

TDCEL = DEWPT-TFREEZE
DEWPT2VPAIR = SATVPFRZ * EXP( (17.27*TDCEL) / (237.30 + TDCEL) )   ! Vapour Press  (Pa)

END FUNCTION DEWPT2VPAIR
! ----------------------------------------------------------------------
! ----------------------------------------------------------------------
FUNCTION SPHM2RELHM(SPHM, PRESS, TAIR)
! Compute specific humidity from dewpoint temp with respect to water
! ---- This is done with respect to water ONLY ----
!
! All units are SI standard - i.e. Kelvin and pascals
! Based on Tetens' formula (1930)
! VPAIR is the current vapor pressure as it used dewpoint to compute staurated VP
IMPLICIT NONE

REAL(DP), INTENT(IN)         :: SPHM      ! Specific Humidity (g/g)
REAL(DP), INTENT(IN)         :: PRESS     ! Pressure              (Pa)
REAL(DP), INTENT(IN)         :: TAIR      ! Air temp

REAL(DP)                     :: SPHM2RELHM ! Dewpoint Temp (K)

REAL(DP)                     :: VPSAT     ! vapour pressure at T  (Pa)
REAL(DP)                     :: TDCEL     ! Dewpt in celcius      (C)
!REAL(DP)                     :: DUM       ! Intermediate

! Units note :              Pa = N m-2 = kg m-1 s-2
! SATVPFRZ=     610.8       ! Saturation water vapour pressure at 273.16K (Pa)

TDCEL = TAIR-TFREEZE
VPSAT = SATVPFRZ * EXP( (17.27*TDCEL) / (237.30 + TDCEL) )       ! Vapour Press      (Pa)
SPHM2RELHM = (SPHM * PRESS)/(VPSAT * (W_RATIO + SPHM*(1.-W_RATIO)))

END FUNCTION SPHM2RELHM
! ----------------------------------------------------------------------
! ----------------------------------------------------------------------
FUNCTION RELHM2SPHM(RELHM, PRESS, TAIR)
! Compute specific humidity from dewpoint temp with respect to water
! ---- This is done with respect to water ONLY ----
!
! All units are SI standard - i.e. Kelvin and pascals
! Based on Tetens' formula (1930)
! VPAIR is the current vapor pressure as it used dewpoint to compute staurated VP
IMPLICIT NONE

REAL(DP), INTENT(IN)         :: RELHM     ! Relative Humidity     (%)
REAL(DP), INTENT(IN)         :: PRESS     ! Pressure              (Pa)
REAL(DP), INTENT(IN)         :: TAIR      ! Air temp

REAL(DP)                     :: RELHM2SPHM ! Specific Humidity (g/g)

REAL(DP)                     :: PVP       ! Partial vapour pressure at T  (Pa)
REAL(DP)                     :: TDCEL     ! Dewpt in celcius      (C)
!REAL(DP)                     :: DUM       ! Intermediate

! Units note :              Pa = N m-2 = kg m-1 s-2
! SATVPFRZ=     610.8       ! Saturation water vapour pressure at 273.16K (Pa)

TDCEL = TAIR-TFREEZE
PVP = RELHM * SATVPFRZ * EXP( (17.27*TDCEL)/(237.30 + TDCEL) ) ! Partial Vapour Press (Pa)
RELHM2SPHM = (PVP * W_RATIO)/(PRESS - (1. - W_RATIO)*PVP)

END FUNCTION RELHM2SPHM

! ----------------------------------------------------------------------
! ----------------------------------------------------------------------
FUNCTION WETBULBTMP(TAIR, RELHM, PRESS)
! Compute wet bulb temperature based on humidity and pressure
IMPLICIT NONE
! input
REAL(DP), INTENT(IN)         :: TAIR      ! Air temp              (K)
REAL(DP), INTENT(IN)         :: RELHM     ! Relative Humidity     (-)
REAL(DP), INTENT(IN)         :: PRESS     ! Pressure              (Pa)
! output
REAL(DP)                     :: WETBULBTMP ! Wet bulb temperature (K)
! locals
REAL(DP)                     :: Tcel           ! Temperature in celcius      (C)
REAL(DP)                     :: PVP            ! Partial vapor pressure (Pa)
REAL(DP)                     :: TWcel          ! Wet bulb temperature in celcius (C)
REAL(DP),PARAMETER           :: k=6.54E-4_DP   ! normalizing factor in wet bulb estimate (C-1)
REAL(DP)                     :: Twet_trial0    ! trial value for wet bulb temperature (C) 
REAL(DP)                     :: Twet_trial1    ! trial value for wet bulb temperature (C) 
REAL(DP)                     :: f0,f1          ! function evaluations (C)
REAL(DP)                     :: df_dT          ! derivative (-)
REAL(DP)                     :: TWinc          ! wet bulb temperature increment (C)
INTEGER(I4B)                 :: iter           ! iterattion index
REAL(DP),PARAMETER           :: Xoff=1.E-5_DP  ! finite difference increment (C)
REAL(DP),PARAMETER           :: Xtol=1.E-8_DP  ! convergence tolerance (C)
INTEGER(I4B)                 :: maxiter=15     ! maximum number of iterations 
! convert temperature to Celcius
Tcel = TAIR-TFREEZE
! compute partial vapor pressure based on temperature (Pa)
PVP = RELHM * SATVPRESS(Tcel)
! define an initial trial value for wetbulb temperature
TWcel = Tcel - 5._dp
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
enddo  ! (iterating)

! return value in K
WETBULBTMP = TWcel + TFREEZE

END FUNCTION WETBULBTMP


! ----------------------------------------------------------------------
! ----------------------------------------------------------------------
FUNCTION SATVPRESS(TCEL)
! Compute saturated vapor pressure (Pa)
! Units note :              Pa = N m-2 = kg m-1 s-2
! SATVPFRZ=     610.8       ! Saturation water vapour pressure at 273.16K (Pa)
IMPLICIT NONE
REAL(DP),INTENT(IN) :: TCEL      ! Temperature (C)
REAL(DP)            :: SATVPRESS ! Saturated vapor pressure (Pa)
SATVPRESS = SATVPFRZ * EXP( (17.27_dp*TCEL)/(237.30_dp + TCEL) ) ! Saturated Vapour Press (Pa)
END FUNCTION SATVPRESS


end module conv_funcs_module

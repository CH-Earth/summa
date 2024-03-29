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

MODULE sunGeomtry_module
USE nrtype
implicit none
private
public::clrsky_rad
contains

! *************************************************************************************************
! public subroutine CLRSKY_RAD: get hourly radiation index
! *************************************************************************************************
 SUBROUTINE CLRSKY_RAD(MONTH,DAY,HOUR,DT,SLOPE,AZI,LAT,HRI,COSZEN)
 ! ----------------------------------------------------------------------------------------
 ! Used to get hourly radiation index
 !
 ! Modification history
 !  - comments added by David Rupp 2006.
 !
 ! Procedure appears similar to Stull(1988) "An Introduction to Boundary Layer
 ! Meteorology " as seen in matlab routine obtained from Joe Kidston <joekidston@yahoo.co.uk>.
 ! Note that equation of time is not used.  Also, solar declination is assumed to stay
 ! constant over the period of one day. What other assumptions are made?  Is this
 ! adequate for the level of precision we require?  Worth reading Stull(1988).
 !
 !  - Modified to integrate over time step up to, but not greater than, 24 hours (D. Rupp, July 2006)
 !
 ! ----------------------------------------------------------------------------------------
 IMPLICIT NONE
 ! Input variables
 INTEGER(I4B), INTENT(IN)                  :: MONTH   ! month as mm integer
 INTEGER(I4B), INTENT(IN)                  :: DAY     ! day of month as dd integer
 real(rkind), INTENT(IN)                      :: HOUR    ! hour of day as real
 real(rkind), INTENT(IN)                      :: DT      ! time step in units of hours
 real(rkind), INTENT(IN)                      :: SLOPE   ! slope of ground surface in degrees
 real(rkind), INTENT(IN)                      :: AZI     ! aspect (azimuth) of ground surface in degrees
 real(rkind), INTENT(IN)                      :: LAT     ! latitude in degrees (negative for southern hemisphere)
 ! Outputs
 real(rkind), INTENT(OUT)                     :: HRI     ! average radiation index over time step DT
 real(rkind), INTENT(OUT)                     :: COSZEN  ! average cosine of the zenith angle over time step DT
 ! Internal
 real(rkind)                                  :: CRAD    ! conversion from degrees to radians
 real(rkind)                                  :: YRAD    ! conversion from year to radians
 real(rkind)                                  :: T       ! time from noon in radians
 real(rkind)                                  :: DELT1   ! time step in radians
 real(rkind)                                  :: SLOPE1  ! slope of ground surface in radians
 real(rkind)                                  :: AZI1    ! aspect (azimuth) of ground surface in radians
 real(rkind)                                  :: LAT1    ! latitude in radians
 real(rkind)                                  :: FJULIAN ! julian date as real
 real(rkind)                                  :: D       ! solar declination
 real(rkind)                                  :: LP      ! latitude adjusted for non-level surface (= LAT1 for level surface)
 real(rkind)                                  :: TD      ! used to calculate sunrise/set
 real(rkind)                                  :: TPI     ! used to calculate sunrise/set
 real(rkind)                                  :: TP      ! used to calculate sunrise/set
 real(rkind)                                  :: DDT     ! used to calculate sunrise/set(= 0 for level surface)
 real(rkind)                                  :: T1      ! first time in time step or sunrise
 real(rkind)                                  :: T2      ! last time in time step or sunset
 real(rkind)                                  :: AUX     ! Auxiliary variable used to check whether the sunset/sunrise time calculation can succeed
 ! ----------------------------------------------------------------------------------------
 ! CONVERSION FACTORS
 !   degrees to radians
 CRAD=PI_D/180.0D0
 !   days-of-year to radians
 YRAD=2.0D0*PI_D/365.0D0
 ! CONVERT TIME TO RADIANS FROM NOON
 T=(HOUR-12.0)*PI_D/12.0D0
 ! Convert time step to radians
 DELT1=DT*PI_D/12.0D0
 ! CONVERT ground slope, ground aspect, and latitude TO RADIANS
 SLOPE1=SLOPE*CRAD  ! tilt angle
 AZI1=AZI*CRAD ! surface-solar Azimuth ??
 LAT1=LAT*CRAD ! latitude
 ! Calculate julian date
 FJULIAN=dble(JULIAN(MONTH,DAY))
 ! Calculate solar declination
 D=CRAD*23.5*SIN((FJULIAN-82.0)*YRAD)
 ! Calculate latitude "adjustment" for ground slope, aspect and latitude (LP = LAT1 for level surface)
 LP=ASIN(SIN(SLOPE1)*COS(AZI1)*COS(LAT1) + COS(SLOPE1)*SIN(LAT1)) ! angle between solar rays and surface (tilted) ??
 ! Calculate time of sunrise/sunset on level surface as radians from noon
 ! Account for high latitude locations, where there might not be a sunrise/sunset time on a given day
 ! In such cases AUX > 1 or AUX < -1. Fix AUX at (-)1 in those cases, to fix sunrise at 00.00 or 24.00 of the current day (instead of some time before/after the current day)
 AUX=-TAN(LAT1)*TAN(D)
 IF(abs(AUX) > 1.) THEN
  TD=ACOS(SIGN(1._rkind, AUX))
 ELSE
  TD=ACOS(AUX)
 END IF
 ! print *, 'Sunrise = ', TD
 ! Calculate time of sunrise/sunset adjusted for inclined ground surface as radians from noon???
 TPI=-TAN(LP)*TAN(D)
 IF(ABS(TPI).LT.1.0) THEN
  TP=ACOS(TPI)
 ELSE IF(TPI.LT.-1.0) THEN ! 24h daylight
  TP=ACOS(-1.0)
 ELSE IF(TPI.GT.1.0) THEN ! 24h dark
  TP=ACOS(1.0)
 ENDIF
 ! Calculate time adjustment for ground slope, aspect and latitude (DDT = 0 for level surface)
 DDT=ATAN(SIN(AZI1)*SIN(SLOPE1)/(COS(SLOPE1)*COS(LAT1)-COS(AZI1)*SIN(SLOPE1)*SIN(LAT1)))
 ! print*, 'ddt = ', ddt
 ! Set beginning time of time step (set to sunrise if before sunrise)
 T1=MAX(T,-TP-DDT,-TD)
 ! Set end time of time step (adjust if after sunset)
 T2=MIN(T+DELT1,TD,TP-DDT)
 ! print *, 'First t1 and t2 = ', t1, t2
 IF(T2.LE.T1) THEN
  HRI=0.0 ! nighttime
 ELSE
 ! Calculate integral of radiation index from T1 to T2 and divide by time step DELTA1
 ! NOTE: this assumes the declination does not change from T1 to T2
  HRI=(SIN(D)*SIN(LP)*(T2-T1)+COS(D)*COS(LP)*(SIN(T2+DDT) &
      -SIN(T1+DDT)))/(COS(SLOPE1)*DELT1)  ! radiation index
 ENDIF
 ! print *, hri
 ! ----------------- for time intervals that extend to following day ----------------------
 ! Check to see of timestep extends to following day
 IF((T+DELT1).GT.PI_D) THEN
  ! Advance julian day by 1
  FJULIAN = FJULIAN + 1
  ! Calculate solar declination
  D=CRAD*23.5*SIN((FJULIAN-82.0)*YRAD)
  ! Calculate time of sunrise/sunset on level surface as radians from noon
  ! Account for high latitude locations, where there might not be a sunrise/sunset time on a given day
  ! In such cases AUX > 1 or AUX < -1. Fix AUX at (-)1 in those cases
  AUX=-TAN(LAT1)*TAN(D)
  IF(abs(AUX) > 1.) THEN
   TD=ACOS(SIGN(1._rkind, AUX))
  ELSE
   TD=ACOS(AUX)
  END IF

  ! print *, 'Sunrise #2 = ', TD, DELT1
  ! Calculate time of sunrise/sunset adjusted for inclined ground surface as radians from noon???
  TPI=-TAN(LP)*TAN(D)
  IF(ABS(TPI).LT.1.0) THEN
   TP=ACOS(TPI)
  ELSE IF(TPI.LT.-1.0) THEN ! 24h daylight
   TP=ACOS(-1.0)
  ELSE IF(TPI.GT.1.0) THEN ! 24h dark
   TP=ACOS(1.0)
  ENDIF
  ! Set beginning time to sunrise
  T1=MAX(-TP-DDT,-TD)
  ! Set end time of time step
  T2=MIN(T+DELT1-2*PI_D,TD,TP-DDT)
  ! print *, 'Second t1 and t2 = ', t1, t2
  IF(T2.LE.T1) THEN
   HRI=HRI ! still nighttime in day 2
  ELSE
   ! Calculate integral of radiation index from T1 to T2 and divide by time step DELTA1
   ! NOTE: this assumes the declination does not change from T1 to T2
   HRI=HRI+(SIN(D)*SIN(LP)*(T2-T1)+COS(D)*COS(LP)*(SIN(T2+DDT) &
           -SIN(T1+DDT)))/(COS(SLOPE1)*DELT1)  ! radiation index
  ENDIF
  ! print *, hri
 ENDIF
 ! ----------------------------------------------------------------------------------------
 ! Calculate cosine of solar zenith angle (= HRI for level surface)
 COSZEN = HRI*COS(SLOPE1)
 ! this is assumed to be an appropriate representative value over the
 ! time step.  It is used for albedo calculations.
 ! ----------------------------------------------------------------------------------------
 CONTAINS


 ! *************************************************************************************************
 ! internal function JULIAN: calculate day of year
 ! *************************************************************************************************
  FUNCTION JULIAN(MONTH,DAY)
  USE nrtype
  IMPLICIT NONE
  ! input
  INTEGER(I4B)                             :: MONTH,DAY  ! month and day
  ! output
  INTEGER(I4B)                             :: JULIAN     ! julian day
  ! internal
  INTEGER(I4B),DIMENSION(12)               :: MADD       ! julian day at start of each month
  ! specify the julian day at the start of each month (-1)
  MADD = (/0,31,59,90,120,151,181,212,243,273,304,334/)
  ! compute the julian day
  JULIAN=DAY+MADD(MONTH)
  RETURN
  END FUNCTION JULIAN
 ! ----------------------------------------------------------------------------------------
 END SUBROUTINE CLRSKY_RAD


end module sunGeomtry_module



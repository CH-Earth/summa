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

module soil_utilsAddPrime_module

! data types
USE nr_type

USE multiconst,only: gravity, & ! acceleration of gravity       (m s-2)
                     Tfreeze, & ! temperature at freezing    (K)
                     LH_fus,  & ! latent heat of fusion      (J kg-1, or m2 s-2)
                     R_wv       ! gas constant for water vapor  (J kg-1 K-1; [J = Pa m3])

! missing values                 
USE globalData,only:realMissing    ! missing real number       

USE soil_utils_module,only:matricHead
USE soil_utils_module,only:dPsi_dTheta
USE soil_utils_module,only:volFracLiq
USE soil_utils_module,only:dTheta_dPsi
USE soil_utils_module,only:satDeficit
USE soil_utils_module,only:log1pAcc
USE soil_utils_module,only:expm1Acc

! privacy
implicit none
private

! routines to make public

public::liquidHeadPrime
public::d2Theta_dPsi2
public::d2Theta_dTk2
contains


! ******************************************************************************************************************************
! public subroutine: compute the liquid water matric potential (and the derivatives w.r.t. total matric potential and temperature)
! ******************************************************************************************************************************
subroutine liquidHeadPrime(&
                      ! input
                      matricHeadTotal                          ,& ! intent(in)    : total water matric potential (m)
                      matricHeadTotalPrime                     ,& ! intent(in)    : total water matric potential time derivative (m s-1)
                      volFracLiq                               ,& ! intent(in)    : volumetric fraction of liquid water (-)
                      volFracIce                               ,& ! intent(in)    : volumetric fraction of ice (-)
                      vGn_alpha,vGn_n,theta_sat,theta_res,vGn_m,& ! intent(in)    : soil parameters
                      dVolTot_dPsi0                            ,& ! intent(in)    : derivative in total volumetric water content w.r.t. matric head (m-1)
                      dTheta_dT                                ,& ! intent(in)    : derivative in volumetric liquid water w.r.t. temperature (K-1)
                      volFracLiqPrime                          ,& ! intent(in)    : volumetric fraction of liquid water time derivative (-)
                      volFracIcePrime                          ,& ! intent(in)    : volumetric fraction of ice time derivative (-)
                      ! output
                      matricHeadLiq                            ,& ! intent(out)   : liquid water matric potential (m)
                      matricHeadLiqPrime                       ,& ! intent(out)   : liquid water matric potential time derivative (m s-1)
                      dPsiLiq_dPsi0                            ,& ! intent(out)   : derivative in the liquid water matric potential w.r.t. the total water matric potential (-)
                      dPsiLiq_dTemp                            ,& ! intent(out)   : derivative in the liquid water matric potential w.r.t. temperature (m K-1)
                      err,message)                                ! intent(out)   : error control
  ! computes the liquid water matric potential (and the derivatives w.r.t. total matric potential and temperature)
 implicit none
  ! input
  real(rkind),intent(in)      :: matricHeadTotal                           ! total water matric potential (m)
  real(rkind),intent(in)      :: matricHeadTotalPrime                      ! total water matric potential time derivative (m s-1)
  real(rkind),intent(in)      :: volFracLiq                                ! volumetric fraction of liquid water (-)
  real(rkind),intent(in)      :: volFracIce                                ! volumetric fraction of ice (-)
  real(rkind),intent(in)      :: vGn_alpha,vGn_n,theta_sat,theta_res,vGn_m ! soil parameters
  real(rkind),intent(in)      :: dVolTot_dPsi0                             ! derivative in total volumetric water content w.r.t. matric head (m-1)
  real(rkind),intent(in)      :: dTheta_dT                                 ! derivative in volumetric liquid water w.r.t. temperature (K-1)
  real(rkind),intent(in)      :: volFracLiqPrime                           ! volumetric fraction of liquid water time derivative ()
  real(rkind),intent(in)      :: volFracIcePrime                           ! volumetric fraction of ice time derivative ()
  ! output
  real(rkind),intent(out)     :: matricHeadLiq                             ! liquid water matric potential (m)
  real(rkind),intent(out)     :: matricHeadLiqPrime                        ! liquid water matric potential time derivative (m s-1)
  real(rkind),intent(out)     :: dPsiLiq_dPsi0                             ! derivative in the liquid water matric potential w.r.t. the total water matric potential (-)
  real(rkind),intent(out)     :: dPsiLiq_dTemp                             ! derivative in the liquid water matric potential w.r.t. temperature (m K-1)
  ! output: error control
  integer(i4b),intent(out)    :: err                                       ! error code
  character(*),intent(out)    :: message                                   ! error message
  ! local
  real(rkind)                 :: xNum,xDen                                 ! temporary variables (numeratir, denominator)
  real(rkind)                 :: effSat                                    ! effective saturation (-)
  real(rkind)                 :: dPsiLiq_dEffSat                           ! derivative in liquid water matric potential w.r.t. effective saturation (m)
  real(rkind)                 :: dEffSat_dTemp                             ! derivative in effective saturation w.r.t. temperature (K-1)
  real(rkind)                 :: effSatPrime                               ! effective saturation time derivative (-)
  real(rkind)                 :: xDef                                      ! saturation deficit of the total water, theta_sat - volFracLiq - volFracIce (-)
  real(rkind)                 :: satTerm                                   ! effSat**(-1/m) - 1 (-)
  ! ------------------------------------------------------------------------------------------------------------------------------
  ! initialize error control
  err=0; message='liquidHeadPrime/'

  ! ** partially frozen soil
  if(volFracIce > epsilon(1._rkind) .and. matricHeadTotal < 0._rkind)then  ! check that ice exists and that the soil is unsaturated
    ! -----
    ! - compute liquid water matric potential...
    ! ------------------------------------------

    ! - compute effective saturation
    ! NOTE: include ice content as part of the solid porosity - major effect of ice is to reduce the pore size; ensure that effSat=1 at saturation
    ! (from Zhao et al., J. Hydrol., 1997: Numerical analysis of simultaneous heat and mass transfer...)
    xNum   = volFracLiq - theta_res
    xDef   = satDeficit(matricHeadTotal,vGn_alpha,theta_sat,theta_res,vGn_n,vGn_m)
    xDen   = xNum + xDef  ! theta_sat - volFracIce - theta_res, without cancellation near saturation
    effSat = xNum/xDen    ! effective saturation

    ! - matric head associated with liquid water, and its derivative w.r.t. effective saturation (m)
    satTerm = 0._rkind
    if(xNum > 0._rkind) satTerm = expm1Acc(-log1pAcc(-xDef/xDen)/vGn_m) ! effSat**(-1/m) - 1 from 1-effSat, smooth as effSat -> 1
    if(satTerm > 0._rkind)then
      matricHeadLiq   = satTerm**(1._rkind/vGn_n)/vGn_alpha
      dPsiLiq_dEffSat = -satTerm**(1._rkind/vGn_n - 1._rkind)*effSat**(-1._rkind/vGn_m - 1._rkind)/(vGn_alpha*vGn_n*vGn_m)
    else
      matricHeadLiq   = matricHead(effSat,vGn_alpha,0._rkind,1._rkind,vGn_n,vGn_m)  ! argument is effective saturation, so theta_res=0 and theta_sat=1
      dPsiLiq_dEffSat = dPsi_dTheta(effSat,vGn_alpha,0._rkind,1._rkind,vGn_n,vGn_m)
    endif
    if (satTerm > 0._rkind .and. volFracLiqPrime.ne.realMissing .and. volFracIcePrime.ne.realMissing)then
      effSatPrime = (volFracLiqPrime * xDen + volFracIcePrime * xNum) / xDen**2_i4b
      matricHeadLiqPrime = dPsiLiq_dEffSat*effSatPrime
    else
      matricHeadLiqPrime = 0._rkind
    endif

    ! derivative in the liquid water matric potential w.r.t. the total water matric potential
    dPsiLiq_dPsi0 = dVolTot_dPsi0*dPsiLiq_dEffSat*xNum/(xDen**2_i4b)
    ! derivative in the liquid water matric potential w.r.t. temperature
    dEffSat_dTemp = -dTheta_dT*xNum/(xDen**2_i4b) + dTheta_dT/xDen
    dPsiLiq_dTemp = dPsiLiq_dEffSat*dEffSat_dTemp

  ! ** unfrozen soil
  else   ! (no ice)
    matricHeadLiq = matricHeadTotal
    matricHeadLiqPrime = matricHeadTotalPrime
    dPsiLiq_dPsi0 = 1._rkind  ! derivative=1 because values are identical
    dPsiLiq_dTemp = 0._rkind  ! derivative=0 because no impact of temperature for unfrozen conditions
  end if  ! (if ice exists)
  
  if(volFracLiqPrime==realMissing .or. volFracIcePrime==realMissing) matricHeadLiqPrime = realMissing

end subroutine liquidHeadPrime

! ******************************************************************************************************************************
! public function d2Theta_dPsi2: compute the second derivative of the soil water characteristic (m-1)
! ******************************************************************************************************************************
function d2Theta_dPsi2(psi,alpha,theta_res,theta_sat,n,m)
  implicit none
  real(rkind),intent(in) :: psi         ! soil water suction (m)
  real(rkind),intent(in) :: alpha       ! scaling parameter (m-1)
  real(rkind),intent(in) :: theta_res   ! residual volumetric water content (-)
  real(rkind),intent(in) :: theta_sat   ! porosity (-)
  real(rkind),intent(in) :: n           ! vGn "n" parameter (-)
  real(rkind),intent(in) :: m           ! vGn "m" parameter (-)
  real(rkind)            :: d2Theta_dPsi2 ! derivative of the soil water characteristic (m-1)
  real(rkind)            :: mult_fcn
  real(rkind)            :: mult_fcnp
  real(rkind)            :: xn1         ! (alpha psi)**(n-1)
  real(rkind)            :: s           ! 1 + (alpha psi)**n
  real(rkind)            :: sm          ! s**(-m)
  real(rkind)            :: dTheta      ! dTheta_dPsi (m-1)
  if(psi<0._rkind)then
    xn1 = (alpha*psi)**(n-1._rkind)
    s   = 1._rkind + xn1*alpha*psi
    sm  = s**(-m)
    mult_fcn  = (-m*n*alpha*xn1) / s
    mult_fcnp = -m*n*alpha*(n-1._rkind)*alpha*(xn1/(alpha*psi)) / s - &
                ( n*alpha*xn1/(s*s) ) * ( -m*n*alpha*xn1 )
    dTheta = (theta_sat-theta_res) * (-m*sm/s) * n*xn1 * alpha ! dTheta_dPsi and volFracLiq inlined to share the powers
    if(abs(dTheta) < epsilon(psi)) dTheta = epsilon(psi)
    d2Theta_dPsi2 = mult_fcn * dTheta + mult_fcnp * (theta_sat - theta_res)*sm
  else
    d2Theta_dPsi2 = 0._rkind
  end if
end function d2Theta_dPsi2


! ******************************************************************************************************************************
! public function d2Theta_dTk2: differentiate the freezing curve w.r.t. temperature
! ******************************************************************************************************************************
function d2Theta_dTk2(Tk,theta_res,theta_sat,alpha,n,m)
  implicit none
  real(rkind),intent(in) :: Tk            ! temperature (K)
  real(rkind),intent(in) :: theta_res     ! residual liquid water content (-)
  real(rkind),intent(in) :: theta_sat     ! porosity (-)
  real(rkind),intent(in) :: alpha         ! vGn scaling parameter (m-1)
  real(rkind),intent(in) :: n             ! vGn "n" parameter (-)
  real(rkind),intent(in) :: m             ! vGn "m" parameter (-)
  real(rkind)            :: d2Theta_dTk2    ! derivative of the freezing curve w.r.t. temperature (K-1)
  ! local variables
  real(rkind)            :: kappa         ! constant (m K-1)
  real(rkind)            :: xtemp         ! alpha*kappa*(Tk-Tfreeze) -- dimensionless variable (used more than once)
  real(rkind)            :: xn1, s, sm1   ! xtemp**(n-1), 1 + xtemp**n, s**(-m-1)
  ! compute kappa (m K-1)
  kappa =  LH_fus/(gravity*Tfreeze)    ! NOTE: J = kg m2 s-2
  ! define a tempory variable that is used more than once (-)
  xtemp = alpha*kappa*(Tk-Tfreeze)
  ! differentiate the freezing curve w.r.t. temperature -- making use of the chain rule
  xn1 = xtemp**(n - 1._rkind)
  s   = 1._rkind + xn1*xtemp
  sm1 = s**(-m - 1._rkind)
  d2Theta_dTk2 = (-alpha*kappa*m*n*alpha*kappa)* (theta_sat - theta_res) * (  (n-1)*(xn1/xtemp) * sm1 &
                                                                          + n*(-m-1)*xn1*xn1 * sm1/s )
end function d2Theta_dTk2

end module soil_utilsAddPrime_module

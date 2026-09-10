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

module soil_utils_module

! data types
USE nr_type

! constants
USE globalData,only: verySmall  ! a small number
USE multiconst,only: gravity, & ! acceleration of gravity       (m s-2)
                     Tfreeze, & ! temperature at freezing       (K)
                     LH_fus,  & ! latent heat of fusion         (J kg-1, or m2 s-2)
                     R_wv       ! gas constant for water vapor  (J kg-1 K-1; [J = Pa m3])

! privacy
implicit none
private
public::iceImpede
public::dIceImpede_dTemp
public::hydCond_psi
public::hydCondMP_liq
public::dHydCond_dPsi
public::volFracLiq
public::matricHead
public::dTheta_dPsi
public::dPsi_dTheta
public::RH_soilair
public::dTheta_dTk
public::crit_soilT
public::liquidHead
public::gammp,gammp_complex
public::LogSumExp
public::SoftArgMax
contains

! ******************************************************************************************************************************
! public subroutine iceImpede: compute the ice impedence factor
! ******************************************************************************************************************************
subroutine iceImpede(volFracIce,f_impede, &            ! input
                    iceImpedeFactor,dIceImpede_dLiq)  ! output
  implicit none
  ! input variables
  real(rkind),intent(in)     :: volFracIce        ! volumetric fraction of ice (-)
  real(rkind),intent(in)     :: f_impede          ! ice impedence parameter (-)
  ! output variables
  real(rkind)                :: iceImpedeFactor   ! ice impedence factor (-)
  real(rkind)                :: dIceImpede_dLiq   ! derivative in ice impedence factor w.r.t. volumetric liquid water content (-)

  ! compute ice impedance factor as a function of volumetric ice content
  iceImpedeFactor = 10._rkind**(-f_impede*volFracIce)
  dIceImpede_dLiq = 0._rkind

end subroutine iceImpede


! ******************************************************************************************************************************
! public subroutine dIceImpede_dTemp: compute the derivative in the ice impedence factor w.r.t. temperature
! ******************************************************************************************************************************
subroutine dIceImpede_dTemp(volFracIce,dTheta_dT,f_impede,dIceImpede_dT)
  implicit none
  ! input variables
  real(rkind),intent(in)     :: volFracIce        ! volumetric fraction of ice (-)
  real(rkind),intent(in)     :: dTheta_dT         ! derivative in volumetric liquid water content w.r.t temperature (K-1)
  real(rkind),intent(in)     :: f_impede          ! ice impedence parameter (-)
  ! output variables
  real(rkind)                :: dIceImpede_dT     ! derivative in the ice impedance factor w.r.t. temperature (K-1)

  dIceImpede_dT = log(10._rkind)*f_impede*(10._rkind**(-f_impede*volFracIce))*dTheta_dT
end subroutine dIceImpede_dTemp


! ******************************************************************************************************************************
! public subroutine: compute the liquid water matric potential (and the derivatives w.r.t. total matric potential and temperature)
! ******************************************************************************************************************************
subroutine liquidHead(&
                     ! input
                     matricHeadTotal                          ,& ! intent(in)    : total water matric potential (m)
                     volFracLiq                               ,& ! intent(in)    : volumetric fraction of liquid water (-)
                     volFracIce                               ,& ! intent(in)    : volumetric fraction of ice (-)
                     vGn_alpha,vGn_n,theta_sat,theta_res,vGn_m,& ! intent(in)    : soil parameters
                     dVolTot_dPsi0                            ,& ! intent(in)    : derivative in total volumetric water content w.r.t. matric head (m-1)
                     dTheta_dT                                ,& ! intent(in)    : derivative in volumetric liquid water w.r.t. temperature (K-1)
                     ! output
                     matricHeadLiq                            ,& ! intent(out)   : liquid water matric potential (m)
                     dPsiLiq_dPsi0                            ,& ! intent(out)   : derivative in the liquid water matric potential w.r.t. the total water matric potential (-)
                     dPsiLiq_dTemp                            ,& ! intent(out)   : derivative in the liquid water matric potential w.r.t. temperature (m K-1)
                     err,message)                                ! intent(out)   : error control
  implicit none
  ! input
  real(rkind),intent(in)            :: matricHeadTotal                           ! total water matric potential (m)
  real(rkind),intent(in)            :: volFracLiq                                ! volumetric fraction of liquid water (-)
  real(rkind),intent(in)            :: volFracIce                                ! volumetric fraction of ice (-)
  real(rkind),intent(in)            :: vGn_alpha,vGn_n,theta_sat,theta_res,vGn_m ! soil parameters
  real(rkind),intent(in)  ,optional :: dVolTot_dPsi0                             ! derivative in total volumetric water content w.r.t. matric head (m-1)
  real(rkind),intent(in)  ,optional :: dTheta_dT                                 ! derivative in volumetric liquid water w.r.t. temperature (K-1)
  ! output
  real(rkind),intent(out)           :: matricHeadLiq                             ! liquid water matric potential (m)
  real(rkind),intent(out) ,optional :: dPsiLiq_dPsi0                             ! derivative in the liquid water matric potential w.r.t. the total water matric potential (-)
  real(rkind),intent(out) ,optional :: dPsiLiq_dTemp                             ! derivative in the liquid water matric potential w.r.t. temperature (m K-1)
  ! output: error control
  integer(i4b),intent(out)          :: err                                       ! error code
  character(*),intent(out)          :: message                                   ! error message
  ! local
  real(rkind)                       :: xNum,xDen                                 ! temporary variables (numeratir, denominator)
  real(rkind)                       :: effSat                                    ! effective saturation (-)
  real(rkind)                       :: dPsiLiq_dEffSat                           ! derivative in liquid water matric potential w.r.t. effective saturation (m)
  real(rkind)                       :: dEffSat_dTemp                             ! derivative in effective saturation w.r.t. temperature (K-1)
   ! ------------------------------------------------------------------------------------------------------------------------------
  ! initialize error control
  err=0; message='liquidHead/'

  ! ** partially frozen soil
  if(volFracIce > epsilon(1._rkind) .and. matricHeadTotal < 0._rkind)then  ! check that ice exists and that the soil is unsaturated

    ! -----
    ! - compute liquid water matric potential...
    ! ------------------------------------------

    ! - compute effective saturation
    ! NOTE: include ice content as part of the solid porosity - major effect of ice is to reduce the pore size; ensure that effSat=1 at saturation
    ! (from Zhao et al., J. Hydrol., 1997: Numerical analysis of simultaneous heat and mass transfer...)
    xNum   = volFracLiq - theta_res
    xDen   = theta_sat - volFracIce - theta_res
    effSat = xNum/xDen          ! effective saturation

    ! - matric head associated with liquid water
    matricHeadLiq = matricHead(effSat,vGn_alpha,0._rkind,1._rkind,vGn_n,vGn_m)  ! argument is effective saturation, so theta_res=0 and theta_sat=1

    ! compute derivative in liquid water matric potential w.r.t. effective saturation (m)
    if(present(dPsiLiq_dPsi0).or.present(dPsiLiq_dTemp))then
      dPsiLiq_dEffSat = dPsi_dTheta(effSat,vGn_alpha,0._rkind,1._rkind,vGn_n,vGn_m)
    endif

    ! -----
    ! - compute derivative in the liquid water matric potential w.r.t. the total water matric potential...
    ! ----------------------------------------------------------------------------------------------------

    ! check if the derivative is desired
    if(present(dPsiLiq_dPsi0))then
      ! (check required input derivative is present)
      if(.not.present(dVolTot_dPsi0))then
        message=trim(message)//'dVolTot_dPsi0 argument is missing'
        err=20; return
      endif
      ! (compute derivative in the liquid water matric potential w.r.t. the total water matric potential)
      dPsiLiq_dPsi0 = dVolTot_dPsi0*dPsiLiq_dEffSat*xNum/(xDen**2_i4b)
    endif  ! if dPsiLiq_dTemp is desired

    ! -----
    ! - compute the derivative in the liquid water matric potential w.r.t. temperature...
    ! -----------------------------------------------------------------------------------

    ! check if the derivative is desired
    if(present(dPsiLiq_dTemp))then
      ! (check required input derivative is present)
      if(.not.present(dTheta_dT))then
        message=trim(message)//'dTheta_dT argument is missing'
        err=20; return
      endif
      ! (compute the derivative in the liquid water matric potential w.r.t. temperature)
      dEffSat_dTemp = -dTheta_dT*xNum/(xDen**2_i4b) + dTheta_dT/xDen
      dPsiLiq_dTemp = dPsiLiq_dEffSat*dEffSat_dTemp
    endif  ! if dPsiLiq_dTemp is desired

  ! ** unfrozen soil
  else   ! (no ice)
    matricHeadLiq = matricHeadTotal
    if(present(dPsiLiq_dPsi0)) dPsiLiq_dPsi0 = 1._rkind  ! derivative=1 because values are identical
    if(present(dPsiLiq_dTemp)) dPsiLiq_dTemp = 0._rkind  ! derivative=0 because no impact of temperature for unfrozen conditions
  end if  ! (if ice exists)

end subroutine liquidHead


! ******************************************************************************************************************************
! public function hydCondMP_liq: compute the hydraulic conductivity of macropores as a function of liquid water content (m s-1)
! ******************************************************************************************************************************
function hydCondMP_liq(volFracLiq,theta_sat,theta_mp,mpExp,satHydCond_ma,satHydCond_mi)
  ! computes hydraulic conductivity given volFracLiq and soil hydraulic parameters
  !  theta_sat, theta_mp, mpExp, satHydCond_ma, and satHydCond_mi
  implicit none
  ! dummies
  real(rkind),intent(in) :: volFracLiq    ! volumetric liquid water content (-)
  real(rkind),intent(in) :: theta_sat     ! soil porosity (-)
  real(rkind),intent(in) :: theta_mp      ! minimum volumetric liquid water content for macropore flow (-)
  real(rkind),intent(in) :: mpExp         ! empirical exponent in macropore flow equation (-)
  real(rkind),intent(in) :: satHydCond_ma ! saturated hydraulic conductivity for macropores (m s-1)
  real(rkind),intent(in) :: satHydCond_mi ! saturated hydraulic conductivity for micropores (m s-1)
  real(rkind)            :: hydCondMP_liq ! hydraulic conductivity (m s-1)
  ! locals
  real(rkind)            :: theta_e     ! effective soil moisture

  if(volFracLiq > theta_mp)then
    theta_e       = (volFracLiq - theta_mp) / (theta_sat - theta_mp)
    hydCondMP_liq = (satHydCond_ma - satHydCond_mi) * (theta_e**mpExp)
  else
    hydCondMP_liq = 0._rkind
  end if
end function hydCondMP_liq


! ******************************************************************************************************************************
! public function hydCond_psi: compute the hydraulic conductivity as a function of matric head (m s-1)
! ******************************************************************************************************************************
function hydCond_psi(psi,k_sat,alpha,n,m)
  implicit none
  ! dummies
  real(rkind),intent(in)     :: psi           ! soil water suction (m)
  real(rkind),intent(in)     :: k_sat         ! saturated hydraulic conductivity (m s-1)
  real(rkind),intent(in)     :: alpha         ! scaling parameter (m-1)
  real(rkind),intent(in)     :: n             ! vGn "n" parameter (-)
  real(rkind),intent(in)     :: m             ! vGn "m" parameter (-)
  real(rkind)                :: hydCond_psi   ! hydraulic conductivity (m s-1)
  ! Smooth transition to k_sat as psi -> 0 from below.
  ! Blend over the interval [-delta,0] using a cubic smoothstep so that value
  ! and first derivative are continuous at the join.
  real(rkind) :: t, s, orig_val

  if (psi >= 0._rkind) then
    hydCond_psi = k_sat
  else
    orig_val = k_sat * ( ( (1._rkind - (psi*alpha)**(n-1._rkind) * (1._rkind + (psi*alpha)**n)**(-m))**2_i4b ) / ( (1._rkind + (psi*alpha)**n)**(m/2._rkind) ) )
    if (psi < -verySmall) then
      hydCond_psi = orig_val
    else ! compute original formula and blend to k_sat using a cubic smoothstep
      t = (psi + verySmall) / verySmall
      s = 3._rkind*t*t - 2._rkind*t*t*t     ! cubic smoothstep: s(0)=0, s(1)=1, continuous first derivative
      hydCond_psi = orig_val*(1._rkind - s) + k_sat*s
    end if
  end if
end function hydCond_psi


! ******************************************************************************************************************************
! public function volFracLiq: compute the volumetric liquid water content as a function of matric head (-)
! ******************************************************************************************************************************
function volFracLiq(psi,alpha,theta_res,theta_sat,n,m)
  implicit none
  real(rkind),intent(in) :: psi         ! soil water suction (m)
  real(rkind),intent(in) :: alpha       ! scaling parameter (m-1)
  real(rkind),intent(in) :: theta_res   ! residual volumetric water content (-)
  real(rkind),intent(in) :: theta_sat   ! porosity (-)
  real(rkind),intent(in) :: n           ! vGn "n" parameter (-)
  real(rkind),intent(in) :: m           ! vGn "m" parameter (-)
  real(rkind)            :: volFracLiq  ! volumetric liquid water content (-)

  if(psi<0._rkind)then
    volFracLiq = theta_res + (theta_sat - theta_res)*(1._rkind + (alpha*psi)**n)**(-m)
  else
    volFracLiq = theta_sat
  end if
end function volFracLiq


! ******************************************************************************************************************************
! public function matricHead: compute the matric head (m) based on the volumetric liquid water content
! ******************************************************************************************************************************
function matricHead(theta,alpha,theta_res,theta_sat,n,m)
  implicit none
  ! dummy variables
  real(rkind),intent(in) :: theta       ! volumetric liquid water content (-)
  real(rkind),intent(in) :: alpha       ! scaling parameter (m-1)
  real(rkind),intent(in) :: theta_res   ! residual volumetric water content (-)
  real(rkind),intent(in) :: theta_sat   ! porosity (-)
  real(rkind),intent(in) :: n           ! vGn "n" parameter (-)
  real(rkind),intent(in) :: m           ! vGn "m" parameter (-)
  real(rkind)            :: matricHead  ! matric head (m)
  ! local variables
  real(rkind)            :: effSat      ! effective saturation (-)
  real(rkind),parameter  :: eps=epsilon(1._rkind) ! a very small number (avoid effective saturation of zero)

  ! compute effective saturation
  effSat = max(eps, (theta - theta_res) / (theta_sat - theta_res))
  ! compute matric head
  if (effSat < 1._rkind .and. effSat > 0._rkind)then
    matricHead = (1._rkind/alpha)*( effSat**(-1._rkind/m) - 1._rkind)**(1._rkind/n)
  else
    matricHead = 0._rkind
  end if
end function matricHead


! ******************************************************************************************************************************
! public function dTheta_dPsi: compute the derivative of the soil water characteristic (m-1)
! ******************************************************************************************************************************
function dTheta_dPsi(psi,alpha,theta_res,theta_sat,n,m)
  implicit none
  real(rkind),intent(in) :: psi         ! soil water suction (m)
  real(rkind),intent(in) :: alpha       ! scaling parameter (m-1)
  real(rkind),intent(in) :: theta_res   ! residual volumetric water content (-)
  real(rkind),intent(in) :: theta_sat   ! porosity (-)
  real(rkind),intent(in) :: n           ! vGn "n" parameter (-)
  real(rkind),intent(in) :: m           ! vGn "m" parameter (-)
  real(rkind)            :: dTheta_dPsi ! derivative of the soil water characteristic (m-1)

  if(psi<=0._rkind)then
    dTheta_dPsi = (theta_sat-theta_res) * &
                 (-m*(1._rkind + (psi*alpha)**n)**(-m-1._rkind)) * n*(psi*alpha)**(n-1._rkind) * alpha
    if(abs(dTheta_dPsi) < epsilon(psi)) dTheta_dPsi = epsilon(psi) ! use to avoid division by zero if divide by dTheta_dPsi
  else
    dTheta_dPsi = epsilon(psi) ! use to avoid division by zero if divide by dTheta_dPsi
  end if
end function dTheta_dPsi


! ******************************************************************************************************************************
! public function dPsi_dTheta: compute the derivative of the soil water characteristic (m)
! ******************************************************************************************************************************
function dPsi_dTheta(volFracLiq,alpha,theta_res,theta_sat,n,m)
  implicit none
  ! dummies
  real(rkind),intent(in) :: volFracLiq  ! volumetric liquid water content (-)
  real(rkind),intent(in) :: alpha       ! scaling parameter (m-1)
  real(rkind),intent(in) :: theta_res   ! residual volumetric water content (-)
  real(rkind),intent(in) :: theta_sat   ! porosity (-)
  real(rkind),intent(in) :: n           ! vGn "n" parameter (-)
  real(rkind),intent(in) :: m           ! vGn "m" parameter (-)
  real(rkind)            :: dPsi_dTheta ! derivative of the soil water characteristic (m)
  ! locals
  real(rkind)            :: y1,d1       ! 1st function and derivative
  real(rkind)            :: y2,d2       ! 2nd function and derivative
  real(rkind)            :: theta_e     ! effective soil moisture
  real(rkind),parameter  :: theta_e_min=0.001_rkind            ! minimum effective soil moisture
  real(rkind),parameter  :: y1_min=10._rkind*epsilon(1._rkind) ! minimum y1 value (to avoid division by zero and complex values)

  ! check if less than saturation
  if(volFracLiq < theta_sat)then
   ! compute effective water content
   theta_e = max(theta_e_min,(volFracLiq - theta_res) / (theta_sat - theta_res))
   ! compute the 1st function and derivative
   y1 = theta_e**(-1._rkind/m) - 1._rkind
   d1 = (-1._rkind/m)*theta_e**(-1._rkind/m - 1._rkind) / (theta_sat - theta_res)
   ! compute the 2nd function and derivative
   ! note: impose a minimum value for y1 to avoid divison by zero and complex values
   !y2 = y1**(1._rkind/n)                         ! original expression
   !d2 = (1._rkind/n)*y1**(1._rkind/n - 1._rkind) ! original expression
   y2 = max(y1_min,y1)**(1._rkind/n)
   d2 = (1._rkind/n)*max(y1_min,y1)**(1._rkind/n - 1._rkind) ! impose a minimum value for y1 to avoid divison by zero and complex values
   ! compute the final function value
   dPsi_dTheta = d1*d2/alpha
  else
   dPsi_dTheta = 0._rkind
  end if
end function dPsi_dTheta


! ******************************************************************************************************************************
! public function dHydCond_dPsi: compute the derivative in hydraulic conductivity w.r.t. matric head (s-1)
! ******************************************************************************************************************************
function dHydCond_dPsi(psi,k_sat,alpha,n,m)
  implicit none
  ! dummies
  real(rkind),intent(in)     :: psi            ! soil water suction (m)
  real(rkind),intent(in)     :: k_sat          ! saturated hydraulic conductivity (m s-1)
  real(rkind),intent(in)     :: alpha          ! scaling parameter (m-1)
  real(rkind),intent(in)     :: n              ! vGn "n" parameter (-)
  real(rkind),intent(in)     :: m              ! vGn "m" parameter (-)
  real(rkind)                :: dHydCond_dPsi  ! derivative in hydraulic conductivity w.r.t. matric head (s-1)
  ! locals for analytical derivatives
  real(rkind)                :: f_x1          ! f(x) for part of the numerator
  real(rkind)                :: f_x2          ! f(x) for part of the numerator
  real(rkind)                :: f_nm          ! f(x) for the numerator
  real(rkind)                :: f_dm          ! f(x) for the denominator
  real(rkind)                :: d_x1          ! df(x)/dpsi for part of the numerator
  real(rkind)                :: d_x2          ! df(x)/dpsi for part of the numerator
  real(rkind)                :: d_nm          ! df(x)/dpsi for the numerator
  real(rkind)                :: d_dm          ! df(x)/dpsi for the denominator
  real(rkind)                :: t,s,orig_val,orig_d,ds_dpsi
  
  if (psi >= 0._rkind) then
    dHydCond_dPsi = 0._rkind
  else
    ! compute the derivative for the numerator (original Van Genuchten form)
    f_x1 = (psi*alpha)**(n - 1._rkind)
    f_x2 = (1._rkind + (psi*alpha)**n)**(-m)
    d_x1 = alpha * (n - 1._rkind)*(psi*alpha)**(n - 2._rkind)
    d_x2 = alpha * n*(psi*alpha)**(n - 1._rkind) * (-m)*(1._rkind + (psi*alpha)**n)**(-m - 1._rkind)
    f_nm = (1._rkind - f_x1*f_x2)**2_i4b
    d_nm = (-d_x1*f_x2 - f_x1*d_x2) * 2._rkind*(1._rkind - f_x1*f_x2)
    ! compute the derivative for the denominator
    f_dm = (1._rkind + (psi*alpha)**n)**(m/2._rkind)
    d_dm = alpha * n*(psi*alpha)**(n - 1._rkind) * (m/2._rkind)*(1._rkind + (psi*alpha)**n)**(m/2._rkind - 1._rkind)
    ! and combine to get the original derivative
    orig_val = k_sat * (f_nm / f_dm)
    orig_d   = k_sat*(d_nm*f_dm - d_dm*f_nm) / (f_dm**2_i4b)
    if(psi < -verySmall) then
      dHydCond_dPsi = orig_d
    else ! compute original formula and blend to k_sat using a cubic smoothstep
      t = (psi + verySmall) / verySmall
      t = max(0._rkind, min(1._rkind, t))
      s = 3._rkind*t*t - 2._rkind*t*t*t
      ds_dpsi = (6._rkind*t*(1._rkind - t)) / verySmall
      dHydCond_dPsi = orig_d*(1._rkind - s) + (k_sat - orig_val)*ds_dpsi
    end if
  end if
end function dHydCond_dPsi


! ******************************************************************************************************************************
! public function RH_soilair: compute relative humidity of air in soil pore space
! ******************************************************************************************************************************
function RH_soilair(matpot,Tk)
  implicit none
  real(rkind),intent(in) :: matpot        ! soil water suction -- matric potential (m)
  real(rkind),intent(in) :: Tk            ! temperature (K)
  real(rkind)            :: RH_soilair    ! relative humidity of air in soil pore space
  ! compute relative humidity (UNITS NOTE: Pa = kg m-1 s-2, so R_wv units = m2 s-2 K-1)
  RH_soilair = exp( (gravity*matpot) / (R_wv*Tk) )
end function RH_soilair


! ******************************************************************************************************************************
! public function crit_soilT: compute the critical temperature above which all water is unfrozen
! ******************************************************************************************************************************
function crit_soilT(psi)
  implicit none
  real(rkind),intent(in) :: psi           ! matric head (m)
  real(rkind)            :: crit_soilT    ! critical soil temperature (K)
  crit_soilT = Tfreeze + min(psi,0._rkind)*gravity*Tfreeze/LH_fus
end function crit_soilT


! ******************************************************************************************************************************
! public function dTheta_dTk: differentiate the freezing curve w.r.t. temperature
! ******************************************************************************************************************************
function dTheta_dTk(Tk,theta_res,theta_sat,alpha,n,m)
  implicit none
  real(rkind),intent(in) :: Tk            ! temperature (K)
  real(rkind),intent(in) :: theta_res     ! residual liquid water content (-)
  real(rkind),intent(in) :: theta_sat     ! porosity (-)
  real(rkind),intent(in) :: alpha         ! vGn scaling parameter (m-1)
  real(rkind),intent(in) :: n             ! vGn "n" parameter (-)
  real(rkind),intent(in) :: m             ! vGn "m" parameter (-)
  real(rkind)            :: dTheta_dTk    ! derivative of the freezing curve w.r.t. temperature (K-1)
  ! local variables
  real(rkind)            :: kappa         ! constant (m K-1)
  real(rkind)            :: xtemp         ! alpha*kappa*(Tk-Tfreeze) -- dimensionless variable (used more than once)
  ! compute kappa (m K-1)
  kappa =  LH_fus/(gravity*Tfreeze)    ! NOTE: J = kg m2 s-2
  ! define a tempory variable that is used more than once (-)
  xtemp = alpha*kappa*(Tk-Tfreeze)
  ! differentiate the freezing curve w.r.t. temperature -- making use of the chain rule
  dTheta_dTk = (alpha*kappa) * n*xtemp**(n - 1._rkind) * (-m)*(1._rkind + xtemp**n)**(-m - 1._rkind) * (theta_sat - theta_res)
end function dTheta_dTk


! ******************************************************************************************************************************
! public function gammp: compute cumulative probability using the Gamma distribution (Gamma CDF)
! ******************************************************************************************************************************
function gammp(a,x)
  implicit none
  ! input
  real(rkind), intent(in) :: a,x
  ! output
  real(rkind)             :: gammp
  ! validation
  if (a < 0._rkind) then
   stop "Error in gammp: a >= 0 required."
  end if
  if (x < 0._rkind) then
   stop "Error in gammp: x >= 0 required."
  end if
  ! computation
  if (x<a+1.0_rkind) then
   gammp=gser(a,x)
  else
   gammp=1.0_rkind-gcf(a,x)
  end if
end function gammp


! ******************************************************************************************************************************
! private function gcf: continued fraction development of the incomplete Gamma function
! ******************************************************************************************************************************
function gcf(a,x,gln)
  implicit none
  ! input
  real(rkind),           intent(in)  :: a,x
  ! output
  real(rkind), optional, intent(out) :: gln
  real(rkind)                        :: gcf
  ! local variables
  integer(i4b), parameter :: ITMAX=100
  real(rkind),  parameter :: EPS=epsilon(x),FPMIN=tiny(x)/EPS
  integer(i4b)            :: i
  real(rkind)             :: an,b,c,d,del,h
  if (x == 0.0_rkind) then
   gcf=1.0_rkind
   return
  end if
  b=x+1.0_rkind-a
  c=1.0_rkind/FPMIN
  d=1.0_rkind/b
  h=d
  do i=1,ITMAX
   an=-i*(i-a)
   b=b+2.0_rkind
   d=an*d+b
   if (abs(d) < FPMIN) d=FPMIN
   c=b+an/c
   if (abs(c) < FPMIN) c=FPMIN
   d=1.0_rkind/d
   del=d*c
   h=h*del
   if (abs(del-1.0_rkind) <= EPS) exit
  end do
  if (i > ITMAX) stop 'a too large, ITMAX too small in gcf'
  if (present(gln)) then
   gln=log_gamma(a)
   gcf=exp(-x+a*log(x)-gln)*h
  else
   gcf=exp(-x+a*log(x)-log_gamma(a))*h
  end if
end function gcf

! ******************************************************************************************************************************
! private function gser: series development of the incomplete Gamma function
! ******************************************************************************************************************************
function gser(a,x,gln)
  implicit none
  ! input
  real(rkind),           intent(in)  :: a,x
  ! output
  real(rkind), optional, intent(out) :: gln
  real(rkind)                        :: gser
  ! local variables
  integer(i4b), parameter :: ITMAX=100
  real(rkind),  parameter :: EPS=epsilon(x)
  integer(i4b)            :: n
  real(rkind)             :: ap,del,summ
  if (x == 0.0_rkind) then
   gser=0.0_rkind
   return
  end if
  ap=a
  summ=1.0_rkind/a
  del=summ
  do n=1,ITMAX
   ap=ap+1.0_rkind
   del=del*x/ap
   summ=summ+del
   if (abs(del) < abs(summ)*EPS) exit
  end do
  if (n > ITMAX) stop 'a too large, ITMAX too small in gser'
  if (present(gln)) then
   gln=log_gamma(a)
   gser=summ*exp(-x+a*log(x)-gln) 
  else
   gser=summ*exp(-x+a*log(x)-log_gamma(a)) 
  end if
end function gser

! ******************************************************************************************************************************
! public function gammp_complex: regularized lower incomplete gamma function (complex output)
! ******************************************************************************************************************************
! Note: input parameters are real but output may have non-zero imaginary parts
function gammp_complex(a,x)
  implicit none
  ! input
  real(rkind), intent(in) :: a,x
  ! output
  complex(rkind)          :: gammp_complex
  ! validation
  if (a < 0._rkind) then
   stop "Error in gammp_complex: a >= 0 required."
  end if
  ! computation
  if (x<a+1.0_rkind) then
   gammp_complex=gser_complex(a,x)
  else
   gammp_complex=1.0_rkind-gcf_complex(a,x)
  end if
end function gammp_complex

! ******************************************************************************************************************************
! private function gcf_complex: continued fraction development of the incomplete Gamma function (complex output)
! ******************************************************************************************************************************
function gcf_complex(a,x,gln)
  implicit none
  ! input
  real(rkind),           intent(in)  :: a,x
  ! output
  real(rkind), optional, intent(out) :: gln
  complex(rkind)                     :: gcf_complex
  ! local variables
  integer(i4b),          parameter   :: ITMAX=100
  real(rkind),           parameter   :: EPS=epsilon(x),FPMIN=tiny(x)/EPS
  integer(i4b)                       :: i
  real(rkind)                        :: an,b,c,d,del,h
  if (x == 0.0_rkind) then
   gcf_complex=1.0_rkind
   return
  end if
  b=x+1.0_rkind-a
  c=1.0_rkind/FPMIN
  d=1.0_rkind/b
  h=d
  do i=1,ITMAX
   an=-i*(i-a)
   b=b+2.0_rkind
   d=an*d+b
   if (abs(d) < FPMIN) d=FPMIN
   c=b+an/c
   if (abs(c) < FPMIN) c=FPMIN
   d=1.0_rkind/d
   del=d*c
   h=h*del
   if (abs(del-1.0_rkind) <= EPS) exit
  end do
  if (i > ITMAX) stop 'a too large, ITMAX too small in gcf'
  if (present(gln)) then
   gln=log_gamma(a)
   gcf_complex=exp(-x-gln)*cmplx(x,0._rkind,rkind)**a*h       ! allows x<0
  else
   gcf_complex=exp(-x-log_gamma(a))*cmplx(x,0._rkind,rkind)**a*h ! allows x<0
  end if
end function gcf_complex


! ******************************************************************************************************************************
! private function gser_complex: series development of the incomplete Gamma function (complex output)
! ******************************************************************************************************************************
function gser_complex(a,x,gln)
  implicit none
  ! input
  real(rkind),           intent(in)  :: a,x
  ! output
  real(rkind), optional, intent(out) :: gln
  complex(rkind)                     :: gser_complex
  ! local variables
  integer(i4b),          parameter   :: ITMAX=100
  real(rkind),           parameter   :: EPS=epsilon(x)
  integer(i4b)                       :: n
  real(rkind)                        :: ap,del,summ
  if (x == 0.0_rkind) then
   gser_complex=(0.0_rkind,0.0_rkind)
   return
  end if
  ap=a
  summ=1.0_rkind/a
  del=summ
  do n=1,ITMAX
   ap=ap+1.0_rkind
   del=del*x/ap
   summ=summ+del
   if (abs(del) < abs(summ)*EPS) exit
  end do
  if (n > ITMAX) stop 'a too large, ITMAX too small in gser'
  if (present(gln)) then
   gln=log_gamma(a)
   gser_complex=summ*exp(-x-gln)*cmplx(x,0._rkind,rkind)**a       ! allows x<0
  else
   gser_complex=summ*exp(-x-log_gamma(a))*cmplx(x,0._rkind,rkind)**a ! allows x<0
  end if
end function gser_complex


! ******************************************************************************************************************************
! public function LogSumExp: LSE (or RealSoftMax) function used for smooth approximations of max or min functions
! ******************************************************************************************************************************
function LogSumExp(alpha,x,err) result(LSE)
  use, intrinsic :: ieee_arithmetic,only:ieee_value,ieee_is_normal,ieee_quiet_nan
  use, intrinsic :: iso_fortran_env,only:real128
  implicit none
  ! input
  real(rkind),intent(in) :: alpha ! smoothness parameter (LSE --> max as alpha --> +Inf, LSE --> min as alpha --> -Inf)
  real(rkind),intent(in) :: x(:)  ! vector of input values
  ! output
  real(rkind)              :: LSE ! LogSumExp value 
  integer(i4b),intent(out) :: err ! error code
  ! local variables
  real(real128),allocatable :: x_qp(:) ! quadruple precision x vector
  real(real128) :: x_star   ! quadruple precision shift value for numerical stability
  real(real128) :: alpha_qp ! quadruple precision alpha
  real(real128) :: LSE_qp   ! quadruple precision LSE value 

  err = 0_i4b ! initialize error code

  ! validation of input parameters
  if (alpha == 0._rkind) then
   err = 20_i4b ! positive error code to indicate failure
   LSE = ieee_value(0._rkind,ieee_quiet_nan) ! assign NaN return value
   return
  end if

  ! use quadruple precision variables to prevent over/underflow
  alpha_qp = real(alpha,real128)
  allocate(x_qp(size(x)))
  x_qp     = real(x,real128)

  ! shift value to improve numerical stability
  x_star = maxval(abs(x_qp))

  LSE_qp= x_star + log(sum(exp(alpha_qp*(x_qp-x_star))))/alpha_qp
  LSE=real(LSE_qp,rkind)
  
  ! check if value is normal (not NaN, -Infinity, or +Infinity)
  ! note: mainly to account for overflow/underflow that may occur in extreme cases
  if (ieee_is_normal(LSE)) then ! return if value is not NaN or infinity
    return
  else                          ! revert to analytic max/min function as a failsafe (accurate but not smoothed)
    if (alpha < 0._rkind) then  ! min
      LSE = minval(x)
    else                        ! max (alpha cannot be zero)
      LSE = maxval(x)
    end if
  end if

end function LogSumExp


! ******************************************************************************************************************************
! public function SoftArgMax: SoftArgMax (aliases: softmax, normalized exponential) function for smooth approximations to argument max or min
! ******************************************************************************************************************************
! Note: Can be used to evaluate the derivatives of LogSumExp
! dLogSumExp(alpha,x)_dx(i) = SoftArgMax(alpha,x)
function SoftArgMax(alpha,x) result(SAM)
  use, intrinsic :: ieee_arithmetic,only:ieee_is_normal
  use, intrinsic :: iso_fortran_env,only:real128
  implicit none
  ! input
  real(rkind),intent(in) :: alpha ! smoothness parameter (SAM --> arg max as alpha --> +Inf, SAM --> arg min as alpha --> -Inf)
  real(rkind),intent(in) :: x(:) ! vector of input values
  ! output
  real(rkind),allocatable  :: SAM(:) ! SoftArgMax value 
  ! local variables
  real(real128) :: alpha_qp ! quadruple precision alpha
  real(real128) :: x_star   ! quadruple precision shift value for numerical stability
  real(real128),allocatable :: x_qp(:)   ! quadruple precision x vector
  real(real128),allocatable :: SAM_qp(:) ! quadruple precision SAM value 

  ! use quadruple precision variables to prevent over/underflow
  alpha_qp = real(alpha,real128)
  allocate(x_qp(size(x)))
  x_qp     = real(x,real128)

  ! shift value to improve numerical stability
  x_star = maxval(abs(x_qp))

  allocate(SAM_qp(size(x)))
  SAM_qp = exp(alpha_qp*(x_qp-x_star)) / sum(exp(alpha_qp*(x_qp-x_star)))
  SAM = real(SAM_qp,rkind)

  ! check if all values are normal (not NaN, -Infinity, or +Infinity)
  ! note: mainly to account for overflow/underflow that may occur in extreme cases
  if (all(ieee_is_normal(SAM))) then ! return if value is not NaN or infinity
    return
  else                          ! revert to analytic arg max/min function in one-hot representation as a failsafe (accurate but not smoothed)
    SAM(:) = 0._rkind
    if (alpha < 0._rkind) then  ! arg min
      SAM(minloc(x)) = 1._rkind
    else                        ! arg max
      SAM(maxloc(x)) = 1._rkind 
    end if
  end if
end function SoftArgMax

end module soil_utils_module

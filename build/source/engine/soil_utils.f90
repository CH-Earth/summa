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
USE nrtype

USE multiconst,only: gravity, & ! acceleration of gravity       (m s-2)
                     Tfreeze, & ! temperature at freezing    (K)
                     LH_fus,  & ! latent heat of fusion      (J kg-1, or m2 s-2)
                     R_wv       ! gas constant for water vapor  (J kg-1 K-1; [J = Pa m3])

! privacy
implicit none
private

! routines to make public
public::iceImpede
public::dIceImpede_dTemp
public::hydCond_psi
public::hydCond_liq
public::hydCondMP_liq
public::dHydCond_dPsi
public::dHydCond_dLiq
public::volFracLiq
public::matricHead
public::dTheta_dPsi
public::dPsi_dTheta
public::dPsi_dTheta2
public::RH_soilair
public::dTheta_dTk
public::crit_soilT
public::liquidHead
public::gammp

! constant parameters
real(rk),parameter     :: valueMissing=-9999._rk    ! missing value parameter
real(rk),parameter     :: verySmall=epsilon(1.0_rk) ! a very small number (used to avoid divide by zero)
real(rk),parameter     :: dx=-1.e-12_rk             ! finite difference increment
contains


 ! ******************************************************************************************************************************
 ! public subroutine iceImpede: compute the ice impedence factor
 ! ******************************************************************************************************************************
 subroutine iceImpede(volFracIce,f_impede, &            ! input
                      iceImpedeFactor,dIceImpede_dLiq)  ! output
 ! computes the ice impedence factor (separate function, as used multiple times)
 implicit none
 ! input variables
 real(rk),intent(in)     :: volFracIce        ! volumetric fraction of ice (-)
 real(rk),intent(in)     :: f_impede          ! ice impedence parameter (-)
 ! output variables
 real(rk)                :: iceImpedeFactor   ! ice impedence factor (-)
 real(rk)                :: dIceImpede_dLiq   ! derivative in ice impedence factor w.r.t. volumetric liquid water content (-)
 ! compute ice impedance factor as a function of volumetric ice content
 iceImpedeFactor = 10._rk**(-f_impede*volFracIce)
 dIceImpede_dLiq = 0._rk

 end subroutine iceImpede


 ! ******************************************************************************************************************************
 ! public subroutine dIceImpede_dTemp: compute the derivative in the ice impedence factor w.r.t. temperature
 ! ******************************************************************************************************************************
 subroutine dIceImpede_dTemp(volFracIce,dTheta_dT,f_impede,dIceImpede_dT)
 ! computes the derivative in the ice impedance factor w.r.t. temperature
 implicit none
 ! input variables
 real(rk),intent(in)     :: volFracIce        ! volumetric fraction of ice (-)
 real(rk),intent(in)     :: dTheta_dT         ! derivative in volumetric liquid water content w.r.t temperature (K-1)
 real(rk),intent(in)     :: f_impede          ! ice impedence parameter (-)
 ! output variables
 real(rk)                :: dIceImpede_dT     ! derivative in the ice impedance factor w.r.t. temperature (K-1)
 ! --
 dIceImpede_dT = log(10._rk)*f_impede*(10._rk**(-f_impede*volFracIce))*dTheta_dT
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
                       dVolTot_dPsi0                            ,& ! intent(in)    : derivative in the soil water characteristic (m-1)
                       dTheta_dT                                ,& ! intent(in)    : derivative in volumetric total water w.r.t. temperature (K-1)
                       ! output
                       matricHeadLiq                            ,& ! intent(out)   : liquid water matric potential (m)
                       dPsiLiq_dPsi0                            ,& ! intent(out)   : derivative in the liquid water matric potential w.r.t. the total water matric potential (-)
                       dPsiLiq_dTemp                            ,& ! intent(out)   : derivative in the liquid water matric potential w.r.t. temperature (m K-1)
                       err,message)                                ! intent(out)   : error control
 ! computes the liquid water matric potential (and the derivatives w.r.t. total matric potential and temperature)
 implicit none
 ! input
 real(rk),intent(in)            :: matricHeadTotal                           ! total water matric potential (m)
 real(rk),intent(in)            :: volFracLiq                                ! volumetric fraction of liquid water (-)
 real(rk),intent(in)            :: volFracIce                                ! volumetric fraction of ice (-)
 real(rk),intent(in)            :: vGn_alpha,vGn_n,theta_sat,theta_res,vGn_m ! soil parameters
 real(rk),intent(in)  ,optional :: dVolTot_dPsi0                             ! derivative in the soil water characteristic (m-1)
 real(rk),intent(in)  ,optional :: dTheta_dT                                 ! derivative in volumetric total water w.r.t. temperature (K-1)
 ! output
 real(rk),intent(out)           :: matricHeadLiq                             ! liquid water matric potential (m)
 real(rk),intent(out) ,optional :: dPsiLiq_dPsi0                             ! derivative in the liquid water matric potential w.r.t. the total water matric potential (-)
 real(rk),intent(out) ,optional :: dPsiLiq_dTemp                             ! derivative in the liquid water matric potential w.r.t. temperature (m K-1)
 ! output: error control
 integer(i4b),intent(out)       :: err                                       ! error code
 character(*),intent(out)       :: message                                   ! error message
 ! local
 real(rk)                       :: xNum,xDen                                 ! temporary variables (numeratir, denominator)
 real(rk)                       :: effSat                                    ! effective saturation (-)
 real(rk)                       :: dPsiLiq_dEffSat                           ! derivative in liquid water matric potential w.r.t. effective saturation (m)
 real(rk)                       :: dEffSat_dTemp                             ! derivative in effective saturation w.r.t. temperature (K-1)
 ! ------------------------------------------------------------------------------------------------------------------------------
 ! initialize error control
 err=0; message='liquidHead/'

 ! ** partially frozen soil
 if(volFracIce > verySmall .and. matricHeadTotal < 0._rk)then  ! check that ice exists and that the soil is unsaturated

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
  matricHeadLiq = matricHead(effSat,vGn_alpha,0._rk,1._rk,vGn_n,vGn_m)  ! argument is effective saturation, so theta_res=0 and theta_sat=1

  ! compute derivative in liquid water matric potential w.r.t. effective saturation (m)
  if(present(dPsiLiq_dPsi0).or.present(dPsiLiq_dTemp))then
   dPsiLiq_dEffSat = dPsi_dTheta(effSat,vGn_alpha,0._rk,1._rk,vGn_n,vGn_m)
  endif

  ! -----
  ! - compute derivative in the liquid water matric potential w.r.t. the total water matric potential...
  ! ----------------------------------------------------------------------------------------------------

  ! check if the derivative is desired
  if(present(dPsiLiq_dTemp))then

   ! (check required input derivative is present)
   if(.not.present(dVolTot_dPsi0))then
    message=trim(message)//'dVolTot_dPsi0 argument is missing'
    err=20; return
   endif

   ! (compute derivative in the liquid water matric potential w.r.t. the total water matric potential)
   dPsiLiq_dPsi0 = dVolTot_dPsi0*dPsiLiq_dEffSat*xNum/(xDen**2._rk)

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
   dEffSat_dTemp = -dTheta_dT*xNum/(xDen**2._rk) + dTheta_dT/xDen
   dPsiLiq_dTemp = dPsiLiq_dEffSat*dEffSat_dTemp

  endif  ! if dPsiLiq_dTemp is desired

 ! ** unfrozen soil
 else   ! (no ice)
  matricHeadLiq = matricHeadTotal
  if(present(dPsiLiq_dTemp)) dPsiLiq_dPsi0 = 1._rk  ! derivative=1 because values are identical
  if(present(dPsiLiq_dTemp)) dPsiLiq_dTemp = 0._rk  ! derivative=0 because no impact of temperature for unfrozen conditions
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
 real(rk),intent(in) :: volFracLiq    ! volumetric liquid water content (-)
 real(rk),intent(in) :: theta_sat     ! soil porosity (-)
 real(rk),intent(in) :: theta_mp      ! minimum volumetric liquid water content for macropore flow (-)
 real(rk),intent(in) :: mpExp         ! empirical exponent in macropore flow equation (-)
 real(rk),intent(in) :: satHydCond_ma ! saturated hydraulic conductivity for macropores (m s-1)
 real(rk),intent(in) :: satHydCond_mi ! saturated hydraulic conductivity for micropores (m s-1)
 real(rk)            :: hydCondMP_liq ! hydraulic conductivity (m s-1)
 ! locals
 real(rk)            :: theta_e     ! effective soil moisture
 if(volFracLiq > theta_mp)then
  theta_e       = (volFracLiq - theta_mp) / (theta_sat - theta_mp)
  hydCondMP_liq = (satHydCond_ma - satHydCond_mi) * (theta_e**mpExp)
 else
  hydCondMP_liq = 0._rk
 end if
 !write(*,'(a,4(f9.3,1x),2(e20.10))') 'in soil_utils: theta_mp, theta_sat, volFracLiq, hydCondMP_liq, satHydCond_ma, satHydCond_mi = ', &
 !                                                    theta_mp, theta_sat, volFracLiq, hydCondMP_liq, satHydCond_ma, satHydCond_mi
 end function hydCondMP_liq


 ! ******************************************************************************************************************************
 ! public function hydCond_psi: compute the hydraulic conductivity as a function of matric head (m s-1)
 ! ******************************************************************************************************************************
 function hydCond_psi(psi,k_sat,alpha,n,m)
 ! computes hydraulic conductivity given psi and soil hydraulic parameters k_sat, alpha, n, and m
 implicit none
 ! dummies
 real(rk),intent(in)     :: psi           ! soil water suction (m)
 real(rk),intent(in)     :: k_sat         ! saturated hydraulic conductivity (m s-1)
 real(rk),intent(in)     :: alpha         ! scaling parameter (m-1)
 real(rk),intent(in)     :: n             ! vGn "n" parameter (-)
 real(rk),intent(in)     :: m             ! vGn "m" parameter (-)
 real(rk)                :: hydCond_psi   ! hydraulic conductivity (m s-1)
 if(psi<0._rk)then
  hydCond_psi = k_sat * &
                ( ( (1._rk - (psi*alpha)**(n-1._rk) * (1._rk + (psi*alpha)**n)**(-m))**2._rk ) &
                 / ( (1._rk + (psi*alpha)**n)**(m/2._rk) ) )
 else
  hydCond_psi = k_sat
 end if
 end function hydCond_psi


 ! ******************************************************************************************************************************
 ! public function hydCond_liq: compute the hydraulic conductivity as a function of volumetric liquid water content (m s-1)
 ! ******************************************************************************************************************************
 function hydCond_liq(volFracLiq,k_sat,theta_res,theta_sat,m)
 ! computes hydraulic conductivity given volFracLiq and soil hydraulic parameters k_sat, theta_sat, theta_res, and m
 implicit none
 ! dummies
 real(rk),intent(in) :: volFracLiq  ! volumetric liquid water content (-)
 real(rk),intent(in) :: k_sat       ! saturated hydraulic conductivity (m s-1)
 real(rk),intent(in) :: theta_res   ! residual volumetric liquid water content (-)
 real(rk),intent(in) :: theta_sat   ! soil porosity (-)
 real(rk),intent(in) :: m           ! vGn "m" parameter (-)
 real(rk)            :: hydCond_liq ! hydraulic conductivity (m s-1)
 ! locals
 real(rk)            :: theta_e     ! effective soil moisture
 if(volFracLiq < theta_sat)then
  theta_e = (volFracLiq - theta_res) / (theta_sat - theta_res)
  hydCond_liq = k_sat*theta_e**(1._rk/2._rk) * (1._rk - (1._rk - theta_e**(1._rk/m) )**m)**2._rk
 else
  hydCond_liq = k_sat
 end if
 end function hydCond_liq


 ! ******************************************************************************************************************************
 ! public function volFracLiq: compute the volumetric liquid water content (-)
 ! ******************************************************************************************************************************
 function volFracLiq(psi,alpha,theta_res,theta_sat,n,m)
 ! computes the volumetric liquid water content given psi and soil hydraulic parameters theta_res, theta_sat, alpha, n, and m
 implicit none
 real(rk),intent(in) :: psi         ! soil water suction (m)
 real(rk),intent(in) :: alpha       ! scaling parameter (m-1)
 real(rk),intent(in) :: theta_res   ! residual volumetric water content (-)
 real(rk),intent(in) :: theta_sat   ! porosity (-)
 real(rk),intent(in) :: n           ! vGn "n" parameter (-)
 real(rk),intent(in) :: m           ! vGn "m" parameter (-)
 real(rk)            :: volFracLiq  ! volumetric liquid water content (-)
 if(psi<0._rk)then
  volFracLiq = theta_res + (theta_sat - theta_res)*(1._rk + (alpha*psi)**n)**(-m)
 else
  volFracLiq = theta_sat
 end if
 end function volFracLiq


 ! ******************************************************************************************************************************
 ! public function matricHead: compute the matric head (m) based on the volumetric liquid water content
 ! ******************************************************************************************************************************
 function matricHead(theta,alpha,theta_res,theta_sat,n,m)
 ! computes the volumetric liquid water content given psi and soil hydraulic parameters theta_res, theta_sat, alpha, n, and m
 implicit none
 ! dummy variables
 real(rk),intent(in) :: theta       ! volumetric liquid water content (-)
 real(rk),intent(in) :: alpha       ! scaling parameter (m-1)
 real(rk),intent(in) :: theta_res   ! residual volumetric water content (-)
 real(rk),intent(in) :: theta_sat   ! porosity (-)
 real(rk),intent(in) :: n           ! vGn "n" parameter (-)
 real(rk),intent(in) :: m           ! vGn "m" parameter (-)
 real(rk)            :: matricHead  ! matric head (m)
 ! local variables
 real(rk)            :: effSat      ! effective saturation (-)
 real(rk),parameter  :: verySmall=epsilon(1._rk)  ! a very small number (avoid effective saturation of zero)
 ! compute effective saturation
 effSat = max(verySmall, (theta - theta_res) / (theta_sat - theta_res))
 ! compute matric head
 if (effSat < 1._rk .and. effSat > 0._rk)then
  matricHead = (1._rk/alpha)*( effSat**(-1._rk/m) - 1._rk)**(1._rk/n)
 else
  matricHead = 0._rk
 end if
 end function matricHead


 ! ******************************************************************************************************************************
 ! public function dTheta_dPsi: compute the derivative of the soil water characteristic (m-1)
 ! ******************************************************************************************************************************
 function dTheta_dPsi(psi,alpha,theta_res,theta_sat,n,m)
 implicit none
 real(rk),intent(in) :: psi         ! soil water suction (m)
 real(rk),intent(in) :: alpha       ! scaling parameter (m-1)
 real(rk),intent(in) :: theta_res   ! residual volumetric water content (-)
 real(rk),intent(in) :: theta_sat   ! porosity (-)
 real(rk),intent(in) :: n           ! vGn "n" parameter (-)
 real(rk),intent(in) :: m           ! vGn "m" parameter (-)
 real(rk)            :: dTheta_dPsi ! derivative of the soil water characteristic (m-1)
 if(psi<=0._rk)then
  dTheta_dPsi = (theta_sat-theta_res) * &
     (-m*(1._rk + (psi*alpha)**n)**(-m-1._rk)) * n*(psi*alpha)**(n-1._rk) * alpha
  if(abs(dTheta_dPsi) < epsilon(psi)) dTheta_dPsi = epsilon(psi)
 else
  dTheta_dPsi = epsilon(psi)
 end if
 end function dTheta_dPsi


 ! ******************************************************************************************************************************
 ! public function dPsi_dTheta: compute the derivative of the soil water characteristic (m-1)
 ! ******************************************************************************************************************************
 function dPsi_dTheta(volFracLiq,alpha,theta_res,theta_sat,n,m)
 implicit none
 ! dummies
 real(rk),intent(in) :: volFracLiq  ! volumetric liquid water content (-)
 real(rk),intent(in) :: alpha       ! scaling parameter (m-1)
 real(rk),intent(in) :: theta_res   ! residual volumetric water content (-)
 real(rk),intent(in) :: theta_sat   ! porosity (-)
 real(rk),intent(in) :: n           ! vGn "n" parameter (-)
 real(rk),intent(in) :: m           ! vGn "m" parameter (-)
 real(rk)            :: dPsi_dTheta ! derivative of the soil water characteristic (m)
 ! locals
 real(rk)            :: y1,d1       ! 1st function and derivative
 real(rk)            :: y2,d2       ! 2nd function and derivative
 real(rk)            :: theta_e     ! effective soil moisture
 ! check if less than saturation
 if(volFracLiq < theta_sat)then
  ! compute effective water content
  theta_e = max(0.001,(volFracLiq - theta_res) / (theta_sat - theta_res))
  ! compute the 1st function and derivative
  y1 = theta_e**(-1._rk/m) - 1._rk
  d1 = (-1._rk/m)*theta_e**(-1._rk/m - 1._rk) / (theta_sat - theta_res)
  ! compute the 2nd function and derivative
  y2 = y1**(1._rk/n)
  d2 = (1._rk/n)*y1**(1._rk/n - 1._rk)
  ! compute the final function value
  dPsi_dTheta = d1*d2/alpha
 else
  dPsi_dTheta = 0._rk
 end if
 end function dPsi_dTheta


 ! ******************************************************************************************************************************
 ! public function dPsi_dTheta2: compute the derivative of dPsi_dTheta (m-1)
 ! ******************************************************************************************************************************
 function dPsi_dTheta2(volFracLiq,alpha,theta_res,theta_sat,n,m,lTangent)
 implicit none
 ! dummies
 real(rk),intent(in)     :: volFracLiq   ! volumetric liquid water content (-)
 real(rk),intent(in)     :: alpha        ! scaling parameter (m-1)
 real(rk),intent(in)     :: theta_res    ! residual volumetric water content (-)
 real(rk),intent(in)     :: theta_sat    ! porosity (-)
 real(rk),intent(in)     :: n            ! vGn "n" parameter (-)
 real(rk),intent(in)     :: m            ! vGn "m" parameter (-)
 logical(lgt),intent(in) :: lTangent     ! method used to compute derivative (.true. = analytical)
 real(rk)                :: dPsi_dTheta2 ! derivative of the soil water characteristic (m)
 ! locals for analytical derivatives
 real(rk)                :: xx           ! temporary variable
 real(rk)                :: y1,d1        ! 1st function and derivative
 real(rk)                :: y2,d2        ! 2nd function and derivative
 real(rk)                :: theta_e      ! effective soil moisture
 ! locals for numerical derivative
 real(rk)                :: func0,func1  ! function evaluations
 ! check if less than saturation
 if(volFracLiq < theta_sat)then
  ! ***** compute analytical derivatives
  if(lTangent)then
   ! compute the effective saturation
   theta_e = (volFracLiq - theta_res) / (theta_sat - theta_res)
   ! get the first function and derivative
   y1 = (-1._rk/m)*theta_e**(-1._rk/m - 1._rk) / (theta_sat - theta_res)
   d1 = ( (m + 1._rk) / (m**2._rk * (theta_sat - theta_res)**2._rk) ) * theta_e**(-1._rk/m - 2._rk)
   ! get the second function and derivative
   xx = theta_e**(-1._rk/m) - 1._rk
   y2 = (1._rk/n)*xx**(1._rk/n - 1._rk)
   d2 = ( -(1._rk - n)/((theta_sat - theta_res)*m*n**2._rk) ) * xx**(1._rk/n - 2._rk) * theta_e**(-1._rk/m - 1._rk)
   ! return the derivative
   dPsi_dTheta2 = (d1*y2 + y1*d2)/alpha
  ! ***** compute numerical derivatives
  else
   func0 = dPsi_dTheta(volFracLiq,   alpha,theta_res,theta_sat,n,m)
   func1 = dPsi_dTheta(volFracLiq+dx,alpha,theta_res,theta_sat,n,m)
   dPsi_dTheta2 = (func1 - func0)/dx
  end if
 ! (case where volumetric liquid water content exceeds porosity)
 else
  dPsi_dTheta2 = 0._rk
 end if
 end function dPsi_dTheta2


 ! ******************************************************************************************************************************
 ! public function dHydCond_dPsi: compute the derivative in hydraulic conductivity w.r.t. matric head (s-1)
 ! ******************************************************************************************************************************
 function dHydCond_dPsi(psi,k_sat,alpha,n,m,lTangent)
 ! computes the derivative in hydraulic conductivity w.r.t matric head,
 !  given psi and soil hydraulic parameters k_sat, alpha, n, and m
 implicit none
 ! dummies
 real(rk),intent(in)     :: psi        ! soil water suction (m)
 real(rk),intent(in)     :: k_sat      ! saturated hydraulic conductivity (m s-1)
 real(rk),intent(in)     :: alpha      ! scaling parameter (m-1)
 real(rk),intent(in)     :: n          ! vGn "n" parameter (-)
 real(rk),intent(in)     :: m          ! vGn "m" parameter (-)
 logical(lgt),intent(in) :: lTangent   ! method used to compute derivative (.true. = analytical)
 real(rk)                :: dHydCond_dPsi  ! derivative in hydraulic conductivity w.r.t. matric head (s-1)
 ! locals for analytical derivatives
 real(rk)                :: f_x1          ! f(x) for part of the numerator
 real(rk)                :: f_x2          ! f(x) for part of the numerator
 real(rk)                :: f_nm          ! f(x) for the numerator
 real(rk)                :: f_dm          ! f(x) for the denominator
 real(rk)                :: d_x1          ! df(x)/dpsi for part of the numerator
 real(rk)                :: d_x2          ! df(x)/dpsi for part of the numerator
 real(rk)                :: d_nm          ! df(x)/dpsi for the numerator
 real(rk)                :: d_dm          ! df(x)/dpsi for the denominator
 ! locals for numerical derivatives
 real(rk)                :: hydCond0   ! hydraulic condictivity value for base case
 real(rk)                :: hydCond1   ! hydraulic condictivity value for perturbed case
 ! derivative is zero if saturated
 if(psi<0._rk)then
  ! ***** compute analytical derivatives
  if(lTangent)then
   ! compute the derivative for the numerator
   f_x1 = (psi*alpha)**(n - 1._rk)
   f_x2 = (1._rk + (psi*alpha)**n)**(-m)
   d_x1 = alpha * (n - 1._rk)*(psi*alpha)**(n - 2._rk)
   d_x2 = alpha * n*(psi*alpha)**(n - 1._rk) * (-m)*(1._rk + (psi*alpha)**n)**(-m - 1._rk)
   f_nm = (1._rk - f_x1*f_x2)**2._rk
   d_nm = (-d_x1*f_x2 - f_x1*d_x2) * 2._rk*(1._rk - f_x1*f_x2)
   ! compute the derivative for the denominator
   f_dm = (1._rk + (psi*alpha)**n)**(m/2._rk)
   d_dm = alpha * n*(psi*alpha)**(n - 1._rk) * (m/2._rk)*(1._rk + (psi*alpha)**n)**(m/2._rk - 1._rk)
   ! and combine
   dHydCond_dPsi = k_sat*(d_nm*f_dm - d_dm*f_nm) / (f_dm**2._rk)
  else
   ! ***** compute numerical derivatives
   hydcond0  = hydCond_psi(psi,   k_sat,alpha,n,m)
   hydcond1  = hydCond_psi(psi+dx,k_sat,alpha,n,m)
   dHydCond_dPsi = (hydcond1 - hydcond0)/dx
  end if
 else
  dHydCond_dPsi = 0._rk
 end if
 end function dHydCond_dPsi


 ! ******************************************************************************************************************************
 ! public function dHydCond_dLiq: compute the derivative in hydraulic conductivity w.r.t. volumetric liquid water content (m s-1)
 ! ******************************************************************************************************************************
 ! computes the derivative in hydraulic conductivity w.r.t the volumetric fraction of liquid water,
 ! given volFracLiq and soil hydraulic parameters k_sat, theta_sat, theta_res, and m
 ! ******************************************************************************************************************************
 function dHydCond_dLiq(volFracLiq,k_sat,theta_res,theta_sat,m,lTangent)
 implicit none
 ! dummies
 real(rk),intent(in)     :: volFracLiq ! volumetric fraction of liquid water (-)
 real(rk),intent(in)     :: k_sat      ! saturated hydraulic conductivity (m s-1)
 real(rk),intent(in)     :: theta_res  ! soil residual volumetric water content (-)
 real(rk),intent(in)     :: theta_sat  ! soil porosity (-)
 real(rk),intent(in)     :: m          ! vGn "m" parameter (-)
 logical(lgt),intent(in) :: lTangent   ! method used to compute derivative (.true. = analytical)
 real(rk)                :: dHydCond_dLiq  ! derivative in hydraulic conductivity w.r.t. matric head (s-1)
 ! locals for analytical derivatives
 real(rk)                :: theta_e  ! effective soil moisture
 real(rk)                :: f1       ! f(x) for the first function
 real(rk)                :: d1       ! df(x)/dLiq for the first function
 real(rk)                :: x1,x2    ! f(x) for different parts of the second function
 real(rk)                :: p1,p2,p3 ! df(x)/dLiq for different parts of the second function
 real(rk)                :: f2       ! f(x) for the second function
 real(rk)                :: d2       ! df(x)/dLiq for the second function
 ! locals for numerical derivatives
 real(rk)                :: hydCond0 ! hydraulic condictivity value for base case
 real(rk)                :: hydCond1 ! hydraulic condictivity value for perturbed case
 ! derivative is zero if super-saturated
 if(volFracLiq < theta_sat)then
  ! ***** compute analytical derivatives
  if(lTangent)then
   ! compute the effective saturation
   theta_e = (volFracLiq - theta_res) / (theta_sat - theta_res)
   ! compute the function and derivative of the first fuction
   f1 = k_sat*theta_e**0.5_rk
   d1 = k_sat*0.5_rk*theta_e**(-0.5_rk) / (theta_sat - theta_res)
   ! compute the function and derivative of the second function
   ! (first part)
   x1 = 1._rk - theta_e**(1._rk/m)
   p1 = (-1._rk/m)*theta_e**(1._rk/m - 1._rk) / (theta_sat - theta_res)   ! differentiate (1.d - theta_e**(1.d/m)
   ! (second part)
   x2 = x1**m
   p2 = m*x1**(m - 1._rk)
   ! (final)
   f2 = (1._rk - x2)**2._rk
   p3 = -2._rk*(1._rk - x2)
   ! (combine)
   d2 = p1*p2*p3
   ! pull it all together
   dHydCond_dLiq = (d1*f2 + d2*f1)
  else
   ! ***** compute numerical derivatives
   hydcond0 = hydCond_liq(volFracLiq,   k_sat,theta_res,theta_sat,m)
   hydcond1 = hydCond_liq(volFracLiq+dx,k_sat,theta_res,theta_sat,m)
   dHydCond_dLiq = (hydcond1 - hydcond0)/dx
  end if
 else
  dHydCond_dLiq = 0._rk
 end if
 end function dHydCond_dLiq


 ! ******************************************************************************************************************************
 ! public function RH_soilair: compute relative humidity of air in soil pore space
 ! ******************************************************************************************************************************
 function RH_soilair(matpot,Tk)
 implicit none
 real(rk),intent(in) :: matpot        ! soil water suction -- matric potential (m)
 real(rk),intent(in) :: Tk            ! temperature (K)
 real(rk)            :: RH_soilair    ! relative humidity of air in soil pore space
 ! compute relative humidity (UNITS NOTE: Pa = kg m-1 s-2, so R_wv units = m2 s-2 K-1)
 RH_soilair = exp( (gravity*matpot) / (R_wv*Tk) )
 end function RH_soilair


 ! ******************************************************************************************************************************
 ! public function crit_soilT: compute the critical temperature above which all water is unfrozen
 ! ******************************************************************************************************************************
 function crit_soilT(psi)
 implicit none
 real(rk),intent(in) :: psi           ! matric head (m)
 real(rk)            :: crit_soilT    ! critical soil temperature (K)
 crit_soilT = Tfreeze + min(psi,0._rk)*gravity*Tfreeze/LH_fus
 end function crit_soilT


 ! ******************************************************************************************************************************
 ! public function dTheta_dTk: differentiate the freezing curve w.r.t. temperature
 ! ******************************************************************************************************************************
 function dTheta_dTk(Tk,theta_res,theta_sat,alpha,n,m)
 implicit none
 real(rk),intent(in) :: Tk            ! temperature (K)
 real(rk),intent(in) :: theta_res     ! residual liquid water content (-)
 real(rk),intent(in) :: theta_sat     ! porosity (-)
 real(rk),intent(in) :: alpha         ! vGn scaling parameter (m-1)
 real(rk),intent(in) :: n             ! vGn "n" parameter (-)
 real(rk),intent(in) :: m             ! vGn "m" parameter (-)
 real(rk)            :: dTheta_dTk    ! derivative of the freezing curve w.r.t. temperature (K-1)
 ! local variables
 real(rk)            :: kappa         ! constant (m K-1)
 real(rk)            :: xtemp         ! alpha*kappa*(Tk-Tfreeze) -- dimensionless variable (used more than once)
 ! compute kappa (m K-1)
 kappa =  LH_fus/(gravity*Tfreeze)    ! NOTE: J = kg m2 s-2
 ! define a tempory variable that is used more than once (-)
 xtemp = alpha*kappa*(Tk-Tfreeze)
 ! differentiate the freezing curve w.r.t. temperature -- making use of the chain rule
 dTheta_dTk = (alpha*kappa) * n*xtemp**(n - 1._rk) * (-m)*(1._rk + xtemp**n)**(-m - 1._rk) * (theta_sat - theta_res)
 end function dTheta_dTk


 ! ******************************************************************************************************************************
 ! public function gammp: compute cumulative probability using the Gamma distribution
 ! ******************************************************************************************************************************
 FUNCTION gammp(a,x)
 IMPLICIT NONE
 real(rk), INTENT(IN) :: a,x
 real(rk) :: gammp
 if (x<a+1.0_rk) then
  gammp=gser(a,x)
 else
  gammp=1.0_rk-gcf(a,x)
 end if
 END FUNCTION gammp


 ! ******************************************************************************************************************************
 ! private function gcf: continued fraction development of the incomplete Gamma function
 ! ******************************************************************************************************************************
 FUNCTION gcf(a,x,gln)
 IMPLICIT NONE
 real(rk), INTENT(IN) :: a,x
 real(rk), OPTIONAL, INTENT(OUT) :: gln
 real(rk) :: gcf
 INTEGER(I4B), PARAMETER :: ITMAX=100
 real(rk), PARAMETER :: EPS=epsilon(x),FPMIN=tiny(x)/EPS
 INTEGER(I4B) :: i
 real(rk) :: an,b,c,d,del,h
 if (x == 0.0) then
  gcf=1.0
  RETURN
 end if
 b=x+1.0_rk-a
 c=1.0_rk/FPMIN
 d=1.0_rk/b
 h=d
 do i=1,ITMAX
  an=-i*(i-a)
  b=b+2.0_rk
  d=an*d+b
  if (abs(d) < FPMIN) d=FPMIN
  c=b+an/c
  if (abs(c) < FPMIN) c=FPMIN
  d=1.0_rk/d
  del=d*c
  h=h*del
  if (abs(del-1.0_rk) <= EPS) exit
 end do
 if (i > ITMAX) stop 'a too large, ITMAX too small in gcf'
 if (present(gln)) then
  gln=gammln(a)
  gcf=exp(-x+a*log(x)-gln)*h
 else
  gcf=exp(-x+a*log(x)-gammln(a))*h
 end if
 END FUNCTION gcf


 ! ******************************************************************************************************************************
 ! private function gser: series development of the incomplete Gamma function
 ! ******************************************************************************************************************************
 FUNCTION gser(a,x,gln)
 IMPLICIT NONE
 real(rk), INTENT(IN) :: a,x
 real(rk), OPTIONAL, INTENT(OUT) :: gln
 real(rk) :: gser
 INTEGER(I4B), PARAMETER :: ITMAX=100
 real(rk), PARAMETER :: EPS=epsilon(x)
 INTEGER(I4B) :: n
 real(rk) :: ap,del,summ
 if (x == 0.0) then
  gser=0.0
  RETURN
 end if
 ap=a
 summ=1.0_rk/a
 del=summ
 do n=1,ITMAX
  ap=ap+1.0_rk
  del=del*x/ap
  summ=summ+del
  if (abs(del) < abs(summ)*EPS) exit
 end do
 if (n > ITMAX) stop 'a too large, ITMAX too small in gser'
 if (present(gln)) then
  gln=gammln(a)
  gser=summ*exp(-x+a*log(x)-gln)
 else
  gser=summ*exp(-x+a*log(x)-gammln(a))
 end if
 END FUNCTION gser


 ! ******************************************************************************************************************************
 ! private function gammln: gamma function
 ! ******************************************************************************************************************************
 FUNCTION gammln(xx)
 USE nr_utility_module,only:arth  ! use to build vectors with regular increments
 IMPLICIT NONE
 real(rk), INTENT(IN) :: xx
 real(rk) :: gammln
 real(rk) :: tmp,x
 real(rk) :: stp = 2.5066282746310005_rk
 real(rk), DIMENSION(6) :: coef = (/76.18009172947146_rk,&
  -86.50532032941677_rk,24.01409824083091_rk,&
  -1.231739572450155_rk,0.1208650973866179e-2_rk,&
  -0.5395239384953e-5_rk/)
 if(xx <= 0._rk) stop 'xx > 0 in gammln'
 x=xx
 tmp=x+5.5_rk
 tmp=(x+0.5_rk)*log(tmp)-tmp
 gammln=tmp+log(stp*(1.000000000190015_rk+&
  sum(coef(:)/arth(x+1.0_rk,1.0_rk,size(coef))))/x)
 END FUNCTION gammln


end module soil_utils_module

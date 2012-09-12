module soil_utils_module
USE nrtype
implicit none
private
public::iceImpede
public::satArea_matric
public::satArea_liquid
public::hydCond_psi
public::hydCond_liq
public::darcyFlux_matric
public::darcyFlux_liquid
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
contains

 ! ***********************************************************************************************************
 ! new function: compute the saturated area based on matric head
 ! ***********************************************************************************************************
 function satArea_matric(matricHead,volFracIce,theta_res,theta_sat,alpha,n,m,vic_bpar)
 implicit none
 ! dummy variables
 real(dp),intent(in) :: matricHead  ! matricHead (m)
 real(dp),intent(in) :: volFracIce  ! volumetric fraction of ice (-)
 real(dp),intent(in) :: theta_res   ! residual volumetric liquid water content (-)
 real(dp),intent(in) :: theta_sat   ! soil porosity (-)
 real(dp),intent(in) :: alpha       ! scaling parameter (m-1)
 real(dp),intent(in) :: n           ! vGn "n" parameter (-)
 real(dp),intent(in) :: m           ! vGn "m" parameter (-)
 real(dp),intent(in) :: vic_bpar    ! VIC "b" parameter (dimensionless surface runoff exponent)
 real(dp)            :: satArea_matric
 ! local variables
 real(dp)            :: vFracLiq    ! volumetric fraction of liquid water (-)
 real(dp)            :: fracCap     ! fractional storage capacity (-)
 ! compute volumetric fraction of liquid water
 vFracLiq = volFracLiq(matricHead,alpha,theta_res,theta_sat,n,m)
 ! compute fraction of pore space filled with liquid water and ice
 fracCap  = min((vFracLiq + volFracIce)/theta_sat, 1._dp)  ! fraction of pore space filled with liquid water and ice
 ! compute the saturated area (-)
 satArea_matric = 1._dp - (1._dp - fracCap)**vic_bpar      ! saturated fraction (-)
 end function satArea_matric


 ! ***********************************************************************************************************
 ! new function: compute the saturated area based on volumetric liquid water content
 ! ***********************************************************************************************************
 function satArea_liquid(volFracLiq,volFracIce,theta_sat,vic_bpar)
 implicit none
 ! dummy variables
 real(dp),intent(in) :: volFracLiq  ! volumetric fraction of liquid water (-)
 real(dp),intent(in) :: volFracIce  ! volumetric fraction of ice (-)
 real(dp),intent(in) :: theta_sat   ! soil porosity (-)
 real(dp),intent(in) :: vic_bpar    ! VIC "b" parameter (dimensionless surface runoff exponent)
 real(dp)            :: satArea_liquid
 ! local variables
 real(dp)            :: fracCap     ! fractional storage capacity (-)
 ! compute fraction of pore space filled with liquid water and ice
 fracCap    = min((volFracLiq + volFracIce)/theta_sat, 1._dp)  ! fraction of pore space filled with liquid water and ice
 ! compute the saturated area (-)
 satArea_liquid = 1._dp - (1._dp - fracCap)**vic_bpar          ! saturated fraction (-)
 end function satArea_liquid


 ! ***********************************************************************************************************
 ! new function: compute the ice impedence factor
 ! ***********************************************************************************************************
 function iceImpede(volFracIce,theta_res,theta_sat,f_impede)
 ! computes the ice impedence factor (separate function, as used multiple times)
 implicit none
 ! dummy variables
 real(dp),intent(in) :: volFracIce  ! volumetric fraction of ice
 real(dp),intent(in) :: theta_res   ! residual volumetric liquid water content (-)
 real(dp),intent(in) :: theta_sat   ! soil porosity (-)
 real(dp),intent(in) :: f_impede    ! ice impedence parameter (-)
 real(dp)            :: iceImpede   ! ice impedence factor (-)
 ! local variables
 real(dp)            :: relIce      ! relative fraction of ice
 ! compute relative fraction of ice
 relIce = volFracIce/(theta_sat - theta_res)
 ! compute the ice impedence factor
 iceImpede = 10._dp**(-f_impede*relIce)
 end function iceImpede


 ! ***********************************************************************************************************
 ! new function: compute the hydraulic conductivity as a function of matric head (m s-1)
 ! ***********************************************************************************************************
 function hydCond_psi(psi,k_sat,alpha,n,m)
 ! computes hydraulic conductivity given psi and soil hydraulic parameters k_sat, alpha, n, and m
 implicit none
 real(dp),intent(in) :: psi      ! soil water suction (m)
 real(dp),intent(in) :: k_sat    ! saturated hydraulic conductivity (m s-1)
 real(dp),intent(in) :: alpha    ! scaling parameter (m-1)
 real(dp),intent(in) :: n        ! vGn "n" parameter (-)
 real(dp),intent(in) :: m        ! vGn "m" parameter (-)
 real(dp)            :: hydCond_psi  ! hydraulic conductivity (m s-1)
 if(psi<0._dp)then
  hydCond_psi = k_sat * &
                ( ( (1._dp - (psi*alpha)**(n-1._dp) * (1._dp + (psi*alpha)**n)**(-m))**2._dp ) &
                 / ( (1._dp + (psi*alpha)**n)**(m/2._dp) ) )
 else
  hydCond_psi = k_sat
 endif
 end function hydCond_psi


 ! ***********************************************************************************************************
 ! new function: compute the hydraulic conductivity as a function of volumetric liquid water content (m s-1)
 ! ***********************************************************************************************************
 function hydCond_liq(volFracLiq,k_sat,theta_res,theta_sat,m)
 ! computes hydraulic conductivity given volFracLiq and soil hydraulic parameters k_sat, theta_sat, theta_res, and m
 implicit none
 ! dummies
 real(dp),intent(in) :: volFracLiq  ! volumetric liquid water content (-)
 real(dp),intent(in) :: k_sat       ! saturated hydraulic conductivity (m s-1)
 real(dp),intent(in) :: theta_res   ! residual volumetric liquid water content (-)
 real(dp),intent(in) :: theta_sat   ! soil porosity (-)
 real(dp),intent(in) :: m           ! vGn "m" parameter (-)
 real(dp)            :: hydCond_liq ! hydraulic conductivity (m s-1)
 ! locals
 real(dp)            :: theta_e     ! effective soil moisture
 if(volFracLiq < theta_sat)then
  theta_e = (volFracLiq - theta_res) / (theta_sat - theta_res)
  hydCond_liq = k_sat*theta_e**(1._dp/2._dp) * (1._dp - (1._dp - theta_e**(1._dp/m) )**m)**2._dp
 else
  hydCond_liq = k_sat
 endif
 end function hydCond_liq


 ! ***********************************************************************************************************
 ! new function: compute the volumetric liquid water content (-)
 ! ***********************************************************************************************************
 function volFracLiq(psi,alpha,theta_res,theta_sat,n,m)
 ! computes the volumetric liquid water content given psi and soil hydraulic parameters theta_res, theta_sat, alpha, n, and m
 implicit none
 real(dp),intent(in) :: psi         ! soil water suction (m)
 real(dp),intent(in) :: alpha       ! scaling parameter (m-1)
 real(dp),intent(in) :: theta_res   ! residual volumetric water content (-)
 real(dp),intent(in) :: theta_sat   ! porosity (-)
 real(dp),intent(in) :: n           ! vGn "n" parameter (-)
 real(dp),intent(in) :: m           ! vGn "m" parameter (-)
 real(dp)            :: volFracLiq  ! volumetric liquid water content (-)
 if(psi<0._dp)then
  volFracLiq = theta_res + (theta_sat - theta_res)*(1._dp + (alpha*psi)**n)**(-m)
 else
  volFracLiq = theta_sat
 endif
 end function volFracLiq


 ! ***********************************************************************************************************
 ! new function: compute the matric head (m) based on the volumetric liquid water content
 ! ***********************************************************************************************************
 function matricHead(theta,alpha,theta_res,theta_sat,n,m)
 ! computes the volumetric liquid water content given psi and soil hydraulic parameters theta_res, theta_sat, alpha, n, and m
 implicit none
 real(dp),intent(in) :: theta       ! volumetric liquid water content (-)
 real(dp),intent(in) :: alpha       ! scaling parameter (m-1)
 real(dp),intent(in) :: theta_res   ! residual volumetric water content (-)
 real(dp),intent(in) :: theta_sat   ! porosity (-)
 real(dp),intent(in) :: n           ! vGn "n" parameter (-)
 real(dp),intent(in) :: m           ! vGn "m" parameter (-)
 real(dp)            :: matricHead  ! matric head (m)
 matricHead = (1._dp/alpha)*( ( (theta - theta_res) / (theta_sat - theta_res) )**(-1._dp/m) - 1._dp)**(1._dp/n)
 end function matricHead


 ! ***********************************************************************************************************
 ! new function: compute the derivative of the soil water characteristic (m-1)
 ! ***********************************************************************************************************
 function dTheta_dPsi(psi,alpha,theta_res,theta_sat,n,m)
 implicit none
 real(dp),intent(in) :: psi         ! soil water suction (m)
 real(dp),intent(in) :: alpha       ! scaling parameter (m-1)
 real(dp),intent(in) :: theta_res   ! residual volumetric water content (-)
 real(dp),intent(in) :: theta_sat   ! porosity (-)
 real(dp),intent(in) :: n           ! vGn "n" parameter (-)
 real(dp),intent(in) :: m           ! vGn "m" parameter (-)
 real(dp)            :: dTheta_dPsi ! derivative of the soil water characteristic (m-1)
 if(psi<=0._dp)then
  dTheta_dPsi = (theta_sat-theta_res) * &
     (-m*(1._dp + (psi*alpha)**n)**(-m-1._dp)) * n*(psi*alpha)**(n-1._dp) * alpha
 else
  dTheta_dPsi = 0._dp
 endif
 end function dTheta_dPsi


 ! ***********************************************************************************************************
 ! new function: compute the derivative of the soil water characteristic (m-1)
 ! ***********************************************************************************************************
 function dPsi_dTheta(volFracLiq,alpha,theta_res,theta_sat,n,m)
 implicit none
 ! dummies
 real(dp),intent(in) :: volFracLiq  ! volumetric liquid water content (-)
 real(dp),intent(in) :: alpha       ! scaling parameter (m-1)
 real(dp),intent(in) :: theta_res   ! residual volumetric water content (-)
 real(dp),intent(in) :: theta_sat   ! porosity (-)
 real(dp),intent(in) :: n           ! vGn "n" parameter (-)
 real(dp),intent(in) :: m           ! vGn "m" parameter (-)
 real(dp)            :: dPsi_dTheta ! derivative of the soil water characteristic (m)
 ! locals
 real(dp)            :: y1,d1       ! 1st function and derivative
 real(dp)            :: y2,d2       ! 2nd function and derivative
 real(dp)            :: theta_e     ! effective soil moisture
 ! check if less than saturation
 if(volFracLiq < theta_sat)then
  ! compute effective water content
  theta_e = (volFracLiq - theta_res) / (theta_sat - theta_res)
  ! compute the 1st function and derivative
  y1 = theta_e**(-1._dp/m) - 1._dp
  d1 = (-1._dp/m)*theta_e**(-1._dp/m - 1._dp) / (theta_sat - theta_res)
  ! compute the 2nd function and derivative
  y2 = y1**(1._dp/n)
  d2 = (1._dp/n)*y1**(1._dp/n - 1._dp)
  ! compute the final function value
  dPsi_dTheta = d1*d2/alpha
 else
  dPsi_dTheta = 0._dp
 endif
 end function dPsi_dTheta

 ! ***********************************************************************************************************
 ! new function: compute the derivative of dPsi_dTheta (m-1)
 ! ***********************************************************************************************************
 function dPsi_dTheta2(volFracLiq,alpha,theta_res,theta_sat,n,m,lTangent)
 implicit none
 ! dummies
 real(dp),intent(in)     :: volFracLiq   ! volumetric liquid water content (-)
 real(dp),intent(in)     :: alpha        ! scaling parameter (m-1)
 real(dp),intent(in)     :: theta_res    ! residual volumetric water content (-)
 real(dp),intent(in)     :: theta_sat    ! porosity (-)
 real(dp),intent(in)     :: n            ! vGn "n" parameter (-)
 real(dp),intent(in)     :: m            ! vGn "m" parameter (-)
 logical(lgt),intent(in) :: lTangent     ! method used to compute derivative (.true. = analytical)
 real(dp)                :: dPsi_dTheta2 ! derivative of the soil water characteristic (m)
 ! locals for analytical derivatives
 real(dp)                :: xx           ! temporary variable
 real(dp)                :: y1,d1        ! 1st function and derivative
 real(dp)                :: y2,d2        ! 2nd function and derivative
 real(dp)                :: theta_e      ! effective soil moisture
 ! locals for numerical derivative
 real(dp)                :: func0,func1  ! function evaluations
 real(dp),parameter      :: dx=1.e-5_dp  ! finite difference increment
 ! check if less than saturation
 if(volFracLiq < theta_sat)then
  ! ***** compute analytical derivatives
  if(lTangent)then
   ! compute the effective saturation
   theta_e = (volFracLiq - theta_res) / (theta_sat - theta_res)
   ! get the first function and derivative
   y1 = (-1._dp/m)*theta_e**(-1._dp/m - 1._dp) / (theta_sat - theta_res)
   d1 = ( (m + 1._dp) / (m**2._dp * (theta_sat - theta_res)**2._dp) ) * theta_e**(-1._dp/m - 2._dp)
   ! get the second function and derivative
   xx = theta_e**(-1._dp/m) - 1._dp
   y2 = (1._dp/n)*xx**(1._dp/n - 1._dp)
   d2 = ( -(1._dp - n)/((theta_sat - theta_res)*m*n**2._dp) ) * xx**(1._dp/n - 2._dp) * theta_e**(-1._dp/m - 1._dp)
   ! return the derivative
   dPsi_dTheta2 = (d1*y2 + y1*d2)/alpha
  ! ***** compute numerical derivatives
  else
   func0 = dPsi_dTheta(volFracLiq,   alpha,theta_res,theta_sat,n,m)
   func1 = dPsi_dTheta(volFracLiq+dx,alpha,theta_res,theta_sat,n,m)
   dPsi_dTheta2 = (func1 - func0)/dx
  endif
 ! (case where volumetric liquid water content exceeds porosity)
 else
  dPsi_dTheta2 = 0._dp
 endif
 end function dPsi_dTheta2


 ! ***********************************************************************************************************
 ! new function: compute the derivative in hydraulic conductivity w.r.t. matric head (s-1)
 ! ***********************************************************************************************************
 function dHydCond_dPsi(psi,k_sat,alpha,n,m,iceImpede,lTangent)
 ! computes the derivative in hydraulic conductivity w.r.t matric head,
 !  given psi and soil hydraulic parameters k_sat, alpha, n, and m
 implicit none
 ! dummies
 real(dp),intent(in)     :: psi        ! soil water suction (m)
 real(dp),intent(in)     :: k_sat      ! saturated hydraulic conductivity (m s-1)
 real(dp),intent(in)     :: alpha      ! scaling parameter (m-1)
 real(dp),intent(in)     :: n          ! vGn "n" parameter (-)
 real(dp),intent(in)     :: m          ! vGn "m" parameter (-)
 real(dp),intent(in)     :: iceImpede  ! ice impedence factor (-)
 logical(lgt),intent(in) :: lTangent   ! method used to compute derivative (.true. = analytical)
 real(dp)                :: dHydCond_dPsi  ! derivative in hydraulic conductivity w.r.t. matric head (s-1)
 ! locals for analytical derivatives
 real(dp)                :: f_x1       ! f(x) for part of the numerator
 real(dp)                :: f_x2       ! f(x) for part of the numerator
 real(dp)                :: f_nm       ! f(x) for the numerator
 real(dp)                :: f_dm       ! f(x) for the denominator
 real(dp)                :: d_x1       ! df(x)/dpsi for part of the numerator
 real(dp)                :: d_x2       ! df(x)/dpsi for part of the numerator
 real(dp)                :: d_nm       ! df(x)/dpsi for the numerator
 real(dp)                :: d_dm       ! df(x)/dpsi for the denominator
 ! locals for numerical derivatives
 real(dp),parameter      :: dx=1.e-5_dp  ! finite difference increment
 real(dp)                :: hydCond0   ! hydraulic condictivity value for base case
 real(dp)                :: hydCond1   ! hydraulic condictivity value for perturbed case
 ! derivative is zero if saturated
 if(psi<0._dp)then
  ! ***** compute analytical derivatives
  if(lTangent)then
   ! compute the derivative for the numerator
   f_x1 = (psi*alpha)**(n - 1._dp)
   f_x2 = (1._dp + (psi*alpha)**n)**(-m)
   d_x1 = alpha * (n - 1._dp)*(psi*alpha)**(n - 2._dp)
   d_x2 = alpha * n*(psi*alpha)**(n - 1._dp) * (-m)*(1._dp + (psi*alpha)**n)**(-m - 1._dp)
   f_nm = (1._dp - f_x1*f_x2)**2._dp
   d_nm = (-d_x1*f_x2 - f_x1*d_x2) * 2._dp*(1._dp - f_x1*f_x2)
   ! compute the derivative for the denominator
   f_dm = (1._dp + (psi*alpha)**n)**(m/2._dp)
   d_dm = alpha * n*(psi*alpha)**(n - 1._dp) * (m/2._dp)*(1._dp + (psi*alpha)**n)**(m/2._dp - 1._dp)
   ! and combine
   dHydCond_dPsi = iceImpede*k_sat*(d_nm*f_dm - d_dm*f_nm) / (f_dm**2._dp)
  else
   ! ***** compute numerical derivatives
   hydcond0  = hydCond_psi(psi,   k_sat,alpha,n,m) * iceImpede
   hydcond1  = hydCond_psi(psi+dx,k_sat,alpha,n,m) * iceImpede
   dHydCond_dPsi = (hydcond1 - hydcond0)/dx
  endif
 else
  dHydCond_dPsi = 0._dp
 endif
 end function dHydCond_dPsi


 ! ***********************************************************************************************************
 ! new function: compute the derivative in hydraulic conductivity w.r.t. volumetric liquid water content (m s-1)
 ! ***********************************************************************************************************
 function dHydCond_dLiq(volFracLiq,k_sat,theta_res,theta_sat,m,iceImpede,lTangent)
 ! computes the derivative in hydraulic conductivity w.r.t the volumetric fraction of liquid water,
 !  given volFracLiq and soil hydraulic parameters k_sat, theta_sat, theta_res, and m
 implicit none
 ! dummies
 real(dp),intent(in)     :: volFracLiq ! volumetric fraction of liquid water (-)
 real(dp),intent(in)     :: k_sat      ! saturated hydraulic conductivity (m s-1)
 real(dp),intent(in)     :: theta_res  ! soil residual volumetric water content (-)
 real(dp),intent(in)     :: theta_sat  ! soil porosity (-)
 real(dp),intent(in)     :: m          ! vGn "m" parameter (-)
 real(dp),intent(in)     :: iceImpede  ! ice impedence factor (-)
 logical(lgt),intent(in) :: lTangent   ! method used to compute derivative (.true. = analytical)
 real(dp)                :: dHydCond_dLiq  ! derivative in hydraulic conductivity w.r.t. matric head (s-1)
 ! locals for analytical derivatives
 real(dp)                :: theta_e  ! effective soil moisture
 real(dp)                :: f1       ! f(x) for the first function
 real(dp)                :: d1       ! df(x)/dLiq for the first function
 real(dp)                :: x1,x2    ! f(x) for different parts of the second function
 real(dp)                :: p1,p2,p3 ! df(x)/dLiq for different parts of the second function
 real(dp)                :: f2       ! f(x) for the second function
 real(dp)                :: d2       ! df(x)/dLiq for the second function
 ! locals for numerical derivatives
 real(dp),parameter      :: dx=1.e-5_dp  ! finite difference increment
 real(dp)                :: hydCond0 ! hydraulic condictivity value for base case
 real(dp)                :: hydCond1 ! hydraulic condictivity value for perturbed case

 ! derivative is zero id super-saturated
 if(volFracLiq < theta_sat)then
  ! ***** compute analytical derivatives
  if(lTangent)then
   ! compute the effective saturation
   theta_e = (volFracLiq - theta_res) / (theta_sat - theta_res)
   ! compute the function and derivative of the first fuction
   f1 = k_sat*theta_e**0.5_dp
   d1 = k_sat*0.5_dp*theta_e**(-0.5_dp) / (theta_sat - theta_res)
   ! compute the function and derivative of the second function
   ! (first part)
   x1 = 1._dp - theta_e**(1._dp/m)
   p1 = (-1._dp/m)*theta_e**(1._dp/m - 1._dp) / (theta_sat - theta_res)   ! differentiate (1.d - theta_e**(1.d/m)
   ! (second part)
   x2 = x1**m
   p2 = m*x1**(m - 1._dp)
   ! (final)
   f2 = (1._dp - x2)**2._dp
   p3 = -2._dp*(1._dp - x2)
   ! (combine)
   d2 = p1*p2*p3
   ! pull it all together
   dHydCond_dLiq = (d1*f2 + d2*f1) * iceImpede
  else
   ! ***** compute numerical derivatives
   hydcond0 = hydCond_liq(volFracLiq,   k_sat,theta_res,theta_sat,m) * iceImpede
   hydcond1 = hydCond_liq(volFracLiq+dx,k_sat,theta_res,theta_sat,m) * iceImpede
   dHydCond_dLiq = (hydcond1 - hydcond0)/dx
  endif
 else
  dHydCond_dLiq = 0._dp
 endif
 end function dHydCond_dLiq


 ! ***********************************************************************************************************
 ! new function: compute Darcy's flux as a function of matric head
 ! ***********************************************************************************************************
 function darcyFlux_matric(ice_above,ice_below,psi_above,psi_below,z_above,z_below,theta_res,theta_sat,k_sat,alpha,n,m,f_impede)
 ! computes Darcy's flux (m s-1)
 implicit none
 ! dummy variables
 real(dp),intent(in) :: ice_above   ! volumetric ice content in the layer above (-)
 real(dp),intent(in) :: ice_below   ! volumetric ice content in the layer below (-)
 real(dp),intent(in) :: psi_above   ! matric head in the layer above the interface (m)
 real(dp),intent(in) :: psi_below   ! matric head in the layer below the interface (m)
 real(dp),intent(in) :: z_above     ! height in the layer above the interface (m)
 real(dp),intent(in) :: z_below     ! height in the layer below the interface (m)
 real(dp),intent(in) :: theta_res   ! soil residual volumetric water content (-)
 real(dp),intent(in) :: theta_sat   ! soil porosity (-)
 real(dp),intent(in) :: k_sat       ! saturated hydraulic conductivity (m s-1)
 real(dp),intent(in) :: alpha       ! scaling parameter (m-1)
 real(dp),intent(in) :: n           ! vGn "n" parameter (-)
 real(dp),intent(in) :: m           ! vGn "m" parameter (-)
 real(dp),intent(in) :: f_impede    ! ice impedence factor (-)
 real(dp)            :: darcyFlux_matric   ! the Darcy flux (m s-1)
 ! local variables
 real(dp)            :: imp_above   ! ice impedence factor in the layer above (-)
 real(dp)            :: imp_below   ! ice impedence factor in the layer below (-)
 real(dp)            :: mHC_above   ! hydraulic conductivity in the layer above (m s-1)
 real(dp)            :: mHC_below   ! hydraulic conductivity in the layer below (m s-1)
 real(dp)            :: iHC         ! hydraulic conductivity at layer interface  (m s-1)
 real(dp)            :: dPsi        ! spatial differences in matric head (m)
 real(dp)            :: dz          ! spatial differences in height (m)
 ! compute the ice impedence factor for the layer above and the layer below
 imp_above = iceImpede(ice_above,theta_res,theta_sat,f_impede)
 imp_below = iceImpede(ice_below,theta_res,theta_sat,f_impede)
 ! compute hydraulic conductivity at layer mid-points (m s-1)
 mHC_above = hydCond_psi(psi_above,k_sat,alpha,n,m) * imp_above
 mHC_below = hydCond_psi(psi_below,k_sat,alpha,n,m) * imp_below
 ! compute hydraulic conductivity at layer interface (m s-1)
 iHC   = (mHC_above*mHC_below)**0.5_dp
 ! compute the spatial differences in matric head and height
 dPsi  = psi_below - psi_above ! m
 dz    = z_below - z_above     ! m
 ! compute the flux
 darcyFlux_matric = -iHC*dPsi/dz + iHC
 end function darcyFlux_matric


 ! ***********************************************************************************************************
 ! new function: compute Darcy's flux as a function of volumetric liquid water content
 ! ***********************************************************************************************************
 function darcyFlux_liquid(ice_above,ice_below,liq_above,liq_below,z_above,z_below,k_sat,theta_res,theta_sat,alpha,n,m,f_impede)
 ! computes Darcy's flux (m s-1)
 implicit none
 ! dummy variables
 real(dp),intent(in) :: ice_above   ! volumetric ice content in the layer above (-)
 real(dp),intent(in) :: ice_below   ! volumetric ice content in the layer below (-)
 real(dp),intent(in) :: liq_above   ! matric head in the layer above the interface (m)
 real(dp),intent(in) :: liq_below   ! matric head in the layer beloe the interface (m)
 real(dp),intent(in) :: z_above     ! height in the layer above the interface (m)
 real(dp),intent(in) :: z_below     ! height in the layer below the interface (m)
 real(dp),intent(in) :: k_sat       ! saturated hydraulic conductivity (m s-1)
 real(dp),intent(in) :: theta_res   ! soil residual volumetric water content (-)
 real(dp),intent(in) :: theta_sat   ! soil porosity (-)
 real(dp),intent(in) :: alpha       ! scaling parameter (m-1)
 real(dp),intent(in) :: n           ! vGn "n" parameter (-)
 real(dp),intent(in) :: m           ! vGn "m" parameter (-)
 real(dp),intent(in) :: f_impede    ! ice impedence factor (-)
 real(dp)            :: darcyFlux_liquid   ! the Darcy flux (m s-1)
 ! local variables
 real(dp)            :: imp_above   ! ice impedence factor in the layer above (-)
 real(dp)            :: imp_below   ! ice impedence factor in the layer below (-)
 real(dp)            :: mHC_above   ! hydraulic conductivity in the layer above (m s-1)
 real(dp)            :: mHC_below   ! hydraulic conductivity in the layer below (m s-1)
 real(dp)            :: mD_above    ! hydraulic diffusivity in the layer above (m2 s-1)
 real(dp)            :: mD_below    ! hydraulic diffusivity in the layer below (m2 s-1)
 real(dp)            :: iHC         ! hydraulic conductivity at layer interface  (m s-1)
 real(dp)            :: iD          ! hydraulic diffusivity at layer interface  (m2 s-1)
 real(dp)            :: dLiq        ! spatial differences in volumetric liquid water content (-)
 real(dp)            :: dz          ! spatial differences in height (m)
 ! compute the ice impedence factor for the layer above and the layer below
 imp_above = iceImpede(ice_above,theta_res,theta_sat,f_impede)
 imp_below = iceImpede(ice_below,theta_res,theta_sat,f_impede)
 ! compute hydraulic conductivity at the layer mid-points (m s-1)
 mHC_above = hydCond_liq(liq_above,k_sat,theta_res,theta_sat,m) * imp_above
 mHC_below = hydCond_liq(liq_below,k_sat,theta_res,theta_sat,m) * imp_below
 ! compute hydraulic diffusivity at layer mid-points (m2 s-1)
 mD_above = dPsi_dTheta(liq_above,alpha,theta_res,theta_sat,n,m) * mHC_above
 mD_below = dPsi_dTheta(liq_below,alpha,theta_res,theta_sat,n,m) * mHC_below
 ! compute the hydraulic conductivity and diffusivity at the layer interface
 iHC = (mHC_above*mHC_below)**0.5_dp  ! (m s-1)
 iD  = (mD_above*mD_below)**0.5_dp    ! (m2 s-1)
 ! compute the spatial differences in volumetric liquid water content and height
 dLiq  = liq_below - liq_above ! (-)
 dz    = z_below - z_above     ! (m)
 ! compute the flux
 darcyFlux_liquid = -iD*dLiq/dz + iHC
 end function darcyFlux_liquid


 ! ***********************************************************************************************************
 ! new function: compute relative humidity of air in soil pore space
 ! ***********************************************************************************************************
 function RH_soilair(matpot,Tk)
 USE multiconst,only: gravity, &      ! acceleration of gravity       (m s-2)
                      R_wv            ! gas constant for water vapor  (J kg-1 K-1; [J = Pa m3])
 implicit none
 real(dp),intent(in) :: matpot        ! soil water suction -- matric potential (m)
 real(dp),intent(in) :: Tk            ! temperature (K)
 real(dp)            :: RH_soilair    ! relative humidity of air in soil pore space
 ! compute relative humidity (UNITS NOTE: Pa = kg m-1 s-2, so R_wv units = m2 s-2 K-1)
 RH_soilair = exp( (gravity*matpot) / (R_wv*Tk) )
 end function RH_soilair

 ! ***********************************************************************************************************
 ! new function: compute the critical temperature above which all water is unfrozen
 ! ***********************************************************************************************************
 function crit_soilT(theta,theta_res,theta_sat,alpha,n,m)
 USE multiconst,only: gravity,   &    ! acceleration of gravity    (m s-2)
                      Tfreeze,   &    ! temperature at freezing    (K)
                      LH_fus,    &    ! latent heat of fusion      (J kg-1, or m2 s-2)
                      iden_ice,  &    ! intrinsic density of ice   (kg m-3)
                      iden_water      ! intrinsic density of water (kg m-3)
 implicit none
 ! dummy variables
 real(dp),intent(in) :: theta         ! total soil water content, frozen plus unfrozen (-)
 real(dp),intent(in) :: theta_res     ! residual liquid water content (-)
 real(dp),intent(in) :: theta_sat     ! porosity (-)
 real(dp),intent(in) :: alpha         ! vGn scaling parameter (m-1)
 real(dp),intent(in) :: n             ! vGn "n" parameter (-)
 real(dp),intent(in) :: m             ! vGn "m" parameter (-)
 real(dp)            :: crit_soilT    ! critical soil temperature (K)
 ! local variables
 real(dp)            :: relsat        ! relative saturation (-)
 real(dp)            :: kappa         ! constant (m K-1)
 ! compute kappa (m K-1)
 kappa  = (iden_ice/iden_water)*(LH_fus/(gravity*Tfreeze))  ! NOTE: J = kg m2 s-2
 ! compute relative saturation (-)
 relsat = (min(theta,theta_sat) - theta_res)/(theta_sat - theta_res)
 ! compute the critical temperature above which all water is unfrozen (K)
 !print*,'in soil_utils',Tfreeze,relsat,m,n,alpha,kappa
 crit_soilT = Tfreeze + ((relsat**(-1._dp/m) - 1._dp)**(1._dp/n))/(alpha*kappa)
 end function crit_soilT

 ! ***********************************************************************************************************
 ! new function: differentiate the freezing curve w.r.t. temperature
 ! *********************************************************************************************************** 
 function dTheta_dTk(Tk,theta_res,theta_sat,alpha,n,m)
 USE multiconst,only: gravity,   &    ! acceleration of gravity    (m s-2)
                      Tfreeze,   &    ! temperature at freezing    (K)
                      LH_fus,    &    ! latent heat of fusion      (J kg-1, or m2 s-2)
                      iden_ice,  &    ! intrinsic density of ice   (kg m-3)
                      iden_water      ! intrinsic density of water (kg m-3)
 implicit none
 real(dp),intent(in) :: Tk            ! temperature (K)
 real(dp),intent(in) :: theta_res     ! residual liquid water content (-)
 real(dp),intent(in) :: theta_sat     ! porosity (-)
 real(dp),intent(in) :: alpha         ! vGn scaling parameter (m-1)
 real(dp),intent(in) :: n             ! vGn "n" parameter (-)
 real(dp),intent(in) :: m             ! vGn "m" parameter (-)
 real(dp)            :: dTheta_dTk    ! derivative of the freezing curve w.r.t. temperature (K-1)
 ! local variables
 real(dp)            :: kappa         ! constant (m K-1)
 real(dp)            :: xtemp         ! alpha*kappa*(Tk-Tfreeze) -- dimensionless variable (used more than once)
 ! compute kappa (m K-1)
 kappa = (iden_ice/iden_water)*(LH_fus/(gravity*Tfreeze))  ! NOTE: J = kg m2 s-2
 ! define a tempory variable that is used more than once (-)
 xtemp = alpha*kappa*(Tk-Tfreeze)
 ! differentiate the freezing curve w.r.t. temperature -- making use of the chain rule
 dTheta_dTk = (alpha*kappa) * n*xtemp**(n - 1._dp) * (-m)*(1._dp + xtemp**n)**(-m - 1._dp) * (theta_sat - theta_res)
 end function dTheta_dTk






end module soil_utils_module

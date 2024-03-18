module hyp_2F1_module

! data types
USE nrtype

! privacy
implicit none

private::INF_NORM
private::TANZ
private::LOG1P
private::EXPM1
private::GAMMA_INV
private::GAMMA_RATIO_DIFF_SMALL_EPS
private::GAMMA_INV_DIFF_EPS
private::A_SUM_INIT
private::LOG_A_SUM_INIT
private::B_SUM_INIT_PS_ONE
private::B_SUM_INIT_PS_INFINITY
private::CV_POLY_DER_CALC
private::MIN_N_CALC
private::HYP_PS_ZERO
private::HYP_PS_ONE
private::HYP_PS_INFINITY
private::HYP_PS_COMPLEX_PLANE_REST
private::TEST_2F1
public::HYP_2F1

! constant parameters
real(rkind),parameter :: EPS15=1.e-15_rkind
real(rkind),parameter :: ZERO=0._rkind,ONE=1._rkind,TWO=2._rkind,HALF=0.5_rkind
real(rkind),parameter :: M_PI=3.14159265358979323846_rkind
real(rkind),parameter :: M_PI_2=1.57079632679489661923_rkind
real(rkind),parameter :: M_1_PI=0.31830988618379067154_rkind

contains

 !============== START HYP_2F1 FILE ====================================
 !
 ! Gamma_inv denotes the entire inverse of the Gamma function.
 ! F(z) means 2F1(a,b,c,z) with the a, b, c and z given as inputs 
 ! in the routine.
 !
 ! Elementary functions and standard constants 
 ! are defined in the module.
 ! See N.J.~Higham, ``Accuracy and Stability of Numerical Algorithms'',
 ! SIAM, Philadelphia, 1996 for expm1 implementation.
 ! log1p follows instantly.
 !
 ! 19/04/2012    Modifications by Daniel Sabanes Bove:
 !               - renamed LOG_GAMMA to LOG_GAMMA_FUN to avoid name
 !                 clash with intrinsic function
 ! 11/01/2024    Modifications by Ashley Van Beusekom:
 !               - made one module and removed functions as variable declarations
 !               - made precision rkind dependent
 !               - lowercase for readability
 !----------------------------------------------------------------------
 !
 function INF_NORM(Z)
   complex(rkind),intent(in) :: Z
   real(rkind)  :: INF_NORM
   INF_NORM=MAX(ABS(REAL(Z,rkind)),ABS(AIMAG(Z)))
   return
 end function INF_NORM
 !
 function TANZ(Z)
   complex(rkind),intent(in) :: Z
   complex(rkind) :: TANZ
   TANZ=SIN(Z)/COS(Z)
   return
 end function TANZ
 !
 function LOG1P(Z)
   complex(rkind),intent(in) :: Z
   real(rkind) :: X,XP1,LOG1P_X
   real(rkind) :: Y,YX,YX2,YX2P1,LOG1P_YX2
   real(rkind) :: RE_LOG1P,IM_LOG1P
   complex(rkind) :: LOG1P
   if(INF_NORM(Z).lt.ONE) then
      X = REAL(Z,rkind); XP1 = X+ONE
      if(XP1.eq.ONE) then
         LOG1P_X = X
      else
         LOG1P_X = LOG(XP1)*X/(XP1-ONE)
      endif
      Y = AIMAG(Z)
      YX = Y/XP1; YX2 = YX*YX; YX2P1 = YX2+ONE
      if(YX2P1.eq.ONE) then
         LOG1P_YX2 = YX2
      else
         LOG1P_YX2 = LOG(YX2P1)*YX2/(YX2P1-ONE)
      endif
      RE_LOG1P = LOG1P_X + HALF*LOG1P_YX2
      IM_LOG1P = ATAN2(Y,XP1)
      LOG1P = CMPLX(RE_LOG1P,IM_LOG1P,rkind)
      return
   else
      LOG1P=LOG(ONE+Z)
      return
   endif
 end function LOG1P
 !
 function EXPM1(Z)
   complex(rkind),intent(in) :: Z
   real(rkind) :: X,EXPM1_X,EXP_X,Y,SIN_HALF_Y
   real(rkind) :: RE_EXPM1,IM_EXPM1
   complex(rkind) :: EXPM1
   if(INF_NORM(Z).lt.ONE) then
      X = real(Z,rkind); EXP_X = EXP(X)
      Y = AIMAG(Z); SIN_HALF_Y=SIN(HALF*Y)
      if(EXP_X.eq.ONE) then
         EXPM1_X = X
      else 
         EXPM1_X = (EXP_X-ONE)*X/LOG(EXP_X)
      endif
      RE_EXPM1 = EXPM1_X-TWO*EXP_X*SIN_HALF_Y*SIN_HALF_Y 
      IM_EXPM1 = EXP_X*SIN(Y)
      EXPM1 = CMPLX(RE_EXPM1,IM_EXPM1,rkind)
      return
   else
      EXPM1=EXP(Z)-ONE
      return
   endif
 end function EXPM1
 !
 !----------------------------------------------------------------------
 recursive function LOG_GAMMA_FUN(Z) result(RES)
   !----------------------------------------------------------------------
   ! Logarithm of Gamma[z] and Gamma inverse function
   ! ------------------------------------------------
   !
   ! For log[Gamma[z]],if z is not finite 
   ! or is a negative integer, the program 
   ! returns an error message and stops.
   ! The Lanczos method is used. Precision : ~ 1E-15
   ! The method works for Re[z]>0.5 .
   ! If Re[z]<=0.5, one uses the formula Gamma[z].Gamma[1-z]=Pi/sin(Pi.z)
   ! log[sin(Pi.z)] is calculated with the Kolbig method 
   ! (K.S. Kolbig, Comp. Phys. Comm., Vol. 4, p.221(1972)): 
   ! If z=x+iy and y>=0, log[sin(Pi.z)]=log[sin(Pi.eps)]-i.Pi.n, 
   ! with z=n+eps so 0<=Re[eps]< 1 and n integer.
   ! If y>110, log[sin(Pi.z)]=-i.Pi.z+log[0.5]+i.Pi/2 
   ! numerically so that no overflow can occur.
   ! If z=x+iy and y< 0, log[Gamma(z)]=[log[Gamma(z*)]]*, 
   ! so that one can use the previous formula with z*.
   !
   ! For Gamma inverse, Lanczos method is also used 
   ! with Euler reflection formula.
   ! sin (Pi.z) is calculated as sin (Pi.(z-n)) 
   ! to avoid inaccuracy with z = n + eps 
   ! with n integer and |eps| as small as possible.
   !
   !
   ! Variables:
   ! ----------
   ! x,y: Re[z], Im[z]
   ! log_sqrt_2Pi,log_Pi : log[sqrt(2.Pi)], log(Pi).
   ! sum : Rational function in the Lanczos method
   ! log_Gamma_z : log[Gamma(z)] value.
   ! c : table containing the fifteen coefficients in the expansion 
   ! used in the Lanczos method.
   ! eps,n : z=n+eps so 0<=Re[eps]< 1 and n integer for Log[Gamma].
   !         z=n+eps and n integer 
   !         so |eps| is as small as possible for Gamma_inv.
   ! log_const : log[0.5]+i.Pi/2
   ! g : coefficient used in the Lanczos formula. It is here 607/128.
   ! z,z_m_0p5,z_p_g_m0p5,zm1 : argument of the Gamma function, 
   ! z-0.5, z-0.5+g, z-1 
   ! res: returned value
   !----------------------------------------------------------------------
   implicit none
   complex(rkind),intent(in) :: Z
   integer(i4b)              :: N,I
   real(rkind)               :: X,Y,LOG_SQRT_2PI,G,LOG_PI,M_LN2,C(0:14)
   complex(rkind)            :: GAMMA_SUM,Z_M_0P5,Z_P_G_M0P5,ZM1
   complex(rkind)            :: LOG_CONST,I_PI,EPS,LOG_SIN_PI_Z,RES
   !
   M_LN2=0.69314718055994530942_rkind; X=REAL(Z,rkind); Y=AIMAG(Z)
   if((Z.eq.NINT(X)).and.(X.le.ZERO)) &
        print*,'Z IS NEGATIVE integer IN LOG_GAMMA_FUN'
   if(X.ge.HALF) then
      LOG_SQRT_2PI=0.91893853320467274177_rkind; G=4.7421875_rkind
      Z_M_0P5=Z-HALF; Z_P_G_M0P5=Z_M_0P5+G; ZM1=Z-ONE
      C=(/   0.99999999999999709182_rkind   , 57.156235665862923517_rkind,       &
           -59.597960355475491248_rkind     , 14.136097974741747174_rkind,       &
            -0.49191381609762019978_rkind   ,  0.33994649984811888699e-4_rkind,  &
             0.46523628927048575665e-4_rkind, -0.98374475304879564677e-4_rkind,  &
             0.15808870322491248884e-3_rkind, -0.21026444172410488319e-3_rkind,  &
             0.21743961811521264320e-3_rkind, -0.16431810653676389022e-3_rkind,  &
             0.84418223983852743293e-4_rkind, -0.26190838401581408670e-4_rkind,  &
             0.36899182659531622704e-5_rkind /)
 
      GAMMA_SUM=C(0)
      do I=1,14
         GAMMA_SUM=GAMMA_SUM+C(I)/(ZM1+I)
      enddo
      RES=LOG_SQRT_2PI+LOG(GAMMA_SUM)+Z_M_0P5*LOG(Z_P_G_M0P5) &
           -Z_P_G_M0P5
      return
   else if(Y.ge.ZERO) then
      if(X.lt.NINT(X)) then
         N=NINT(X)-1
      else
         N=NINT(X)
      endif
      LOG_PI=1.1447298858494002_rkind
      LOG_CONST=CMPLX(-M_LN2,M_PI_2,rkind); I_PI=CMPLX(ZERO,M_PI,rkind)
      EPS=Z-N
      if(Y.gt.110._rkind) then
         LOG_SIN_PI_Z=-I_PI*Z+LOG_CONST
      else
         LOG_SIN_PI_Z=LOG(SIN(M_PI*EPS))-I_PI*N
      endif
      RES=LOG_PI-LOG_SIN_PI_Z-LOG_GAMMA_FUN(ONE-Z);
      return
   else
      RES=CONJG(LOG_GAMMA_FUN(CONJG(Z)))
      return
   endif
 end function LOG_GAMMA_FUN
 !
 !----------------------------------------------------------------------
 ! Inverse of the Gamma function [1/Gamma](z)
 ! ------------------------------------------
 ! It is calculated with the Lanczos method for Re[z] >= 0.5 
 ! and is precise up to 10^{-15}.
 ! If Re[z] <= 0.5, one uses the formula 
 ! Gamma[z].Gamma[1-z] = Pi/sin (Pi.z).
 ! sin (Pi.z) is calculated as sin (Pi.(z-n)) to avoid inaccuracy,
 ! with z = n + eps with n integer and |eps| as small as possible.
 ! 
 ! Variables 
 ! ---------
 ! z : argument of the function
 ! x: Re[z]
 ! eps,n : z = n + eps with n integer and |eps| as small as possible.
 ! res: returned value
 !----------------------------------------------------------------------
 recursive function GAMMA_INV(Z) result(RES)
   !--------------------------------------------------------------------
   implicit none
   complex(rkind),intent(in) :: Z
   integer(i4b)              :: N,I
   real(rkind)               :: X,LOG_SQRT_2PI,G,C(0:14)
   complex(rkind)            :: RES,GAMMA_SUM,Z_M_0P5,Z_P_G_M0P5,ZM1,EPS
   !
   X=REAL(Z,rkind)
   if(X.ge.HALF) then
      LOG_SQRT_2PI=0.91893853320467274177_rkind; G=4.7421875_rkind
      Z_M_0P5=Z-HALF; Z_P_G_M0P5=Z_M_0P5+G; ZM1=Z-ONE
      C=(/   0.99999999999999709182_rkind   , 57.156235665862923517_rkind,       &
           -59.597960355475491248_rkind     , 14.136097974741747174_rkind,       &
            -0.49191381609762019978_rkind   ,  0.33994649984811888699e-4_rkind,  &
             0.46523628927048575665e-4_rkind, -0.98374475304879564677e-4_rkind,  &
             0.15808870322491248884e-3_rkind, -0.21026444172410488319e-3_rkind,  &
             0.21743961811521264320e-3_rkind, -0.16431810653676389022e-3_rkind,  &
             0.84418223983852743293e-4_rkind, -0.26190838401581408670e-4_rkind,  &
             0.36899182659531622704e-5_rkind /)
 
      GAMMA_SUM=C(0)
      do I=1,14
         GAMMA_SUM=GAMMA_SUM+C(I)/(ZM1+I);
      enddo
      RES=EXP(Z_P_G_M0P5-Z_M_0P5*LOG(Z_P_G_M0P5)-LOG_SQRT_2PI) &
           /GAMMA_SUM
      return
   else
      X=REAL(Z,rkind); N=NINT(X)
      EPS=Z-N
      if(MOD(N,2).eq.0) then
         RES=SIN(M_PI*EPS)*M_1_PI/GAMMA_INV(ONE-Z)
         return
      else
         RES=-SIN(M_PI*EPS)*M_1_PI/GAMMA_INV(ONE-Z)
         return
      endif
   endif
 end function GAMMA_INV
 !----------------------------------------------------------------------
 !
 ! Calculation of H(z,eps) = [Gamma(z+eps)/Gamma(z) - 1]/eps, with e and
 ! ---------------------------------------------------------------------
 ! z complex so z,z+eps are not negative integers and 0 <= |eps|oo < 0.1
 ! ---------------------------------------------------------------------
 ! The function H(z,eps) = [Gamma(z+eps)/Gamma(z) - 1]/e is calculated 
 ! here with the Lanczos method.
 ! For the Lanczos method, the gamma parameter, denoted as g, 
 ! is 4.7421875 and one uses a sum of 15 numbers with the table c[15], 
 ! so that it is precise up to machine accuracy.
 ! The H(z,eps) function is used in formulas occuring in1-z and 1/z 
 ! transformations (see Comp. Phys. Comm. paper).
 !
 ! One must have z and z+eps not negative integers as otherwise 
 ! it is clearly not defined.
 ! As this function is meant to be precise for small |eps|oo, 
 ! one has to have 0 <= |eps|oo < 0.1 .
 ! Indeed, a direct implementation of H(z,eps) with Gamma_inv or 
 ! log_Gamma for |eps|oo >= 0.1 is numerically stable.
 ! The returned function has full numerical accuracy 
 ! even if |eps|oo is very small.
 !
 ! eps not equal to zero
 ! ---------------------
 ! If Re(z) >= 0.5 or Re(z+eps) >= 0.5, one clearly has Re(z) > 0.4 
 ! and Re(z+eps) > 0.4, 
 ! so that the Lanczos summation can be used for both Gamma(z) 
 ! and Gamma(z+eps).
 ! One then has:
 ! log[Gamma(z+eps)/Gamma(z)] = 
 ! (z-0.5) log1p[eps/(z+g-0.5)] + eps log(z+g-0.5+eps) - eps 
 ! + log1p[-eps \sum_{i=1}^{14} c[i]/((z-1+i)(z-1+i+eps)) 
 ! / (c[0] + \sum_{i=1}^{14} c[i]/(z-1+i))]
 ! H(z,eps) = expm1[log[Gamma(z+eps)/Gamma(z)]]/eps .
 !
 ! If Re(z) < 0.5 and Re(z+eps) < 0.5, 
 ! Euler reflection formula is used for both Gamma(z) and Gamma(z+eps).
 ! One then has: 
 ! H(z+eps,-eps) = [cos(pi.eps) + sin(pi.eps)/tan(pi(z-n))].H(1-z,-eps) 
 ! + (2/eps).sin^2(eps.pi/2) - sin(pi.eps)/(eps.tan(pi.(z-n)))
 ! H(1-z,-eps) is calculated with the Lanczos summation 
 ! as Re(1-z) >= 0.5 and Re(1-z-eps) >= 0.5 .
 ! z-n is used in tan(pi.z) instead of z to avoid inaccuracies 
 ! due the finite number of digits of pi.
 ! H(z,eps) = H(z+eps,-eps)/(1 - eps.H(z+eps,-eps)) 
 ! provides the final result.
 !
 ! eps equal to zero
 ! -----------------
 ! It is obtained with the previous case and eps -> 0 :
 ! If Re(z) >= 0.5, one has:
 ! H(z,eps) = (z-0.5)/(z+g-0.5) + log(z+g-0.5) - 1 -
 ! \sum_{i=1}^{14} c[i]/((z-1+i)^2)/(c[0]+\sum_{i=1}^{14} c[i]/(z-1+i))
 !
 ! If Re(z) < 0.5, one has:
 ! H(z,0) = H(1-z,0) - pi/tan(pi.(z-n))
 !
 ! Variables
 ! ---------
 ! z,eps: input variables of the function H(z,eps)
 ! g,c[15]: double and table of 15 doubles defining the Lanczos sum 
 ! so that it provides the Gamma function 
 ! precise up to machine accuracy.
 ! eps_pz,z_m_0p5,z_pg_m0p5,eps_pz_pg_m0p5,zm1,zm1_p_eps: 
 ! z+eps,z-0.5,z+g-0.5,z+eps+g-0.5,z-1,z-1+eps
 ! x,eps_px: real parts of z and z+eps.
 ! n,m: closest integer ot the real part of z, same for z+eps.
 ! sum_num,sum_den: \sum_{i=1}^{14} c[i]/((z-1+i)(z-1+i+eps)) 
 ! and (c[0] + \sum_{i=1}^{14} c[i]/(z-1+i)). 
 ! They appear respectively as numerator and denominator in formulas.
 ! Pi_eps,term,T1_eps_z: pi.eps, sin (pi.eps)/tan(pi.(z-n)), 
 ! [cos(pi.eps) + sin(pi.eps)/tan(pi(z-n))].H(1-z,-eps)
 ! sin_Pi_2_eps,T2_eps_z,T_eps_z: sin^2(eps.pi/2), 
 ! (2/eps).sin^2(eps.pi/2) - sin(pi.eps)/(eps.tan(pi.(z-n))), 
 ! H(z+eps,-eps)
 ! res: returned value
 !----------------------------------------------------------------------
 recursive function GAMMA_RATIO_DIFF_SMALL_EPS(Z,EPS) result(RES)
   !--------------------------------------------------------------------
   implicit none
   complex(rkind),intent(in) :: Z,EPS
   integer(i4b)              :: N,M,I
   real(rkind)               :: G,X,EPS_PX,C(0:14)
   complex(rkind)            :: RES,SUM_NUM,SUM_DEN
   complex(rkind)            :: EPS_PZ,Z_M_0P5,Z_PG_M0P5,EPS_PZ_PG_M0P5,ZM1
   complex(rkind)            :: CI_ZM1_PI_INV,PI_EPS,TT,T1_EPS_Z,SIN_PI_2_EPS
   complex(rkind)            :: ZM1_P_EPS,T2_EPS_Z,T_EPS_Z
   !
   G=4.74218750_rkind
   if(INF_NORM(EPS).gt.0.1_rkind) &
        print*,'ONE MUST HAVE |EPS|< 0.1 IN GAMMA_RATIO_DIFF_SMALL_EPS'
   EPS_PZ=Z+EPS; Z_M_0P5=Z-HALF; Z_PG_M0P5=Z_M_0P5+G
   EPS_PZ_PG_M0P5=Z_PG_M0P5+EPS; ZM1=Z-ONE; ZM1_P_EPS=ZM1+EPS
   X=REAL(Z,rkind); EPS_PX=REAL(EPS_PZ,rkind); N=NINT(X); M=NINT(EPS_PX)
   if((Z.eq.N).and.(N.le.0)) then
      print*,'Z IS NEGATIVE integer IN GAMMA_RATIO_DIFF_SMALL_EPS'
   endif
   if((EPS_PZ.eq.M).and.(M.le.0)) then
      print*,'Z+EPS IS NEGATIVE integer IN GAMMA_RATIO_DIFF_SMALL_EPS'
   endif
   C=(/   0.99999999999999709182_rkind   , 57.156235665862923517_rkind,       &
        -59.597960355475491248_rkind     , 14.136097974741747174_rkind,       &
         -0.49191381609762019978_rkind   ,  0.33994649984811888699e-4_rkind,  &
          0.46523628927048575665e-4_rkind, -0.98374475304879564677e-4_rkind,  &
          0.15808870322491248884e-3_rkind, -0.21026444172410488319e-3_rkind,  &
          0.21743961811521264320e-3_rkind, -0.16431810653676389022e-3_rkind,  &
          0.84418223983852743293e-4_rkind, -0.26190838401581408670e-4_rkind,  &
          0.36899182659531622704e-5_rkind /)
   if((X.ge.HALF).or.(EPS_PX.ge.HALF)) then
      SUM_NUM=ZERO;SUM_DEN=C(0)
      do I=1,14
         CI_ZM1_PI_INV=C(I)/(ZM1+I)
         SUM_NUM=SUM_NUM+CI_ZM1_PI_INV/(ZM1_P_EPS+I)
         SUM_DEN=SUM_DEN+CI_ZM1_PI_INV
      enddo
      if(EPS.ne.ZERO) then
         RES=EXPM1(Z_M_0P5*LOG1P(EPS/Z_PG_M0P5) &
              +EPS*LOG(EPS_PZ_PG_M0P5)-EPS+LOG1P(-EPS*SUM_NUM/SUM_DEN))&
              /EPS
         return
      else
         RES=Z_M_0P5/Z_PG_M0P5 &
              +LOG(EPS_PZ_PG_M0P5)-ONE-SUM_NUM/SUM_DEN
         return
      endif
   else
      if(EPS.ne.ZERO) then
         PI_EPS=M_PI*EPS
         TT=SIN(PI_EPS)/TANZ(M_PI*(Z-N))
         T1_EPS_Z=(COS(PI_EPS)+TT)*& 
              GAMMA_RATIO_DIFF_SMALL_EPS(ONE-Z,-EPS)
         SIN_PI_2_EPS=SIN(M_PI_2*EPS)
         T2_EPS_Z=(TWO*SIN_PI_2_EPS*SIN_PI_2_EPS-TT)/EPS
         T_EPS_Z=T1_EPS_Z+T2_EPS_Z
         RES=(T_EPS_Z/(ONE-EPS*T_EPS_Z))
         return
      else
         RES=GAMMA_RATIO_DIFF_SMALL_EPS(ONE-Z,-EPS) &
              -M_PI/TANZ(M_PI*(Z-N))
         return
      endif
   endif
 end function GAMMA_RATIO_DIFF_SMALL_EPS
 !
 !----------------------------------------------------------------------
 ! Calculation of G(z,eps) = [Gamma_inv(z) - Gamma_inv(z+eps)]/eps 
 ! ---------------------------------------------------------------
 ! with e and z complex
 !---------------------
 ! The G(z,eps) function is used in formulas occuring in 1-z 
 ! and 1/z transformations (see Comp. Phys. Comm. paper).
 ! Several case have to be considered for its evaluation. 
 ! eps is considered equal to zero 
 ! if z+eps and z are equal numerically.
 !
 ! |eps|oo > 0.1
 ! -------------
 ! A direct evaluation with the values Gamma_inv(z) 
 ! and Gamma_inv(z+eps) is stable and returned.
 !
 ! |eps|oo <= 0.1 with z+eps and z numerically different
 ! -----------------------------------------------------
 ! If z is a negative integer, z+eps is not, 
 ! so that G(z,eps) = -Gamma_inv(z+eps)/eps, 
 ! for which a direct evaluation is precise and returned.
 ! If z+eps is a negative integer, z is not, 
 ! so that G(z,eps) = Gamma_inv(z)/eps, 
 ! for which a direct evaluation is precise and returned.
 ! If both of them are not negative integers, 
 ! one looks for the one of z and z+eps 
 ! which is the closest to a negative integer.
 ! If it is z, one returns H(z,eps).Gamma_inv(z+eps). 
 ! If it is z+eps, one returns H(z+eps,-eps).Gamma_inv(z).
 ! Both values are equal, so that one chooses the one 
 ! which makes the Gamma ratio Gamma(z+eps)/Gamma(z) 
 ! in H(z,eps) the smallest in modulus.
 !
 ! z+eps and z numerically equal
 ! -----------------------------
 ! If z is negative integer, G(z,0) = (-1)^(n+1) n!, 
 ! where z = -n, n integer, which is returned.
 ! If z is not negative integer, one returns H(z,eps).Gamma_inv(z+eps)
 !
 ! Variables
 ! ---------
 ! z,eps: input variables of the function G(z,eps)
 ! eps_pz,x,eps_px: z+eps,real parts of z and z+eps.
 ! n,m: closest integer ot the real part of z, same for z+eps.
 ! fact,k: (-1)^(n+1) n!, returned when z = -n, n integer 
 ! and z and z+eps identical numerically (eps ~ 0). 
 ! It is calculated with integer index k.
 ! is_z_negative_integer,is_eps_pz_negative_integer: 
 ! true if z is a negative integer, false if not, same for z+eps.
 ! z_neg_int_distance, eps_pz_neg_int_distance: 
 ! |z + |n||oo, |z + eps + |m||oo. 
 ! If |z + |n||oo < |z + eps + |m||oo, 
 ! z is closer to the set of negative integers than z+eps.
 ! Gamma_inv(z+eps) is then of moderate modulus 
 ! if Gamma_inv(z) is very small. 
 ! If z ~ n, H(z,eps) ~ -1/eps, 
 ! that so returning 
 ! G(z,eps) = H(z,eps).Gamma_inv(z+eps) here is preferred.
 ! Same for |z + |n||oo > |z + eps + |m||oo with z <-> z+eps.
 !
 !----------------------------------------------------------------------
 function GAMMA_INV_DIFF_EPS(Z,EPS)
   !--------------------------------------------------------------------
   implicit none
   complex(rkind),intent(in) :: Z,EPS
   integer(i4b)              :: M,N,K
   real(rkind)               :: X,EPS_PX,FACT
   real(rkind)               :: Z_NEG_INT_DISTANCE
   real(rkind)               :: EPS_PZ_NEG_INT_DISTANCE
   complex(rkind)            :: GAMMA_INV_DIFF_EPS,EPS_PZ
   logical(lgt)              :: IS_Z_NEG_INT,IS_EPS_PZ_NEG_INT
 
   EPS_PZ=Z+EPS; X=REAL(Z,rkind); EPS_PX=REAL(EPS_PZ,rkind)
   N=NINT(X); M=NINT(EPS_PX)
   IS_Z_NEG_INT=(Z.eq.N).and.(N.le.0)
   IS_EPS_PZ_NEG_INT=(EPS_PZ.eq.M).and.(M.le.0)
   if(INF_NORM(EPS).gt.0.10_rkind) then
      GAMMA_INV_DIFF_EPS = (GAMMA_INV(Z) - GAMMA_INV(EPS_PZ))/EPS
      return
   else if(EPS_PZ.ne.Z) then 
      if(IS_Z_NEG_INT) then
         GAMMA_INV_DIFF_EPS = (-GAMMA_INV(EPS_PZ)/EPS)
         return
      else if(IS_EPS_PZ_NEG_INT) then
         GAMMA_INV_DIFF_EPS = (GAMMA_INV(Z)/EPS)
         return
      else
         Z_NEG_INT_DISTANCE = INF_NORM (Z + ABS (N))
         EPS_PZ_NEG_INT_DISTANCE = INF_NORM (EPS_PZ + ABS (M))
         if(Z_NEG_INT_DISTANCE.lt.EPS_PZ_NEG_INT_DISTANCE) then
            GAMMA_INV_DIFF_EPS= &
                 GAMMA_RATIO_DIFF_SMALL_EPS(Z,EPS)*GAMMA_INV(EPS_PZ)
            return
         else
            GAMMA_INV_DIFF_EPS= &
                 GAMMA_RATIO_DIFF_SMALL_EPS(EPS_PZ,-EPS)*GAMMA_INV(Z)
            return
         endif
      endif
   else if(IS_Z_NEG_INT.and.IS_EPS_PZ_NEG_INT) then
      FACT = -ONE;K=-1
      do while (K.ge.N) 
         FACT=FACT*K
         K=K-1 
      enddo
      GAMMA_INV_DIFF_EPS = FACT
      return
   else
      GAMMA_INV_DIFF_EPS = &
           GAMMA_RATIO_DIFF_SMALL_EPS(Z,EPS)*GAMMA_INV(EPS_PZ)
      return
   endif
 end function GAMMA_INV_DIFF_EPS
 !----------------------------------------------------------------------
 !
 ! Calculation of Gamma_inv(1-m-eps)/eps of the A(z) polynomial in 1-z
 ! -------------------------------------------------------------------
 ! and 1/z transformations
 ! -----------------------
 ! This value occurs in A(z) in 1-z and 1/z transformations 
 ! (see Comp. Phys. Comm. paper) for m > 0.
 ! Both cases of 1-m-eps numerically negative integer 
 ! or not have to be considered
 ! 
 ! 1-eps-m and 1-m numerically different
 ! -------------------------------------
 ! One returns Gamma_inv(1-m-eps)/eps directly 
 ! as its value is accurate.
 ! To calculate Gamma_inv(1-m-eps), 
 ! one uses the value Gamma_inv(1-eps), 
 ! needed in considered transformations,
 ! and one uses the equality 
 ! Gamma_inv(1-m-eps) = Gamma_inv(1-eps) \prod_{i=1}^{m} (1-eps-i) 
 ! for m > 0.
 ! It is trivially demonstrated 
 ! from the equality Gamma(x+1) = x.Gamma(x). 
 ! One Gamma function evaluation is removed this way 
 ! from the calculation.
 ! 
 ! 1-eps-m and 1-m numerically equal
 ! ---------------------------------
 ! This implies that 1-m-eps is negative integer numerically.
 ! Here, eps~0, so that one returns the limit of Gamma_inv(1-m-eps)/eps
 ! for eps -> 0, which is (-1)^m (m-1)!
 !
 ! Variables
 ! ---------
 ! m,eps: variable inputs of the function 
 ! (m,eps) -> Gamma_inv(1-m-eps)/eps
 ! Gamma_inv_one_meps: Gamma_inv(1-eps), 
 ! previously calculated and here recycled 
 ! to quickly calculate Gamma_inv(1-m-eps).
 ! one_meps: 1-eps
 !----------------------------------------------------------------------
 function A_SUM_INIT(M,EPS,GAMMA_INV_ONE_MEPS)
   !--------------------------------------------------------------------
   implicit none
   integer(i4b),intent(in)   :: M
   complex(rkind),intent(in) :: EPS,GAMMA_INV_ONE_MEPS
   integer(i4b)              :: N,I
   real(rkind)               :: FACT
   complex(rkind)            :: A_SUM_INIT,ONE_MEPS
   complex(rkind)            :: GAMMA_INV_ONE_MEPS_MM
   !
   ONE_MEPS = ONE - EPS
   if(ONE_MEPS-M.ne.1-M) then
      GAMMA_INV_ONE_MEPS_MM = GAMMA_INV_ONE_MEPS
      do I=1,M
         GAMMA_INV_ONE_MEPS_MM = GAMMA_INV_ONE_MEPS_MM*(ONE_MEPS-I)
      enddo
      A_SUM_INIT=GAMMA_INV_ONE_MEPS_MM/EPS
      return
   else
      FACT=ONE
      do N=2,M-1
         FACT=FACT*N
      enddo
      if(MOD(M,2).eq.0) then
         A_SUM_INIT=FACT
      else
         A_SUM_INIT=-FACT
      endif
      return
   endif
 end function A_SUM_INIT
 !
 !----------------------------------------------------------------------
 ! Calculation of the log of Gamma_inv(1-m-eps)/eps
 ! ------------------------------------------------
 ! See previous function. 
 ! It is used in case Gamma_inv(1-m-eps)/eps might overflow.
 !
 ! Variables
 ! ---------
 ! m,eps: variable inputs of the function 
 ! (m,eps) -> log[Gamma_inv(1-m-eps)/eps]
 ! one_meps_mm: 1-eps-m
 ! i_Pi: i.Pi
 ! log_fact: logarithm of (-1)^m (m-1)!, 
 ! here defined as log((m-1)!) + i.Pi if m is odd.
 !----------------------------------------------------------------------
 function LOG_A_SUM_INIT(M,EPS)
   !--------------------------------------------------------------------
   implicit none
   integer(i4b),intent(in)   :: M
   complex(rkind),intent(in) :: EPS
   integer(i4b)              :: N
   real(rkind)               :: LOG_FACT
   complex(rkind)            :: ONE_MEPS_MM,LOG_A_SUM_INIT
   !
   ONE_MEPS_MM=ONE-EPS-M
   if(ONE_MEPS_MM.ne.1-M) then
      LOG_A_SUM_INIT=(-LOG_GAMMA_FUN(ONE_MEPS_MM) - LOG(EPS))
      return
   else
      LOG_FACT=ZERO
      do N=2,M-1
         LOG_FACT=LOG_FACT + LOG(DBLE(N))
      enddo
      if(MOD(M,2).eq.0) then
         LOG_A_SUM_INIT=LOG_FACT
      else
         LOG_A_SUM_INIT=CMPLX(LOG_FACT,M_PI,rkind)
      endif
      return
   endif
 end function LOG_A_SUM_INIT
 !----------------------------------------------------------------------
 ! Calculation of the first term of the B(z) power series
 ! ------------------------------------------------------
 ! in the 1-z transformation, divided by (1-z)^m
 ! ----------------------------------------------
 ! In the 1-z transformation, 
 ! the power series B(z) = \sum_{n=0}^{+oo} \beta_n (1-z)^n occurs 
 ! (see Comp. Phys. Comm. paper).
 ! The first term \beta_0, divided by (1-z)^m, is calculated here. 
 ! m is the closest integer to Re(c-a-b) >= 0 and eps = c-a-b-m.
 !
 ! One has to consider |eps|oo > 0.1 and |eps|oo <= 0.1, 
 ! where 1-m-eps and 1-m can be different or equal numerically, 
 ! leading to some changes in this last case.
 !
 ! |eps|oo > 0.1
 ! -------------
 ! One has \beta_0/(1-z)^m = [(a)_m (b)_m Gamma_inv(1-eps) 
 ! Gamma_inv(a+m+eps) Gamma_inv(b+m+eps) Gamma_inv(m+1)
 ! - (1-z)^eps Gamma_inv(a) Gamma_inv(b) Gamma_inv(1+m+eps)]
 ! [Gamma(c)/eps], stable in this regime for a direct evaluation.
 !
 ! The values of Gamma(c), Gamma_inv(a+m+eps) 
 ! and Gamma_inv(b+m+eps) were already calculated and recycled here.
 ! Gamma_inv(m+1) is calculated as 1/(m!).
 !
 ! Gamma_inv(1+m+eps) is calculated from Gamma_inv(1-eps), 
 ! using the equalities:
 ! Gamma_inv(1-m-eps) = Gamma_inv(1-eps) \prod_{i=1}^{m} (1-eps-i), 
 ! where the product is 1 by definition if m = 0,
 ! Gamma_inv(1+m+eps) = (-1)^m sin (pi.eps)
 ! /[pi.(eps+m).Gamma_inv(1-m-eps)] 
 ! from Euler reflection formula, Gamma(x+1) = x.Gamma(x) equality, 
 ! and m+eps no zero.
 ! This scheme is much faster than 
 ! to recalculate Gamma_inv(1+m+eps) directly.
 ! 
 ! |eps|oo <= 0.1
 ! --------------
 ! The \beta_0/(1-z)^m expression is rewritten 
 ! so that it contains no instabilities:
 ! \beta_0/(1-z)^m = Gamma_inv(a+m+eps) Gamma_inv(b+m+eps) 
 ! [(G(1,-eps) Gamma_inv(m+1) + G(m+1,eps))
 ! - Gamma_inv(1+m+eps) (G(a+m,eps) Gamma_inv(b+m+eps) 
 ! + G(b+m,eps) Gamma_inv(a+m)) 
 ! - E(log(1-z),eps) Gamma_inv(a+m) Gamma_inv(b+m) Gamma_inv(1+m+eps)] 
 ! (a)_m (b)_m Gamma(c)
 !
 ! E(log(1-z),eps) is [(1-z)^eps - 1]/eps 
 ! if 1-m-eps and 1-m are different numerically, 
 ! and log(1-z) otherwise (eps ~ 0).
 ! If 1-m-eps and 1-m are equal numerically, 
 ! Gamma_inv(1+m+eps) is numerically equal to Gamma_inv(1+m), 
 ! already calculated as 1/(m!).
 ! See |eps|oo > 0.1 case for data recycling of other values 
 ! or for 1-m-eps and 1-m different numerically.
 !
 !----------------------------------------------------------------------
 ! Variables
 ! ---------
 ! a,b,c,one_minus_z: a,b,c and 1-z parameters and arguments 
 ! of the 2F1(a,b,c,z) function.
 ! m,eps: closest integer to c-a-b, with Re(c-a-b) >= 0 
 ! and eps = c-a-b-m
 ! Gamma_c,Gamma_inv_one_meps,Gamma_inv_eps_pa_pm, Gamma_inv_eps_pb_pm: 
 ! recycled values of Gamma(c), Gamma_inv(1-eps), 
 ! Gamma_inv(a+m+eps) and Gamma_inv(b+m+eps).
 ! inf_norm_eps,phase,a_pm,b_pm,one_meps,Pi_eps,Pi_eps_pm: 
 ! |eps|oo,(-1)^m,a+m,b+m,1-eps,pi.eps,pi.(eps+m)
 ! Gamma_inv_one_meps_mm,Gamma_inv_eps_pm_p1: 
 ! Gamma_inv(1-m-eps) and Gamma_inv(1+m+eps) 
 ! calculated with the recycling scheme.
 ! prod1: (a)_m (b)_m Gamma_inv(1-eps) Gamma_inv(a+m+eps) 
 ! x Gamma_inv(b+m+eps) Gamma_inv(m+1) in |eps|oo > 0.1 case.
 ! prod2: (1-z)^eps Gamma_inv(a) Gamma_inv(b) Gamma_inv(1+m+eps) 
 ! in |eps|oo > 0.1 case.
 ! Gamma_inv_mp1,prod_ab: Gamma_inv(m+1) calculated as 1/(m!) 
 ! and (a)_m (b)_m in |eps|oo <= 0.1 case.
 ! is_eps_non_zero: true if 1-m-eps and 1-m are different numerically,
 ! false if not.
 ! Gamma_inv_a_pm,Gamma_inv_b_pm,z_term: Gamma_inv(a+m),Gamma_inv(b+m),
 ! E(eps,log(1-z))
 ! prod1: Gamma_inv(a+m+eps) Gamma_inv(b+m+eps) 
 ! x [(G(1,-eps) Gamma_inv(m+1) + G(m+1,eps)) in |eps|oo <= 0.1 case.
 ! prod2: Gamma_inv(1+m+eps) (G(a+m,eps) Gamma_inv(b+m+eps) 
 ! + G(b+m,eps) Gamma_inv(a+m))
 ! prod3: E(eps,log(1-z)) Gamma_inv(a+m) Gamma_inv(b+m) 
 ! Gamma_inv(1+m+eps) 
 ! res: returned \beta_0/(1-z)^m value in all cases.
 !----------------------------------------------------------------------
 function B_SUM_INIT_PS_ONE(A,B,GAMMA_C,GAMMA_INV_ONE_MEPS, &
      GAMMA_INV_EPS_PA_PM,GAMMA_INV_EPS_PB_PM,MZP1,M,EPS)
   !--------------------------------------------------------------------
   implicit none
   integer(i4b),intent(in)   :: M
   complex(rkind),intent(in) :: A,B,GAMMA_C,GAMMA_INV_ONE_MEPS
   complex(rkind),intent(in) :: GAMMA_INV_EPS_PA_PM,GAMMA_INV_EPS_PB_PM,MZP1,EPS
   integer(i4b)              :: M_M1,N,I,PHASE
   real(rkind)               :: INF_NORM_EPS,GAMMA_INV_MP1
   complex(rkind)            :: A_PM,B_SUM_INIT_PS_ONE,PI_EPS,GAMMA_INV_ONE_MEPS_MM
   complex(rkind)            :: B_PM,TMP1,TMP2
   complex(rkind)            :: Z_TERM,PROD1,PROD2,PROD3,ONE_MEPS,PI_EPS_PM
   complex(rkind)            :: GAMMA_INV_A_PM,PROD_AB,GAMMA_INV_B_PM
   complex(rkind)            :: GAMMA_INV_EPS_PM_P1
   !       
   INF_NORM_EPS=INF_NORM(EPS); M_M1=M-1; A_PM=A+M; B_PM=B+M
   ONE_MEPS=ONE-EPS; PI_EPS=M_PI*EPS; PI_EPS_PM = M_PI*(EPS+M)
   if(MOD(M,2).eq.0) then
      PHASE = 1
   else
      PHASE = -1
   endif
   GAMMA_INV_ONE_MEPS_MM = GAMMA_INV_ONE_MEPS
   do I=1,M
      GAMMA_INV_ONE_MEPS_MM = GAMMA_INV_ONE_MEPS_MM*(ONE_MEPS - I)
   enddo
   if(INF_NORM_EPS.gt.0.10_rkind) then
      GAMMA_INV_EPS_PM_P1 = PHASE*SIN(PI_EPS) &
           /(PI_EPS_PM*GAMMA_INV_ONE_MEPS_MM)
      PROD1=GAMMA_INV_ONE_MEPS*GAMMA_INV_EPS_PA_PM*GAMMA_INV_EPS_PB_PM
      do N=0,M_M1
         PROD1=PROD1*(A+N)*(B+N)/(N+ONE)
      enddo
      PROD2=GAMMA_INV(A)*GAMMA_INV(B)*GAMMA_INV_EPS_PM_P1*(MZP1**EPS)
      B_SUM_INIT_PS_ONE=GAMMA_C*(PROD1-PROD2)/EPS
      return
   else
      GAMMA_INV_MP1=ONE;PROD_AB=ONE
      do N=0,M_M1
         GAMMA_INV_MP1 = GAMMA_INV_MP1/(N+ONE)
         PROD_AB = PROD_AB*(A+N)*(B+N)
      enddo
      if(ONE_MEPS-M.ne.1-M) then
         Z_TERM=EXPM1(EPS*LOG(MZP1))/EPS
         GAMMA_INV_EPS_PM_P1 = PHASE*SIN(PI_EPS) &
              /(PI_EPS_PM*GAMMA_INV_ONE_MEPS_MM)
      else
         Z_TERM=LOG(MZP1)
         GAMMA_INV_EPS_PM_P1 = GAMMA_INV_MP1
      endif
      GAMMA_INV_A_PM=GAMMA_INV(A_PM);GAMMA_INV_B_PM=GAMMA_INV(B_PM)
      TMP1=ONE; TMP2=M+1;
      PROD1 = GAMMA_INV_EPS_PA_PM*GAMMA_INV_EPS_PB_PM    &
           *(GAMMA_INV_MP1*GAMMA_INV_DIFF_EPS(TMP1,-EPS) &
           +GAMMA_INV_DIFF_EPS(TMP2,EPS))
      PROD2 = GAMMA_INV_EPS_PM_P1 &
           *(GAMMA_INV_EPS_PB_PM*GAMMA_INV_DIFF_EPS(A_PM,EPS) &
           +GAMMA_INV_A_PM*GAMMA_INV_DIFF_EPS(B_PM,EPS))
      PROD3 = GAMMA_INV_A_PM*GAMMA_INV_B_PM*GAMMA_INV_EPS_PM_P1*Z_TERM
      B_SUM_INIT_PS_ONE=GAMMA_C*PROD_AB*(PROD1-PROD2-PROD3)
      return
   endif
 end function B_SUM_INIT_PS_ONE
 !
 !----------------------------------------------------------------------
 ! Calculation of the first term of the B(z) power series 
 ! ------------------------------------------------------
 ! in the 1/z transformation, divided by z^{-m}
 !---------------------------------------------
 ! In the 1/z transformation, the power series 
 ! B(z) = \sum_{n=0}^{+oo} \beta_n z^{-n} occurs 
 ! (see Comp. Phys. Comm. paper).
 ! The first term \beta_0, divided by z^{-m}, is calculated here. 
 ! m is the closest integer to Re(b-a) >= 0 and eps = b-a-m.
 !
 ! One has to consider |eps|oo > 0.1 and |eps|oo <= 0.1, 
 ! where 1-m-eps and 1-m can be different or equal numerically, 
 ! leading to some changes in this last case.
 !
 ! |eps|oo > 0.1
 ! -------------
 ! One has \beta_0/z^{-m} = [(a)_m (1-c+a)_m Gamma_inv(1-eps) 
 ! Gamma_inv(a+m+eps) Gamma_inv(c-a) Gamma_inv(m+1)
 ! - (-z)^{-eps} (1-c+a+eps)_m Gamma_inv(a) Gamma_inv(c-a-eps) 
 ! Gamma_inv(1+m+eps)].[Gamma(c)/eps], 
 ! stable in this regime for a direct evaluation.
 !
 ! The values of Gamma(c), Gamma_inv(c-a) and Gamma_inv(a+m+eps) 
 ! were already calculated and recycled here.
 ! Gamma_inv(m+1) is calculated as 1/(m!). 
 ! Gamma_inv(1+m+eps) is calculated from Gamma_inv(1-eps) 
 ! as in the 1-z transformation routine.
 ! 
 ! |eps|oo <= 0.1
 ! --------------
 ! The \beta_0/z^{-m} expression is rewritten 
 ! so that it contains no instabilities:
 ! \beta_0/z^{-m} = [((1-c+a+eps)_m G(1,-eps) - P(m,eps,1-c+a) 
 ! Gamma_inv(1-eps)) Gamma_inv(c-a) Gamma_inv(a+m+eps) Gamma_inv(m+1)
 ! + (1-c+a+eps)_m [G(m+1,eps) Gamma_inv(c-a) Gamma_inv(a+m+eps) 
 ! - G(a+m,eps) Gamma_inv(c-a) Gamma_inv(m+1+eps)]
 ! - (G(c-a,-eps) - E(log(-z),-eps)) Gamma_inv(m+1+eps) 
 ! Gamma_inv(a+m)]] (a)_m Gamma(c)
 !
 ! Definitions and method are the same 
 ! as in the 1-z transformation routine, except for P(m,eps,1-c+a).
 ! P(m,eps,s) = [(s+eps)_m - (s)_m]/eps 
 ! for eps non zero and has a limit for eps -> 0.
 ! Let n0 be the closest integer to -Re(s) for s complex. 
 ! A stable formula available for eps -> 0 for P(m,eps,s) is:
 ! P(m,eps,s) = (s)_m E(\sum_{n=0}^{m-1} L(1/(s+n),eps),eps) 
 ! if n0 is not in [0:m-1],
 ! P(m,eps,s) = \prod_{n=0, n not equal to n0}^{m-1} (s+eps+n) 
 ! + (s)_m E(\sum_{n=0, n not equal to n0}^{m-1} L(1/(s+n),eps),eps) 
 ! if n0 is in [0:m-1].
 ! L(s,eps) is log1p(s eps)/eps if eps is not zero, 
 ! and L(s,0) = s.
 ! This expression is used in the code.
 !
 ! Variables
 ! ---------
 ! a,b,c,z: a,b,c and z parameters 
 ! and arguments of the 2F1(a,b,c,z) function.
 ! m,eps: closest integer to b-a, with Re(b-a) >= 0 and eps = b-a-m.
 ! Gamma_c,Gamma_inv_cma,Gamma_inv_one_meps,Gamma_inv_eps_pa_pm: 
 ! recycled values of Gamma(c), Gamma_inv(c-a), Gamma_inv(1-eps) 
 ! and Gamma_inv(a+m+eps).
 ! inf_norm_eps,phase,cma,a_mc_p1,a_mc_p1_pm,cma_eps,eps_pa_mc_p1,a_pm: 
 ! |eps|oo,(-1)^m,c-a,1-c+a+m,c-a-eps,1-c+a+eps,a+m
 ! Gamma_inv_cma_meps,one_meps,Pi_eps,Pi_eps_pm: 
 ! Gamma_inv(c-a-eps),1-eps,pi.eps,pi.(eps+m)
 ! Gamma_inv_one_meps_mm,Gamma_inv_eps_pm_p1: Gamma_inv(1-m-eps) 
 ! and Gamma_inv(1+m+eps) calculated with the recycling scheme.
 ! prod1: (a)_m (1-c+a)_m Gamma_inv(1-eps) Gamma_inv(a+m+eps) 
 ! x Gamma_inv(c-a) Gamma_inv(m+1) in |eps|oo > 0.1 case.
 ! prod2: (-z)^{-eps} (1-c+a+eps)_m Gamma_inv(a) 
 ! x Gamma_inv(c-a-eps) Gamma_inv(1+m+eps) in |eps|oo > 0.1 case.
 ! n0: closest integer to -Re(1-c+a)
 ! is_n0_here: true is n0 belongs to [0:m-1], false if not.
 ! is_eps_non_zero: true if 1-m-eps and 1-m are different numerically, 
 ! false if not.
 ! Gamma_inv_mp1,prod_a,prod_a_mc_p1: 
 ! Gamma_inv(m+1) calculated as 1/(m!), 
 ! (a)_m and (1-c+a)_m in |eps|oo <= 0.1 case.
 ! prod_eps_pa_mc_p1_n0: 
 ! \prod_{n=0, n not equal to n0}^{m-1} (1-c+a+eps+n) 
 ! if n0 belongs to [0:m-1], 0.0 if not, in |eps|oo <= 0.1 case.
 ! prod_eps_pa_mc_p1: (1-c+a+eps)_m in |eps|oo <= 0.1 case.
 ! sum: \sum_{n=0, n not equal to n0}^{m-1} L(1/(s+n),eps) if 1-m-eps 
 ! and 1-m are different numerically, 
 ! \sum_{n=0, n not equal to n0}^{m-1} 1/(s+n) if not.
 ! a_pn,a_mc_p1_pn,eps_pa_mc_p1_pn: a+n,1-c+a+n,1-c+a+eps+n values 
 ! used in (a)_m, (1-c+a)_m and (1-c+a+eps)_m evaluations.
 ! sum_term,prod_diff_eps,z_term: 
 ! E(\sum_{n=0, n not equal to n0}^{m-1} L(1/(s+n),eps),eps), 
 ! P(m,eps,1-c+a), -E(-eps,log(-z))
 ! Gamma_inv_a_pm,Gamma_prod1: Gamma_inv(a+m), 
 ! Gamma_inv(c-a).Gamma_inv(a+m+eps)
 ! prod1: ((1-c+a+eps)_m G(1,-eps) 
 ! - P(m,eps,1-c+a) Gamma_inv(1-eps)) Gamma_inv(c-a) 
 ! x Gamma_inv(a+m+eps) Gamma_inv(m+1)
 ! prod_2a: Gamma_inv(c-a).Gamma_inv(a+m+eps).G(m+1,eps)
 ! prod_2b: G(a+m,eps) Gamma_inv(c-a) Gamma_inv(m+1+eps)
 ! prod_2c: (G(c-a,-eps) 
 ! - E(log(-z),-eps)) Gamma_inv(m+1+eps) Gamma_inv(a+m)
 ! prod2: (1-c+a+eps)_m [G(m+1,eps) Gamma_inv(c-a) Gamma_inv(a+m+eps) 
 ! - G(a+m,eps) Gamma_inv(c-a) Gamma_inv(m+1+eps)] 
 ! - (G(c-a,-eps) - E(log(-z),-eps)) 
 ! x Gamma_inv(m+1+eps) Gamma_inv(a+m)]]
 ! res: returned \beta_0/z^{-m} value in all cases.
 !----------------------------------------------------------------------
 function B_SUM_INIT_PS_INFINITY(A,C,GAMMA_C,GAMMA_INV_CMA, &
      GAMMA_INV_ONE_MEPS,GAMMA_INV_EPS_PA_PM,Z,M,EPS)
   !--------------------------------------------------------------------
   implicit none
   integer(i4b),intent(in)   :: M
   complex(rkind),intent(in) :: A,C,GAMMA_C,GAMMA_INV_CMA,Z,EPS
   complex(rkind),intent(in) :: GAMMA_INV_ONE_MEPS,GAMMA_INV_EPS_PA_PM
   integer(i4b)              :: M_M1,I,N,N0,PHASE
   logical(lgt)              :: IS_N0_HERE,IS_EPS_NON_ZERO
   real(rkind)               :: INF_NORM_EPS,NP1,GAMMA_INV_MP1
   complex(rkind)            :: B_SUM_INIT_PS_INFINITY,TMP1
   complex(rkind)            :: CMA,A_MC_P1,A_MC_P1_PM,CMA_MEPS,EPS_PA_MC_P1,A_PM
   complex(rkind)            :: GAMMA_INV_EPS_PM_P1,GAMMA_INV_CMA_MEPS,PI_EPS
   complex(rkind)            :: PROD1,PROD2,A_PN,A_MC_P1_PN,ONE_MEPS
   complex(rkind)            :: PROD_A,PROD_A_MC_P1,PROD_EPS_PA_MC_P1_N0,PI_EPS_PM
   complex(rkind)            :: PROD_EPS_PA_MC_P1,SUM_N0,Z_TERM,SUM_TERM
   complex(rkind)            :: PROD_DIFF_EPS,GAMMA_INV_A_PM,GAMMA_PROD1
   complex(rkind)            :: PROD_2A,PROD_2B,PROD_2C
   complex(rkind)            :: EPS_PA_MC_P1_PN,GAMMA_INV_ONE_MEPS_MM
   !
   INF_NORM_EPS=INF_NORM(EPS); CMA=C-A; A_MC_P1=A-C+ONE
   A_MC_P1_PM=A_MC_P1+M; CMA_MEPS=CMA-EPS; EPS_PA_MC_P1=EPS+A_MC_P1
   A_PM=A+M; M_M1=M-1; ONE_MEPS=ONE-EPS; PI_EPS=M_PI*EPS
   PI_EPS_PM=M_PI*(EPS+M); GAMMA_INV_CMA_MEPS=GAMMA_INV(CMA_MEPS)
   if(MOD(M,2).eq.0) then
      PHASE = 1
   else
      PHASE = -1
   endif
   GAMMA_INV_ONE_MEPS_MM = GAMMA_INV_ONE_MEPS
   do I=1,M
      GAMMA_INV_ONE_MEPS_MM = GAMMA_INV_ONE_MEPS_MM*(ONE_MEPS - I)
   enddo
   if(INF_NORM_EPS.gt.0.1_rkind) then
      GAMMA_INV_EPS_PM_P1 = PHASE*SIN(PI_EPS) &
           /(PI_EPS_PM*GAMMA_INV_ONE_MEPS_MM)
      PROD1 = GAMMA_INV_CMA*GAMMA_INV_EPS_PA_PM*GAMMA_INV_ONE_MEPS
      PROD2 = GAMMA_INV(A)*GAMMA_INV_CMA_MEPS*GAMMA_INV_EPS_PM_P1 &
           *((-Z)**(-EPS))
      do N=0,M_M1
         A_PN=A+N; A_MC_P1_PN=A_MC_P1+N
         EPS_PA_MC_P1_PN=EPS+A_MC_P1_PN;NP1=N+ONE
         PROD1 = PROD1*A_PN*A_MC_P1_PN/NP1
         PROD2 = PROD2*EPS_PA_MC_P1_PN
      enddo
      B_SUM_INIT_PS_INFINITY = GAMMA_C*(PROD1-PROD2)/EPS
      return
   else
      N0=-NINT(REAL(A_MC_P1,rkind))
      IS_EPS_NON_ZERO=ONE_MEPS-M.ne.1-M
      IS_N0_HERE=(N0.ge.0).and.(N0.lt.M)     
      GAMMA_INV_MP1=ONE; PROD_A=ONE; PROD_A_MC_P1=ONE
      PROD_EPS_PA_MC_P1=ONE; SUM_N0=ZERO
      if(IS_N0_HERE) then
         PROD_EPS_PA_MC_P1_N0 = ONE
      else
         PROD_EPS_PA_MC_P1_N0 = ZERO
      endif
      do N=0,M_M1
         A_PN=A+N; A_MC_P1_PN=A_MC_P1+N
         EPS_PA_MC_P1_PN=EPS+A_MC_P1_PN; NP1=N+ONE
         PROD_A = PROD_A*A_PN
         PROD_A_MC_P1 = PROD_A_MC_P1*A_MC_P1_PN
         PROD_EPS_PA_MC_P1 = PROD_EPS_PA_MC_P1*EPS_PA_MC_P1_PN
         GAMMA_INV_MP1 = GAMMA_INV_MP1/NP1
         if(N.ne.N0) then
            if(IS_N0_HERE) then
               PROD_EPS_PA_MC_P1_N0=PROD_EPS_PA_MC_P1_N0 &
                    *EPS_PA_MC_P1_PN
            endif
            if(IS_EPS_NON_ZERO) then
               SUM_N0 = SUM_N0 + LOG1P(EPS/A_MC_P1_PN)
            else
               SUM_N0 = SUM_N0 + ONE/A_MC_P1_PN
            endif
         endif
      enddo
      if(IS_EPS_NON_ZERO) then
         GAMMA_INV_EPS_PM_P1 = PHASE*SIN(PI_EPS) &
              /(PI_EPS_PM*GAMMA_INV_ONE_MEPS_MM)
         SUM_TERM = EXPM1(SUM_N0)/EPS
         Z_TERM = EXPM1(-EPS*LOG(-Z))/EPS
      else
         GAMMA_INV_EPS_PM_P1 = GAMMA_INV_MP1
         SUM_TERM = SUM_N0
         Z_TERM = -LOG(-Z)
      endif
      PROD_DIFF_EPS = PROD_EPS_PA_MC_P1_N0 + PROD_A_MC_P1*SUM_TERM
      GAMMA_INV_A_PM = GAMMA_INV(A_PM)
      GAMMA_PROD1=GAMMA_INV_CMA*GAMMA_INV_EPS_PA_PM
      TMP1=ONE
      PROD1 = GAMMA_PROD1*GAMMA_INV_MP1*(GAMMA_INV_DIFF_EPS(TMP1,-EPS) &
           *PROD_EPS_PA_MC_P1 - GAMMA_INV_ONE_MEPS*PROD_DIFF_EPS)
      TMP1=M+1
      PROD_2A = GAMMA_PROD1*GAMMA_INV_DIFF_EPS(TMP1,EPS) 
      PROD_2B = GAMMA_INV_CMA*GAMMA_INV_EPS_PM_P1  &
           *GAMMA_INV_DIFF_EPS(A_PM,EPS)
      PROD_2C = GAMMA_INV_EPS_PM_P1*GAMMA_INV_A_PM &
           *(GAMMA_INV_DIFF_EPS(CMA,-EPS) + GAMMA_INV_CMA_MEPS*Z_TERM)
      PROD2 = PROD_EPS_PA_MC_P1*(PROD_2A - PROD_2B - PROD_2C)
      B_SUM_INIT_PS_INFINITY = GAMMA_C*PROD_A*(PROD1+PROD2)
      return
   endif
 end function B_SUM_INIT_PS_INFINITY
 !
 !----------------------------------------------------------------------
 ! Calculation of the derivative of the polynomial P(X) 
 ! ----------------------------------------------------
 ! testing power series convergence
 ! --------------------------------
 ! P(X) = |z(a+X)(b+X)|^2 - |(c+X)(X+1)|^2 
 !      = \sum_{i=0}^{4} c[i] X^{i}, for |z| < 1.
 ! It is positive when the power series term modulus increases 
 ! and negative when it decreases, 
 ! so that its derivative provides information on its convergence 
 ! (see Comp. Phys. Comm. paper).
 ! Its derivative components cv_poly_der_tab[i] = (i+1) c[i+1] 
 ! for i in [0:3] 
 ! so that P'(X) = \sum_{i=0}^{3} cv_poly_der_tab[i] X^{i} 
 ! are calculated.
 !
 ! Variables:
 ! ----------
 ! a,b,c,z: a,b,c and z parameters and arguments 
 ! of the 2F1(a,b,c,z) function.
 ! cv_poly_der_tab[3]: table of four doubles 
 ! containing the P'(X) components.
 ! mod_a2,mod_b2,mod_c2,mod_z2,R_a,Re_b,Re_c: |a|^2, |b|^2, |c|^2, 
 ! |z|^2, Re(a), Re(b), Re(c), with which P(X) can be expressed.
 !----------------------------------------------------------------------
 subroutine CV_POLY_DER_TAB_CALC(A,B,C,Z,CV_POLY_DER_TAB)
   !--------------------------------------------------------------------
   implicit none
   complex(rkind),intent(in) :: A,B,C,Z
   real(rkind),intent(out)   :: CV_POLY_DER_TAB(0:3)
   real(rkind)               :: MOD_A2,MOD_B2,MOD_C2,MOD_Z2
   real(rkind)               :: RE_A,RE_B,RE_C,IM_A,IM_B,IM_C,RE_Z,IM_Z
   !
   RE_A=REAL(A,rkind); IM_A=AIMAG(A); MOD_A2=RE_A*RE_A+IM_A*IM_A
   RE_B=REAL(B,rkind); IM_B=AIMAG(B); MOD_B2=RE_B*RE_B+IM_B*IM_B
   RE_C=REAL(C,rkind); IM_C=AIMAG(C); MOD_C2=RE_C*RE_C+IM_C*IM_C
   RE_Z=REAL(Z,rkind); IM_Z=AIMAG(Z); MOD_Z2=RE_Z*RE_Z+IM_Z*IM_Z
   CV_POLY_DER_TAB(0)=TWO*((RE_A*MOD_B2+RE_B*MOD_A2)*MOD_Z2-RE_C-MOD_C2)
   CV_POLY_DER_TAB(1)=TWO*((MOD_A2+MOD_B2+4._rkind*RE_A*RE_B)*MOD_Z2 &
        -ONE-4._rkind*RE_C-MOD_C2)
   CV_POLY_DER_TAB(2)=6._rkind*((RE_A+RE_B)*MOD_Z2-RE_C-ONE)
   CV_POLY_DER_TAB(3)=4._rkind*(MOD_Z2-ONE)
 end subroutine CV_POLY_DER_TAB_CALC
 !
 !----------------------------------------------------------------------
 ! Calculation of the derivative of the polynomial P(X) 
 ! ----------------------------------------------------
 ! testing power series convergence at one x value
 ! -----------------------------------------------
 ! P'(x) is calculated for a real x. 
 ! See P'(X) components calculation routine for definitions.
 !----------------------------------------------------------------------
 function CV_POLY_DER_CALC(CV_POLY_DER_TAB,X)
   !--------------------------------------------------------------------
   implicit none
   real(rkind),intent(in) :: X
   real(rkind),intent(in) :: CV_POLY_DER_TAB(0:3)
   real(rkind)            :: CV_POLY_DER_CALC
   !
   CV_POLY_DER_CALC=CV_POLY_DER_TAB(0)+X*(CV_POLY_DER_TAB(1) &
        +X*(CV_POLY_DER_TAB(2)+X*CV_POLY_DER_TAB(3)))
   return
 end function CV_POLY_DER_CALC
 !
 !----------------------------------------------------------------------
 ! Calculation of an integer after which false convergence cannot occur
 ! --------------------------------------------------------------------
 ! See cv_poly_der_tab_calc routine for definitions.
 ! If P'(x) < 0 and P''(x) < 0 for x > xc, it will be so for all x > xc 
 ! as P(x) -> -oo for x -> +oo 
 ! and P(x) can have at most one maximum for x > xc. 
 ! It means that the 2F1 power series term modulus will increase 
 ! or decrease to 0 for n > nc, 
 ! with nc the smallest positive integer larger than xc.
 !
 ! If P'(X) = C0 + C1.X + C2.X^2 + C3.X^3, 
 ! the discriminant of P''(X) is Delta = C2^2 - 3 C1 C3.
 !
 ! If Delta > 0, P''(X) has two different real roots 
 ! and its largest root is -(C2 + sqrt(Delta))/(3 C3), 
 ! because C3 = 4(|z|^2 - 1) < 0.
 ! One can take xc = -(C2 + sqrt(Delta))/(3 C3) 
 ! and one returns its associated nc integer.
 !
 ! If Delta <= 0, P''(X) has at most one real root, 
 ! so that P'(X) has only one root and then P(X) only one maximum.
 ! In this case, one can choose xc = nc = 0, which is returned.
 !
 ! Variables
 ! ---------
 ! cv_poly_der_tab: table of four doubles 
 ! containing the P'(X) coefficients
 ! C1,C2,three_C3: cv_poly_der_tab[1], cv_poly_der_tab[2] 
 ! and 3.0*cv_poly_der_tab[3], so that P''(X) = C1 + 2.C2.x + three_C3.x^2
 ! Delta: discriminant of P''(X), equal to C2^2 - 3 C1 C3.
 ! largest_root: if Delta > 0, 
 ! P''(X) largest real root equal to -(C2 + sqrt(Delta))/(3 C3).
 !----------------------------------------------------------------------
 function MIN_N_CALC(CV_POLY_DER_TAB)
   !--------------------------------------------------------------------
   implicit none
   real(rkind),intent(in) :: CV_POLY_DER_TAB(0:3)
   integer(i4b)           :: MIN_N_CALC
   real(rkind)            :: C1,C2,THREE_C3,DELTA,LARGEST_ROOT
   !
   C1=CV_POLY_DER_TAB(1); C2=CV_POLY_DER_TAB(2)
   THREE_C3=3._rkind*CV_POLY_DER_TAB(3); DELTA = C2*C2 - THREE_C3*C1
   if(DELTA.le.ZERO) then
      MIN_N_CALC = 0
      return
   else
      LARGEST_ROOT = -(C2 + SQRT (DELTA))/THREE_C3
      MIN_N_CALC = MAX(CEILING(LARGEST_ROOT),0)
      return
   endif
 end function MIN_N_CALC
 !
 !----------------------------------------------------------------------
 ! Calculation of the 2F1 power series converging for |z| < 1
 ! ----------------------------------------------------------
 ! One has 2F1(a,b,c,z) 
 ! = \sum_{n = 0}^{+oo} (a)_n (b)_n / ((c)_n n!) z^n,
 ! so that 2F1(a,b,c,z) = \sum_{n = 0}^{+oo} t[n] z^n, 
 ! with t[0] = 1 and t[n+1] = (a+n)(b+n)/((c+n)(n+1)) t[n] for n >= 0.
 ! If a or b are negative integers, 
 ! F(z) is a polynomial of degree -a or -b, evaluated directly.
 ! If not, one uses the test of convergence |t[n] z^n|oo < 1E-15 
 ! to truncate the series after it was checked 
 ! that false convergence cannot occur.
 ! Variables:
 ! ----------
 ! a,b,c,z: a,b,c and z parameters and arguments 
 ! of the 2F1(a,b,c,z) function. One must have here |z| < 1.
 ! term,sum: term of the 2F1 power series equal to t[n] z^n, 
 ! truncated sum at given n of the 2F1 power series.
 ! na,nb: absolute values of the closest integers to Re(a) and Re(b). 
 ! a = -na or b = -nb means one is in the polynomial case.
 ! cv_poly_der_tab: coefficients of the derivative 
 ! of the polynomial P(X) = |z(a+X)(b+X)|^2 - |(c+X)(X+1)|^2
 ! min_n: smallest integer after which false convergence cannot occur. 
 ! It is calculated in min_n_calc.
 ! possible_false_cv: always true if n < min_n. 
 ! If n >= min_n, it is true if P'(n) > 0. 
 ! If n >= min_n and P'(n) < 0, 
 ! it becomes false and remains as such for the rest of the calculation. 
 ! One can then check if |t[n] z^n|oo < 1E-15 to truncate the series.
 !----------------------------------------------------------------------
 function HYP_PS_ZERO(A,B,C,Z)
   !--------------------------------------------------------------------
   implicit none
   complex(rkind),intent(in) :: A,B,C,Z
   integer(i4b)              :: N,NA,NB,MIN_N
   complex(rkind)            :: HYP_PS_ZERO,TERM
   logical(lgt)              :: POSSIBLE_FALSE_CV
   real(rkind)               :: CV_POLY_DER_TAB(0:3)
   !
   NA = ABS(NINT(REAL(A,rkind)))
   NB = ABS(NINT(REAL(B,rkind)))
   TERM=ONE; HYP_PS_ZERO=ONE  
   if(A.eq.(-NA)) then
      do N=0,NA-1
         TERM = TERM*Z*(A+N)*(B+N)/((N+ONE)*(C+N))
         HYP_PS_ZERO = HYP_PS_ZERO + TERM
      enddo
      return
   else if(B.eq.(-NB)) then
      do N=0,NB-1
         TERM = TERM*Z*(A+N)*(B+N)/((N+ONE)*(C+N))
         HYP_PS_ZERO = HYP_PS_ZERO + TERM
      enddo
      return
   else
      call CV_POLY_DER_TAB_CALC(A,B,C,Z,CV_POLY_DER_TAB)
      POSSIBLE_FALSE_CV=.TRUE.
      MIN_N=MIN_N_CALC(CV_POLY_DER_TAB);N=0
      do while(POSSIBLE_FALSE_CV.or.(INF_NORM(TERM).gt.EPS15))
         TERM = TERM*Z*(A+N)*(B+N)/((N+ONE)*(C+N))
         HYP_PS_ZERO = HYP_PS_ZERO + TERM
         if(POSSIBLE_FALSE_CV.and.(N.gt.MIN_N)) then
            POSSIBLE_FALSE_CV = &
                 (CV_POLY_DER_CALC(CV_POLY_DER_TAB,DBLE(N)).gt.ZERO)
         endif
         N=N+1 
      enddo
      return
   endif
 end function HYP_PS_ZERO
 !
 !----------------------------------------------------------------------
 ! Calculation of the 2F1 power series 
 ! -----------------------------------
 ! converging with the 1-z transformation
 ! --------------------------------------
 ! The formula for F(z) in the 1-z transformation holds:
 ! F(z) = (-1)^m (pi.eps)/sin (pi.eps) [A(z) + B(z)] 
 ! for eps not equal to zero, F(z) = (-1)^m [A(z) + B(z)] for eps = 0
 ! where m = |Re(c-a-b)], eps = c-a-b-m, 
 ! A(z) = \sum_{n=0}^{m-1} alpha[n] (1-z)^n, 
 ! B(z) = \sum_{n=0}^{+oo} beta[n] (1-z)^n, and:
 !
 ! alpha[0] = [Gamma_inv(1-m-eps)/eps] Gamma_inv(a+m+eps) 
 !          x Gamma_inv(b+m+eps) Gamma(c)
 ! [Gamma_inv(1-m-eps)/eps] is calculated in A_sum_init. 
 ! alpha[0] is calculated with log[Gamma] 
 ! if the previous expression might overflow, 
 ! and its imaginary part removed if a, b and c are real.
 ! alpha[n+1] = (a+n)(b+n)/[(n+1)(1-m-eps+n)] alpha[n], n in [0:m-2].
 !
 ! beta[0] is defined in B_sum_init_PS_one function comments.
 ! gamma[0] = Gamma(c) (a)_m (b)_m (1-z)^m Gamma_inv(a+m+eps) 
 !          x Gamma_inv(b+m+eps) Gamma_inv(m+1) Gamma_inv(1-eps)
 !
 ! beta[n+1] = (a+m+n+eps)(b+m+n+eps)/[(m+n+1+eps)(n+1)] beta[n]
 ! + [(a+m+n)(b+m+n)/(m+n+1) - (a+m+n) - (b+m+n) - eps 
 ! + (a+m+n+eps)(b+m+n+eps)/(n+1)]
 !             x gamma[n]/[(n+m+1+eps)(n+1+eps)], n >= 0.
 ! gamma[n+1] = (a+m+n)(b+m+n)/[(m+n+1)(n+1-eps)] gamma[n], n >= 0.
 !
 ! B(z) converges <=> |1-z| < 1
 ! The test of convergence is |beta[n] (1-z)^n|oo < 1E-15 |beta[0]|oo
 ! for n large enough so that false convergence cannot occur.
 !
 ! Variables
 ! ---------
 ! a,b,c,one_minus_z: a,b,c parameters 
 ! and 1-z from z argument of 2F1(a,b,c,z)
 ! m,phase,m_p1,eps,eps_pm,eps_pm_p1,
 ! a_pm,b_pm,one_meps,one_meps_pm: 
 ! |Re(c-a-b)], (-1)^m, m+1, c-a-b-m, 
 ! eps+m, eps+m+1, a+m, b+m, 1-eps, 1-eps-m
 ! eps_pa,eps_pb,eps_pa_pm,eps_pb_pm,Pi_eps,Gamma_c: 
 ! eps+a, eps+b, eps+a+m, eps+b+m, pi.eps, Gamma(c)
 ! Gamma_inv_eps_pa_pm,Gamma_inv_eps_pb_pm,Gamma_prod: 
 ! Gamma_inv(eps+a+m), Gamma_inv(eps+b+m), 
 ! Gamma(c).Gamma_inv(eps+a+m).Gamma_inv(eps+b+m)
 ! Gamma_inv_one_meps,A_first_term,A_sum,A_term: 
 ! Gamma_inv(1-eps), alpha[0], A(z), alpha[n] (1-z)^n
 ! pow_mzp1_m,B_first_term,prod_B,ratio: (1-z)^m, beta[0], 
 ! (a)_m (b)_m (1-z)^m, (a+n)(b+n)/(n+1) for n in [0:m-2].
 ! B_extra_term,B_term,B_sum,B_prec: 
 ! gamma[n], beta[n] (1-z)^n, B(z), 1E-15 |beta[0|oo
 ! cv_poly1_der_tab,cv_poly2_der_tab: P1'(X) and P2'(X) coefficients 
 ! of the potentials derivatives of P1(X) and P2(X) 
 ! defined in cv_poly_der_tab_calc with parameters 
 ! a1 = a, b1 = b, c1 = 1-m-eps, z1 = 1-z 
 ! and a2 = eps+b+m, b2 = eps+a+m,c2 = eps+m+1, z2 = 1-z.
 ! min_n: smallest integer after which false convergence cannot occur. 
 ! It is calculated in min_n_calc with both P1'(X) and P2'(X), 
 ! so one takes the largest integer coming from both calculations.
 ! possible_false_cv: always true if n < min_n. 
 ! If n >= min_n, it is true if P1'(n) > 0 or P2'(n) > 0. 
 ! If n >= min_n and P1'(n) < 0 and P2'(n) < 0, 
 ! it becomes false and remains as such for the rest of the calculation.
 ! One can then check if |beta[n] z^n|oo < 1E-15 to truncate the series.
 ! n,n_pm_p1,n_p1,a_pm_pn,b_pm_pn,eps_pm_p1_pn,n_p1_meps,eps_pa_pm_pn,
 ! eps_pb_pm_pn,eps_pm_pn: index of power series, n+m+1, n+1, 
 ! a+m+n, b+m+n, eps+m+n+1, n+1-eps, eps+a+m+n, eps+b+m+n, eps+m+n,
 ! prod1,prod2,prod3: (eps+a+m+n)(eps+b+m+n), 
 ! (eps+m+1+n)(n+1), (a+m+n)(b+m+n)
 !----------------------------------------------------------------------
 function HYP_PS_ONE(A,B,C,MZP1)
   !--------------------------------------------------------------------
   implicit none
   complex(rkind),intent(in) :: A,B,C,MZP1
   integer(i4b)              :: N,M,PHASE,M_M2,MIN_N,M_P1
   real(rkind)               :: B_PREC,N_P1,N_PM_P1 
   complex(rkind)            :: HYP_PS_ONE,EPS,EPS_PM,EPS_PM_P1,A_PM
   complex(rkind)            :: B_PM,ONE_MEPS_MM,EPS_PA,EPS_PB,PI_EPS,GAMMA_PROD
   complex(rkind)            :: EPS_PA_PM,EPS_PB_PM, A_SUM,A_TERM,ONE_MEPS
   complex(rkind)            :: B_EXTRA_TERM,B_TERM,B_SUM,GAMMA_C,RATIO
   complex(rkind)            :: A_PM_PN,B_PM_PN,EPS_PM_P1_PN,N_P1_MEPS
   complex(rkind)            :: PROD1,PROD2,PROD3
   complex(rkind)            :: EPS_PA_PM_PN,EPS_PB_PM_PN,EPS_PM_PN,PROD_B,POW_MZP1_M
   complex(rkind)            :: GAMMA_INV_EPS_PA_PM,GAMMA_INV_EPS_PB_PM
   complex(rkind)            :: GAMMA_INV_ONE_MEPS
   logical(lgt)              :: POSSIBLE_FALSE_CV
   real(rkind)               :: CV_POLY1_DER_TAB(0:3),CV_POLY2_DER_TAB(0:3)
   !
   M=NINT(REAL(C-A-B,rkind)); M_M2=M-2; M_P1=M+1
   if(MOD(M,2).eq.0) then
      PHASE=1
   else
      PHASE=-1
   endif
   EPS=C-A-B-M; EPS_PM=EPS+M; EPS_PM_P1=EPS_PM+ONE; A_PM=A+M;B_PM=B+M
   ONE_MEPS=ONE-EPS; ONE_MEPS_MM=ONE_MEPS-M; EPS_PA=EPS+A; EPS_PB=EPS+B 
   PI_EPS=M_PI*EPS; EPS_PA_PM=EPS_PA+M; EPS_PB_PM=EPS_PB+M
   GAMMA_C=ONE/GAMMA_INV(C)
   GAMMA_INV_EPS_PA_PM=GAMMA_INV(EPS_PA_PM)
   GAMMA_INV_EPS_PB_PM=GAMMA_INV(EPS_PB_PM)
   GAMMA_PROD=GAMMA_C*GAMMA_INV_EPS_PA_PM*GAMMA_INV_EPS_PB_PM
   GAMMA_INV_ONE_MEPS=GAMMA_INV(ONE_MEPS)
   if(M.eq.0) then
      A_TERM=ZERO
   else if(INF_NORM(ONE_MEPS_MM &
        *(LOG(ONE + ABS(ONE_MEPS_MM))-ONE)).lt.300.0d0) then
      A_TERM=GAMMA_PROD*A_SUM_INIT(M,EPS,GAMMA_INV_ONE_MEPS)
   else
      A_TERM=EXP(LOG_GAMMA_FUN(C)-LOG_GAMMA_FUN(EPS_PA_PM)&
           -LOG_GAMMA_FUN(EPS_PB_PM)+LOG_A_SUM_INIT(M,EPS))
      if((AIMAG(A).eq.ZERO).and.(AIMAG(B).eq.ZERO)&
           .and.(AIMAG(C).eq.ZERO)) then
         A_TERM=REAL(A_TERM,rkind)
      endif
   endif
   A_SUM=A_TERM
   POW_MZP1_M = MZP1**M
   B_TERM=B_SUM_INIT_PS_ONE(A,B,GAMMA_C,GAMMA_INV_ONE_MEPS, &
        GAMMA_INV_EPS_PA_PM,GAMMA_INV_EPS_PB_PM,MZP1,M,EPS)*POW_MZP1_M
   PROD_B=POW_MZP1_M
   do N=0,M_M2
      RATIO=(A+N)*(B+N)/(N+ONE)
      A_TERM=A_TERM*MZP1*RATIO/(N+ONE_MEPS_MM)
      A_SUM=A_SUM+A_TERM
      PROD_B = PROD_B*RATIO
   enddo
   if(M.gt.0) then
      PROD_B = PROD_B*(A+M-ONE)*(B+M-ONE)/DBLE(M)
   endif
   B_EXTRA_TERM = PROD_B*GAMMA_PROD*GAMMA_INV_ONE_MEPS; B_SUM=B_TERM
   B_PREC=EPS15*INF_NORM(B_TERM)
   call CV_POLY_DER_TAB_CALC(A,B,ONE_MEPS_MM,MZP1,CV_POLY1_DER_TAB)
   call CV_POLY_DER_TAB_CALC(EPS_PB_PM,EPS_PA_PM,EPS_PM_P1,MZP1, &
        CV_POLY2_DER_TAB)
   MIN_N=MAX(MIN_N_CALC(CV_POLY1_DER_TAB),MIN_N_CALC(CV_POLY2_DER_TAB))
   POSSIBLE_FALSE_CV=.TRUE.; N=0
   do while(POSSIBLE_FALSE_CV.or.(INF_NORM(B_TERM).gt.B_PREC))
      N_PM_P1=N+M_P1; N_P1=N+ONE; A_PM_PN=A_PM+N; B_PM_PN=B_PM+N
      EPS_PM_P1_PN=EPS_PM_P1+N; N_P1_MEPS=ONE_MEPS+N
      EPS_PM_PN=EPS_PM+N; EPS_PA_PM_PN=EPS_PA_PM+N 
      EPS_PB_PM_PN=EPS_PB_PM+N
      PROD1=EPS_PA_PM_PN*EPS_PB_PM_PN
      PROD2=EPS_PM_P1_PN*N_P1
      PROD3=A_PM_PN*B_PM_PN
      B_TERM = MZP1*(B_TERM*PROD1/PROD2+B_EXTRA_TERM*(PROD3/N_PM_P1 &
           -A_PM_PN-B_PM_PN-EPS+PROD1/N_P1)/(EPS_PM_P1_PN*N_P1_MEPS))
      B_SUM=B_SUM+B_TERM
      B_EXTRA_TERM=B_EXTRA_TERM*MZP1*PROD3/(N_PM_P1*N_P1_MEPS)
      if(POSSIBLE_FALSE_CV.and.(N.gt.MIN_N)) then
         POSSIBLE_FALSE_CV = &
              (CV_POLY_DER_CALC(CV_POLY1_DER_TAB,DBLE(N)).gt.ZERO).or. &
              (CV_POLY_DER_CALC(CV_POLY2_DER_TAB,DBLE(N)).gt.ZERO)
      endif
      N=N+1
   enddo
   if(EPS.eq.ZERO) then
      HYP_PS_ONE=PHASE*(A_SUM+B_SUM)
      return
   else
      HYP_PS_ONE=PHASE*(A_SUM+B_SUM)*PI_EPS/SIN(PI_EPS)
      return
   endif
 end function HYP_PS_ONE
 !
 !----------------------------------------------------------------------
 ! Calculation of the 2F1 power series 
 ! -----------------------------------
 ! converging with the 1/z transformation
 ! --------------------------------------
 ! The formula for F(z) in the 1/z transformation holds:
 ! F(z) = (-1)^m (pi.eps)/sin (pi.eps) [A(z) + B(z)] 
 ! for eps not equal to zero, 
 ! F(z) = (-1)^m [A(z) + B(z)] for eps = 0
 ! where m = |Re(b-a)], eps = b-a-m, 
 ! A(z) = \sum_{n=0}^{m-1} alpha[n] z^{-n}, 
 ! B(z) = \sum_{n=0}^{+oo} beta[n] z^{-n}, and:
 !
 ! alpha[0] = [Gamma_inv(1-m-eps)/eps] Gamma_inv(c-a) 
 !          x Gamma_inv(a+m+eps) Gamma(c)
 ! [Gamma_inv(1-m-eps)/eps] is calculated in A_sum_init. 
 ! alpha[0] is calculated with log[Gamma] 
 ! if the previous expression might overflow, 
 ! and its imaginary part removed if a, b and c are real.
 ! alpha[n+1] = (a+n)(1-c+a+n)/[(n+1)(1-m-eps+n)] alpha[n], 
 ! n in [0:m-2].
 !
 ! beta[0] is defined in B_sum_init_PS_infinity function comments.
 ! gamma[0] = Gamma(c) (a)_m (1-c+a)_m z^{-m} Gamma_inv(a+m+eps) 
 !          x Gamma_inv(c-a) Gamma_inv(m+1) Gamma_inv(1-eps)
 !
 ! beta[n+1] = (a+m+n+eps)(1-c+a+m+n+eps)/[(m+n+1+eps)(n+1)] beta[n] 
 ! + [(a+m+n)(1-c+a+m+n)/(m+n+1) - (a+m+n) - (1-c+a+m+n) 
 ! - eps + (a+m+n+eps)(1-c+a+m+n+eps)/(n+1)]
 ! x gamma[n]/[(n+m+1+eps)(n+1+eps)], n >= 0.
 ! gamma[n+1] = (a+m+n)(b+m+n)/[(m+n+1)(n+1-eps)] gamma[n], n >= 0.
 !
 ! B(z) converges <=> |z| > 1
 ! The test of convergence is |beta[n] z^{-n}|oo < 1E-15 |beta[0]|oo
 ! for n large enough so that false convergence cannot occur.
 !
 ! Variables
 ! ---------
 ! a,b,c,z: a,b,c parameters and z argument of 2F1(a,b,c,z)
 ! m,phase,m_p1,eps,a_mc_p1,one_meps,
 ! one_meps_pm,a_pm,a_mc_p1_pm,cma: |Re(b-a)], (-1)^m, m+1, b-a-m, 
 ! 1-c+a, 1-eps, 1-eps-m, a+m, 1-c+a+m, c-a
 ! eps_pa,eps_pm_p1,eps_pa_mc_p1_pm,Pi_eps,eps_pa_pm,eps_pm,Gamma_c: 
 ! eps+a, eps+m+1, eps+1-c+a+m, pi.eps, eps+a+m, eps+m, Gamma(c)
 ! Gamma_inv_eps_pa_pm,Gamma_inv_cma,z_inv,pow_mz_ma,
 ! Gamma_inv_one_meps,Gamma_prod: Gamma_inv(eps+a+m), Gamma_inv(c-a), 
 ! 1/z, (-z)^(-a), Gamma_inv(1-eps), 
 ! Gamma(c) Gamma_inv(c-a) Gamma_inv(eps+a+m)
 ! A_first_term,A_sum,A_term: alpha[0], A(z), alpha[n] z^{-n}
 ! pow_z_inv_m,B_first_term,prod_B,ratio: z^{-m}, beta[0], 
 ! (a)_m (1-c+a)_m z^{-m}, (a+n)(1-c+a+n)/(n+1) for n in [0:m-2].
 ! B_extra_term,B_term,B_sum,B_prec: 
 ! gamma[n], beta[n] z^{-n}, B(z), 1E-15 |beta[0|oo
 ! cv_poly1_der_tab,cv_poly2_der_tab: P1'(X) and P2'(X) coefficients 
 ! of the potentials derivatives of P1(X) and P2(X) 
 ! defined in cv_poly_der_tab_calc 
 ! with parameters a1 = a, b1 = 1-c+a, c1 = 1-m-eps, z1 = 1/z 
 ! and a2 = b, b2 = eps+1-c+a+m,c2 = eps+m+1, z2 = 1/z.
 ! min_n: smallest integer after which false convergence cannot occur. 
 !        It is calculated in min_n_calc with both P1'(X) and P2'(X), 
 ! so one takes the largest integer coming from both calculations.
 ! possible_false_cv: always true if n < min_n. If n >= min_n, 
 ! it is true if P1'(n) > 0 or P2'(n) > 0. 
 ! If n >= min_n and P1'(n) < 0 and P2'(n) < 0, 
 ! it becomes false and remains as such for the rest of the calculation. 
 ! One can then check if |beta[n] z^n|oo < 1E-15 to truncate the series.
 ! n,n_pm_p1,n_p1,a_pm_pn,a_mc_p1_pm_pn,eps_pm_p1_pn,n_p1_meps,
 ! eps_pa_pm_pn,eps_pa_mc_p1_pm_pn,eps_pm_pn: 
 ! index of power series, n+m+1, n+1, a+m+n, 1-c+a+m+n, eps+m+n+1,
 ! n+1-eps, eps+a+m+n, eps+1-c+a+m+n, eps+m+n,
 ! prod1,prod2,prod3: (eps+a+m+n)(eps+1-c+a+m+n),
 ! (eps+m+1+n)(n+1), (a+m+n)(1-c+a+m+n)
 !----------------------------------------------------------------------
 function HYP_PS_INFINITY(A,B,C,Z)
   !--------------------------------------------------------------------
   implicit none
   complex(rkind),intent(in) :: A,B,C,Z
   integer(i4b)              :: N,M,PHASE,M_M2,MIN_N,M_P1
   real(rkind)               :: B_PREC,N_P1,N_PM_P1
   complex(rkind)            :: POW_Z_INV_M,HYP_PS_INFINITY,Z_INV,RATIO
   complex(rkind)            :: EPS,A_MC_P1,ONE_MEPS,ONE_MEPS_MM,A_PM,A_MC_P1_PM
   complex(rkind)            :: CMA,EPS_PA,EPS_PM_P1,EPS_PA_MC_P1_PM,PI_EPS
   complex(rkind)            :: EPS_PA_PM,EPS_PM,GAMMA_C,GAMMA_INV_CMA,POW_MZ_MA
   complex(rkind)            :: A_SUM,A_TERM
   complex(rkind)            :: GAMMA_INV_EPS_PA_PM,GAMMA_INV_ONE_MEPS
   complex(rkind)            :: PROD_B,B_EXTRA_TERM,B_TERM,B_SUM,PROD1
   complex(rkind)            :: A_PM_PN,A_MC_P1_PM_PN,EPS_PM_P1_PN,N_P1_MEPS
   complex(rkind)            :: PROD2,PROD3,GAMMA_PROD
   complex(rkind)            :: EPS_PA_PM_PN,EPS_PA_MC_P1_PM_PN,EPS_PM_PN
   logical(lgt)              :: POSSIBLE_FALSE_CV
   real(rkind)               :: CV_POLY1_DER_TAB(0:3),CV_POLY2_DER_TAB(0:3)
   !
   M=NINT(REAL(B-A,rkind)); M_M2=M-2;M_P1=M+1
   if(MOD(M,2).eq.0) then
      PHASE=1
   else
      PHASE=-1
   endif
   EPS=B-A-M; A_MC_P1=ONE-C+A; ONE_MEPS=ONE-EPS; ONE_MEPS_MM=ONE_MEPS-M
   A_PM=A+M; A_MC_P1_PM=A_MC_P1+M; CMA=C-A; EPS_PA=EPS+A
   EPS_PM=EPS+M; EPS_PM_P1=EPS_PM+ONE; EPS_PA_MC_P1_PM=EPS+A_MC_P1_PM
   PI_EPS=M_PI*EPS; EPS_PA_PM=EPS_PA+M
   GAMMA_C=ONE/GAMMA_INV(C); GAMMA_INV_EPS_PA_PM = GAMMA_INV(EPS_PA_PM)
   GAMMA_INV_ONE_MEPS = GAMMA_INV(ONE_MEPS)
   GAMMA_INV_CMA=GAMMA_INV(CMA); Z_INV=ONE/Z;POW_MZ_MA=(-Z)**(-A)
   GAMMA_PROD=GAMMA_C*GAMMA_INV_CMA*GAMMA_INV_EPS_PA_PM
   if(M.eq.0) then
      A_TERM=ZERO
   else if(INF_NORM(ONE_MEPS_MM &
        *(LOG(ONE + ABS(ONE_MEPS_MM))-ONE)).lt.300._rkind) then
      A_TERM=GAMMA_PROD*A_SUM_INIT(M,EPS,GAMMA_INV_ONE_MEPS)
   else
      A_TERM=EXP(LOG_GAMMA_FUN(C)-LOG_GAMMA_FUN(CMA)-LOG_GAMMA_FUN(B) &
           + LOG_A_SUM_INIT(M,EPS))
      if((AIMAG(A).eq.ZERO).and.(AIMAG(B).eq.ZERO).and.     &
           (AIMAG(C).eq.ZERO)) then
         A_TERM=REAL(A_TERM,rkind)
      endif
   endif
   A_SUM=A_TERM
   POW_Z_INV_M=Z_INV**M
   B_TERM=B_SUM_INIT_PS_INFINITY(A,C,GAMMA_C,GAMMA_INV_CMA, &
        GAMMA_INV_ONE_MEPS,GAMMA_INV_EPS_PA_PM,Z,M,EPS)*POW_Z_INV_M
   PROD_B=POW_Z_INV_M
   do N=0,M_M2
      RATIO=(A+N)*(A_MC_P1+N)/(N+ONE)
      A_TERM = A_TERM*Z_INV*RATIO/(N+ONE_MEPS_MM)
      A_SUM = A_SUM+A_TERM
      PROD_B = PROD_B*RATIO
   enddo
   if (M.gt.0) then
      PROD_B=PROD_B*(A+M-ONE)*(A_MC_P1+M-ONE)/DBLE(M)
   endif
   B_EXTRA_TERM = PROD_B*GAMMA_PROD*GAMMA_INV_ONE_MEPS
   B_SUM=B_TERM
   B_PREC=EPS15*INF_NORM(B_TERM)
   call CV_POLY_DER_TAB_CALC(A,A_MC_P1,ONE_MEPS_MM,Z_INV, &
        CV_POLY1_DER_TAB)
   call CV_POLY_DER_TAB_CALC(B,EPS_PA_MC_P1_PM,EPS_PM_P1, &
        Z_INV,CV_POLY2_DER_TAB)
   MIN_N=MAX(MIN_N_CALC(CV_POLY1_DER_TAB),MIN_N_CALC(CV_POLY2_DER_TAB))
   POSSIBLE_FALSE_CV=.TRUE.; N=0
   do while(POSSIBLE_FALSE_CV.or.(INF_NORM(B_TERM).gt.B_PREC))
      N_PM_P1=N+M_P1; N_P1=N+ONE; A_PM_PN=A_PM+N
      A_MC_P1_PM_PN=A_MC_P1_PM+N; EPS_PM_P1_PN=EPS_PM_P1+N
      N_P1_MEPS=N_P1-EPS; EPS_PA_PM_PN=EPS_PA_PM+N
      EPS_PA_MC_P1_PM_PN=EPS_PA_MC_P1_PM+N; EPS_PM_PN=EPS_PM+N
      PROD1=EPS_PA_PM_PN*EPS_PA_MC_P1_PM_PN; PROD2=EPS_PM_P1_PN*N_P1
      PROD3=A_PM_PN*A_MC_P1_PM_PN
      B_TERM = Z_INV*(B_TERM*PROD1/PROD2+B_EXTRA_TERM*(PROD3/N_PM_P1 &
           -A_PM_PN-A_MC_P1_PM_PN-EPS+PROD1/N_P1)                    &
           /(EPS_PM_P1_PN*N_P1_MEPS))
      B_SUM=B_SUM+B_TERM
      B_EXTRA_TERM=B_EXTRA_TERM*Z_INV*PROD3/(N_PM_P1*N_P1_MEPS)
      if(POSSIBLE_FALSE_CV.and.(N.gt.MIN_N)) then
         POSSIBLE_FALSE_CV = (CV_POLY_DER_CALC( &
              CV_POLY1_DER_TAB,DBLE(N)).gt.ZERO).or.(&
              CV_POLY_DER_CALC(CV_POLY2_DER_TAB,DBLE(N)).gt.ZERO)
      endif
      N=N+1
   enddo
   if(EPS.eq.ZERO) then
      HYP_PS_INFINITY=PHASE*POW_MZ_MA*(A_SUM+B_SUM)
      return
   else
      HYP_PS_INFINITY=PHASE*POW_MZ_MA*(A_SUM+B_SUM)*PI_EPS &
           /SIN(PI_EPS)
      return
   endif
 end function HYP_PS_INFINITY
 !
 !----------------------------------------------------------------------
 ! Calculation of F(z) in transformation theory missing zones 
 ! ----------------------------------------------------------
 ! of the complex plane with a Taylor series
 ! -----------------------------------------
 ! If z is close to exp(+/- i.pi/3), no transformation in 1-z, z, 
 ! z/(z-1) or combination of them can transform z in a complex number 
 ! of modulus smaller than a given Rmax < 1 .
 ! Rmax is a radius for which one considers power series summation 
 ! for |z| > Rmax is too slow to be processed. One takes Rmax = 0.9 .
 ! Nevertheless, for Rmax = 0.9, 
 ! these zones are small enough to be handled 
 ! with a Taylor series expansion around a point z0 close to z 
 ! where transformation theory can be used to calculate F(z).
 ! One then chooses z0 to be 0.9 z/|z| if |z| < 1, and 1.1 z/|z| 
 ! if |z| > 1, 
 ! so that hyp_PS_zero or hyp_PS_infinity can be used 
 ! (see comments of these functions above).
 ! For this z0, F(z) = \sum_{n=0}^{+oo} q[n] (z-z0)^n, with:
 ! q[0] = F(z0), q[1] = F'(z0) = (a b/c) 2F1(a+1,b+1,c+1,z0)
 ! q[n+2] = [q[n+1] (n (2 z0 - 1) - c + (a+b+c+1) z0) 
 ! + q[n] (a+n)(b+n)/(n+1)]/(z0(1-z0)(n+2))
 ! As |z-z0| < 0.1, it converges with around 15 terms, 
 ! so that no instability can occur for moderate a, b and c.
 ! Convergence is tested 
 ! with |q[n] (z-z0)^n|oo + |q[n+1] (z-z0)^{n+1}|oo. 
 ! Series is truncated when this test is smaller 
 ! than 1E-15 (|q[0]|oo + |q[1] (z-z0)|oo).
 ! No false convergence can happen here 
 ! as q[n] behaves smoothly for n -> +oo.
 !
 ! Variables
 ! ---------
 ! a,b,c,z: a,b,c parameters and z argument of 2F1(a,b,c,z)
 ! abs_z,is_abs_z_small: |z|, true if |z| < 1 and false if not.
 ! z0,zc_z0_ratio,z0_term1,z0_term2: 0.9 z/|z| if |z| < 1, 
 ! and 1.1 z/|z| if |z| > 1, (z-z0)/(z0 (1-z0)), 
 ! 2 z0 - 1, c - (a+b+c+1) z0
 ! hyp_PS_z0,dhyp_PS_z0,prec: F(z0), F'(z0) calculated with 2F1 
 ! as F'(z0) = (a b/c) 2F1(a+1,b+1,c+1,z0), 
 ! precision demanded for series truncation 
 ! equal to 1E-15 (|q[0]|oo + |q[1] (z-z0)|oo).
 ! n,an,anp1,anp2,sum: index of the series, q[n] (z-z0)^n, 
 ! q[n+1] (z-z0)^{n+1}, q[n+2] (z-z0)^{n+2}, 
 ! truncated sum of the power series.
 !----------------------------------------------------------------------
 function HYP_PS_COMPLEX_PLANE_REST(A,B,C,Z)
   !--------------------------------------------------------------------
   implicit none
   complex(rkind),intent(in) :: A,B,C,Z
   integer(i4b)              :: N
   real(rkind)               :: ABS_Z,PREC
   complex(rkind)            :: HYP_PS_COMPLEX_PLANE_REST
   complex(rkind)            :: Z0,ZC,ZC_Z0_RATIO,Z0_TERM1,Z0_TERM2
   complex(rkind)            :: HYP_PS_Z0,DHYP_PS_Z0,AN,ANP1,ANP2
   !
   ABS_Z=ABS(Z)
   if(ABS_Z.lt.ONE) then
      Z0=0.9_rkind*Z/ABS_Z; ZC=Z-Z0; ZC_Z0_RATIO=ZC/(Z0*(ONE-Z0))
      Z0_TERM1=TWO*Z0 - ONE; Z0_TERM2=C-(A+B+ONE)*Z0
      HYP_PS_Z0=HYP_PS_ZERO(A,B,C,Z0)
      DHYP_PS_Z0=HYP_PS_ZERO(A+ONE,B+ONE,C+ONE,Z0)*A*B/C 
   else
      Z0=1.1_rkind*Z/ABS_Z; ZC=Z-Z0; ZC_Z0_RATIO=ZC/(Z0*(ONE-Z0))
      Z0_TERM1=TWO*Z0 - ONE; Z0_TERM2=C-(A+B+ONE)*Z0
      HYP_PS_Z0=HYP_PS_INFINITY(A,B,C,Z0)
      DHYP_PS_Z0=HYP_PS_INFINITY(A+ONE,B+ONE,C+ONE,Z0)*A*B/C 
   endif
   AN=HYP_PS_Z0;ANP1=ZC*DHYP_PS_Z0;HYP_PS_COMPLEX_PLANE_REST=AN+ANP1
   PREC=EPS15*(INF_NORM(AN)+INF_NORM(ANP1)); N=0
   do while(INF_NORM(AN).gt.PREC)
      ANP2=ZC_Z0_RATIO*(ANP1*(N*Z0_TERM1-Z0_TERM2)+AN*ZC*(A+N)*(B+N) &
           /(N+ONE))/(N+TWO)
      HYP_PS_COMPLEX_PLANE_REST = HYP_PS_COMPLEX_PLANE_REST + ANP2
      N=N+1
      AN=ANP1
      ANP1=ANP2
   enddo
   return
 end function HYP_PS_COMPLEX_PLANE_REST
 !
 !----------------------------------------------------------------------
 ! Calculation of F(z) for arbitrary z using previous routines
 ! -----------------------------------------------------------
 ! Firstly, it is checked if a,b and c are negative integers.
 ! If neither a nor b is negative integer but c is, 
 ! F(z) is undefined so that the program stops with an error message.
 ! If a and c are negative integers with c < a, 
 ! or b and c are negative integers with b < a, 
 ! or c is not negative integer integer but a or b is, 
 ! one is in the polynomial case.
 ! In this case, if |z| < |z/(z-1)| or z = 1, 
 ! hyp_PS_zero is used directly, as then |z| <= 2 
 ! and no instability arises with hyp_PS_zero 
 ! as long the degree of the polynomial is small (<= 10 typically).
 ! If not, one uses the transformation 
 ! F(z) = (1-z)^{-a} 2F1(a,c-b,c,z/(z-1)) if a is negative integer 
 ! or F(z) = (1-z)^{-b} 2F1(b,c-a,c,z/(z-1)) if b is negative integer 
 ! along with hyp_PS_zero.
 ! Indeed, 2F1(a,c-b,c,X) is a polynomial if a is negative integer, 
 ! and so is 2F1(b,c-a,c,X) if b is negative integer, 
 ! so that one has here |z/(z-1)| <= 2 
 ! and the stability of the method is the same 
 ! as for the |z| < |z/(z-1)| case.
 ! If one is in the non-polynomial case, one checks if z >= 1. 
 ! If it is, one is the cut of F(z) 
 ! so that z is replaced by z - 10^{-307}i.
 ! Then, using F(z) = 2F1(b,a,c,z) 
 ! and F(z) = (1-z)^{c-a-b} 2F1(c-a,c-b,c,z), 
 ! one replaces a,b,c parameters by combinations of them 
 ! so that Re(b-a) >= 0 and Re(c-a-b) >= 0.
 ! Exchanging a and b does not change convergence properties, 
 ! while having Re(c-a-b) >= 0 accelerates it 
 ! (In hyp_PS_zero, t[n] z^n ~ z^n/(n^{c-a-b}) for n -> +oo).
 ! If |1-z| < 1E-5, one uses hyp_PS_one 
 ! as the vicinity of the singular point z = 1 is treated properly.
 ! After that, one compares |z| and |z/(z-1)| 
 ! to R in {0.5,0.6,0.7,0.8,0.9}. 
 ! If one of them is smaller than R, 
 ! one uses hyp_PS_zero without transformation
 ! or with the transformation F(z) = (1-z)^{-a} 2F1(a,c-b,c,z/(z-1)).
 ! Then, if both of them are larger than 0.9, 
 ! one compares |1/z|, |(z-1)/z|, |1-z| and |1/(1-z)| 
 ! to R in {0.5,0.6,0.7,0.8,0.9}. 
 ! If one of them is found smaller than R, 
 ! with the condition that |c-b|oo < 5 for (z-1)/z transformation, 
 ! |a,b,c|oo < 5 for |1-z| transformation 
 ! and |a,c-b,c|oo < 5 for |1/(1-z)| transformation,
 ! the corresponding transformation is used. 
 ! If none of them was smaller than 0.9, 
 ! one is in the missing zones of transformation theory 
 ! so that the Taylor series of hyp_PS_complex_plane_rest is used.
 !
 ! Variables
 ! ---------
 ! a,b,c,z: a,b,c parameters and z argument of 2F1(a,b,c,z)
 ! Re_a,Re_b,Re_c,na,nb,nc,is_a_neg_int,is_b_neg_int,is_c_neg_int: 
 ! real parts of a,b,c, closest integers to a,b,c, 
 ! true if a,b,c is negative integers and false if not.
 ! zm1,z_over_zm1,z_shift: z-1, z/(z-1), z - 10^{-307}i in case z >= 1.
 ! ab_condition, cab_condition: true if Re(b-a) >= 0 and false if not, 
 ! true if Re(c-a-b) >= 0 and false if not.
 ! abs_zm1,abz_z,abs_z_inv,abs_z_over_zm1,abs_zm1_inv,abs_zm1_over_z: 
 ! |z-1|, |z|, |1/z|, |z/(z-1)|, |1/(z-1)|, |(z-1)/z|
 ! are_ac_small: true if |a|oo < 5 and |c|oo < 5, false if not.
 ! is_cmb_small: true if |c-b|oo < 5, false if not.
 ! are_abc_small: true if |a|oo < 5, |b|oo < 5 and |c|oo < 5, 
 ! false if not.
 ! are_a_cmb_c_small: true if |a|oo < 5, |c-b|oo < 5 and |c|oo < 5, 
 ! false if not.
 ! R_tab,R: table of radii {0.5,0.6,0.7,0.8,0.9}, one of these radii.
 ! res: returned result
 !----------------------------------------------------------------------
 recursive function HYP_2F1(A,B,C,Z) result(RES)
   !--------------------------------------------------------------------
   implicit none
   complex(rkind),intent(in) :: A,B,C,Z
   integer(i4b)              :: NA,NB,NC,I
   real(rkind)               :: RE_A,RE_B,RE_C,ABS_Z,ABS_ZM1,ABS_Z_OVER_ZM1
   real(rkind)               :: ABS_ZM1_OVER_Z,ABS_ZM1_INV,R_TABLE(1:5),R,ABS_Z_INV
   complex(rkind)            :: RES,Z_SHIFT,Z_OVER_ZM1,ZM1
   logical(lgt)              :: IS_A_NEG_INT,IS_B_NEG_INT,IS_C_NEG_INT
   logical(lgt)              :: AB_CONDITION,CAB_CONDITION,ARE_A_CMB_C_SMALL
   logical(lgt)              :: IS_CMB_SMALL,ARE_AC_SMALL,ARE_ABC_SMALL
   !
   RE_A=REAL(A,rkind); RE_B=REAL(B,rkind); RE_C=REAL(C,rkind);
   NA=NINT(RE_A); NB=NINT(RE_B); NC=NINT(RE_C);
   IS_A_NEG_INT=A.eq.NA.and.NA.le.0
   IS_B_NEG_INT=B.eq.NB.and.NB.le.0
   IS_C_NEG_INT=C.eq.NC.and.NC.le.0
   ZM1=Z-ONE
   if(IS_C_NEG_INT) then
      ABS_Z=ABS(Z); Z_OVER_ZM1 = Z/ZM1
      ABS_Z_OVER_ZM1=ABS(Z_OVER_ZM1)
      if(IS_A_NEG_INT.and.(NC.lt.NA)) then
         if((Z.eq.ONE).or.(ABS_Z.lt.ABS_Z_OVER_ZM1)) then
            RES=HYP_PS_ZERO(A,B,C,Z)
            return
         else
            RES=((-ZM1)**(-A))*HYP_PS_ZERO(A,C-B,C,Z_OVER_ZM1)
            return
         endif
      else if(IS_B_NEG_INT.and.(NC.lt.NB)) then
         if((Z.eq.ONE).or.(ABS_Z.lt.ABS_Z_OVER_ZM1)) then
            RES=HYP_PS_ZERO(A,B,C,Z)
            return
         else
            RES=((-ZM1)**(-B))*HYP_PS_ZERO(B,C-A,C,Z_OVER_ZM1)
            return
         endif
      else
         print*,'2F1 UNDEFINED'
      endif
   endif
   if(IS_A_NEG_INT) then
      ABS_Z=ABS(Z); Z_OVER_ZM1 = Z/ZM1
      ABS_Z_OVER_ZM1=ABS(Z_OVER_ZM1)
      if((Z.eq.ONE).or.(ABS_Z.lt.ABS_Z_OVER_ZM1)) then
         RES=HYP_PS_ZERO(A,B,C,Z)
         return
      else
         RES=((-ZM1)**(-A))*HYP_PS_ZERO(A,C-B,C,Z_OVER_ZM1)
         return
      endif
   else if(IS_B_NEG_INT) then
      ABS_Z=ABS(Z); Z_OVER_ZM1 = Z/ZM1
      ABS_Z_OVER_ZM1=ABS(Z_OVER_ZM1)
      if((Z.eq.ONE).or.(ABS_Z.lt.ABS_Z_OVER_ZM1)) then
         RES=HYP_PS_ZERO(A,B,C,Z)
         return
      else
         RES=((-ZM1)**(-B))*HYP_PS_ZERO(B,C-A,C,Z_OVER_ZM1)
         return
      endif
   endif
   if((REAL(Z,rkind).ge.ONE).and.(AIMAG(Z).eq.ZERO)) then
      Z_SHIFT=CMPLX(REAL(Z,rkind),-1.e-307_rkind,rkind)
      RES=HYP_2F1(A,B,C,Z_SHIFT)
      return
   endif
   AB_CONDITION = (RE_B.ge.RE_A)
   CAB_CONDITION = (RE_C.ge.RE_A + RE_B)
   if ((.NOT.AB_CONDITION).or.(.NOT.CAB_CONDITION)) then
      if ((.NOT.AB_CONDITION).and.(CAB_CONDITION)) then
         RES=HYP_2F1(B,A,C,Z)
         return
      else if((.NOT.CAB_CONDITION).and.(AB_CONDITION)) then
         RES=((-ZM1)**(C-A-B))*HYP_2F1(C-B,C-A,C,Z)
         return
      else
         RES=((-ZM1)**(C-A-B))*HYP_2F1(C-A,C-B,C,Z)
         return 
      endif
   endif
   ABS_ZM1=ABS(ZM1)
   if(ABS_ZM1.lt.1.e-5_rkind) then 
      RES=HYP_PS_ONE(A,B,C,-ZM1)
      return
   endif
   ABS_Z=ABS(Z); ABS_Z_OVER_ZM1=ABS_Z/ABS_ZM1; ABS_Z_INV=ONE/ABS_Z
   ABS_ZM1_OVER_Z=ONE/ABS_Z_OVER_ZM1; ABS_ZM1_INV=ONE/ABS_ZM1
   IS_CMB_SMALL = INF_NORM(C-B).lt.5._rkind; 
   ARE_AC_SMALL = (INF_NORM(A).lt.5._rkind).and.(INF_NORM(C).lt.5._rkind)
   ARE_ABC_SMALL = ARE_AC_SMALL.and.(INF_NORM(B).lt.5._rkind)
   ARE_A_CMB_C_SMALL = ARE_AC_SMALL.and.IS_CMB_SMALL
   R_TABLE=(/0.5_rkind,0.6_rkind,0.7_rkind,0.8_rkind,0.9_rkind/)
   do I=1,5
      R=R_TABLE(I)
      if(ABS_Z.le.R) then 
         RES=HYP_PS_ZERO(A,B,C,Z)
         return
      endif
      if(IS_CMB_SMALL.and.(ABS_Z_OVER_ZM1.le.R)) then
         RES=((-ZM1)**(-A))*HYP_PS_ZERO(A,C-B,C,Z/ZM1)
         return
      endif
   enddo
   do I=1,5
      R=R_TABLE(I)
      if(ABS_Z_INV.le.R) then 
         RES=HYP_PS_INFINITY(A,B,C,Z)
         return 
      endif
      if(IS_CMB_SMALL.and.(ABS_ZM1_OVER_Z.le.R)) then 
         RES=((-ZM1)**(-A))*HYP_PS_INFINITY(A,C-B,C,Z/ZM1)
         return
      endif
      if(ARE_ABC_SMALL.and.(ABS_ZM1.le.R)) then 
         RES=HYP_PS_ONE(A,B,C,-ZM1)
         return
      endif
      if(ARE_A_CMB_C_SMALL.and.(ABS_ZM1_INV.le.R)) then 
         RES=((-ZM1)**(-A))*HYP_PS_ONE(A,C-B,C,-ONE/ZM1)
         return
      endif
   enddo
   RES=HYP_PS_COMPLEX_PLANE_REST(A,B,C,Z)
   return
 end function HYP_2F1
 !
 !----------------------------------------------------------------------
 ! Test of 2F1 numerical accuracy 
 ! ------------------------------
 ! using hypergeometric differential equation
 ! ------------------------------------------
 ! If z = 0, F(z) = 1 so that this value is trivially tested.
 ! To test otherwise if the value of F(z) is accurate, 
 ! one uses the fact that 
 ! z(z-1) F''(z) + (c - (a+b+1) z) F'(z) - a b F(z) = 0.
 ! If z is not equal to one, a relative precision test is provided 
 ! by |F''(z) + [(c - (a+b+1) z) F'(z) - a b F(z)]/[z(z-1)]|oo
 ! /(|F(z)|oo + F'(z)|oo + |F''(z)|oo).
 ! If z is equal to one, one uses |(c - (a+b+1)) F'(z) - a b F(z)|oo
 ! /(|F(z)|oo + F'(z)|oo + 1E-307).
 ! F'(z) and F''(z) are calculated using equalities 
 ! F'(z) = (a b/c) 2F1(a+1,b+1,c+1,z) 
 ! and F'(z) = ((a+1)(b+1)/(c+1)) (a b/c) 2F1(a+2,b+2,c+2,z).
 !
 ! Variables
 ! ---------
 ! a,b,c,z: a,b,c parameters and z argument of 2F1(a,b,c,z)
 ! F,dF,d2F: F(z), F'(z) and F''(z) calculated with hyp_2F1 
 ! using F'(z) = (a b/c) 2F1(a+1,b+1,c+1,z) 
 ! and F'(z) = ((a+1)(b+1)/(c+1)) (a b/c) 2F1(a+2,b+2,c+2,z).
 !----------------------------------------------------------------------
 function TEST_2F1(A,B,C,Z,F)
   !--------------------------------------------------------------------
   implicit none
   complex(rkind),intent(in) :: A,B,C,Z
   real(rkind)               :: TEST_2F1
   complex(rkind)            :: F,DF,D2F  
   !
   if(Z.eq.ZERO) then
      TEST_2F1=INF_NORM(F-ONE)
      return
   else if(Z.eq.ONE) then
      DF = HYP_2F1(A+ONE,B+ONE,C+ONE,Z)*A*B/C
      TEST_2F1=INF_NORM((C-(A+B+ONE))*DF-A*B*F) &
           /(INF_NORM (F)+INF_NORM(DF)+1.e-307_rkind)
      return
   else
      DF = HYP_2F1(A+ONE,B+ONE,C+ONE,Z)*A*B/C
      D2F = HYP_2F1(A+TWO,B+TWO,C+TWO,Z)*A*(A+ONE)*B*(B+ONE) &
           /(C*(C+ONE))
      TEST_2F1=INF_NORM(D2F+((C-(A+B+ONE)*Z)*DF-A*B*F)/(Z*(ONE-Z))) &
           /(INF_NORM(F)+INF_NORM(DF)+INF_NORM(D2F))
      return
   endif
 end function TEST_2F1
 !============== END HYP_2F1 FILE ======================================

end module hyp_2F1_module

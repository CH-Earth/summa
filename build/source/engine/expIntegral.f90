module expIntegral_module
USE nrtype
implicit none
private
public::expint
contains

 ! Numerical recipes routines removed; use code from UEB-Veg

 ! ****************** EXPONENTIAL INTEGRAL FUNCTION *****************************************
 ! From UEB-Veg
 ! Computes the exponential integral function for the given value
 FUNCTION EXPINT (LAI)
 real(summa_prec) LAI
 real(summa_prec) EXPINT
 real(summa_prec) a0,a1,a2,a3,a4,a5,b1,b2,b3,b4
 real(summa_prec),parameter :: verySmall=tiny(1.0_summa_prec)     ! a very small number
 IF (LAI < verySmall)THEN
  EXPINT=1._summa_prec

 ELSEIF (LAI.LE.1.0) THEN
  a0=-.57721566_summa_prec
  a1=.99999193_summa_prec
  a2=-.24991055_summa_prec
  a3=.05519968_summa_prec
  a4=-.00976004_summa_prec
  a5=.00107857_summa_prec

  EXPINT = a0+a1*LAI+a2*LAI**2+a3*LAI**3+a4*LAI**4+a5*LAI**5 - log(LAI)

 ELSE
  a1=8.5733287401_summa_prec
  a2=18.0590169730_summa_prec
  a3=8.6347637343_summa_prec
  a4=.2677737343_summa_prec
  b1=9.5733223454_summa_prec
  b2=25.6329561486_summa_prec
  b3=21.0996530827_summa_prec
  b4=3.9584969228_summa_prec

  EXPINT=(LAI**4+a1*LAI**3+a2*LAI**2+a3*LAI+a4)/ &
      ((LAI**4+b1*LAI**3+b2*LAI**2+b3*LAI+b4)*LAI*exp(LAI))

 END IF
 RETURN
 END FUNCTION EXPINT

END MODULE expIntegral_module

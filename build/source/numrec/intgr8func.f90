module integr8func_module
USE nrtype
implicit none
private
public::qromb
contains

 ! ---------------------------------------------------------------------------------------------------------------------------
 ! public subroutine: integrate function "func" from a to b
 ! ---------------------------------------------------------------------------------------------------------------------------
 FUNCTION qromb(func,a,b)
 IMPLICIT NONE
 REAL(DP), INTENT(IN) :: a,b
 REAL(DP) :: qromb
 INTERFACE
  FUNCTION func(x)
  USE nrtype
  REAL(DP), DIMENSION(:), INTENT(IN) :: x
  REAL(DP), DIMENSION(size(x)) :: func
  END FUNCTION func
 END INTERFACE
 INTEGER(I4B), PARAMETER :: JMAX=20,JMAXP=JMAX+1,K=5,KM=K-1
 REAL(DP), PARAMETER :: EPS=1.0e-10_dp
 REAL(DP), DIMENSION(JMAXP) :: h,s
 REAL(DP) :: dqromb
 INTEGER(I4B) :: j
 h(1)=1.0_dp
 do j=1,JMAX
  call trapzd(func,a,b,s(j),j)
  if (j >= K) then
   call polint(h(j-KM:j),s(j-KM:j),0.0_dp,qromb,dqromb)
   if (abs(dqromb) <= EPS*abs(qromb)) RETURN
  end if
  s(j+1)=s(j)
  h(j+1)=0.25_dp*h(j)
 end do
 stop 'qromb: too many steps'
 END FUNCTION qromb


 ! ---------------------------------------------------------------------------------------------------------------------------
 ! private subroutine: use the trapezoidal rule to integrate function "func" from a to b, using 2**(n-2) interior points
 ! ---------------------------------------------------------------------------------------------------------------------------
 SUBROUTINE trapzd(func,a,b,s,n)
 USE nr_utility_module, ONLY : arth
 IMPLICIT NONE
 REAL(DP), INTENT(IN) :: a,b
 REAL(DP), INTENT(INOUT) :: s
 INTEGER(I4B), INTENT(IN) :: n
 INTERFACE
  FUNCTION func(x)
  USE nrtype
  REAL(DP), DIMENSION(:), INTENT(IN) :: x
  REAL(DP), DIMENSION(size(x)) :: func
  END FUNCTION func
 END INTERFACE
 REAL(DP) :: del,fsum
 INTEGER(I4B) :: it
 if (n == 1) then
  s=0.5_dp*(b-a)*sum(func( (/ a,b /) ))
 else
  it=2**(n-2)
  del=(b-a)/it
  fsum=sum(func(arth(a+0.5_dp*del,del,it)))
  s=0.5_dp*(s+del*fsum)
 end if
 END SUBROUTINE trapzd


 ! ---------------------------------------------------------------------------------------------------------------------------
 ! private subroutine: use a polynomial to interpolate between the points used in the trapezoidal function
 ! ---------------------------------------------------------------------------------------------------------------------------
 SUBROUTINE polint(xa,ya,x,y,dy)
 IMPLICIT NONE
 REAL(DP), DIMENSION(:), INTENT(IN) :: xa,ya
 REAL(DP), INTENT(IN) :: x
 REAL(DP), INTENT(OUT) :: y,dy
 INTEGER(I4B) :: m,n,ns
 REAL(DP), DIMENSION(size(xa)) :: c,d,den,ho
 n=size(xa); if(n/=size(ya)) stop 'polint: size mis-match'
 c=ya
 d=ya
 ho=xa-x
 ns=iminloc(abs(x-xa))
 y=ya(ns)
 ns=ns-1
 do m=1,n-1
  den(1:n-m)=ho(1:n-m)-ho(1+m:n)
  if (any(den(1:n-m) == 0.0)) stop 'polint: calculation failure'
  den(1:n-m)=(c(2:n-m+1)-d(1:n-m))/den(1:n-m)
  d(1:n-m)=ho(1+m:n)*den(1:n-m)
  c(1:n-m)=ho(1:n-m)*den(1:n-m)
  if (2*ns < n-m) then
   dy=c(ns+1)
  else
   dy=d(ns)
   ns=ns-1
  end if
  y=y+dy
 end do
 END SUBROUTINE polint


 ! ---------------------------------------------------------------------------------------------------------------------------
 ! private function: utility to estimate the location of the minimum value as a single integer, rather than one-element vector
 ! ---------------------------------------------------------------------------------------------------------------------------
 FUNCTION iminloc(arr)
 REAL(DP), DIMENSION(:), INTENT(IN) :: arr
 INTEGER(I4B), DIMENSION(1) :: imin
 INTEGER(I4B) :: iminloc
 imin=minloc(arr(:))
 iminloc=imin(1)
 END FUNCTION iminloc


end module integr8func_module

module expIntegral_module
USE nrtype
implicit none
private
public::expIntegral
contains

 SUBROUTINE expIntegral(n,x,expint,err,message)
 USE nr_utility_module, ONLY : arth
 IMPLICIT NONE
 ! dummy variables
 INTEGER(I4B), INTENT(IN) :: n
 REAL(DP), INTENT(IN) :: x
 REAL(DP),INTENT(OUT) :: expint
 integer(i4b),intent(out)      :: err
 character(*),intent(out)      :: message
 ! local variables
 INTEGER(I4B), PARAMETER :: MAXIT=100
 REAL(DP), PARAMETER :: EPS=epsilon(x),BIG=huge(x)*EPS
 INTEGER(I4B) :: i,nm1
 REAL(DP) :: a,b,c,d,del,fact,h
 ! initialize error control
 err=0; message="expint/"
 ! check that the size of the vectors match
 if(.not. ((n >= 0) .and. (x >= 0.0_dp) .and. (x > 0.0_dp .or. n > 1)))then
  message=trim(message)//'incorrect arguments'
  err=20; return
 endif
 
 ! n=0
 if (n == 0) then
  expint=exp(-x)/x
  RETURN
 end if
 
 ! x=0
 nm1=n-1
 if (x == 0.0) then
  expint=1.0_dp/nm1
 
 ! continued fraction
 else if (x > 1.0) then
  b=x+n
  c=BIG
  d=1.0_dp/b
  h=d
  do i=1,MAXIT
   a=-i*(nm1+i)
   b=b+2.0_dp
   d=1.0_dp/(a*d+b)
   c=b+a/c
   del=c*d
   h=h*del
   if (abs(del-1.0_dp) <= EPS) exit
  end do
  if (i > MAXIT)then
   message=trim(message)//'continued fraction failed (exceeded maximum number of iterations)'
   err=20; return
  endif
  expint=h*exp(-x)
 
 ! series
 else
  if (nm1 /= 0) then
   expint=1.0_dp/nm1
  else
   expint=-log(x)-EULER
  end if
  fact=1.0
  do i=1,MAXIT
   fact=-fact*x/i
   if (i /= nm1) then
    del=-fact/(i-nm1)
   else
    del=fact*(-log(x)-EULER+sum(1.0_dp/arth(1,1,nm1)))
   end if
   expint=expint+del
   if (abs(del) < abs(expint)*EPS) exit
  end do
  if (i > MAXIT)then
   message=trim(message)//'series failed (exceeded maximum number of iterations)'
   err=20; return
  endif
 end if

 END SUBROUTINE expIntegral

END MODULE expIntegral_module

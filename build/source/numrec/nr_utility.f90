module nr_utility_module
USE nrtype
! contains functions that should really be part of the fortran standard, but are not
implicit none
INTEGER(I4B), PARAMETER :: NPAR_ARTH=16,NPAR2_ARTH=8
INTERFACE outerdiff
 MODULE PROCEDURE outerdiff_r,outerdiff_d,outerdiff_i
END INTERFACE
INTERFACE arth
 MODULE PROCEDURE arth_r, arth_d, arth_i
END INTERFACE
! (everything private unless otherwise specifed)
private
! matrix intrinsics
public::get_diag
public::put_diag
public::unit_matrix
public::vabs
public::lower_triangle
public::outerdiff
public::outerprod
public::ifirstloc
! the error functiom
public::erf
! build vectors of regularly spaced numbers
public::arth
contains

 ! *************************************************************************************************
 ! *************************************************************************************************
 ! (A) matrix intrinsics
 ! *************************************************************************************************
 ! *************************************************************************************************

 ! *************************************************************************************************
 ! (A1) extract the diagonal vector from a matrix
 ! *************************************************************************************************
 FUNCTION get_diag(mat)
 REAL(DP), DIMENSION(:,:), INTENT(IN) :: mat
 REAL(DP), DIMENSION(size(mat,1)) :: get_diag
 INTEGER(I4B) :: j
 do j=1,size(mat,1)
  get_diag(j)=mat(j,j)
 end do
 END FUNCTION get_diag

 ! *************************************************************************************************
 ! (A2) put a diagonal vector in a matrix
 ! *************************************************************************************************
 SUBROUTINE put_diag(diagv,mat,err,message)
 REAL(DP), DIMENSION(:), INTENT(IN) :: diagv
 REAL(DP), DIMENSION(:,:), INTENT(INOUT) :: mat
 INTEGER(I4B) :: j
 integer(i4b),intent(out)      :: err                    ! error code
 character(*),intent(out)      :: message                ! error message
 err=0; message="f-put_diag/"
 if(size(diagv)/=min(size(mat,1),size(mat,2)))then
  err=10; message=trim(message)//'/sizeMismatch'; return
 endif
 do j=1,size(diagv)
  mat(j,j)=diagv(j)
 end do
 END SUBROUTINE put_diag
 
 ! *************************************************************************************************
 ! (A3) construct a unit (identity) matrix
 ! *************************************************************************************************
 SUBROUTINE unit_matrix(mat)
 REAL(DP), DIMENSION(:,:), INTENT(OUT) :: mat
 INTEGER(I4B) :: i,n
 n=min(size(mat,1),size(mat,2))
 mat(:,:)=0.0_sp
 do i=1,n
  mat(i,i)=1.0_sp
 end do
 END SUBROUTINE unit_matrix

 ! *************************************************************************************************
 ! (A4) compute length of vector in L2 norm
 ! *************************************************************************************************
 FUNCTION vabs(v)
 REAL(DP), DIMENSION(:), INTENT(IN) :: v
 REAL(DP) :: vabs
 vabs=sqrt(dot_product(v,v))
 END FUNCTION vabs

 ! *************************************************************************************************
 ! (A5) get a logical mask for the lower triangle
 ! *************************************************************************************************
 FUNCTION lower_triangle(j,k,extra)
 INTEGER(I4B), INTENT(IN) :: j,k
 INTEGER(I4B), OPTIONAL, INTENT(IN) :: extra
 LOGICAL(LGT), DIMENSION(j,k) :: lower_triangle
 INTEGER(I4B) :: n
 n=0
 if (present(extra)) n=extra
 lower_triangle=(outerdiff(arth_i(1,1,j),arth_i(1,1,k)) > -n)
 END FUNCTION lower_triangle

 ! *************************************************************************************************
 ! (A6) compute outer difference of two vectors
 ! *************************************************************************************************
 FUNCTION outerdiff_r(a,b)
 REAL(SP), DIMENSION(:), INTENT(IN) :: a,b
 REAL(SP), DIMENSION(size(a),size(b)) :: outerdiff_r
 outerdiff_r = spread(a,dim=2,ncopies=size(b)) - &
  spread(b,dim=1,ncopies=size(a))
 END FUNCTION outerdiff_r
 ! -------------------------------------------------------------------------------------------------
 FUNCTION outerdiff_d(a,b)
 REAL(DP), DIMENSION(:), INTENT(IN) :: a,b
 REAL(DP), DIMENSION(size(a),size(b)) :: outerdiff_d
 outerdiff_d = spread(a,dim=2,ncopies=size(b)) - &
  spread(b,dim=1,ncopies=size(a))
 END FUNCTION outerdiff_d
 ! -------------------------------------------------------------------------------------------------
 FUNCTION outerdiff_i(a,b)
 INTEGER(I4B), DIMENSION(:), INTENT(IN) :: a,b
 INTEGER(I4B), DIMENSION(size(a),size(b)) :: outerdiff_i
 outerdiff_i = spread(a,dim=2,ncopies=size(b)) - &
  spread(b,dim=1,ncopies=size(a))
 END FUNCTION outerdiff_i
 ! *************************************************************************************************
 ! (A7) compute outer product of two vectors
 ! *************************************************************************************************
 FUNCTION outerprod(a,b)
 REAL(DP), DIMENSION(:), INTENT(IN) :: a,b
 REAL(DP), DIMENSION(size(a),size(b)) :: outerprod
 outerprod = spread(a,dim=2,ncopies=size(b)) * &
  spread(b,dim=1,ncopies=size(a))
 END FUNCTION outerprod

 ! *************************************************************************************************
 ! (A8) find the location of the first .true. element in a vector, returned as an integer
 ! *************************************************************************************************
 FUNCTION ifirstloc(mask)
 LOGICAL(LGT), DIMENSION(:), INTENT(IN) :: mask
 INTEGER(I4B) :: ifirstloc
 INTEGER(I4B), DIMENSION(1) :: loc
 loc=maxloc(merge(1,0,mask))
 ifirstloc=loc(1)
 if (.not. mask(ifirstloc)) ifirstloc=size(mask)+1
 END FUNCTION ifirstloc


 ! *************************************************************************************************
 ! *************************************************************************************************
 ! (B) the error function
 ! *************************************************************************************************
 ! *************************************************************************************************

 ! *************************************************************************************************
 ! (B1) the error function
 ! *************************************************************************************************
 FUNCTION erf(x)
 ! Purpose: returns the error function erf(x)
 IMPLICIT NONE
 REAL(DP), INTENT(IN) :: x
 REAL(DP) :: erf
 erf=gammp(0.5_dp,x**2)
 if (x < 0.0) erf=-erf
 END FUNCTION erf

 ! *************************************************************************************************
 ! (B2) the incomplete Gamma function
 ! *************************************************************************************************
 FUNCTION gammp(a,x)
 ! Purpose: jacket for incomplete gamma function P(a,x) evaluation
 ! Selects from two algorithms:
 ! a. series representation,scheme g_ser;
 ! b. continued fraction representation,scheme gcf;
 ! Returns large negative number for both illegal input and
 ! non-converged approximations
 IMPLICIT NONE
 REAL(DP), INTENT(IN) :: a,x
 REAL(DP) :: gammp
 if(a<=0._dp.or.x<0._dp)then ! early return
  gammp=-huge(kind(dp)); return
 elseif(x<a+1.0_dp) then     ! use the series representations
  gammp=gser(a,x) 
 else                        ! use continued fraction    
  gammp=1.0_dp-gcf(a,x)      ! and take its complement
 end if
 END FUNCTION gammp

 ! *************************************************************************************************
 ! (B3) returns the incomplete Gamma function P(a,x) evaluated by its series representation
 ! *************************************************************************************************
 FUNCTION gser(a,x)
 IMPLICIT NONE
 REAL(DP), INTENT(IN) :: a,x
 REAL(DP) :: gser
 INTEGER(I4B), PARAMETER :: ITMAX=100
 REAL(DP), PARAMETER :: EPS=epsilon(x)
 INTEGER(I4B) :: n
 REAL(DP) :: ap,del,summ,lng_a
 if (x == 0.0) then
  gser=0.0
  RETURN
 end if
 ap=a
 summ=1.0_dp/a
 del=summ
 do n=1,ITMAX
  ap=ap+1.0_dp
  del=del*x/ap
  summ=summ+del
  if (abs(del) < abs(summ)*EPS) exit
 end do
 if(n<=ITMAX)then
  lng_a = gammln(a)
  if(lng_a<0.99_dp*(-huge(kind(dp))))then
   gser=-huge(kind(dp))
  endif
  gser=summ*exp(-x+a*log(x)-lng_a)
 else
  gser=-huge(kind(dp))
 endif
 END FUNCTION gser

 ! *************************************************************************************************
 ! (B4) returns the incomplete Gamma function P(a,x) evaluated by its continued fraction representation
 ! *************************************************************************************************
 FUNCTION gcf(a,x)
 IMPLICIT NONE
 REAL(DP), INTENT(IN) :: a,x
 REAL(DP) :: gcf
 INTEGER(I4B), PARAMETER :: ITMAX=100
 REAL(DP), PARAMETER :: EPS=epsilon(x),FPMIN=tiny(x)/EPS
 INTEGER(I4B) :: i
 REAL(DP) :: an,b,c,d,del,h,lng_a
 if (x == 0.0) then
  gcf=1.0
  RETURN
 end if
 b=x+1.0_dp-a
 c=1.0_dp/FPMIN
 d=1.0_dp/b
 h=d
 do i=1,ITMAX
  an=-i*(i-a)
  b=b+2.0_dp
  d=an*d+b
  if (abs(d) < FPMIN) d=FPMIN
  c=b+an/c
  if (abs(c) < FPMIN) c=FPMIN
  d=1.0_dp/d
  del=d*c
  h=h*del
  if (abs(del-1.0_dp) <= EPS) exit
 end do
 if(i<=ITMAX)then
  lng_a = gammln(a)
  if(lng_a<0.99_dp*(-huge(kind(dp))))then
   gcf=-huge(kind(dp))
  endif
  gcf=exp(-x+a*log(x)-lng_a)*h
 else
  gcf=-huge(kind(dp))
 endif
 END FUNCTION gcf

 ! *************************************************************************************************
 ! (B5) evaluates the natural logarithm of the Gamma function
 ! *************************************************************************************************
 FUNCTION gammln(xx)
 IMPLICIT NONE
 REAL(DP), INTENT(IN) :: xx
 REAL(DP) :: gammln
 REAL(DP) :: tmp,x
 REAL(DP) :: stp = 2.5066282746310005_dp
 REAL(DP), DIMENSION(6) :: coef = (/76.18009172947146_dp,&
     -86.50532032941677_dp,24.01409824083091_dp,&
     -1.231739572450155_dp,0.1208650973866179e-2_dp,&
     -0.5395239384953e-5_dp/)
 if(xx<=0._dp)then
  gammln=-huge(kind(dp)); return
 endif
 x=xx
 tmp=x+5.5_dp
 tmp=(x+0.5_dp)*log(tmp)-tmp
 gammln=tmp+log(stp*(1.000000000190015_dp+&
     sum(coef(:)/arth(x+1.0_dp,1.0_dp,size(coef))))/x)
 END FUNCTION gammln

 
 ! *************************************************************************************************
 ! *************************************************************************************************
 ! (C) build a vector of regularly spaced numbers
 ! *************************************************************************************************
 ! *************************************************************************************************

 ! *************************************************************************************************
 ! (C1) the arth function, used to build a vector of regularly spaced numbers 
 ! *************************************************************************************************
 FUNCTION arth_r(first,increment,n)
 implicit none
 REAL(SP), INTENT(IN) :: first,increment
 INTEGER(I4B), INTENT(IN) :: n
 REAL(SP), DIMENSION(n) :: arth_r
 INTEGER(I4B) :: k,k2
 REAL(SP) :: temp
 if (n > 0) arth_r(1)=first
 if (n <= NPAR_ARTH) then
  do k=2,n
   arth_r(k)=arth_r(k-1)+increment
  end do
 else
  do k=2,NPAR2_ARTH
   arth_r(k)=arth_r(k-1)+increment
  end do
  temp=increment*NPAR2_ARTH
  k=NPAR2_ARTH
  do
   if (k >= n) exit
   k2=k+k
   arth_r(k+1:min(k2,n))=temp+arth_r(1:min(k,n-k))
   temp=temp+temp
   k=k2
  end do
 end if
 END FUNCTION arth_r
 ! ------------------------------------------------------------------------------------------------
 FUNCTION arth_d(first,increment,n)
 implicit none
 REAL(DP), INTENT(IN) :: first,increment
 INTEGER(I4B), INTENT(IN) :: n
 REAL(DP), DIMENSION(n) :: arth_d
 INTEGER(I4B) :: k,k2
 REAL(DP) :: temp
 if (n > 0) arth_d(1)=first
 if (n <= NPAR_ARTH) then
  do k=2,n
   arth_d(k)=arth_d(k-1)+increment
  end do
 else
  do k=2,NPAR2_ARTH
   arth_d(k)=arth_d(k-1)+increment
  end do
  temp=increment*NPAR2_ARTH
  k=NPAR2_ARTH
  do
   if (k >= n) exit
   k2=k+k
   arth_d(k+1:min(k2,n))=temp+arth_d(1:min(k,n-k))
   temp=temp+temp
   k=k2
  end do
 end if
 END FUNCTION arth_d
 ! ------------------------------------------------------------------------------------------------
 FUNCTION arth_i(first,increment,n)
 implicit none
 INTEGER(I4B), INTENT(IN) :: first,increment,n
 INTEGER(I4B), DIMENSION(n) :: arth_i
 INTEGER(I4B) :: k,k2,temp
 if (n > 0) arth_i(1)=first
 if (n <= NPAR_ARTH) then
  do k=2,n
   arth_i(k)=arth_i(k-1)+increment
  end do
 else
  do k=2,NPAR2_ARTH
   arth_i(k)=arth_i(k-1)+increment
  end do
  temp=increment*NPAR2_ARTH
  k=NPAR2_ARTH
  do
   if (k >= n) exit
   k2=k+k
   arth_i(k+1:min(k2,n))=temp+arth_i(1:min(k,n-k))
   temp=temp+temp
   k=k2
  end do
 end if
 END FUNCTION arth_i

end module nr_utility_module

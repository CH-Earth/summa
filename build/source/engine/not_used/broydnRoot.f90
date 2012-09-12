module broydnRoot_module
USE nrtype
implicit none
private
public::broydn
contains

 ! *************************************************************
 ! new subroutine: broydn
 ! *************************************************************
 SUBROUTINE broydn(x,check,err,message)
 USE nr_utility_module, ONLY : get_diag,put_diag,unit_matrix,vabs,&
                               lower_triangle,outerprod
 USE data_struc,only:xboundLower,xboundUpper ! apply state constraints
 USE fminln
 IMPLICIT NONE
 REAL(DP), DIMENSION(:), INTENT(INOUT) :: x
 LOGICAL(LGT), INTENT(OUT) :: check
 INTEGER(I4B), PARAMETER :: MAXITS=200
 REAL(DP), PARAMETER :: EPS=epsilon(x),TOLF=1.0e-4_dp,TOLMIN=1.0e-6_dp,&
     TOLX=EPS,STPMX=100.0_dp
 INTEGER(I4B) :: i,its,k,n
 REAL(DP) :: f,fold,stpmax
 REAL(DP), DIMENSION(size(x)), TARGET :: fvec
 REAL(DP), DIMENSION(size(x)) :: c,d,fvcold,g,p,s,t,w,xold
 REAL(DP), DIMENSION(size(x),size(x)) :: qt,r
 LOGICAL :: restrt,sing
 integer(i4b),intent(out)      :: err
 character(*),intent(out)      :: message
 character(len=128)            :: cmessage ! error for downwind routine
 err=0; message='broydn/'
 ! MPC change: ensure x is in bounds
 where(x < xboundLower) x=xboundLower
 where(x > xboundUpper) x=xboundUpper
 ! end MPC change
 fmin_fvecp=>fvec
 n=size(x)
 f=fmin(x)
 if (maxval(abs(fvec(:))) < 0.01_dp*TOLF) then
  check=.false.
  RETURN
 end if
 stpmax=STPMX*max(vabs(x(:)),real(n,dp))
 restrt=.true.
 do its=1,MAXITS
  if (restrt) then
   call fDiffJacbn(x,fvec,r,err,cmessage)
   if(err/=0)then; message=trim(message)//trim(cmessage); return; endif
   call qrdcmp(r,c,d,sing,err,cmessage)
   if(err/=0)then; message=trim(message)//trim(cmessage); return; endif
   if(sing)then; err=10; message=trim(message)//'singular Jacobian in broydn'; return; endif
   call unit_matrix(qt)
   do k=1,n-1
    if (c(k) /= 0.0) then
     qt(k:n,:)=qt(k:n,:)-outerprod(r(k:n,k),&
      matmul(r(k:n,k),qt(k:n,:)))/c(k)
    end if
   end do
   where (lower_triangle(n,n)) r(:,:)=0.0
   call put_diag(d(:),r(:,:),err,cmessage)
   if(err/=0)then; message=trim(message)//trim(cmessage); return; endif
  else
   s(:)=x(:)-xold(:)
   do i=1,n
    t(i)=dot_product(r(i,i:n),s(i:n))
   end do
   w(:)=fvec(:)-fvcold(:)-matmul(t(:),qt(:,:))
   where (abs(w(:)) < EPS*(abs(fvec(:))+abs(fvcold(:)))) &
    w(:)=0.0
   if (any(w(:) /= 0.0)) then
    t(:)=matmul(qt(:,:),w(:))
    s(:)=s(:)/dot_product(s,s)
    call qrupdt(r,qt,t,s,err,cmessage)
    if(err/=0)then; message=trim(message)//trim(cmessage); return; endif
    d(:)=get_diag(r(:,:))
    if (any(d(:) == 0.0)) then
     err=10; message=trim(message)//'r singular in broydn'; return
    endif
   end if
  end if
  p(:)=-matmul(qt(:,:),fvec(:))
  do i=1,n
   g(i)=-dot_product(r(1:i,i),p(1:i))
  end do
  xold(:)=x(:)
  fvcold(:)=fvec(:)
  fold=f
  call rsolv(r,d,p,err,cmessage)
  if(err/=0)then; message=trim(message)//trim(cmessage); return; endif
  call lnsrch(xold,fold,g,p,x,f,stpmax,check,fmin,err,cmessage)
  if(err/=0)then; message=trim(message)//trim(cmessage); return; endif
  if (maxval(abs(fvec(:))) < TOLF) then
   check=.false.
   RETURN
  end if
  if (check) then
   if (restrt .or. maxval(abs(g(:))*max(abs(x(:)), &
    1.0_dp)/max(f,0.5_dp*n)) < TOLMIN) RETURN
   restrt=.true.
  else
   restrt=.false.
   if (maxval((abs(x(:)-xold(:)))/max(abs(x(:)), &
    1.0_dp)) < TOLX) RETURN
  end if
 end do
 err=-50; message=trim(message)//'MAXITS exceeded in broydn'
 END SUBROUTINE broydn

 ! ************************************************************************************************
 ! new subroutine: compute finite-difference approximation of the Jacobian matrix
 ! ************************************************************************************************
 SUBROUTINE fDiffJacbn(xtry,fvec,jmat,err,message)
 USE funcvector_module,only:funcv
 IMPLICIT NONE
 REAL(DP), DIMENSION(:), INTENT(IN)    :: fvec
 REAL(DP), DIMENSION(:), INTENT(INOUT) :: xtry
 REAL(DP), DIMENSION(:,:), INTENT(OUT) :: jmat
 REAL(DP), PARAMETER :: EPS=1.0e-4_dp
 INTEGER(I4B) :: j,n
 REAL(DP), DIMENSION(size(xtry)) :: ftest,xsav,xph,h
 integer(i4b),intent(out)      :: err                    ! error code
 character(*),intent(out)      :: message                ! error message
 err=0; message="fDiffJacbn/"
 ! check that the size of the vectors match
 if ( all((/size(xtry),size(fvec),size(jmat,1),size(jmat,2)/) == size(xtry)) ) then
  n=size(xtry)
 else
  err=20; message=trim(message)//"sizeMismatch"; return
 endif
 ! start procedure
 xsav=xtry
 h=EPS*abs(xsav)
 where (h == 0.0) h=EPS
 xph=xsav+h
 h=xph-xsav
 do j=1,n
  xtry(j)=xph(j)
  ftest = funcv(xtry)
  jmat(:,j)=(ftest(:)-fvec(:))/h(j)
  xtry(j)=xsav(j)
 end do
 end subroutine fDiffJacbn

 ! *************************************************************
 ! new subroutine: QR decomposition
 ! *************************************************************
 SUBROUTINE qrdcmp(a,c,d,sing,err,message)
 USE nr_utility_module, ONLY : outerprod,vabs
 IMPLICIT NONE
 REAL(DP), DIMENSION(:,:), INTENT(INOUT) :: a
 REAL(DP), DIMENSION(:), INTENT(OUT) :: c,d
 LOGICAL(LGT), INTENT(OUT) :: sing
 INTEGER(I4B) :: k,n
 REAL(DP) :: scale,sigma
 integer(i4b),intent(out)      :: err                    ! error code
 character(*),intent(out)      :: message                ! error message
 err=0; message="qrdcmp/"
 ! check that the size of the vectors match
 if ( all((/size(a,1),size(a,2),size(c),size(d)/) == size(a,1)) ) then
  n=size(a,1)
 else
  err=20; message=trim(message)//"sizeMismatch"; return
 endif
 ! start procedure
 sing=.false.
 do k=1,n-1
  scale=maxval(abs(a(k:n,k)))
  if (scale == 0.0) then
   sing=.true.
   c(k)=0.0
   d(k)=0.0
  else
   a(k:n,k)=a(k:n,k)/scale
   sigma=sign(vabs(a(k:n,k)),a(k,k))
   a(k,k)=a(k,k)+sigma
   c(k)=sigma*a(k,k)
   d(k)=-scale*sigma
   a(k:n,k+1:n)=a(k:n,k+1:n)-outerprod(a(k:n,k),&
    matmul(a(k:n,k),a(k:n,k+1:n)))/c(k)
  end if
 end do
 d(n)=a(n,n)
 if (d(n) == 0.0) sing=.true.
 END SUBROUTINE qrdcmp

 ! *************************************************************
 ! new subroutine: QR update
 ! *************************************************************
 SUBROUTINE qrupdt(r,qt,u,v,err,message)
 USE nr_utility_module, ONLY : ifirstloc
 IMPLICIT NONE
 REAL(DP), DIMENSION(:,:), INTENT(INOUT) :: r,qt
 REAL(DP), DIMENSION(:), INTENT(INOUT) :: u
 REAL(DP), DIMENSION(:), INTENT(IN) :: v
 INTEGER(I4B) :: i,k,n
 integer(i4b),intent(out)      :: err                    ! error code
 character(*),intent(out)      :: message                ! error message
 character(len=128)            :: cmessage               ! error message of downwind routine
 err=0; message="qrupdt/"
 ! check that the size of the vectors match
 if ( all((/size(r,1),size(r,2),size(qt,1),size(qt,2),size(u),size(v)/) == size(r,1)) ) then
  n=size(r,1)
 else
  err=20; message=trim(message)//"sizeMismatch"; return
 endif
 ! start procedure
 k=n+1-ifirstloc(u(n:1:-1) /= 0.0)
 if (k < 1) k=1
 do i=k-1,1,-1
  call rotate(r,qt,i,u(i),-u(i+1),err,cmessage)
  if(err/=0)then; message=trim(message)//trim(cmessage); return; endif
  u(i)=pythag(u(i),u(i+1))
 end do
 r(1,:)=r(1,:)+u(1)*v
 do i=1,k-1
  call rotate(r,qt,i,r(i,i),-r(i+1,i),err,cmessage)
  if(err/=0)then; message=trim(message)//trim(cmessage); return; endif
 end do
 END SUBROUTINE qrupdt

 ! *************************************************************
 ! new subroutine: rotate
 ! *************************************************************
 SUBROUTINE rotate(r,qt,i,a,b,err,message)
 IMPLICIT NONE
 REAL(DP), DIMENSION(:,:), TARGET, INTENT(INOUT) :: r,qt
 INTEGER(I4B), INTENT(IN) :: i
 REAL(DP), INTENT(IN) :: a,b
 REAL(DP), DIMENSION(size(r,1)) :: temp
 INTEGER(I4B) :: n
 REAL(DP) :: c,fact,s
 integer(i4b),intent(out)      :: err                    ! error code
 character(*),intent(out)      :: message                ! error message
 err=0; message="rotate/"
 ! check that the size of the vectors match
 if ( all((/size(r,1),size(r,2),size(qt,1),size(qt,2)/) == size(r,1)) ) then
  n=size(r,1)
 else
  err=20; message=trim(message)//"sizeMismatch"; return
 endif
 ! start procedure
 if (a == 0.0) then
  c=0.0
  s=sign(1.0_dp,b)
 else if (abs(a) > abs(b)) then
  fact=b/a
  c=sign(1.0_dp/sqrt(1.0_dp+fact**2),a)
  s=fact*c
 else
  fact=a/b
  s=sign(1.0_dp/sqrt(1.0_dp+fact**2),b)
  c=fact*s
 end if
 temp(i:n)=r(i,i:n)
 r(i,i:n)=c*temp(i:n)-s*r(i+1,i:n)
 r(i+1,i:n)=s*temp(i:n)+c*r(i+1,i:n)
 temp=qt(i,:)
 qt(i,:)=c*temp-s*qt(i+1,:)
 qt(i+1,:)=s*temp+c*qt(i+1,:)
 END SUBROUTINE rotate

 ! *************************************************************
 ! new subroutine: pythag
 ! *************************************************************
 FUNCTION pythag(a,b)
 IMPLICIT NONE
 REAL(DP), INTENT(IN) :: a,b
 REAL(DP) :: pythag
 REAL(DP) :: absa,absb
 absa=abs(a)
 absb=abs(b)
 if (absa > absb) then
  pythag=absa*sqrt(1.0_dp+(absb/absa)**2)
 else
  if (absb == 0.0) then
   pythag=0.0
  else
   pythag=absb*sqrt(1.0_dp+(absa/absb)**2)
  end if
 end if
 END FUNCTION pythag

 ! *************************************************************
 ! new subroutine: rsolv
 ! *************************************************************
 SUBROUTINE rsolv(a,d,b,err,message)
 IMPLICIT NONE
 REAL(DP), DIMENSION(:,:), INTENT(IN) :: a
 REAL(DP), DIMENSION(:), INTENT(IN) :: d
 REAL(DP), DIMENSION(:), INTENT(INOUT) :: b
 INTEGER(I4B) :: i,n
 integer(i4b),intent(out)      :: err                    ! error code
 character(*),intent(out)      :: message                ! error message
 err=0; message="f-rsolv/"
 ! check that the size of the vectors match
 if ( all((/size(a,1),size(a,2),size(b),size(d)/) == size(a,1)) ) then
  n=size(a,1)
 else
  err=20; message=trim(message)//"sizeMismatch"; return
 endif
 ! start procedure
 b(n)=b(n)/d(n)
 do i=n-1,1,-1
  b(i)=(b(i)-dot_product(a(i,i+1:n),b(i+1:n)))/d(i)
 end do
 END SUBROUTINE rsolv

 ! *************************************************************
 ! new subroutine: lnsrch
 ! *************************************************************
 SUBROUTINE lnsrch(xold,fold,g,p,x,f,stpmax,check,func,err,message)
 USE nr_utility_module, ONLY : vabs
 USE data_struc,only:xboundLower,xboundUpper ! apply state constraints
 IMPLICIT NONE
 REAL(DP), DIMENSION(:), INTENT(IN) :: xold,g
 REAL(DP), DIMENSION(:), INTENT(INOUT) :: p
 REAL(DP), INTENT(IN) :: fold,stpmax
 REAL(DP), DIMENSION(:), INTENT(OUT) :: x
 REAL(DP), INTENT(OUT) :: f
 LOGICAL(LGT), INTENT(OUT) :: check
 INTERFACE
  FUNCTION func(x)
  USE nrtype
  IMPLICIT NONE
  REAL(DP) :: func
  REAL(DP), DIMENSION(:), INTENT(IN) :: x
  END FUNCTION func
 END INTERFACE
 REAL(DP), PARAMETER :: ALF=1.0e-4_dp,TOLX=epsilon(x)
 INTEGER(I4B) :: ndum
 REAL(DP) :: a,alam,alam2,alamin,b,disc,f2,fold2,pabs,rhs1,rhs2,slope,&
  tmplam
 LOGICAL(LGT) :: first  ! used to check if in the first loop
 integer(i4b),intent(out)      :: err                    ! error code
 character(*),intent(out)      :: message                ! error message
 err=0; message="lnsrch/"
 ! check that the size of the vectors match
 if ( all((/size(g),size(p),size(x),size(xold)/) == size(g)) ) then
  ndum=size(g)
 else
  err=20; message=trim(message)//"sizeMismatch"; return
 endif
 ! start procedure
 check=.false.
 first=.true.
 pabs=vabs(p(:))
 if (pabs > stpmax) p(:)=p(:)*stpmax/pabs
 slope=dot_product(g,p)
 alamin=TOLX/maxval(abs(p(:))/max(abs(xold(:)),1.0_dp))
 alam=1.0
 do
  x(:)=xold(:)+alam*p(:)
  ! MPC change: ensure x is in bounds
  do
   if(count(x < xboundLower)>0 .or. count(x > xboundUpper)>0)then
    alam=0.1_dp*alam
    if (alam < alamin) then
     x(:)=xold(:)
     check=.true.
     RETURN
    else
     x(:)=xold(:)+alam*p(:)
    endif
   else
    exit
   endif
  end do
  ! end MPC change 
  f=func(x)
  if (alam < alamin) then
   x(:)=xold(:)
   check=.true.
   RETURN
  else if (f <= fold+ALF*alam*slope) then
   RETURN
  else
   if (first) then
    tmplam=-slope/(2.0_dp*(f-fold-slope))
   else
    rhs1=f-fold-alam*slope
    rhs2=f2-fold2-alam2*slope
    a=(rhs1/alam**2-rhs2/alam2**2)/(alam-alam2)
    b=(-alam2*rhs1/alam**2+alam*rhs2/alam2**2)/&
     (alam-alam2)
    if (a == 0.0) then
     tmplam=-slope/(2.0_dp*b)
    else
     disc=b*b-3.0_dp*a*slope
     if (disc < 0.0) then
      err=30; message=trim(message)//'roundoff problem in lnsrch'; return
     endif
     tmplam=(-b+sqrt(disc))/(3.0_dp*a)
    end if
    if (tmplam > 0.5_dp*alam) tmplam=0.5_dp*alam
   end if
  end if
  alam2=alam
  f2=f
  fold2=fold
  alam=max(tmplam,0.1_dp*alam)
  first=.false.
 end do
 END SUBROUTINE lnsrch

end module broydnRoot_module

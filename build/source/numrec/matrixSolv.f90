module matrixSolv_module
USE nrtype
implicit none
private
public::matrixSolv
contains

 subroutine matrixSolv(jMat,rVec,xInc,err,message)
 ! used to solve the matrix xInc*Jmat-1 = -rVec
 implicit none
 ! input
 real(dp),dimension(:,:),intent(in)        :: jMat     ! Jacobian matrix
 real(dp),dimension(:),intent(in)          :: rVec     ! residual vector
 ! output
 real(dp),dimension(:),intent(out)         :: xInc     ! iteration increment
 integer(i4b),intent(out)                  :: err      ! error code
 character(*),intent(out)                  :: message  ! error message
 ! local variables
 character(len=256)                        :: cmessage ! error message of downstream routine
 real(dp),dimension(size(rVec),size(rVec)) :: a        ! copy of jMat
 integer(i4b),dimension(size(rVec))        :: indx     ! records row permutation affected by partial pivoting
 real(dp)                                  :: d        ! inticates if number row interchanges is even/odd
 ! initialize error control
 err=0; message='matrixSolv/'
 ! get a copy of jMat (a is destroyed in ludcmp)
 a = jMat
 ! decompose the matrix
 call ludcmp(a,indx,d,err,cmessage)
 if(err/=0)then; message=trim(message)//trim(cmessage); return; endif
 ! solve the equations
 xInc = rVec ! set Xinc to the right-hand-side vector (which is over-written on output)
 call lubksb(a,indx,xInc,err,cmessage)
 if(err/=0)then; message=trim(message)//trim(cmessage); return; endif
 end subroutine matrixSolv




 ! ==========================================================================
 ! * PRIVATE SUBROUTINES ****************************************************
 ! ==========================================================================

 SUBROUTINE ludcmp(a,indx,d,err,message)
 IMPLICIT NONE
 ! input/output
 REAL(DP), DIMENSION(:,:), INTENT(INOUT) :: a
 INTEGER(I4B), DIMENSION(:), INTENT(OUT) :: indx
 REAL(DP), INTENT(OUT) :: d
 integer(i4b),intent(out)           :: err
 character(*),intent(out)           :: message
 ! local variables
 REAL(DP), DIMENSION(size(a,1)) :: vv
 REAL(DP), PARAMETER :: TINY=1.0e-20_dp
 INTEGER(I4B) :: j,n,imax
 ! initialize error control
 err=0; message='ludcmp/'
 n = size(indx)
 if(size(a,1) /= n .or. size(a,2) /= n)then
  message=trim(message)//'mismatch in size of matrices'
  err=20; return
 endif
 d=1.0_dp
 vv=maxval(abs(a),dim=2)
 if (any(vv == 0.0))then
  message=trim(message)//'singular matrix'
  err=20; return
 endif
 vv=1.0_dp/vv
 do j=1,n
  imax=(j-1)+imaxloc(vv(j:n)*abs(a(j:n,j)))
  if (j /= imax) then
   call swap(a(imax,:),a(j,:))
   d=-d
   vv(imax)=vv(j)
  end if
  indx(j)=imax
  if (a(j,j) == 0.0) a(j,j)=TINY
  a(j+1:n,j)=a(j+1:n,j)/a(j,j)
  a(j+1:n,j+1:n)=a(j+1:n,j+1:n)-outerprod(a(j+1:n,j),a(j,j+1:n))
 end do
 END SUBROUTINE ludcmp

 SUBROUTINE lubksb(a,indx,b,err,message)
 IMPLICIT NONE
 ! input/output
 REAL(DP), DIMENSION(:,:), INTENT(IN) :: a
 INTEGER(I4B), DIMENSION(:), INTENT(IN) :: indx
 REAL(DP), DIMENSION(:), INTENT(INOUT) :: b
 integer(i4b),intent(out)           :: err
 character(*),intent(out)           :: message
 ! local
 INTEGER(I4B) :: i,n,ii,ll
 REAL(DP) :: summ
 ! initialize error control
 err=0; message='lubksb/'
 n = size(indx)
 if(size(a,1) /= n .or. size(a,2) /= n)then
  message=trim(message)//'mismatch in size of matrices'
  err=20; return
 endif
 ii=0
 do i=1,n
  ll=indx(i)
  summ=b(ll)
  b(ll)=b(i)
  if (ii /= 0) then
   summ=summ-dot_product(a(i,ii:i-1),b(ii:i-1))
  else if (summ /= 0.0) then
   ii=i
  end if
  b(i)=summ
 end do
 do i=n,1,-1
  b(i) = (b(i)-dot_product(a(i,i+1:n),b(i+1:n)))/a(i,i)
 end do
 END SUBROUTINE lubksb

 SUBROUTINE swap(a,b)
 REAL(DP), DIMENSION(:), INTENT(INOUT) :: a,b
 REAL(DP), DIMENSION(SIZE(a)) :: dum
 dum=a
 a=b
 b=dum
 END SUBROUTINE swap

 FUNCTION outerprod(a,b)
 REAL(DP), DIMENSION(:), INTENT(IN) :: a,b
 REAL(DP), DIMENSION(size(a),size(b)) :: outerprod
 outerprod = spread(a,dim=2,ncopies=size(b)) * &
  spread(b,dim=1,ncopies=size(a))
 END FUNCTION outerprod

 FUNCTION imaxloc(arr)
 REAL(DP), DIMENSION(:), INTENT(IN) :: arr
 INTEGER(I4B) :: imaxloc
 INTEGER(I4B), DIMENSION(1) :: imax
 imax=maxloc(arr(:))
 imaxloc=imax(1)
 END FUNCTION imaxloc

end module matrixSolv_module

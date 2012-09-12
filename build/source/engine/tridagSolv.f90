module tridagSolv_module
USE nrtype
implicit none
private
public::tridag
contains

 ! *************************************************************
 ! new subroutine: tridag
 ! *************************************************************
 SUBROUTINE tridag(a,b,c,r,u,err,message)
 ! solve the tridiagonal system of equations
 USE nrtype
 IMPLICIT NONE
 REAL(DP), DIMENSION(:), INTENT(IN) :: a,b,c,r
 REAL(DP), DIMENSION(:), INTENT(OUT) :: u
 REAL(DP), DIMENSION(size(b)) :: gam
 INTEGER(I4B) :: n,j
 REAL(DP) :: bet
 integer(i4b),intent(out)      :: err                    ! error code
 character(*),intent(out)      :: message                ! error message
 err=0; message="tridag/"
 ! check that the size of the vectors match
 if ( all((/size(a)+1,size(b),size(c)+1/) == size(r)) ) then
  n=size(r)
 else
  err=20; message=trim(message)//"sizeMismatch"; return
 endif
 ! start procedure
 bet=b(1)
 if (bet == 0.0) then; err=30; message=trim(message)//'Error at code stage 1'; return; endif
 u(1)=r(1)/bet
 do j=2,n
  gam(j)=c(j-1)/bet
  bet=b(j)-a(j-1)*gam(j)
  if (bet == 0.0) then; err=30; message=trim(message)//'Error at code stage 2'; return; endif
  u(j)=(r(j)-a(j-1)*u(j-1))/bet
 end do
 do j=n-1,1,-1
  u(j)=u(j)-gam(j+1)*u(j+1)
 end do
 END SUBROUTINE tridag

end module tridagSolv_module

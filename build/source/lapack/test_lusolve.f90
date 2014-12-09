program test_lusolve
! Use modules
USE nrtype
USE luSolv_module,only:ludcmp,lubksb
! used to test the lusolve functions in lapack
implicit none
integer(i4b), parameter  :: n=3      ! size of the matrix
real(dp),dimension(n,n)  :: aa,a     ! the input array
real(dp),dimension(1,n)  :: bb,bm    ! the rhs of a linear sustem ax=b (matrix)
real(dp),dimension(n)    :: b        ! the rhs of a linear sustem ax=b (vector)
real(dp)                 :: d        ! indicates if number row interchanges is even/odd
integer(i4b)             :: indx(n)  ! records row permutation affected by partial pivoting
integer(i4b)             :: err      ! error code
character(len=256)       :: cmessage ! error message of downstream routine

! define the A matrix
aa(1,1:3) = (/ 2.0,  1.0,  1.0/)
aa(2,1:3) = (/ 4.0, -6.0,  0.0/)
aa(3,1:3) = (/-2.0,  7.0,  2.0/)

! define the rhs vector
bb(1,1:3) = (/ 3.0, -8.0, 10.0/)

! set a and b (overwritten in ludcmp and lubksb)
a = aa
b = bb(1,1:3)

! decompose the matrix
call ludcmp(a,indx,d,err,cmessage)
if(err/=0)then; print*, trim(cmessage); stop; endif

! solve the equations
call lubksb(a,indx,b,err,cmessage)
if(err/=0)then; print*, trim(cmessage); stop; endif
write(*,'(a,1x,3(f9.3,1x),a)') 'numrec: b = ', b, '; should be (1,2,-1)'

! set a and b (overwritten in lapack)
a  = aa
bm = bb

! use lapack
call dgesv(n,1,a,n,indx,bm,n,err); b = bm(1,1:3)
write(*,'(a,1x,3(f9.3,1x),a)') 'lapack: b = ', b, '; should be (1,2,-1)'

end program test_lusolve

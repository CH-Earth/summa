! SUMMA - Structure for Unifying Multiple Modeling Alternatives
! Copyright (C) 2014-2020 NCAR/RAL; University of Saskatchewan; University of Washington
!
! This file is part of SUMMA
!
! For more information see: http://www.ral.ucar.edu/projects/summa
!
! This program is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! This program is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License
! along with this program.  If not, see <http://www.gnu.org/licenses/>.

module scratch_module
! Sizes the per-thread scratch arrays of the solver and flux routines.
! Callers declare them allocatable, save and !$omp threadprivate, and call sizeScratch on entry;
! contents are undefined on entry. Sizes change between HRUs and state subsets, so a real array that
! changes size parks its buffer in a per-thread pool and takes back a parked buffer of the new size.
! real(qp) arrays use the real(rkind) procedures, which assumes qp is the same kind as rkind (nr_type).
USE nr_type
implicit none
private
public :: sizeScratch
! pool of parked buffers, bounded by the number of values it holds
integer(i4b),parameter :: maxParked=65536             ! most values held in the pool per thread (512 KB)
integer(i4b),save      :: nParked=0                   ! values held in the pool
! parked rank-1 buffers, one stack per lower bound (0 or 1) and upper bound
integer(i4b),parameter :: maxPool1=256                ! largest upper bound kept in the pool
type :: buffer1
 real(rkind),allocatable :: a(:)
end type buffer1
type :: stack1
 type(buffer1),allocatable :: b(:)
 integer(i4b)              :: n=0                     ! number of parked buffers
end type stack1
type(stack1),save :: pool1(0:1,0:maxPool1)
! parked rank-2 buffers, any shape
integer(i4b),parameter :: nPool2=16                   ! number of rank-2 buffers kept in the pool
type :: buffer2
 real(rkind),allocatable :: a(:,:)
end type buffer2
type(buffer2),save :: pool2(nPool2)
!$omp threadprivate(nParked, pool1, pool2)
interface sizeScratch
 module procedure sizeScratch1
 module procedure sizeScratch2
 module procedure sizeScratchL
 module procedure sizeScratchI
end interface sizeScratch
contains

 ! give each rank-1 array present the bounds lb:ub
 subroutine sizeScratch1(lb,ub,a1,a2,a3,a4,a5,a6,a7,a8,a9,a10,a11,a12,a13,a14,a15,a16,a17,a18,a19,a20)
 integer(i4b),intent(in)                                :: lb,ub
 real(rkind),allocatable,intent(inout)                  :: a1(:)
 real(rkind),allocatable,intent(inout),optional         :: a2(:),a3(:),a4(:),a5(:),a6(:),a7(:),a8(:),a9(:),a10(:),a11(:)
 real(rkind),allocatable,intent(inout),optional         :: a12(:),a13(:),a14(:),a15(:),a16(:),a17(:),a18(:),a19(:),a20(:)
 call fit(a1)
 if(present(a2)) call fit(a2);   if(present(a3)) call fit(a3);   if(present(a4)) call fit(a4);   if(present(a5)) call fit(a5)
 if(present(a6)) call fit(a6);   if(present(a7)) call fit(a7);   if(present(a8)) call fit(a8);   if(present(a9)) call fit(a9)
 if(present(a10)) call fit(a10); if(present(a11)) call fit(a11); if(present(a12)) call fit(a12); if(present(a13)) call fit(a13)
 if(present(a14)) call fit(a14); if(present(a15)) call fit(a15); if(present(a16)) call fit(a16); if(present(a17)) call fit(a17)
 if(present(a18)) call fit(a18); if(present(a19)) call fit(a19); if(present(a20)) call fit(a20)
 contains
  subroutine fit(a)
  real(rkind),allocatable,intent(inout) :: a(:)
  if(allocated(a))then
   if(lbound(a,1)==lb .and. ubound(a,1)==ub) return
   call park1(a)
  end if
  if(pooled1(lb,ub))then
   if(pool1(lb,ub)%n>0)then
    call move_alloc(pool1(lb,ub)%b(pool1(lb,ub)%n)%a, a); pool1(lb,ub)%n = pool1(lb,ub)%n - 1
    nParked = nParked - size(a); return
   end if
  end if
  allocate(a(lb:ub))
  end subroutine fit
 end subroutine sizeScratch1

 ! true if rank-1 buffers with bounds lb:ub are kept in the pool
 logical(lgt) function pooled1(lb,ub)
 integer(i4b),intent(in) :: lb,ub
 pooled1 = (lb==0 .or. lb==1) .and. ub>=0 .and. ub<=maxPool1
 end function pooled1

 ! move a rank-1 buffer onto the stack for its bounds, or free it if those bounds are not pooled or the pool is full
 subroutine park1(a)
 real(rkind),allocatable,intent(inout) :: a(:)
 type(buffer1),allocatable             :: grown(:)
 integer(i4b)                          :: lb,ub,i
 lb = lbound(a,1); ub = ubound(a,1)
 if(.not.pooled1(lb,ub) .or. nParked+size(a)>maxParked)then; deallocate(a); return; end if
 associate(st => pool1(lb,ub))
  if(.not.allocated(st%b)) allocate(st%b(8))
  if(st%n==size(st%b))then
   allocate(grown(2*st%n))
   do i=1,st%n; call move_alloc(st%b(i)%a, grown(i)%a); end do
   call move_alloc(grown, st%b)
  end if
  st%n = st%n + 1; nParked = nParked + size(a)
  call move_alloc(a, st%b(st%n)%a)
 end associate
 end subroutine park1

 ! give each rank-2 array present the shape (n1,n2)
 subroutine sizeScratch2(n1,n2,a1,a2,a3,a4)
 integer(i4b),intent(in)                        :: n1,n2
 real(rkind),allocatable,intent(inout)          :: a1(:,:)
 real(rkind),allocatable,intent(inout),optional :: a2(:,:),a3(:,:),a4(:,:)
 call fit(a1)
 if(present(a2)) call fit(a2); if(present(a3)) call fit(a3); if(present(a4)) call fit(a4)
 contains
  subroutine fit(a)
  real(rkind),allocatable,intent(inout) :: a(:,:)
  integer(i4b)                          :: i
  if(allocated(a))then
   if(size(a,1)==n1 .and. size(a,2)==n2) return
   if(nParked+size(a)<=maxParked)then
    do i=1,nPool2
     if(.not.allocated(pool2(i)%a))then; nParked = nParked + size(a); call move_alloc(a, pool2(i)%a); exit; end if
    end do
   end if
   if(allocated(a)) deallocate(a)
  end if
  do i=1,nPool2
   if(.not.allocated(pool2(i)%a)) cycle
   if(size(pool2(i)%a,1)==n1 .and. size(pool2(i)%a,2)==n2)then
    call move_alloc(pool2(i)%a, a); nParked = nParked - size(a); return
   end if
  end do
  allocate(a(n1,n2))
  end subroutine fit
 end subroutine sizeScratch2

 ! give each logical array present the bounds lb:ub
 subroutine sizeScratchL(lb,ub,a1,a2)
 integer(i4b),intent(in)                         :: lb,ub
 logical(lgt),allocatable,intent(inout)          :: a1(:)
 logical(lgt),allocatable,intent(inout),optional :: a2(:)
 call fit(a1)
 if(present(a2)) call fit(a2)
 contains
  subroutine fit(a)
  logical(lgt),allocatable,intent(inout) :: a(:)
  if(allocated(a))then
   if(lbound(a,1)==lb .and. ubound(a,1)==ub) return
   deallocate(a)
  end if
  allocate(a(lb:ub))
  end subroutine fit
 end subroutine sizeScratchL
 ! give each integer array present the bounds lb:ub
 subroutine sizeScratchI(lb,ub,a1,a2)
 integer(i4b),intent(in)                         :: lb,ub
 integer(i4b),allocatable,intent(inout)          :: a1(:)
 integer(i4b),allocatable,intent(inout),optional :: a2(:)
 call fit(a1)
 if(present(a2)) call fit(a2)
 contains
  subroutine fit(a)
  integer(i4b),allocatable,intent(inout) :: a(:)
  if(allocated(a))then
   if(lbound(a,1)==lb .and. ubound(a,1)==ub) return
   deallocate(a)
  end if
  allocate(a(lb:ub))
  end subroutine fit
 end subroutine sizeScratchI

end module scratch_module

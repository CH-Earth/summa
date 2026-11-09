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
! an array is reallocated only when its bounds change, so its contents are undefined on entry.
USE nr_type
implicit none
private
public :: sizeScratch
interface sizeScratch
 module procedure sizeScratch1
 module procedure sizeScratch2
 module procedure sizeScratchL
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
   deallocate(a)
  end if
  allocate(a(lb:ub))
  end subroutine fit
 end subroutine sizeScratch1

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
  if(allocated(a))then
   if(size(a,1)==n1 .and. size(a,2)==n2) return
   deallocate(a)
  end if
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

end module scratch_module

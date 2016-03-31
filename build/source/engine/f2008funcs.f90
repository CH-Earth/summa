! SUMMA - Structure for Unifying Multiple Modeling Alternatives
! Copyright (C) 2014-2015 NCAR/RAL
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

module f2008funcs_module
USE nrtype
implicit none
private
public::cloneStruc

! define generic interface
interface cloneStruc
 module procedure cloneStruc_rv, cloneStruc_iv
end interface cloneStruc

contains

 ! ************************************************************************************************
 ! public subroutine cloneStruc_rv: clone a data structure (real vector)
 ! ************************************************************************************************
 subroutine cloneStruc_rv(dataVec,source,mold,err,message)
 implicit none
 ! input-output: data vector for allocation/population
 real(dp),intent(inout),allocatable     :: dataVec(:)     ! data vector
 ! optional input
 class(*),intent(in),optional           :: source(:)      ! dataVec = shape of source + elements of source
 class(*),intent(in),optional           :: mold(:)        ! dataVec = shape of mold
 ! error control
 integer(i4b),intent(out)               :: err            ! error code
 character(*),intent(out)               :: message        ! error message
 ! ----------------------------------------------------------------------------------------------------------------------------------
 ! local variables
 integer(i4b)                           :: lowerBound     ! lower bound of the data vector
 integer(i4b)                           :: upperBound     ! upper bound of the data vector
 character(LEN=256)                     :: cmessage       ! error message from downwind routine
 ! -----------------------------------------------------------------------------------------------------------------------------------
 ! initialize errors
 err=0; message="cloneStruc_rv/"

 ! check that source and mold are present
 if(.not.present(source) .and. .not.present(mold))then
  message=trim(message)//'expect to receive either optional argument "source" or "mold" (neither given)'
  err=20; return
 endif

 ! check that source and mold are not both present
 if(present(source) .and. present(mold))then
  message=trim(message)//'expect to receive either optional argument "source" or "mold" (both given)'
  err=20; return
 endif

 ! get the bounds of the source or the mold vector
 if(present(source)) call getVecBounds(source,lowerBound,upperBound,err,cmessage)
 if(present(mold))   call getVecBounds(mold,  lowerBound,upperBound,err,cmessage)
 if(err/=0)then; message=trim(message)//trim(cmessage); return; endif

 ! reallocate spcae
 if(allocated(dataVec)) deallocate(dataVec)
 allocate(dataVec(lowerBound:upperBound),stat=err)
 if(err/=0)then; err=20; message=trim(message)//'unable to allocate space for the data vector'; return; endif

 ! copy data
 if(present(source))then
  select type(source)
   type is(real(dp)); dataVec(lowerBound:upperBound) = source(lowerBound:upperBound)
   class default; err=20; message=trim(message)//'expect source to be of type real(dp)'; return
  end select
 endif  ! if source is present

 end subroutine cloneStruc_rv

 ! ************************************************************************************************
 ! public subroutine cloneStruc_iv: clone a data structure (integer vector)
 ! ************************************************************************************************
 subroutine cloneStruc_iv(dataVec,source,mold,err,message)
 implicit none
 ! input-output: data vector for allocation/population
 integer(i4b),intent(inout),allocatable :: dataVec(:)     ! data vector
 ! optional input
 class(*),intent(in),optional           :: source(:)      ! dataVec = shape of source + elements of source
 class(*),intent(in),optional           :: mold(:)        ! dataVec = shape of mold
 ! error control
 integer(i4b),intent(out)               :: err            ! error code
 character(*),intent(out)               :: message        ! error message
 ! -----------------------------------------------------------------------------------------------------------------------------------
 ! local variables
 integer(i4b)                           :: lowerBound     ! lower bound of the data vector
 integer(i4b)                           :: upperBound     ! upper bound of the data vector
 character(LEN=256)                     :: cmessage       ! error message from downwind routine
 ! -----------------------------------------------------------------------------------------------------------------------------------
 ! initialize errors
 err=0; message="cloneStruc_iv/"

 ! check that source and mold are present
 if(.not.present(source) .and. .not.present(mold))then
  message=trim(message)//'expect to receive either optional argument "source" or "mold" (neither given)'
  err=20; return
 endif

 ! check that source and mold are not both present
 if(present(source) .and. present(mold))then
  message=trim(message)//'expect to receive either optional argument "source" or "mold" (both given)'
  err=20; return
 endif

 ! get the bounds of the source or the mold vector
 if(present(source)) call getVecBounds(source,lowerBound,upperBound,err,cmessage)
 if(present(mold))   call getVecBounds(mold,  lowerBound,upperBound,err,cmessage)
 if(err/=0)then; message=trim(message)//trim(cmessage); return; endif

 ! reallocate spcae
 if(allocated(dataVec)) deallocate(dataVec)
 allocate(dataVec(lowerBound:upperBound),stat=err)
 if(err/=0)then; err=20; message=trim(message)//'unable to allocate space for the data vector'; return; endif

 ! copy data
 if(present(source))then
  select type(source)
   type is(integer(i4b)); dataVec(lowerBound:upperBound) = source(lowerBound:upperBound)
   class default; err=20; message=trim(message)//'expect source to be of type integer(i4b)'; return
  end select
 endif  ! if source is present

 end subroutine cloneStruc_iv

 ! ************************************************************************************************
 ! private subroutine getVecSize: get the size of a data vector
 ! ************************************************************************************************
 subroutine getVecBounds(dataVec,lowerBound,upperBound,err,message)
 implicit none
 ! dummy variables
 class(*),intent(in)       :: dataVec(:)     ! data vector
 integer(i4b),intent(out)  :: lowerBound     ! lower bound of the data vector
 integer(i4b),intent(out)  :: upperBound     ! upper bound of the data vector
 integer(i4b),intent(out)  :: err            ! error code
 character(*),intent(out)  :: message        ! error message
 ! local variables
 integer(i4b),dimension(1) :: lowerBoundVec  ! lower bound of the data vector
 integer(i4b),dimension(1) :: upperBoundVec  ! upper bound of the data vector
 ! initialize errors
 err=0; message="getVecBounds/"
 ! get the size of the data vector
 select type(dataVec)
  type is(real(dp)    ); lowerBoundVec= lBound(dataVec); upperBoundVec= uBound(dataVec)
  type is(integer(i4b)); lowerBoundVec= lBound(dataVec); upperBoundVec= uBound(dataVec)
  class default; err=20; message=trim(message)//'unable to identify data type'; return
 end select
 ! return bounds
 lowerBound = lowerBoundVec(1)
 upperBound = upperBoundVec(1)
 end subroutine getVecBounds

end module f2008funcs_module

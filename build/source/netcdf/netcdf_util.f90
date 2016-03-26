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

module netcdf_util_module
USE nrtype
USE netcdf
implicit none
private
public::nc_file_open
public::netcdf_err
contains


 ! *********************************************************************************************************
 ! public subroutine file_open: open file
 ! *********************************************************************************************************
 subroutine nc_file_open(infile,mode,ncid,err,message)
 implicit none
 ! declare dummy variables
 character(*),intent(in)              :: infile      ! filename
 integer(i4b),intent(in)              :: mode        ! file open mode
 integer(i4b),intent(out)             :: ncid        ! file unit
 integer(i4b),intent(out)             :: err         ! error code
 character(*),intent(out)             :: message     ! error message
 ! declare local variables
 logical(lgt)                         :: xist        ! .TRUE. if the file exists

 ! initialize errors
 err=0; message="netcdf-file_open/"

 ! check if the file exists
 inquire(file=trim(infile),exist=xist) ! Check for existence of file
 if(.not.xist)then
   message=trim(message)//"FileNotFound[file='"//trim(infile)//"']"
   err=10; return
 endif

 ! open file
 err=nf90_open(infile, mode, ncid) 
 if(err/=nf90_noerr) then
   message=trim(message)//"OpenError['"//trim(infile)//"']"//trim(nf90_strerror(err))
   err=20; return
 endif

 end subroutine nc_file_open


! ***********************************************************************************************
! check the status of netCDF file operation and return error message 
! ***********************************************************************************************
 subroutine netcdf_err(err,message)
  ! used to handle errors for NetCDF calls
  use netcdf
  implicit none
  ! declare dummies
  integer(i4b), intent(inout)   :: err
  character(*), intent(inout)   :: message
  ! start procedure here
  if (err/=nf90_noerr) then
   message=trim(message)//"["//trim(nf90_strerror(err))//"]"
   err=200
  endif
 end subroutine netcdf_err

end module netcdf_util_module

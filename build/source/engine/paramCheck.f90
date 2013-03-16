module paramCheck_module
USE nrtype
implicit none
private
public::paramCheck
contains

 ! ************************************************************************************************
 ! (1) new subroutine: check consistency of model parameters
 ! ************************************************************************************************
 subroutine paramCheck(err,message)
 ! FUSE data structures
 USE data_struc,only:mpar_data    ! data structures for model parameters
 USE var_lookup,only:iLookPARAM   ! named variables for elements of the data structures
 implicit none
 ! define output
 integer(i4b),intent(out)             :: err         ! error code
 character(*),intent(out)             :: message     ! error message
 ! Start procedure here
 err=0; message="paramCheck/"

 ! check that the snow layer bounds are OK
 if(mpar_data%var(iLookPARAM%zmin)/mpar_data%var(iLookPARAM%zmax) > 0.25_dp)then
  message=trim(message)//'zmax must be at least 4 times larger than zmin'
  err=20; return
 endif

 ! check that the maximum transpiration limit is within bounds
 if(mpar_data%var(iLookPARAM%critSoilTranspire)>mpar_data%var(iLookPARAM%theta_sat) .or. &
    mpar_data%var(iLookPARAM%critSoilTranspire)<mpar_data%var(iLookPARAM%theta_res))then
  message=trim(message)//'critSoilTranspire parameter is out of range '// &
                         '[NOTE: if overwriting Noah-MP soil table values in paramTrial, must overwrite all soil parameters]'
  err=20; return
 endif

 ! check that the soil wilting point is within bounds
 if(mpar_data%var(iLookPARAM%critSoilWilting)>mpar_data%var(iLookPARAM%theta_sat) .or. &
    mpar_data%var(iLookPARAM%critSoilWilting)<mpar_data%var(iLookPARAM%theta_res))then
  message=trim(message)//'critSoilWilting parameter is out of range '// &
                         '[NOTE: if overwriting Noah-MP soil table values in paramTrial, must overwrite all soil parameters]'
  err=20; return
 endif

 end subroutine paramCheck

end module paramCheck_module

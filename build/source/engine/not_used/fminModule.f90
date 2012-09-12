module fminln
USE nrtype
implicit none
real(dp),pointer :: fmin_fvecp(:)  ! points to vector of function values
contains

 ! ************************************************************************************************
 ! new function: compute the dot-product of the vector of function values
 ! ************************************************************************************************
 FUNCTION fmin(x)
 USE funcvector_module,only:funcv
 IMPLICIT NONE
 REAL(DP), DIMENSION(:), INTENT(IN) :: x
 REAL(DP) :: fmin
 fmin_fvecp=funcv(x)
 fmin=0.5_dp*dot_product(fmin_fvecp,fmin_fvecp)
 END FUNCTION fmin

end module fminln

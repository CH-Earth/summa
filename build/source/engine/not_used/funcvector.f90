module funcvector_module
USE nrtype
implicit none
private
public::funcv
contains

 ! ************************************************************************************************
 ! new function: compute vector of function evaluations for a given x vector
 ! ************************************************************************************************
 function funcv(x_try)
 ! function = (x_old + dXdt_old*dt/2 + dXdt_try*dt/2) - x_try
 USE data_struc,only:xboundLower,xboundUpper             ! apply constraints
 USE data_struc,only:xold_save,dXdt_xold,dtsub_save      ! provide access to x_old and dt
 USE cmput_dXdt_module,only:cmput_dXdt                   ! compute dX_dt at a given value of X
 implicit none
 ! dummy variables
 real(dp),dimension(:),intent(in)          :: x_try      ! trial x-vector
 real(dp),dimension(size(x_try))           :: funcv      ! returned function vector
 ! local variables
 real(dp),dimension(size(x_try))           :: dX_dt      ! temporal derivatives
 integer(i4b)                              :: err        ! error code
 character(len=256)                        :: cmessage   ! error message for downwind routine

 ! compute dX_dt
 !print*,'new function evaluation '
 !print*, 'x_try = ', x_try
 call cmput_dXdt(x_try,dX_dt,err,cmessage)
 if(err/=0)then; print*, 'FORTRAN STOP: in funcv, '//trim(cmessage); stop; endif

 ! compute function evaluation
 funcv = (xold_save + dXdt_xold*(dtsub_save/2._dp) + dX_dt*(dtsub_save/2._dp) ) - x_try

 end function funcv

end module funcvector_module

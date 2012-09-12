module odeIntegr8_module
USE nrtype
implicit none
private
public::odeIntegr8
contains

 ! ************************************************************************************************
 ! new subroutine: compute value of the model state vector at the end of the time step
 ! ************************************************************************************************
 subroutine odeIntegr8(dt,x_old,x_new,err,message)
 ! used to compute the value of the model state vector at the end of the time step
 ! (just a wrapper for different solution techniques)
 USE data_struc,only:model_decisions        ! model decision structure
 USE var_lookup,only:iLookDECISIONS         ! named variables for elements of the decision structure
 USE iEulerSolv_module,only:iEulerSolv      ! provide access to implicit Euler solver
 implicit none
 ! dummy variables
 real(dp),intent(in)                       :: dt         ! time step (seconds)
 real(dp),intent(in)                       :: x_old(:)   ! state vector at the start of the time step
 real(dp),intent(out)                      :: x_new(:)   ! state vector at the end of the time step
 integer(i4b),intent(out)                  :: err        ! error code
 character(*),intent(out)                  :: message    ! error message
 ! local variables
 character(LEN=256)                        :: cmessage   ! error message of downwind routine 
 ! initialize error control
 err=0; message="odeIntegr8/"

 ! "solve" using selected solution technique
 select case(trim(model_decisions(iLookDECISIONS%num_method)%decision))
  case('implicit')
   call iEulerSolv(dt,x_old,x_new,err,cmessage)
   if(err/=0)then; err=10; message=trim(message)//trim(cmessage); return; endif
  case default
   err=10; message=trim(message)//"case '"//trim(model_decisions(iLookDECISIONS%num_method)%decision)//"' not implemented yet"; return
 end select

 end subroutine odeIntegr8

end module odeIntegr8_module

module cmput_dXdt_module
USE nrtype
implicit none
private
public::cmput_dXdt
contains

 ! ************************************************************************************************
 ! new subroutine: compute dX_dt at a given value of X
 ! ************************************************************************************************
 subroutine cmput_dXdt(x_try,dX_dt,err,message)
 USE MapStr2vec_module,only:vector2str                   ! module to put the trial X vector into data structures
 USE ConvE2Temp_module,only:layer_Temp                   ! module to compute temperature based on enthalpy
 USE cmput_dvar_module,only:cmput_dvar                   ! module to compute diagnostic/ancillary variables
 USE energyflux_module,only:energyflux                   ! module to compute the temporal derivative in layer energy
 USE liquidflux_module,only:liquidflux                   ! module to compute the temporal derivative in liquid water
 USE build_dXdt_module,only:build_dXdt                   ! module to compute the temporal derivative in all model states
 ! used to compute dX_dt at a given value of X
 implicit none
 ! dummy variables
 real(dp),intent(in)                       :: x_try(:)   ! trial x-vector
 real(dp),intent(out)                      :: dX_dt(:)   ! derivative in X
 integer(i4b),intent(out)                  :: err        ! error code
 character(*),intent(out)                  :: message    ! error message
 ! local variables
 character(LEN=256)                        :: cmessage   ! error message of downwind routine 

 ! initialize error control
 err=0; message="cmput_dXdt/"

 ! put the trial X vector into data structures
 call vector2str(x_try,err,cmessage)
 if(err/=0)then; err=10; message=trim(message)//trim(cmessage); return; endif
 ! compute temperature based on enthalpy
 call layer_Temp(err,cmessage)
 if(err/=0)then; err=10; message=trim(message)//trim(cmessage); return; endif
 ! compute diagnostic/ancillary variables
 call cmput_dvar(err,cmessage)
 if(err/=0)then; err=10; message=trim(message)//trim(cmessage); return; endif


 ! compute the time change in layer energy
 call energyflux(err,cmessage)
 if(err/=0)then; err=10; message=trim(message)//trim(cmessage); return; endif

 ! compute the time change in liquid water
 call liquidflux(err,cmessage)
 if(err/=0)then; err=10; message=trim(message)//trim(cmessage); return; endif
 

 ! assemble model state equations
 call build_dXdt(dX_dt,err,cmessage)


 !print*,"x_try = ", x_try
 !print*,"dX_dt = ", dX_dt
 !stop 'FORTRAN STOP: computing dX/dt'

 end subroutine cmput_dXdt

end module cmput_dXdt_module

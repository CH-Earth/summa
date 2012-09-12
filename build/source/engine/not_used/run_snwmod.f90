module run_snwmod_module
USE nrtype
implicit none
private
public::run_snwmod
contains

 ! ************************************************************************************************
 ! new subroutine: run the snow model
 ! ************************************************************************************************
 subroutine run_snwmod(err,message)
 ! used to run the snow model for one timestep
 USE data_struc,only:forcFileInfo                         ! extract time step of forcing data
 USE data_struc,only:xold_save,dXdt_xold,dtsub_save       ! save initial state+deriv and time step in the data_struc module
 USE MapStr2vec_module,only:str2vector,vector2str         ! map a structure to a vector
 USE stateLimit_module,only:stateLimit                    ! identify reasonable bounds for model state variables
 USE odeIntegr8_module,only:odeIntegr8                    ! integrate the ODE forward in time
 implicit none
 ! define output
 integer(i4b),intent(out)             :: err              ! error code
 character(*),intent(out)             :: message          ! error message
 ! define local variables
 real(dp),pointer                     :: dt=>null()       ! length of time step (seconds)
 real(dp),pointer                     :: x_old(:)=>null() ! initial x-vector
 real(dp),pointer                     :: x_new(:)=>null() ! new x-vector
 integer(i4b)                         :: nX               ! number of elements in the X vector
 character(len=256)                   :: cmessage         ! error message from downwind routine

 ! ************************************************************************************************
 ! (1) initialize -- get x-vector and allocate space...
 ! ************************************************************************************************
 ! Start procedure here
 err=0; message="f-fuse/run_snwmod/"
 ! get the x-vector at the start of the time step
 call str2vector(x_old,err,cmessage)
 if(err/=0)then; message=trim(message)//trim(cmessage); return; endif
 ! get the limits of the x-vector, and store in the data structure module
 call stateLimit(err,cmessage)
 if(err/=0)then; message=trim(message)//trim(cmessage); return; endif
 ! get the size of the x-vector
 nX = size(x_old)
 ! allocate space for the x-vectors
 allocate(dt,dtsub_save,xold_save(nX),dXdt_xold(nX),x_new(nX),stat=err)
 if(err/=0)then;err=30;message=trim(message)//"problemAllocate"; return; endif
 ! get the length of the time step (seconds)
 dt = forcFileInfo%data_step

 ! **********************************************************************************************
 ! (2) run the model for one time step
 ! **********************************************************************************************
 call odeIntegr8(dt,x_old,x_new,err,cmessage)
 if(err/=0)then; message=trim(message)//trim(cmessage); return; endif
 !print*, 'x_old = ', x_old
 !print*, 'x_new = ', x_new
 !stop 'FORTRAN STOP: after odeIntegr8'

 ! **********************************************************************************************
 ! (3) put the new vector in the data structures
 ! **********************************************************************************************
 call vector2str(x_new,err,cmessage)
 if(err/=0)then; message=trim(message)//trim(cmessage); return; endif

 ! **********************************************************************************************
 ! deallocate space
 deallocate(dt,dtsub_save,x_old,xold_save,dXdt_xold,x_new,stat=err)
 if(err/=0)then;err=40;message=trim(message)//"problemDeallocate"; return; endif
 end subroutine run_snwmod

end module run_snwmod_module

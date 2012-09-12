module iEulerSolv_module
USE nrtype
implicit none
private
public::iEulerSolv
contains

 ! ************************************************************************************************
 ! new subroutine: compute value of the model state vector at the end of the time step
 ! ************************************************************************************************
 subroutine iEulerSolv(dt,x_old,x_new,err,message)
 ! use implicit Euler method to compute the value of the model state vector at the end of the time step
 USE data_struc,only:xold_save,dXdt_xold,dtsub_save      ! save initial state+deriv and time step in the data_struc module
 USE cmput_dXdt_module,only:cmput_dXdt                   ! compute dX_dt at a given value of X
 USE funcvector_module,only:funcv                        ! provide access to the function call
 USE broydnRoot_module,only:broydn                       ! provide access to the Broyden module
 implicit none
 ! dummy variables
 real(dp),intent(in)                       :: dt         ! time step (seconds)
 real(dp),intent(in)                       :: x_old(:)   ! state vector at the start of the time step
 real(dp),intent(out)                      :: x_new(:)   ! state vector at the end of the time step
 integer(i4b),intent(out)                  :: err        ! error code
 character(*),intent(out)                  :: message    ! error message
 ! local variables
 character(LEN=256)                        :: cmessage           ! error message of downwind routine 
 real(dp),dimension(size(x_old))           :: x0,dXdt0           ! x-vector and derivative at start of a sub-step
 real(dp),dimension(size(x_old))           :: x1,dXdt1           ! x-vector and derivative at end of a sub-step
 real(dp),parameter                        :: atol=10._dp       ! absolute tolerance
 real(dp),parameter                        :: safety=0.85_dp     ! safety factor
 real(dp),parameter                        :: eps=epsilon(x_old) ! avoid divide by zero
 real(dp),parameter                        :: r_min=0.1_dp       ! minimum multiplier
 real(dp),parameter                        :: r_max=4.0_dp       ! maximum multiplier
 real(dp)                                  :: dt_sub             ! length of time step (seconds)
 real(dp)                                  :: dtfrac             ! fraction of time step completed
 real(dp)                                  :: state_err(1)       ! maximum absoute error in model state vector
 
 logical(lgt)                              :: restrt             ! .true. if computing new sub-step
 logical(lgt)                              :: bcheck             ! checks for local minimum in Broyden's method
 ! initialize error control
 err=0; message="iEulerSolv/"

 ! initialize the length of the time step
 dt_sub = dt          ! length of time step (seconds)
 dtfrac = 0._dp       ! fraction of time step completed

 ! initialize x0 (state at the start of a sub-step)
 x0 = x_old ! (x_old is subroutine input)

 ! initialize the re-start flag
 restrt = .true.
 
 ! loop through time step
 do   ! (exit once completed time step)

  ! compute dX_dt at x0
  if(restrt)then
   call cmput_dXdt(x0,dXdt0,err,cmessage)
   if(err/=0)then; err=20; message=trim(message)//trim(cmessage); return; endif
  endif

  ! get trial x-vector (initialize with state vector at the start of the sub-step)
  x1 = x0

  ! save information in the data_struc module
  xold_save  = x0
  dXdt_xold  = dXdt0
  dtsub_save = dt_sub

  ! use Broyden's method to estimate x1
  call broydn(x1,bcheck,err,cmessage)
  if(err>0)then; err=20; message=trim(message)//trim(cmessage); return; endif 

  ! ***** check for non-convergence
  if(bcheck.or.err<0)then
   dt_sub = dt_sub*r_min
   restrt = .false.
   cycle
  ! ***** convergence
  else
   ! compute derivative at the end of the sub-step
   call cmput_dXdt(x1,dXdt1,err,cmessage)
   if(err/=0)then; err=20; message=trim(message)//trim(cmessage); return; endif 
   ! compute difference in successive derivatives
   state_err = maxval(abs( (dt_sub/2._dp)*(dXdt1 - dXdt0) ) )
  endif  ! (if converged)
  
  ! print progress
  print*, dt_sub,dt_sub+dtfrac,state_err

  ! ***** check if satisfied tolerance
  if(state_err(1) < atol)then
   ! increment fluxes
   ! (to-do)
   ! set end of sub-step to start of new sub-step
   x0 = x1
   ! increment fraction of sub-step completed (and exit if finished)
   dtfrac = dtfrac + dt_sub
   if(dtfrac >= dt) exit ! completed time step
   ! compute length of new sub-step
   dt_sub = dt_sub * min(safety*sqrt( atol / max(state_err(1),eps) ), r_max)
   ! truncate sub-step if greater than time step remaining
   if(dt_sub > dt-dtfrac) dt_sub = dt-dtfrac
   ! set flag to start of new sub-step
   restrt = .true.

  ! ***** if not satisfied tolerance
  else
   ! re-compute sub-step, but with a smaller time increment
   dt_sub = dt_sub * max(safety*sqrt( atol / max(state_err(1),eps) ), r_min)
   restrt = .false.
  endif  ! (if satisfied tolerance)

 end do  ! (looping through sub-steps)

 ! save output (x_new)
 x_new = x1

 !print*, 'iEulerSolv; x_old = ', x_old
 !print*, 'iEulerSolv; x_new = ', x_new
 !pause

 end subroutine iEulerSolv

end module iEulerSolv_module

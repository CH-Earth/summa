module surfAlbedo_module
USE nrtype
implicit none
private
public::surfAlbedo
contains

 ! ************************************************************************************************
 ! new subroutine: compute surface albedo
 ! ************************************************************************************************
 subroutine surfAlbedo(dt,&          ! input: time step (seconds)
                       err,message)  ! output: error control
 USE data_struc,only:mpar_data,mvar_data,model_decisions     ! data structures
 USE var_lookup,only:iLookPARAM,iLookMVAR,iLookDECISIONS     ! named variables for structure elements
 USE mDecisions_module,only:funcSnowAge,BATSlike             ! named variables for albedo options
 ! compute the surface albedo
 implicit none
 ! dummy variables
 real(dp),intent(in)           :: dt                         ! time step (seconds)
 integer(i4b),intent(out)      :: err                        ! error code
 character(*),intent(out)      :: message                    ! error message
 ! local pointers to model parameters
 real(dp),pointer              :: snw_crit                   ! critical mass necessary for albedo refreshment (kg m-2)
 real(dp),pointer              :: alb_fresh                  ! fresh snow albedo (-)
 real(dp),pointer              :: alb_dry                    ! minimum snow albedo during winter (-)
 real(dp),pointer              :: alb_wet                    ! minimum snow albedo during spring (-)
 real(dp),pointer              :: alb_decay                  ! temporal decay factor for snow albedo (s-1) 
 real(dp),pointer              :: alb_scale                  ! temporal albedo scaling factor (s)
 real(dp),pointer              :: soot_load                  ! albedo decay associated with soot load (-)
 ! local pointers to model forcing data
 real(dp),pointer              :: snowfall                   ! snowfall (kg m-2 s-1) 
 ! local pointers to model state variables
 real(dp),pointer              :: surfaceAlbedo              ! surface albedo (-) 
 real(dp),pointer              :: surfaceTemp                ! temperature of the top layer (K)
 real(dp),pointer              :: surfaceVolFracLiq          ! volumetric fraction of liquid water the top snow layer (-)
 ! local variables
 real(dp),parameter            :: liqWinter=0.001            ! volumetric fraction of liquid water that defines "winter" or "spring" conditions
 real(dp)                      :: alb_min                    ! minimum surface albedo (-)
 real(dp)                      :: dAlb_dt                    ! change in albedo with time (s-1)
 ! initialize error control
 err=0; message="surfAlbedo/"

 ! assign pointers to the surface albedo
 surfaceAlbedo => mvar_data%var(iLookMVAR%scalarAlbedo)%dat(1)

 ! assign pointers to the snowfall
 snowfall => mvar_data%var(iLookMVAR%scalarSnowfall)%dat(1)

 ! ***** compute increase in snow albedo
 if(snowfall > 0._dp)then
  ! assign pointers to albedo refreshment parameters
  snw_crit   => mpar_data%var(iLookPARAM%snw_crit)           ! critical mass necessary for albedo refreshment (kg m-2)  
  alb_fresh  => mpar_data%var(iLookPARAM%alb_fresh)          ! fresh snow albedo (-)
  ! compute increase in snow albedo with time (s-1)
  dAlb_dt = min(snowfall/snw_crit, (alb_fresh - surfaceAlbedo)/dt)

 ! ***** compute decrease in snow albedo
 else
  
  ! identify the albedo decay method
  select case(model_decisions(iLookDECISIONS%alb_method)%iDecision)

   ! method 1: albedo decay a function of snow age
   case(funcSnowAge)
    ! assign pointers to model parameters
    alb_dry   => mpar_data%var(iLookPARAM%alb_dry)           ! minimum snow albedo during winter (-)
    alb_wet   => mpar_data%var(iLookPARAM%alb_wet)           ! minimum snow albedo during spring (-)
    alb_decay => mpar_data%var(iLookPARAM%alb_decay)         ! temporal decay factor for snow albedo (s-1)
    ! assign pointers to the volumetric fraction of liquid water in the upper snow layer
    surfaceVolFracLiq => mvar_data%var(iLookMVAR%mLayerVolFracLiq)%dat(1)
    ! compute the minimum snow albedo (-)
    if(surfaceVolFracLiq < liqWinter)then; alb_min = alb_dry ! winter conditions
    else; alb_min = alb_wet ! spring conditions
    endif
    ! compute decrease in snow albedo with time (s-1)
    dAlb_dt = -alb_decay*(surfaceAlbedo - alb_min)

   ! method 2: albedo decay based on a BATS-like approach, with destructive metamorphism + soot content
   case(BATSlike)
    err=10; message=trim(message)//'"batslike" method not implemented yet'; return

   ! check for unknown albedo decay method
   case default
    err=10; message=trim(message)//"unknownOption"; return

  end select  ! (identifying the albedo decay method)

 endif ! (if snowing: switch between albedo increase or decay)

 ! update albedo (use explicit solution)
 surfaceAlbedo = surfaceAlbedo + dAlb_dt*dt
 
 end subroutine surfAlbedo

end module surfAlbedo_module

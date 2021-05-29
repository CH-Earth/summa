module computSnowDepth_module

! data types
USE nrtype

! physical constants
USE multiconst,only:&
                    Tfreeze,      & ! temperature at freezing              (K)
                    LH_fus,       & ! latent heat of fusion                (J kg-1)
                    LH_sub,       & ! latent heat of sublimation           (J kg-1)
                    iden_ice,     & ! intrinsic density of ice             (kg m-3)
                    iden_water      ! intrinsic density of liquid water    (kg m-3)

! data types
USE data_types,only:&
                    var_i,               & ! x%var(:)                (i4b)
                    var_d,               & ! x%var(:)                (dp)
                    var_ilength,         & ! x%var(:)%dat            (i4b)
                    var_dlength,         & ! x%var(:)%dat            (dp)
                    zLookup                ! x%z(:)%var(:)%lookup(:) (dp)

! named variables for parent structures
USE var_lookup,only:iLookDECISIONS         ! named variables for elements of the decision structure
USE var_lookup,only:iLookPROG              ! named variables for structure elements
USE var_lookup,only:iLookDIAG              ! named variables for structure elements
USE var_lookup,only:iLookFLUX              ! named variables for structure elements
USE var_lookup,only:iLookPARAM             ! named variables for structure elements
USE var_lookup,only:iLookINDEX             ! named variables for structure elements
USE globalData,only:iname_snow             ! named variables for snow
USE globalData,only:iname_soil             ! named variables for soil


! privacy
implicit none
private
public::computSnowDepth

real(dp),parameter     :: verySmall=1.e-6_dp   ! used as an additive constant to check if substantial difference among real numbers

contains

 ! ************************************************************************************************
 ! public subroutine coupled_em: run the coupled energy-mass model for one timestep
 ! ************************************************************************************************
 subroutine computSnowDepth(&
 							dt_sub,					&
 							nSnow,					& ! intent(in)
 							mLayerVolFracLiq,   	& ! intent(inout)
 							mLayerVolFracIce,		& ! intent(inout)
 							mLayerTemp,				& ! intent(in)
 							mpar_data,				& ! intent(in)
 							flux_data,				& ! intent(in)
 							diag_data,				& ! intent(in)
 					   		! output
 					   		mLayerDepth,			& ! intent(out)
 					   		scalarSnowDepth,		& ! intent(out)
 					   		scalarSWE,				& ! intetn(out)
                       		! error control
                       		err,message)         ! intent(out):   error control
                       		
 USE snwDensify_module,only:snwDensify      ! snow densification (compaction and cavitation)

 implicit none
  real(qp),intent(in)				   :: dt_sub
  integer(i4b),intent(in)              :: nSnow                  ! number of snow layers
  real(dp),intent(inout)			   :: mLayerVolFracLiq(:)
  real(dp),intent(inout)			   :: mLayerVolFracIce(:)
  real(dp),intent(in)				   :: mLayerTemp(:)
  type(var_dlength),intent(in)         :: mpar_data              ! model parameters
  type(var_dlength),intent(in)		   :: flux_data              ! model fluxes for a local HRU
  type(var_dlength),intent(in)         :: diag_data              ! diagnostic variables for a local HRU
  real(dp),intent(out)				   :: mLayerDepth(:)
  real(dp),intent(out)				   :: scalarSnowDepth
  real(dp),intent(out)				   :: scalarSWE
  integer(i4b),intent(out)             :: err                    ! error code
  character(*),intent(out)             :: message                ! error message
 
 ! local variables
 character(len=256)                   :: cmessage               ! error message
 integer(i4b)                         :: iSnow                  ! index of snow layers
 real(dp)                             :: massLiquid             ! mass liquid water (kg m-2)
 
 
 associate(&
   scalarSnowSublimation   => flux_data%var(iLookFLUX%scalarSnowSublimation)%dat(1)   & ! sublimation from the snow surface (kg m-2 s-1)
   )

  ! * compute change in ice content of the top snow layer due to sublimation...
  ! ---------------------------------------------------------------------------
  ! NOTE: this is done BEFORE densification
  if(nSnow > 0)then ! snow layers exist

   ! try to remove ice from the top layer
   iSnow=1

   ! save the mass of liquid water (kg m-2)
   massLiquid = mLayerDepth(iSnow)*mLayerVolFracLiq(iSnow)*iden_water

   ! add/remove the depth of snow gained/lost by frost/sublimation (m)
   ! NOTE: assume constant density
   mLayerDepth(iSnow) = mLayerDepth(iSnow) + dt_sub*scalarSnowSublimation/(mLayerVolFracIce(iSnow)*iden_ice)

   ! update the volumetric fraction of liquid water   
   mLayerVolFracLiq(iSnow) = massLiquid / (mLayerDepth(iSnow)*iden_water)

  ! no snow
  else

   ! no snow: check that sublimation is zero
   if(abs(scalarSnowSublimation) > verySmall)then
    message=trim(message)//'sublimation of snow has been computed when no snow exists'
    err=20; return
   end if

  end if  ! (if snow layers exist)

  ! *** account for compaction and cavitation in the snowpack...
  ! ------------------------------------------------------------
  if(nSnow>0)then
   call snwDensify(&
                   ! intent(in): variables
                   dt_sub,                                                  & ! intent(in): time step (s)
                   nSnow,                 									& ! intent(in): number of snow layers
                   mLayerTemp(1:nSnow),       								& ! intent(in): temperature of each layer (K)
                   diag_data%var(iLookDIAG%mLayerMeltFreeze)%dat(1:nSnow), 	& ! intent(in): volumetric melt in each layer (kg m-3)
                   ! intent(in): parameters
                   mpar_data%var(iLookPARAM%densScalGrowth)%dat(1),         & ! intent(in): density scaling factor for grain growth (kg-1 m3)
                   mpar_data%var(iLookPARAM%tempScalGrowth)%dat(1),         & ! intent(in): temperature scaling factor for grain growth (K-1)
                   mpar_data%var(iLookPARAM%grainGrowthRate)%dat(1),        & ! intent(in): rate of grain growth (s-1)
                   mpar_data%var(iLookPARAM%densScalOvrbdn)%dat(1),         & ! intent(in): density scaling factor for overburden pressure (kg-1 m3)
                   mpar_data%var(iLookPARAM%tempScalOvrbdn)%dat(1),         & ! intent(in): temperature scaling factor for overburden pressure (K-1)
                   mpar_data%var(iLookPARAM%baseViscosity)%dat(1),          & ! intent(in): viscosity coefficient at T=T_frz and snow density=0 (kg m-2 s)
                   ! intent(inout): state variables
                   mLayerDepth(1:nSnow),      								& ! intent(inout): depth of each layer (m)
                   mLayerVolFracLiq(1:nSnow), 								& ! intent(inout):  volumetric fraction of liquid water after itertations (-)
                   mLayerVolFracIce(1:nSnow), 								& ! intent(inout):  volumetric fraction of ice after itertations (-)
                   ! output: error control
                   err,cmessage)                     ! intent(out): error control
   if(err/=0)then; err=55; message=trim(message)//trim(cmessage); return; end if
  end if  ! if snow layers exist

  ! recompute snow depth and SWE
  if(nSnow > 0)then
   scalarSnowDepth = sum( mLayerDepth(1:nSnow) )
   scalarSWE       = sum( (mLayerVolFracLiq(1:nSnow)*iden_water + mLayerVolFracIce(1:nSnow)*iden_ice) * mLayerDepth(1:nSnow) )
  end if
  
  end associate
  
 end subroutine computSnowDepth
  
  
end module computSnowDepth_module


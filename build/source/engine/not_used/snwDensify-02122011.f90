module snwDensify_module
USE nrtype
implicit none
private
public::snwDensify
contains

 ! ************************************************************************************************
 ! new subroutine: compute change in temperature over the time step
 ! ************************************************************************************************
 subroutine snwDensify(dt,            & ! time step (seconds)
                       err,message)     ! error control
 USE multiconst,only:&
                     iden_ice,   &    ! intrinsic density of ice (kg m-3)
                     iden_water       ! intrinsic density of liquid water (kg m-3)
 USE data_struc,only:mpar_meta,forc_meta,mvar_meta,indx_meta                    ! metadata
 USE data_struc,only:mpar_data,forc_data,mvar_data,indx_data,ix_soil,ix_snow    ! data structures
 USE var_lookup,only:iLookPARAM,iLookFORCE,iLookMVAR,iLookINDEX                 ! named variables for structure elements
 ! compute change in temperature over the time step
 implicit none
 real(dp),intent(in)                 :: dt                       ! time step (seconds)
 integer(i4b),intent(out)            :: err                      ! error code
 character(*),intent(out)            :: message                  ! error message
 ! local pointers to model parameters (new snow density)
 real(dp),pointer                    :: scal_phi                 ! density scaling factor for destructive metamorphism (kg m-3)
 real(dp),pointer                    :: metam_m2                 ! temporal scaling factor for destructive metamorphism (s-1)
 ! local pointers to model state variables (surface layer only)
 real(dp),pointer                    :: mLayerDepth(:)           ! depth of each layer (m)
 real(dp),pointer                    :: mLayerVolFracIce(:)      ! volumetric fraction of ice in each layer (-)
 real(dp),pointer                    :: mLayerVolFracLiq(:)      ! volumetric fraction of liquid water in each layer (-)
 ! local pointers to diagnostic variables
 real(dp),pointer                    :: mLayerMeltFreeze(:)      ! melt/freeze in each layer (kg m-3)
 ! local pointers to model index variables
 integer(i4b),pointer                :: nLayers                  ! number of layers
 integer(i4b),pointer                :: layerType(:)             ! type of the layer (ix_soil or ix_snow)
 ! define local variables
 character(LEN=256)                  :: cmessage                 ! error message of downwind routine
 real(dp),parameter                  :: snwden_min=100._dp       ! minimum snow density for reducing metamorphism rate (kg m-3)
 real(dp),parameter                  :: dt_toler=0.1_dp          ! fraction of compaction allowed in a time step (-)
 integer(i4b)                        :: iSnow                    ! index of snow layers
 integer(i4b)                        :: nSnow                    ! number of snow layers
 integer(i4b)                        :: nSoil                    ! number of soil layers
 real(dp)                            :: phi                      ! multiplier in the densification algorithm (-)
 real(dp)                            :: CR_metamorph             ! compaction rate for metamorphism (s-1)
 real(dp)                            :: CR_snowmelt              ! compaction rate for snowmelt (s-1)
 real(dp)                            :: CR_total                 ! total compaction rate (s-1)
 real(dp)                            :: massIceOld               ! mass of ice in the snow layer (kg m-2)
 real(dp)                            :: massLiqOld               ! mass of liquid water in the snow layer (kg m-2)
 real(dp)                            :: tempDepth                ! temporary layer depth (m)
 real(dp)                            :: tempvar                  ! temporary variable
 real(dp)                            :: delDepthMetamorph        ! change in depth from metamorphism
 real(dp)                            :: delDepthSnowmelt         ! change in depth due to snowmelt
 real(dp)                            :: volFracMelt              ! volumetric fraction of snowmelt
 real(dp)                            :: CR_snowmelt1 
 ! initialize error control
 err=0; message="snwDensify/"
 ! assign pointers to model parameters (new snow density)
 scal_phi         => mpar_data%var(iLookPARAM%scal_phi)                ! density scaling factor for destructive metamorphism (kg m-3)
 metam_m2         => mpar_data%var(iLookPARAM%metam_m2)                ! temporal scaling factor for destructive metamorphism (s-1)
 ! assign pointers to model state variables -- all surface layer only
 mLayerDepth      => mvar_data%var(iLookMVAR%mLayerDepth)%dat          ! layer depth (m)
 mLayerVolFracIce => mvar_data%var(iLookMVAR%mLayerVolFracIce)%dat     ! volumetric fraction of ice in each layer  (-)
 mLayerVolFracLiq => mvar_data%var(iLookMVAR%mLayerVolFracLiq)%dat     ! volumetric fraction of liquid water in each layer  (-)
 ! assign pointers to diagnostic variables
 mLayerMeltFreeze => mvar_data%var(iLookMVAR%mLayerMeltFreeze)%dat     ! melt/freeze in each layer (kg m-3)
 ! assign local pointers to the model index structures
 nLayers          => indx_data%var(iLookINDEX%nLayers)%dat(1)          ! number of layers
 layerType        => indx_data%var(iLookINDEX%layerType)%dat           ! layer type (ix_soil or ix_snow)

 ! identify the number of snow and soil layers
 nSnow = count(layerType==ix_snow)
 nSoil = count(layerType==ix_soil)

 ! loop through snow layers
 do iSnow=1,nSnow
  ! save mass of liquid water and ice (mass does not change)
  if(iSnow==nSnow) print*, 'depth old = ', mLayerDepth(iSnow)
  if(iSnow==nSnow) print*, 'VolFracIce old = ', mLayerVolFracIce(iSnow)
  massIceOld = iden_ice*mLayerVolFracIce(iSnow)*mLayerDepth(iSnow)   ! (kg m-2)
  massLiqOld = iden_water*mLayerVolFracLiq(iSnow)*mLayerDepth(iSnow) ! (kg m-2)
  if(iSnow==nSnow) print*, ' massIceOld = ', massIceOld
  ! compute the rate of compaction associated with metamoprphism (s-1)
  if(mLayerVolFracIce(iSnow) <snwden_min) phi=1._dp
  if(mLayerVolFracIce(iSnow)>=snwden_min) phi=exp(-(mLayerVolFracIce(iSnow)*iden_ice - snwden_min)/scal_phi)
  CR_metamorph = phi*metam_m2  ! (should be much less than dt)
  ! compute the compaction rate associated with snow melt (s-1)
  




  if(iSnow==nSnow) print*, 'mLayerMeltFreeze(iSnow)/iden_ice = ', mLayerMeltFreeze(iSnow)/iden_ice 
  CR_snowmelt  = (max(0._dp,-mLayerMeltFreeze(iSnow))/iden_ice)/dt  ! (should be much less than dt)
  ! compute the total compaction rate
  CR_total = CR_metamorph + CR_snowmelt
  ! check that the compaction rate is not too large
  if(CR_total > dt_toler*dt)then; err=20; message=trim(message)//'total compaction rate is too large'; return; endif
  ! test
  if(iSnow==nSnow)then
   tempvar = -CR_snowmelt*dt/(mLayerVolFracIce(iSnow) + CR_snowmelt*dt)
   tempDepth = mLayerDepth(iSnow) + tempvar*mLayerDepth(iSnow)
   print*, 'ratio = ', tempvar, mLayerVolFracIce(iSnow) + CR_snowmelt*dt
   print*, 'tempDepth = ', tempDepth
   print*, 'test VolFracIce new = ', massIceOld/(tempDepth*iden_ice)
  endif
  ! compute the change in depth due to metamorphism
  delDepthMetamorph  = mLayerDepth(iSnow)/(1._dp + CR_metamorph*dt) - mLayerDepth(iSnow)
  ! compute the change in depth associated with snowmelt
  volFracMelt = min(mLayerMeltFreeze(iSnow),0._dp)/iden_ice
  CR_snowmelt1 = (1._dp/dt) * volFracMelt/(mLayerVolFracIce(iSnow) - volFracMelt)
  delDepthSnowmelt = dt*CR_snowmelt1*mLayerDepth(iSnow)
  if(iSnow==nSnow)then
   print*, 'mLayerDepth(iSnow) + delDepthSnowmelt = ', mLayerDepth(iSnow) + delDepthSnowmelt
   tempDepth = mLayerDepth(iSnow)/(1._dp + CR_metamorph*dt) + dt*CR_snowmelt1*mLayerDepth(iSnow)
   print*, 'final check...', massIceOld/(tempDepth*iden_ice)
  endif
  mLayerDepth(iSnow)      = mLayerDepth(iSnow) / (1._dp + (CR_metamorph + CR_snowmelt)*dt)
  mLayerVolFracIce(iSnow) = massIceOld/(mLayerDepth(iSnow)*iden_ice)
  mLayerVolFracLiq(iSnow) = massLiqOld/(mLayerDepth(iSnow)*iden_water)
  if(iSnow==nSnow)then
   print*, 'massIceNew = ', iden_ice*mLayerVolFracIce(iSnow)*mLayerDepth(iSnow)
   print*, 'VolFracIce new = ', mLayerVolFracIce(iSnow)
   print*, 'depth new = ', mLayerDepth(iSnow)
   print*, 'CR_metamorph, CR_snowmelt = ', CR_metamorph, CR_snowmelt
  endif
 end do  ! looping through snow layers

 end subroutine snwDensify

end module snwDensify_module

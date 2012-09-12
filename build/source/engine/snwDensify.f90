module snwDensify_module
USE nrtype
implicit none
private
public::snwDensify
contains

 ! ************************************************************************************************
 ! new subroutine: compute change in snow density over the time step
 ! ************************************************************************************************
 subroutine snwDensify(dt,                  &   ! input:  time step (seconds)
                       mLayerVolFracLiqIter,&   ! input:  volumetric fraction of liquid water after itertations (-)
                       mLayerVolFracIceIter,&   ! input:  volumetric fraction of ice after itertations (-)
                       mLayerVolFracLiqNew, &   ! output: volumetric fraction of liquid water after densification (-)
                       mLayerVolFracIceNew, &   ! output: volumetric fraction of ice after densification (-)
                       err,message)             ! output: error control
 USE multiconst,only:&
                     Tfreeze,    &     ! freezing point of pure water (K)
                     iden_ice,   &     ! intrinsic density of ice (kg m-3)
                     iden_water        ! intrinsic density of liquid water (kg m-3)
 USE var_derive_module,only:calcHeight ! module to calculate height at layer interfaces and layer mid-point
 USE data_struc,only:mpar_meta,forc_meta,mvar_meta,indx_meta                    ! metadata
 USE data_struc,only:mpar_data,forc_data,mvar_data,indx_data,ix_soil,ix_snow    ! data structures
 USE var_lookup,only:iLookPARAM,iLookFORCE,iLookMVAR,iLookINDEX                 ! named variables for structure elements
 ! compute change in snow density over the time step
 implicit none
 real(dp),intent(in)                 :: dt                       ! time step (seconds)
 real(dp),intent(in)                 :: mLayerVolFracLiqIter(:)  ! volumetric fraction of liquid water after iterations (-)
 real(dp),intent(in)                 :: mLayerVolFracIceIter(:)  ! volumetric fraction of ice after iterations (-)
 real(dp),intent(out)                :: mLayerVolFracLiqNew(:)   ! new volumetric fraction of liquid water (-)
 real(dp),intent(out)                :: mLayerVolFracIceNew(:)   ! new volumetric fraction of ice (-)
 integer(i4b),intent(out)            :: err                      ! error code
 character(*),intent(out)            :: message                  ! error message
 ! local pointers to model parameters
 real(dp),pointer                    :: densScalGrowth           ! density scaling factor for grain growth (kg-1 m3)
 real(dp),pointer                    :: tempScalGrowth           ! temperature scaling factor for grain growth (K-1)
 real(dp),pointer                    :: grainGrowthRate          ! rate of grain growth (s-1)
 real(dp),pointer                    :: densScalOvrbdn           ! density scaling factor for overburden pressure (kg-1 m3)
 real(dp),pointer                    :: tempScalOvrbdn           ! temperature scaling factor for overburden pressure (K-1)
 real(dp),pointer                    :: base_visc                ! viscosity coefficient at T=T_frz and snow density=0 (kg m-2 s)
 ! local pointers to model state variables
 real(dp),pointer                    :: mLayerDepth(:)           ! depth of each layer (m)
 real(dp),pointer                    :: mLayerTemp(:)            ! temperature of each layer (m)
 real(dp),pointer                    :: mLayerVolFracIce(:)      ! volumetric fraction of ice in each layer (-)
 real(dp),pointer                    :: mLayerVolFracLiq(:)      ! volumetric fraction of liquid water in each layer (-)
 ! local pointers to model index variables
 integer(i4b),pointer                :: nLayers                  ! number of layers
 integer(i4b),pointer                :: layerType(:)             ! type of the layer (ix_soil or ix_snow)
 ! define local variables
 character(len=256)                  :: cmessage                 ! error message for downwind routine
 real(dp),parameter                  :: dt_toler=0.1_dp          ! fraction of compaction allowed in a time step (-)
 integer(i4b)                        :: iSnow                    ! index of snow layers
 integer(i4b)                        :: nSnow                    ! number of snow layers
 integer(i4b)                        :: nSoil                    ! number of soil layers
 real(dp)                            :: chi1,chi2,chi3,chi4,chi5 ! multipliers in the densification algorithm (-)
 real(dp)                            :: fracMelt                 ! fraction of ice that melts in a time step (-)
 real(dp)                            :: halfWeight               ! half of the weight of the current snow layer (kg m-2)
 real(dp)                            :: weightSnow               ! total weight of snow above the current snow layer (kg m-2)
 real(dp)                            :: CR_grainGrowth           ! compaction rate for grain growth (s-1)
 real(dp)                            :: CR_ovrvdnPress           ! compaction rate associated with over-burden pressure (s-1)
 real(dp)                            :: CR_metamorph             ! compaction rate for metamorphism (s-1)
 real(dp)                            :: CR_snowmelt              ! compaction rate for snowmelt (s-1)
 real(dp)                            :: massIceOld               ! mass of ice in the snow layer (kg m-2)
 real(dp)                            :: massLiqOld               ! mass of liquid water in the snow layer (kg m-2)
 real(dp),parameter                  :: snwden_min=100._dp       ! minimum snow density for reducing metamorphism rate (kg m-3)
 real(dp),parameter                  :: snwDensityMax=550._dp    ! maximum snow density for collapse under melt (kg m-3)
 real(dp),parameter                  :: wetSnowThresh=0.01_dp    ! threshold to discriminate between "wet" and "dry" snow

 ! initialize error control
 err=0; message="snwDensify/"
 ! assign pointers to model parameters (new snow density)
 densScalGrowth   => mpar_data%var(iLookPARAM%densScalGrowth)          ! density scaling factor for grain growth (kg-1 m3)
 tempScalGrowth   => mpar_data%var(iLookPARAM%tempScalGrowth)          ! temperature scaling factor for grain growth (K-1)
 grainGrowthRate  => mpar_data%var(iLookPARAM%grainGrowthRate)         ! rate of grain growth (s-1)
 densScalOvrbdn   => mpar_data%var(iLookPARAM%densScalOvrbdn)          ! density scaling factor for overburden pressure (kg-1 m3)
 tempScalOvrbdn   => mpar_data%var(iLookPARAM%tempScalOvrbdn)          ! temperature scaling factor for overburden pressure (K-1)
 base_visc        => mpar_data%var(iLookPARAM%base_visc)               ! viscosity coefficient at T=T_frz and snow density=0 (kg m-2 s)
 ! assign pointers to model state variables
 mLayerDepth      => mvar_data%var(iLookMVAR%mLayerDepth)%dat          ! depth of each layer (m)
 mLayerTemp       => mvar_data%var(iLookMVAR%mLayerTemp)%dat           ! temperature of each layer (K)
 mLayerVolFracIce => mvar_data%var(iLookMVAR%mLayerVolFracIce)%dat     ! volumetric fraction of ice in each layer  (-)
 mLayerVolFracLiq => mvar_data%var(iLookMVAR%mLayerVolFracLiq)%dat     ! volumetric fraction of liquid water in each layer  (-)
 ! assign local pointers to the model index structures
 nLayers          => indx_data%var(iLookINDEX%nLayers)%dat(1)          ! number of layers
 layerType        => indx_data%var(iLookINDEX%layerType)%dat           ! layer type (ix_soil or ix_snow)

 ! identify the number of snow and soil layers
 nSnow = count(layerType==ix_snow)
 nSoil = count(layerType==ix_soil)

 ! NOTE: still need to process the case of "snow without a layer"
 if(nSnow==0)return

 ! initialize the weight of snow above each layer (kg m-2)
 weightSnow = 0._dp

 ! loop through snow layers
 do iSnow=1,nSnow
  ! save mass of liquid water and ice (mass does not change)
  massIceOld = iden_ice*mLayerVolFracIceIter(iSnow)*mLayerDepth(iSnow)   ! (kg m-2)
  massLiqOld = iden_water*mLayerVolFracLiqIter(iSnow)*mLayerDepth(iSnow) ! (kg m-2)
  ! *** compute the compaction associated with grain growth (s-1)
  ! compute the base rate of grain growth (-)
  if(mLayerVolFracIceIter(iSnow)*iden_ice <snwden_min) chi1=1._dp
  if(mLayerVolFracIceIter(iSnow)*iden_ice>=snwden_min) chi1=exp(-densScalGrowth*(mLayerVolFracIceIter(iSnow)*iden_ice - snwden_min))
  ! compute the reduction of grain growth under colder snow temperatures (-)
  chi2 = exp(-tempScalGrowth*(Tfreeze - mLayerTemp(iSnow)))
  ! compute the acceleration of grain growth in the presence of liquid water (-)
  if(mLayerVolFracLiq(iSnow) > wetSnowThresh)then; chi3=2._dp  ! snow is "wet"
  else; chi3=1._dp; endif                                      ! snow is "dry"
  ! compute the compaction associated with grain growth (s-1)
  CR_grainGrowth = grainGrowthRate*chi1*chi2*chi3
  ! **** compute the compaction associated with over-burden pressure (s-1)
  ! compute the weight imposed on the current layer (kg m-2)
  halfWeight = (massIceOld + massLiqOld)/2._dp  ! there is some over-burden pressure from the layer itself
  weightSnow = weightSnow + halfweight          ! add half of the weight from the current layer
  ! compute the increase in compaction under colder snow temperatures (-)
  chi4 = exp(-tempScalOvrbdn*(Tfreeze - mLayerTemp(iSnow)))
  ! compute the increase in compaction under low density snow (-)
  chi5 = exp(-densScalOvrbdn*mLayerVolFracIceIter(iSnow)*iden_ice)
  ! compute the compaction associated with over-burden pressure (s-1)
  CR_ovrvdnPress = (weightSnow/base_visc)*chi4*chi5
  ! update the snow weight with the halfWeight not yet used
  weightSnow = weightSnow + halfweight          ! add half of the weight from the current layer
  ! *** compute the compaction rate associated with snow melt (s-1)
  if(iden_ice*mLayerVolFracIceIter(iSnow) < snwDensityMax)then ! only collapse layers if below a critical density
   fracMelt    = min(mLayerVolFracIceIter(iSnow) - mLayerVolFracIce(iSnow), 0._dp)/mLayerVolFracIce(iSnow)
   CR_snowmelt = fracMelt/dt
  else
   CR_snowmelt = 0._dp
  endif
  ! compute the total compaction rate associated with metamorphism
  CR_metamorph = CR_grainGrowth + CR_ovrvdnPress
  ! update state variables -- note snowmelt is implicit, so can be updated directly
  mLayerDepth(iSnow)         = mLayerDepth(iSnow)/(1._dp + CR_metamorph*dt) + dt*CR_snowmelt*mLayerDepth(iSnow)
  mLayerVolFracIceNew(iSnow) = massIceOld/(mLayerDepth(iSnow)*iden_ice)
  mLayerVolFracLiqNew(iSnow) = massLiqOld/(mLayerDepth(iSnow)*iden_water)
 end do  ! looping through snow layers

 ! update coordinate variables
 call calcHeight(err,cmessage)
 if(err/=0)then; err=20; message=trim(message)//trim(cmessage); return; endif

 ! check for low/high snow density
 if(any(mLayerVolFracIceNew(1:nSnow)*iden_ice < 50._dp) .or. &
    any(mLayerVolFracIceNew(1:nSnow)*iden_ice > 900._dp))then
  err=20; message=trim(message)//'unreasonable value for snow density'
  return
 endif

 end subroutine snwDensify

end module snwDensify_module

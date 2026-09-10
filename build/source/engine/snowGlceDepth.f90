module snowGlceDepth_module

! data types
USE nr_type
USE data_types,only:&
                    var_ilength,   & ! x%var(:)%dat            (i4b)
                    var_dlength,   & ! x%var(:)%dat            (rkind)
                    zLookup          ! x%z(:)%var(:)%lookup(:) (rkind)

! constants
USE multiconst,only:&
                    Tfreeze,      & ! freezing temperature (K)
                    iden_ice,     & ! intrinsic density of ice             (kg m-3)
                    iden_air,     & ! intrinsic density of air (kg m-3)
                    iden_water      ! intrinsic density of liquid water    (kg m-3)
USE globalData,only:verySmall       ! a small number
USE globalData,only:verySmaller     ! a smaller number than verySmall

! named variables for parent structures
USE var_lookup,only:iLookDECISIONS   ! named variables for elements of the decision structure
USE var_lookup,only:iLookPROG        ! named variables for structure elements
USE var_lookup,only:iLookDIAG        ! named variables for structure elements
USE var_lookup,only:iLookFLUX        ! named variables for structure elements
USE var_lookup,only:iLookPARAM       ! named variables for structure elements
USE var_lookup,only:iLookINDEX       ! named variables for structure elements

! privacy
implicit none
private
public::snowGlceDepth
contains

! ************************************************************************************************
! public subroutine snowGlceDepth: compute snow depth for one sub timestep
! ************************************************************************************************
subroutine snowGlceDepth(&
                         dt_sub,                   & ! intent(in):    time step (s)
                         nSnow,                    & ! intent(in):    number of snow layers
                         nLake,                    & ! intent(in):    number of lake layers
                         nSoil,                    & ! intent(in):    number of soil layers
                         nGlce,                    & ! intent(in):    number of glacier ice layers
                         noThetaChange,            & ! intent(in):    number of layers with no change in total water content (bottom layers)
                         scalarGroundSublimation,  & ! intent(in):    scalar sublimation of snow/ice (kg m-2)
                         mLayerVolFracLiq,         & ! intent(inout): volumetric fraction of liquid water
                         mLayerVolFracIce,         & ! intent(inout): volumetric fraction of ice
                         mLayerTemp,               & ! intent(in):    temperature of each layer (K)
                         mLayerMeltFreeze,         & ! intent(in):    volumetric melt in each layer (kg m-3)
                         mpar_data,                & ! intent(in):    model parameters
                         ! output
                         glceReduceLiq,            & ! intent(out):   liquid water squeezed out of glacier layers (m s-1)
                         tooMuchSublim,            & ! intent(out):   flag to denote that there was too much sublimation in a given time step
                         mLayerDepth,              & ! intent(inout): depth of each layer (m)
                         ! error control
                         err,message)                ! intent(out):   error control
  ! -----------------------------------------------------------------------------------------------------------------------------------------
  implicit none
  real(qp),intent(in)                  :: dt_sub                   ! time step (s)
  integer(i4b),intent(in)              :: nSnow                    ! number of snow layers
  integer(i4b),intent(in)              :: nLake                    ! number of lake layers
  integer(i4b),intent(in)              :: nSoil                    ! number of soil layers
  integer(i4b),intent(in)              :: nGlce                    ! number of glacier ice layers
  integer(i4b),intent(in)              :: noThetaChange            ! number of layers with no change in total water content (bottom layers)
  real(rkind),intent(in)               :: scalarGroundSublimation  ! scalar sublimation of snow/ice (kg m-2 s-1)
  real(rkind),intent(inout)            :: mLayerVolFracLiq(:)      ! volumetric fraction of liquid water
  real(rkind),intent(inout)            :: mLayerVolFracIce(:)      ! volumetric fraction of ice
  real(rkind),intent(in)               :: mLayerTemp(:)            ! temperature of each layer (K)
  real(rkind),intent(in)               :: mLayerMeltFreeze(:)      ! volumetric melt in each layer (kg m-3)
  type(var_dlength),intent(in)         :: mpar_data                ! model parameters
  real(rkind),intent(out)              :: glceReduceLiq            ! liquid water squeezed out of glacier layers (m s-1)
  logical(lgt)                         :: tooMuchSublim            ! flag to denote that there was too much sublimation in a given time step
  real(rkind),intent(inout)            :: mLayerDepth(:)           ! depth of each layer (m)
  integer(i4b),intent(out)             :: err                      ! error code
  character(*),intent(out)             :: message                  ! error message
  ! local variables
  integer(i4b)                         :: nLayers                  ! total number of layers
  character(len=256)                   :: cmessage                 ! error message
  real(rkind)                          :: massLiquid               ! mass liquid water (kg m-2)
  ! -----------------------------------------------------------------------------------------------------------------------------------------
  ! initialize error control
  err=0; message="snowGlceDepth/"
  nLayers = nSnow + nLake + nSoil + nGlce ! total number of layers

  ! *** compute change in ice content of the top snow or glacier ice layer due to sublimation
  tooMuchSublim = .false.  ! initialize too much sublimation (merge layers) to false
  if(nSnow>0 .or. (nLake>0 .and. mLayerTemp(1)<Tfreeze) .or. (nGlce>0 .and. nSoil==0) )then ! snow or ice layers exist on top
    massLiquid = mLayerDepth(1)*mLayerVolFracLiq(1)*iden_water ! save the mass of liquid water (kg m-2)

    ! add/remove the depth of snow/ice gained/lost by frost/sublimation (m)
    mLayerDepth(1) = mLayerDepth(1) + dt_sub*scalarGroundSublimation/(mLayerVolFracIce(1)*iden_ice) ! assume constant density
    if(mLayerDepth(1) < verySmall)then ! check that we did not remove the entire layer
      tooMuchSublim = .true.
      return
    endif
    mLayerVolFracLiq(1) = massLiquid / (mLayerDepth(1)*iden_water) ! update the volumetric fraction of liquid water
  else ! no snow or ice
    if(abs(scalarGroundSublimation) > verySmall)then ! check that sublimation is zero
      message=trim(message)//'sublimation of snow and ice has been computed when no snow or ice exists'
      err=20; return
    end if
  end if  ! (if snow or ice layers exist on top)

  ! *** account for compaction and cavitation in the snowpack or in firn...
  if(nSnow>0)then
    call snowDensify(&
                    ! intent(in): variables
                    dt_sub,                                            & ! intent(in):    time step (s)
                    nSnow,                                             & ! intent(in):    number of snow layers
                    mLayerTemp(1:nSnow),                               & ! intent(in):    temperature of each layer (K)
                    mLayerMeltFreeze(1:nSnow),                         & ! intent(in):    volumetric melt in each layer (kg m-3)
                    ! intent(in): parameters
                    mpar_data%var(iLookPARAM%densScalGrowth)%dat(1),   & ! intent(in):    density scaling factor for grain growth (kg-1 m3)
                    mpar_data%var(iLookPARAM%tempScalGrowth)%dat(1),   & ! intent(in):    temperature scaling factor for grain growth (K-1)
                    mpar_data%var(iLookPARAM%grainGrowthRate)%dat(1),  & ! intent(in):    rate of grain growth (s-1)
                    mpar_data%var(iLookPARAM%densScalOvrbdn)%dat(1),   & ! intent(in):    density scaling factor for overburden pressure (kg-1 m3)
                    mpar_data%var(iLookPARAM%tempScalOvrbdn)%dat(1),   & ! intent(in):    temperature scaling factor for overburden pressure (K-1)
                    mpar_data%var(iLookPARAM%baseViscosity)%dat(1),    & ! intent(in):    viscosity coefficient at T=T_frz and snow density=0 (kg m-2 s)
                    ! intent(inout): state variables
                    mLayerDepth(1:nSnow),                              & ! intent(inout): depth of each layer (m)
                    mLayerVolFracLiq(1:nSnow),                         & ! intent(inout): volumetric fraction of liquid water after iterations (-)
                    mLayerVolFracIce(1:nSnow),                         & ! intent(inout): volumetric fraction of ice after iterations (-)
                    ! output: error control
                    err,cmessage)                     ! intent(out): error control
    if(err/=0)then; err=55; message=trim(message)//trim(cmessage); return; end if
  end if ! if snow layers exist

  ! *** reduce ice depth if melt occurred
  if(nGlce>0)then
    call glceReduce(&
                    ! intent(in): variables
                    dt_sub,                                                      & ! intent(in):    time step (s)
                    nGlce-noThetaChange,                                         & ! intent(in):    number of glacier ice layers to reduce
                    mLayerMeltFreeze(nSnow+nLake+nSoil+1:nLayers-noThetaChange), & ! intent(in):    volumetric melt in each layer (kg m-3)
                    ! intent(inout): state variables
                    mLayerDepth(nSnow+nLake+nSoil+1:nLayers-noThetaChange),      & ! intent(inout): depth of each layer (m)
                    mLayerVolFracLiq(nSnow+nLake+nSoil+1:nLayers-noThetaChange), & ! intent(inout): volumetric fraction of liquid water (-)
                    mLayerVolFracIce(nSnow+nLake+nSoil+1:nLayers-noThetaChange), & ! intent(inout): volumetric fraction of ice (-)
                    glceReduceLiq,                                               & ! intent(out):   liquid water squeezed out of glacier layers (m s-1)
                    ! output: error control
                    tooMuchSublim,                                               & ! intent(inout): flag to denote that there was too much melt in a given time step
                    err,cmessage)                                                  ! intent(out):   error controls
    if(err/=0)then; err=55; message=trim(message)//trim(cmessage); return; end if
  else
    glceReduceLiq = 0._rkind
  endif ! if glacier ice layers exist

end subroutine snowGlceDepth


! ************************************************************************************************
! private subroutine snowDensify: compute change in snow or firn density over the time step
! ************************************************************************************************
subroutine snowDensify(&
                      ! intent(in): variables
                      dt,                             & ! intent(in): time step (s)
                      nSnow,                          & ! intent(in): number of snow layers
                      mLayerTemp,                     & ! intent(in): temperature of each layer (K)
                      mLayerMeltFreeze,               & ! intent(in): volumnetric melt in each layer (kg m-3)
                      ! intent(in): parameters
                      densScalGrowth,                 & ! intent(in): density scaling factor for grain growth (kg-1 m3)
                      tempScalGrowth,                 & ! intent(in): temperature scaling factor for grain growth (K-1)
                      grainGrowthRate,                & ! intent(in): rate of grain growth (s-1)
                      densScalOvrbdn,                 & ! intent(in): density scaling factor for overburden pressure (kg-1 m3)
                      tempScalOvrbdn,                 & ! intent(in): temperature scaling factor for overburden pressure (K-1)
                      baseViscosity,                  & ! intent(in): viscosity coefficient at T=T_frz and snow density=0 (kg m-2 s)
                      ! intent(inout): state variables
                      mLayerDepth,                    & ! intent(inout): depth of each layer (m)
                      mLayerVolFracLiqNew,            & ! intent(inout):  volumetric fraction of liquid water after iterations (-)
                      mLayerVolFracIceNew,            & ! intent(inout):  volumetric fraction of ice after iterations (-)
                      ! output: error control
                      err,message)                      ! intent(out): error control
  ! -----------------------------------------------------------------------------------------------------------------------------------------
  implicit none
  ! intent(in): variables
  real(rkind),intent(in)              :: dt                       ! time step (seconds)
  integer(i4b),intent(in)             :: nSnow                    ! number of snow layers
  real(rkind),intent(in)              :: mLayerTemp(:)            ! temperature of each snow layer after iterations (K)
  real(rkind),intent(in)              :: mLayerMeltFreeze(:)      ! volumetric melt in each layer (kg m-3)
  ! intent(in): parameters
  real(rkind),intent(in)              :: densScalGrowth           ! density scaling factor for grain growth (kg-1 m3)
  real(rkind),intent(in)              :: tempScalGrowth           ! temperature scaling factor for grain growth (K-1)
  real(rkind),intent(in)              :: grainGrowthRate          ! rate of grain growth (s-1)
  real(rkind),intent(in)              :: densScalOvrbdn           ! density scaling factor for overburden pressure (kg-1 m3)
  real(rkind),intent(in)              :: tempScalOvrbdn           ! temperature scaling factor for overburden pressure (K-1)
  real(rkind),intent(in)              :: baseViscosity            ! viscosity coefficient at T=T_frz and snow density=0 (kg m-2 s)
  ! intent(inout): state variables
  real(rkind),intent(inout)           :: mLayerDepth(:)           ! depth of each layer (m)
  real(rkind),intent(inout)           :: mLayerVolFracLiqNew(:)   ! volumetric fraction of liquid water in each snow layer after iterations (-)
  real(rkind),intent(inout)           :: mLayerVolFracIceNew(:)   ! volumetric fraction of ice in each snow layer after iterations (-)
  ! intent(out): error control
  integer(i4b),intent(out)            :: err                      ! error code
  character(*),intent(out)            :: message                  ! error message
  ! -----------------------------------------------------------------------------------------------------------------------------------------
  ! define local variables
  integer(i4b)                        :: iSnow                       ! index of snow layers
  real(rkind)                         :: chi1,chi2,chi3,chi4,chi5    ! multipliers in the densification algorithm (-)
  real(rkind)                         :: halfWeight                  ! half of the weight of the current snow layer (kg m-2)
  real(rkind)                         :: weightSnow                  ! total weight of snow above the current snow layer (kg m-2)
  real(rkind)                         :: CR_grainGrowth              ! compaction rate for grain growth (s-1)
  real(rkind)                         :: CR_ovrbdnPress              ! compaction rate associated with over-burden pressure (s-1)
  real(rkind)                         :: CR_metamorph                ! compaction rate for metamorphism (s-1)
  real(rkind)                         :: massIceOld                  ! mass of ice in the snow layer (kg m-2)
  real(rkind)                         :: massLiqOld                  ! mass of liquid water in the snow layer (kg m-2)
  real(rkind)                         :: scalarDepthNew              ! updated layer depth (m)
  real(rkind)                         :: scalarDepthMin              ! minimum layer depth (m)
  real(rkind)                         :: volFracIceLoss              ! volumetric fraction of ice lost due to melt (-)
  real(rkind), dimension(nSnow)       :: mLayerVolFracAirNew         ! volumetric fraction of air in each layer after compaction (-)
  real(rkind),parameter               :: snwden_min=100._rkind       ! minimum snow density for reducing metamorphism rate (kg m-3)
  real(rkind),parameter               :: snwDensityMax=550._rkind    ! maximum snow density for collapse under melt (kg m-3)
  real(rkind),parameter               :: wetSnowThresh=0.01_rkind    ! threshold to discriminate between "wet" and "dry" snow
  real(rkind),parameter               :: minLayerDensity=40._rkind   ! minimum snow density allowed for any layer (kg m-3)
  ! -----------------------------------------------------------------------------------------------------------------------------------------
  ! initialize error control
  err=0; message="snowDensify/"

  weightSnow = 0._rkind ! initialize the weight of snow above each layer (kg m-2)

  ! loop through snow layers
  do iSnow=1,nSnow
    ! save mass of liquid water and ice (mass does not change)
    massIceOld = iden_ice*mLayerVolFracIceNew(iSnow)*mLayerDepth(iSnow)   ! (kg m-2)
    massLiqOld = iden_water*mLayerVolFracLiqNew(iSnow)*mLayerDepth(iSnow) ! (kg m-2)

    ! *** compute the compaction associated with grain growth (s-1)
    ! base rate of grain growth
    if(mLayerVolFracIceNew(iSnow)*iden_ice <snwden_min)then ! low density snow
      chi1=1._rkind
    else ! high density snow
      chi1=exp(-densScalGrowth*(mLayerVolFracIceNew(iSnow)*iden_ice - snwden_min))
    end if
    ! reduction of grain growth under colder snow temperatures (-)
    chi2 = exp(-tempScalGrowth*(Tfreeze - mLayerTemp(iSnow)))
    ! acceleration of grain growth in the presence of liquid water (-)
    if(mLayerVolFracLiqNew(iSnow) > wetSnowThresh)then ! snow is "wet"
      chi3=2._rkind
    else ! snow is "dry"
      chi3=1._rkind
    end if
    CR_grainGrowth = grainGrowthRate*chi1*chi2*chi3 ! compaction associated with grain growth (s-1)

    ! **** compute the compaction associated with over-burden pressure (s-1)
    halfWeight = (massIceOld + massLiqOld)/2._rkind  ! half of the weight from the current layer
    weightSnow = weightSnow + halfweight             ! there is some over-burden pressure from the layer itself, add half
    chi4 = exp(-tempScalOvrbdn*(Tfreeze - mLayerTemp(iSnow)))       ! increase in compaction under colder snow temperatures (-)
    chi5 = exp(-densScalOvrbdn*mLayerVolFracIceNew(iSnow)*iden_ice) ! increase in compaction under low density snow (-)
    CR_ovrbdnPress = (weightSnow/baseViscosity)*chi4*chi5 ! compaction associated with over-burden pressure (s-1)
    weightSnow = weightSnow + halfweight                  ! update the snow weight, add half of the weight from the current layer

    ! *** compute the compaction rate associated with snowmelt (s-1)
    ! NOTE: loss of ice due to snowmelt is implicit, so can be updated directly
    if(iden_ice*mLayerVolFracIceNew(iSnow) < snwDensityMax)then ! only collapse layers if below a critical density
      volFracIceLoss = max(0._rkind,mLayerMeltFreeze(iSnow)/iden_ice)  ! volumetric fraction of ice lost due to melt (-)
      scalarDepthNew = mLayerDepth(iSnow) * mLayerVolFracIceNew(iSnow)/(mLayerVolFracIceNew(iSnow) + volFracIceLoss) ! adjust snow depth to account for cavitation
    else ! do not allow collapse of snow layers
      scalarDepthNew = mLayerDepth(iSnow)
    end if

    ! **** compute the total compaction rate associated with metamorphism
    CR_metamorph = CR_grainGrowth + CR_ovrbdnPress
    ! update depth due to metamorphism (implicit solution)
    scalarDepthNew = scalarDepthNew/(1._rkind + CR_metamorph*dt)
    ! ensure that the new depth is in line with the maximum amount of compaction that can occur given the masses of ice and liquid in the layer
    scalarDepthMin = (massIceOld / iden_ice) + (massLiqOld / iden_water)
    mLayerDepth(iSnow) = max(scalarDepthMin, scalarDepthNew)

    ! update volumetric ice and liquid water content
    mLayerVolFracIceNew(iSnow) = massIceOld/(mLayerDepth(iSnow)*iden_ice)
    mLayerVolFracLiqNew(iSnow) = massLiqOld/(mLayerDepth(iSnow)*iden_water)
    mLayerVolFracAirNew(iSnow) = 1.0_rkind - mLayerVolFracIceNew(iSnow) - mLayerVolFracLiqNew(iSnow)

  end do  ! looping through snow layers

  ! check depth
  if(any(mLayerDepth(1:nSnow) < 0._rkind))then
    do iSnow=1,nSnow
      write(*,'(a,1x,i4,1x,4(f12.5,1x))') 'iSnow, mLayerDepth(iSnow)', iSnow, mLayerDepth(iSnow)
    end do
    message=trim(message)//'unreasonable value for snow depth'
    err=20; return
  end if

  ! check for low/high snow density
  if(any(mLayerVolFracIceNew(1:nSnow)*iden_ice + mLayerVolFracLiqNew(1:nSnow)*iden_water + mLayerVolFracAirNew(1:nSnow)*iden_air < minLayerDensity) .or. &
     any(mLayerVolFracIceNew(1:nSnow) + mLayerVolFracLiqNew(1:nSnow) + mLayerVolFracAirNew(1:nSnow) > 1._rkind))then
    do iSnow=1,nSnow
      write(*,*) 'iSnow, volFracIce, density = ', iSnow, mLayerVolFracIceNew(iSnow),  mLayerVolFracIceNew(iSnow)*iden_ice
    end do
    message=trim(message)//'unreasonable value for snow density'
    err=20; return
  end if

end subroutine snowDensify


! ************************************************************************************************
! private subroutine glceReduce: compute change ice depth over the time step
! ************************************************************************************************
subroutine glceReduce(&
                      ! intent(in): variables
                      dt,                             & ! intent(in):    time step (s)
                      nGlce,                          & ! intent(in):    number of glacier ice layers
                      mLayerMeltFreeze,               & ! intent(in):    volumetric melt in each layer (kg m-3)
                       ! intent(inout): state variables
                      mLayerDepth,                    & ! intent(inout): depth of each layer (m)
                      mLayerVolFracLiqNew,            & ! intent(inout): volumetric fraction of liquid water (-)
                      mLayerVolFracIceNew,            & ! intent(inout): volumetric fraction of ice (-)
                      glceReduceLiq,                  & ! intent(out):   glacier liquid water squeezed out of the layers (kg m-2)
                      ! output: error control
                      tooMuchMelt,                    & ! intent(inout): flag to denote that there was too much melt in a given time step
                      err,message)                      ! intent(out):   error control
  ! -----------------------------------------------------------------------------------------------------------------------------------------
  ! compute change in snow density over the time step
  implicit none
  ! intent(in): variables
  real(rkind),intent(in)              :: dt                       ! time step (seconds)
  integer(i4b),intent(in)             :: nGlce                    ! number of glacier ice layers
  real(rkind),intent(in)              :: mLayerMeltFreeze(:)      ! volumetric melt in each layer (kg m-3)
  ! intent(inout): state variables
  real(rkind),intent(inout)           :: mLayerDepth(:)           ! depth of each layer (m)
  real(rkind),intent(inout)           :: mLayerVolFracLiqNew(:)   ! volumetric fraction of liquid water in each snow layer after iterations (-)
  real(rkind),intent(inout)           :: mLayerVolFracIceNew(:)   ! volumetric fraction of ice in each snow layer after iterations (-)
  real(rkind),intent(out)             :: glceReduceLiq            ! glacier liquid water squeezed out of the layers (kg m-2)
  ! intent(out): error control
  logical(lgt),intent(inout)          :: tooMuchMelt              ! flag to denote that there was too much melt in a given time step
  integer(i4b),intent(out)            :: err                      ! error code
  character(*),intent(out)            :: message                  ! error message
  ! -----------------------------------------------------------------------------------------------------------------------------------------
  ! define local variables
  integer(i4b)                        :: iGlce                       ! index of glce layers
  real(rkind)                         :: massIceOld                  ! mass of ice in the snow layer (kg m-2)
  real(rkind)                         :: massLiqOld                  ! mass of liquid water in the snow layer (kg m-2)
  real(rkind)                         :: scalarDepthNew              ! updated layer depth (m)
  real(rkind)                         :: scalarDepthMin              ! minimum layer depth (m)
  real(rkind)                         :: volFracIceChange            ! volumetric fraction of ice lost/gained due to melt and sublimation (-)
  real(rkind)                         :: massLiqRetained             ! retained liquid water in the shrunken layer (kg m-2)
  real(rkind)                         :: layerReduceLiq              ! liquid water squeezed from current layer (kg m-2)
  real(rkind), dimension(nGlce)       :: mLayerVolFracAirNew         ! volumetric fraction of air in each layer after compaction (-)
  real(rkind),parameter               :: minLayerDensity=730._rkind  ! minimum glce density allowed for any layer (kg m-3) (roughly 80% ice)
  ! -----------------------------------------------------------------------------------------------------------------------------------------
  ! initialize error control
  err=0; message="glceReduce/"
  glceReduceLiq = 0._rkind ! initialize the amount of liquid water squeezed out of the glacier ice layers (kg m-2)

  ! loop through glce layers
  do iGlce=1,nGlce
    ! save mass of liquid water and ice (mass does not change)
    massIceOld = iden_ice*mLayerVolFracIceNew(iGlce)*mLayerDepth(iGlce)   ! (kg m-2)
    massLiqOld = iden_water*mLayerVolFracLiqNew(iGlce)*mLayerDepth(iGlce) ! (kg m-2)

    ! adjust depth in proportion to the amount of ice lost/gained in the layer (m), keeping ice fraction constant
    volFracIceChange = mLayerMeltFreeze(iGlce)/iden_ice  ! volumetric fraction of ice lost/gained due to melt/refreeze (-)
    scalarDepthNew = massIceOld/((mLayerVolFracIceNew(iGlce)+volFracIceChange)*iden_ice)
    mLayerDepth(iGlce) = scalarDepthNew

    ! Keep liquid water fraction in the layer constant; excess liquid is squeezed out
    ! Note, if it got colder the depth icreases and thus adds a bit of water to the layer, but will be corrected when the layer is reduced to the active layer depth
    massLiqRetained = iden_water*mLayerVolFracLiqNew(iGlce)*mLayerDepth(iGlce) ! kg m-2
    layerReduceLiq = max(0._rkind, massLiqOld - massLiqRetained) ! only take away liquid water, don't add any
    glceReduceLiq = (glceReduceLiq + layerReduceLiq)/(iden_water*dt) ! convert to m s-1

    ! check that we did not remove the entire layer
    if(mLayerDepth(iGlce) < verySmall)then
      tooMuchMelt=.true.
      return
    endif

    ! update volumetric ice and liquid water content, should be the same but for consistency redo
    mLayerVolFracIceNew(iGlce) = massIceOld/(mLayerDepth(iGlce)*iden_ice) ! mLayerVolFracIceNew(iGlce)+volFracIceChange
    mLayerVolFracLiqNew(iGlce) = massLiqRetained/(mLayerDepth(iGlce)*iden_water)
    mLayerVolFracAirNew(iGlce) = 1.0_rkind - mLayerVolFracIceNew(iGlce) - mLayerVolFracLiqNew(iGlce)

  end do  ! looping through glce layers

  ! check for low/high glce density, could happen from initial condition consistency recalculation of liquid if start very near the freezing point
  if(any(mLayerVolFracIceNew(1:nGlce)*iden_ice + mLayerVolFracLiqNew(1:nGlce)*iden_water + mLayerVolFracAirNew(1:nGlce)*iden_air < minLayerDensity) .or. &
     any(mLayerVolFracIceNew(1:nGlce) + mLayerVolFracLiqNew(1:nGlce) + mLayerVolFracAirNew(1:nGlce) > 1._rkind))then
    do iGlce=1,nGlce
      write(*,*) 'iGlce, volFracIce, density = ', iGlce, mLayerVolFracIceNew(iGlce),  mLayerVolFracIceNew(iGlce)*iden_ice
    end do
    message=trim(message)//'unreasonable value for glce density, likely due to initial condition near Tfreeze'
    err=20; return
  end if

end subroutine glceReduce

end module snowGlceDepth_module

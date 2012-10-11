module groundwatr_module
USE nrtype
implicit none
private
public::groundwatr
public::waterTableHeight
contains

 ! ************************************************************************************************
 ! new subroutine: compute change in snow density over the time step
 ! ************************************************************************************************
 subroutine groundwatr(&
                       ! input
                       dt,                                     & ! input:  time step (seconds) 
                       scalarAquiferRcharge,                   & ! input:  aquifer recharge (m s-1)
                       scalarAquiferTranspire,                 & ! input:  aquifer transpiration (m s-1)
                       mLayerVolFracLiqIter,                   & ! input:  volumetric fraction of liquid water after itertations (-)
                       mLayerVolFracIceIter,                   & ! input:  volumetric fraction of ice after itertations (-)
                       iLayerHeight,                           & ! input:  height of each interface (m)
                       mLayerDepth,                            & ! input:  depth of each soil layer (m)
                       theta_sat,                              & ! input:  soil porosity (-)
                       specificYield,                          & ! input:  fraction of water volume drained by gravity in an unconfined aquifer (-)
                       k_soil,                                 & ! input:  hydraulic conductivity (m s-1) 
                       kAnisotropic,                           & ! input:  anisotropy factor for lateral hydraulic conductivity (-)
                       zScale_TOPMODEL,                        & ! input:  scale factor for TOPMODEL-ish baseflow parameterization (m)
                       ! input-output
                       scalarAquiferStorage,                   & ! input-output: aquifer storage (m)
                       ! output
                       scalarWaterTableDepth,                  & ! output: water table depth at the end of the time step (m)
                       scalarAquiferBaseflow,                  & ! output: baseflow from the aquifer (m s-1)
                       err,message)                              ! output: error control
 ! compute change in aquifer storage over the time step
 implicit none
 ! input variables
 real(dp),intent(in)                 :: dt                       ! time step (seconds)
 real(dp),intent(in)                 :: scalarAquiferRcharge     ! aquifer recharge averaged over the time step (m s-1)
 real(dp),intent(in)                 :: scalarAquiferTranspire   ! aquifer transpiration averaged over the time step (m s-1)
 real(dp),intent(in)                 :: mLayerVolFracLiqIter(:)  ! volumetric fraction of liquid water after itertations (-)
 real(dp),intent(in)                 :: mLayerVolFracIceIter(:)  ! volumetric fraction of ice after itertations (-)
 real(dp),intent(in)                 :: iLayerHeight(0:)         ! height at each layer interface (m)
 real(dp),intent(in)                 :: mLayerDepth(:)           ! depth of each soil layer (m)
 real(dp),intent(in)                 :: theta_sat                ! total porosity (-)
 real(dp),intent(in)                 :: specificYield            ! fraction of water volume drained by gravity in an unconfined aquifer (-)
 real(dp),intent(in)                 :: k_soil                   ! hydraulic conductivity (m s-1) 
 real(dp),intent(in)                 :: kAnisotropic             ! anisotropy factor for lateral hydraulic conductivity (-)
 real(dp),intent(in)                 :: zScale_TOPMODEL          ! scale factor for TOPMODEL-ish baseflow parameterization (m)
 ! input-output variables
 real(dp),intent(inout)              :: scalarAquiferStorage     ! aquifer storage (m)
 ! output variables
 real(dp),intent(out)                :: scalarWaterTableDepth    ! water table depth at the end of the time step (m)
 real(dp),intent(out)                :: scalarAquiferBaseflow    ! baseflow from the aquifer (m s-1)
 integer(i4b),intent(out)            :: err                      ! error code
 character(*),intent(out)            :: message                  ! error message
 ! algorithmic control parameters
 real(dp),pointer                    :: wimplicit                ! weight assigned to start-of-step fluxes (-)
 integer(i4b),parameter              :: maxiter=10               ! maximum number of iterations
 real(dp),parameter                  :: incTol= 1.e-6_dp         ! convergence tolerance for the iteration increment (m)
 real(dp),parameter                  :: resTol= 1.e-6_dp         ! convergence tolerance for the residual (m)
 ! define local variables
 character(len=256)                  :: cmessage                 ! error message for downwind routine
 integer(i4b)                        :: iter                     ! iteration index
 integer(i4b)                        :: iSoil                    ! index of each soil layer
 integer(i4b)                        :: nSoil                    ! number of soil layers
 integer(i4b)                        :: nUnsat                   ! number of unsaturated layers
 real(dp)                            :: baseflowMax              ! maximum baseflow rate from the aquifer (m s-1)
 real(dp),dimension(0:size(mLayerVolFracLiqIter)) :: waterReq2FillPore  ! water required to fill pore space (m)
 real(dp)                            :: drainablePorosity        ! drainable porosity (-)
 real(dp)                            :: aquiferStorageTrial      ! trial value for aquifer storage (m)
 real(dp)                            :: zWaterTrial              ! depth of the water table (m)
 real(dp)                            :: dBaseflow_dStorage       ! derivative in the baseflow term w.r.t. aquifer storage (s-1)
 real(dp)                            :: netFlux                  ! recharge minus baseflow (m s-1)
 real(dp)                            :: residual                 ! residual in aquifer storage (m)
 real(dp)                            :: dStorage                 ! iteration increment for aquifer storage (m)

 ! initialize error control
 err=0; message="groundwatr/"

 ! compute the maximum baseflow rate
 baseflowMax = kAnisotropic*k_soil

 ! initialize the aquifer storage (m)
 aquiferStorageTrial = scalarAquiferStorage

 ! compute the water required to fill pore space up to the level of each layer interface (m)
 call waterReq2FillLayer(theta_sat,            & ! intent(in):  total porosity (-)
                         mLayerDepth,          & ! intent(in):  depth of each soil layer (m)
                         mLayerVolFracLiqIter, & ! intent(in):  volumetric liquid water content in each soil layer (-)
                         mLayerVolFracIceIter, & ! intent(in):  volumetric ice content in each soil layer (-)
                         waterReq2FillPore,    & ! intent(out): water required to fill pore space at each layer interface (m)
                         err,cmessage)           ! intent(out): error control
 if(err/=0)then; message=trim(message)//trim(cmessage); return; endif


 ! *****
 ! * start iterations...
 ! *********************
 do iter=1,maxiter

  ! -----
  ! compute the depth to the water table and the drainable porosity...
  ! ---------------------------------------
  ! compute the number of unsaturated layers, the drainable porosity, and the depth to the water table
  call waterTablePosition(mLayerVolFracLiqIter, & ! intent(in):  volumetric liquid water content in each soil layer (-)
                          mLayerVolFracIceIter, & ! intent(in):  volumetric ice content in each soil layer (-)
                          aquiferStorageTrial,  & ! intent(in):  aquifer storage (m)
                          waterReq2FillPore,    & ! intent(in):  water required to fill pore space up to the level of each layer interface (m(
                          iLayerHeight,         & ! intent(in):  height of each interface (m)
                          mLayerDepth,          & ! intent(in):  depth of each soil layer (m)
                          theta_sat,            & ! intent(in):  soil porosity (-)
                          specificYield,        & ! intent(in):  fraction of water volume drained by gravity in an unconfined aquifer (-)
                          nUnsat,               & ! intent(out): number of unsaturated layers
                          drainablePorosity,    & ! intent(out): drainable porosity (-)
                          zWaterTrial,          & ! intent(out): water table depth (m)
                          err,cmessage)           ! intent(out): error control
  if(err/=0)then; message=trim(message)//trim(cmessage); return; endif

  ! -----
  ! compute iteration increment...
  ! ------------------------------
  ! calculate the baseflow term (m s-1)
  scalarAquiferBaseflow = baseflowMax*exp(-zWaterTrial/zScale_TOPMODEL)
  ! compute the net flux (m s-1)
  netFlux = (scalarAquiferRcharge + scalarAquiferTranspire) - scalarAquiferBaseflow
  ! compute the residual in aquifer storage (m)
  residual = aquiferStorageTrial - (scalarAquiferStorage + netFlux*dt)
  ! compute the derivative in baseflow w.r.t. storage (s-1)
  dBaseflow_dStorage = scalarAquiferBaseflow/(zScale_TOPMODEL*drainablePorosity)
  ! compute the iteration increment (m)
  dStorage = -residual/(1._dp + dBaseflow_dStorage*dt)
  write(*,'(a)')            'aquiferStorageTrial, scalarWaterTableDepth, scalarAquiferBaseflow, scalarAquiferRcharge, dBaseflow_dStorage, dStorage = '
  write(*,'(10(e20.10,1x))') aquiferStorageTrial, scalarWaterTableDepth, scalarAquiferBaseflow, scalarAquiferRcharge, dBaseflow_dStorage, dStorage
 
  ! check for convergence and update states
  if(abs(dStorage) < incTol .and. abs(residual) < resTol)then
   ! update aquifer storage (m)
   scalarAquiferStorage  = aquiferStorageTrial + dStorage
   ! update the depth to the water table (m)
   call waterTablePosition(mLayerVolFracLiqIter, & ! intent(in):  volumetric liquid water content in each soil layer (-)
                           mLayerVolFracIceIter, & ! intent(in):  volumetric ice content in each soil layer (-)
                           scalarAquiferStorage, & ! intent(in):  aquifer storage (m)
                           waterReq2FillPore,    & ! intent(in):  water required to fill pore space up to the level of each layer interface (m)
                           iLayerHeight,         & ! intent(in):  height of each interface (m)
                           mLayerDepth,          & ! intent(in):  depth of each soil layer (m)
                           theta_sat,            & ! intent(in):  soil porosity (-)
                           specificYield,        & ! intent(in):  fraction of water volume drained by gravity in an unconfined aquifer (-)
                           nUnsat,               & ! intent(out): number of unsaturated layers
                           drainablePorosity,    & ! intent(out): drainable porosity (-)
                           scalarWaterTableDepth,& ! intent(out): water table depth (m)
                           err,cmessage)           ! intent(out): error control
   if(err/=0)then; message=trim(message)//trim(cmessage); return; endif
   exit
  ! ***** non-convergence -- update the aquifer storage and water table depth
  else
   ! ** update the aquifer storage (m) 
   aquiferStorageTrial = aquiferStorageTrial + dStorage
  endif  ! ***** non-convergence -- update
  ! check for non-convergence
  if(iter==maxiter)then; err=20; message=trim(message)//'failed to converge'; return; endif
 end do  ! (loop through iterations)

 ! clean-up fluxes
 scalarAquiferBaseflow = scalarAquiferBaseflow + dBaseflow_dStorage*dStorage

 end subroutine groundwatr

 ! ---------------------------------------------------------------------------------------------------------------------------
 ! new subroutine: compute the water table depth (m)
 ! ---------------------------------------------------------------------------------------------------------------------------
 subroutine waterTableHeight(mLayerVolFracLiqIter,      & ! intent(in):  volumetric liquid water content in each soil layer (-)
                             mLayerVolFracIceIter,      & ! intent(in):  volumetric ice content in each soil layer (-)
                             aquiferStorageTrial,       & ! intent(in):  aquifer storage (m)
                             iLayerHeight,              & ! intent(in):  height of each interface (m)
                             mLayerDepth,               & ! intent(in):  depth of each soil layer (m)
                             theta_sat,                 & ! intent(in):  soil porosity (-)
                             specificYield,             & ! intent(in):  fraction of water volume drained by gravity in an unconfined aquifer (-)
                             waterTableDepth,           & ! intent(out): water table depth (m)
                             err,message)                 ! intent(out): error control
 implicit none
 ! input
 real(dp),intent(in)        :: mLayerVolFracLiqIter(:) ! volumetric liquid water content in each soil layer (-)
 real(dp),intent(in)        :: mLayerVolFracIceIter(:) ! volumetric ice content in each soil layer (-)
 real(dp),intent(in)        :: aquiferStorageTrial     ! aquifer storage (m)
 real(dp),intent(in)        :: iLayerHeight(0:)        ! height at each layer interface (m)
 real(dp),intent(in)        :: mLayerDepth(:)          ! depth of each soil layer (m)
 real(dp),intent(in)        :: theta_sat               ! total porosity (-)
 real(dp),intent(in)        :: specificYield           ! fraction of water volume drained by gravity in an unconfined aquifer (-)
 ! output
 real(dp),intent(out)       :: waterTableDepth       ! water table depth (m)
 integer(i4b),intent(out)   :: err                   ! error code
 character(*),intent(out)   :: message               ! error message
 ! internal
 character(len=256)         :: cmessage              ! error message for downwind routine
 integer(i4b)               :: nSoil                 ! number of soil layers
 integer(i4b)               :: nUnsat                ! number of unsaturated layers
 real(dp)                   :: drainablePorosity     ! drainable porosity (-)
 real(dp),dimension(0:size(mLayerVolFracLiqIter)) :: waterReq2FillPore  ! water required to fill pore space (m)
 ! initialize error control
 err=0; message='waterTableHeight/'

 ! ***** case where water table is below the soil column
 if(aquiferStorageTrial < 0._dp)then
  nSoil             = size(mLayerVolFracLiqIter)
  waterTableDepth   = iLayerHeight(nSoil) - aquiferStorageTrial/specificYield
 ! ***** water table within the soil column
 else
  ! compute the water required to fill pore space up to the level of each layer interface (m)
  call waterReq2FillLayer(theta_sat,            & ! intent(in):  total porosity (-)
                          mLayerDepth,          & ! intent(in):  depth of each soil layer (m)
                          mLayerVolFracLiqIter, & ! intent(in):  volumetric liquid water content in each soil layer (-)
                          mLayerVolFracIceIter, & ! intent(in):  volumetric ice content in each soil layer (-)
                          waterReq2FillPore,    & ! intent(out): water required to fill pore space at each layer interface (m)
                          err,cmessage)           ! intent(out): error control
  if(err/=0)then; message=trim(message)//trim(cmessage); return; endif

  ! compute the number of unsaturated layers, the drainable porosity, and the depth to the water table
  call waterTablePosition(mLayerVolFracLiqIter, & ! intent(in):  volumetric liquid water content in each soil layer (-)
                          mLayerVolFracIceIter, & ! intent(in):  volumetric ice content in each soil layer (-)
                          aquiferStorageTrial,  & ! intent(in):  aquifer storage (m)
                          waterReq2FillPore,    & ! intent(in):  water required to fill pore space up to the level of each layer interface (m(
                          iLayerHeight,         & ! intent(in):  height of each interface (m)
                          mLayerDepth,          & ! intent(in):  depth of each soil layer (m)
                          theta_sat,            & ! intent(in):  soil porosity (-)
                          specificYield,        & ! intent(in):  fraction of water volume drained by gravity in an unconfined aquifer (-)
                          nUnsat,               & ! intent(out): number of unsaturated layers
                          drainablePorosity,    & ! intent(out): drainable porosity (-)
                          waterTableDepth,      & ! intent(out): water table depth (m)
                          err,cmessage)           ! intent(out): error control
  if(err/=0)then; message=trim(message)//trim(cmessage); return; endif
 endif  ! (if water table is within the soil column)

 end subroutine waterTableHeight


 ! ---------------------------------------------------------------------------------------------------------------------------
 ! private subroutine: compute the water required to fill pore space at each layer interface (m)
 ! ---------------------------------------------------------------------------------------------------------------------------
 subroutine waterReq2FillLayer(theta_sat,            & ! intent(in):  total porosity (-)
                               mLayerDepth,          & ! intent(in):  depth of each soil layer (m)
                               mLayerVolFracLiqIter, & ! intent(in):  volumetric liquid water content in each soil layer (-)
                               mLayerVolFracIceIter, & ! intent(in):  volumetric ice content in each soil layer (-)
                               waterReq2FillPore,    & ! intent(out): water required to fill pore space at each layer interface (m)
                               err,message)            ! intent(out): error control
 implicit none
 ! input
 real(dp),intent(in)        :: theta_sat               ! total porosity (-)
 real(dp),intent(in)        :: mLayerDepth(:)          ! depth of each soil layer (m)
 real(dp),intent(in)        :: mLayerVolFracLiqIter(:) ! volumetric liquid water content in each soil layer (-)
 real(dp),intent(in)        :: mLayerVolFracIceIter(:) ! volumetric ice content in each soil layer (-)
 ! output
 real(dp),intent(out)       :: waterReq2FillPore(0:)   ! water required to fill pore space at each layer interface (m)
 integer(i4b),intent(out)   :: err                     ! error code
 character(*),intent(out)   :: message                 ! error message
 ! internal
 integer(i4b)               :: iSoil                   ! index of soil layers
 integer(i4b)               :: nSoil                   ! number of soil layers
 real(dp)                   :: drainablePorosity       ! drainable porosity (-)
 ! initialize error control
 err=0; message='waterReq2FillLayer/'

 ! get the number of soil layers
 nSoil = size(mLayerVolFracLiqIter)

 ! initialize the water required to fill pore space (m)
 waterReq2FillPore(nSoil) = 0._dp ! (no water required to fill pore space at the bottom interface)

 ! compute the water required to fill pore space at each layer interface (m)
 do iSoil=nSoil,1,-1  ! (loop through soil layers, starting at the bottom)
  drainablePorosity = theta_sat - mLayerVolFracLiqIter(iSoil) - mLayerVolFracIceIter(iSoil)
  waterReq2FillPore(iSoil-1) = waterReq2FillPore(iSoil) + drainablePorosity*mLayerDepth(iSoil)
 end do

 end subroutine waterReq2FillLayer


 ! ---------------------------------------------------------------------------------------------------------------------------
 ! private subroutine: given aquifer storage,
 !                     compute the number of unsaturated layers, the drainable porosity, and the depth to the water table
 ! ---------------------------------------------------------------------------------------------------------------------------
 subroutine waterTablePosition(mLayerVolFracLiqIter, & ! intent(in):  volumetric liquid water content in each soil layer (-)
                               mLayerVolFracIceIter, & ! intent(in):  volumetric ice content in each soil layer (-)
                               aquiferStorageTrial,  & ! intent(in):  aquifer storage (m)
                               waterReq2FillPore,    & ! intent(in):  water required to fill pore space up to the level of each layer interface (m)
                               iLayerHeight,         & ! intent(in):  height of each interface (m)
                               mLayerDepth,          & ! intent(in):  depth of each soil layer (m)
                               theta_sat,            & ! intent(in):  soil porosity (-)
                               specificYield,        & ! intent(in):  fraction of water volume drained by gravity in an unconfined aquifer (-)
                               nUnsat,               & ! intent(out): number of unsaturated layers
                               drainablePorosity,    & ! intent(out): drainable porosity (-)
                               waterTableDepth,      & ! intent(out): water table depth (m)
                               err,message)            ! intent(out): error control
 implicit none
 ! input
 real(dp),intent(in)        :: mLayerVolFracLiqIter(:) ! volumetric liquid water content in each soil layer (-)
 real(dp),intent(in)        :: mLayerVolFracIceIter(:) ! volumetric ice content in each soil layer (-)
 real(dp),intent(in)        :: aquiferStorageTrial     ! aquifer storage (m)
 real(dp),intent(in)        :: waterReq2FillPore(0:)   ! water required to fill pore space up to the level of each layer interface (m)
 real(dp),intent(in)        :: iLayerHeight(0:)        ! height of each interface (m)
 real(dp),intent(in)        :: mLayerDepth(:)          ! depth of each soil layer (m)
 real(dp),intent(in)        :: theta_sat               ! soil porosity (-)
 real(dp),intent(in)        :: specificYield           ! fraction of water volume drained by gravity in an unconfined aquifer (-)
 ! output
 integer(i4b),intent(out)   :: nUnsat                  ! number of unsaturated layers
 real(dp),intent(out)       :: drainablePorosity       ! drainable porosity (-)
 real(dp),intent(out)       :: waterTableDepth         ! water table depth (m)
 integer(i4b),intent(out)   :: err                     ! error code
 character(*),intent(out)   :: message                 ! error message
 ! local variables
 integer(i4b)               :: nSoil                   ! number of soil layers
 real(dp)                   :: fracFilled              ! fraction layer filled with water (-)
 ! initialize error control 
 err=0; message='waterTablePosition/'
 ! get the number of soil layers
 nSoil = size(mLayerDepth)
 ! ***** water table below the soil column (note: waterReq2FillPore(nSoil) = zero)
 if(aquiferStorageTrial < waterReq2FillPore(nSoil))then
  nUnsat            = nSoil
  drainablePorosity = specificYield
  waterTableDepth   = iLayerHeight(nSoil) - aquiferStorageTrial/specificYield
 ! ***** water table within the soil column
 else
  ! check that the water table is below the surface
  if(aquiferStorageTrial > waterReq2FillPore(0))then; err=20; message=trim(message)//'entire soil column is saturated'; return; endif
  ! get the number of unsaturated layers
  nUnsat = count(aquiferStorageTrial < waterReq2FillPore(1:nSoil))
  if(nUnsat >= nSoil)then; err=20; message=trim(message)//'expect some soil layers to be saturated'; return; endif
  ! assign drainable porsity for the layer that contains the water table
  drainablePorosity = theta_sat - mLayerVolFracLiqIter(nUnsat+1) - mLayerVolFracIceIter(nUnsat+1)
  ! compute depth to the water table
  fracFilled      = (aquiferStorageTrial - waterReq2FillPore(nUnsat+1)) / (waterReq2FillPore(nUnsat) - waterReq2FillPore(nUnsat+1))  ! fraction layer filled with water
  waterTableDepth = iLayerHeight(nUnsat+1) + fracFilled*mLayerDepth(nUnsat+1)
 endif
 end subroutine waterTablePosition


end module groundwatr_module

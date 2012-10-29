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
                       dt,                                              & ! input:  time step (seconds) 
                       ! input (state variables)
                       scalarAquiferStorage,                            & ! input: aquifer storage at the start of the time step (m)
                       scalarAquiferStorageTrial,                       & ! input: trial value of aquifer storage (m)
                       ! input (fluxes)
                       scalarAquiferRcharge,                            & ! input: aquifer recharge (m s-1)
                       scalarAquiferTranspire,                          & ! input: aquifer transpiration (m s-1)
                       ! input (available depth)
                       mLayerVolFracLiqIter,                            & ! input: volumetric fraction of liquid water after itertations (-)
                       mLayerVolFracIceIter,                            & ! input: volumetric fraction of ice after itertations (-)
                       iLayerHeight,                                    & ! input: height of each interface (m), NOTE: include height of soil (iLayer=nSnow)
                       mLayerDepth,                                     & ! input: depth of each soil layer (m)
                       ! input (diagnostic variables)
                       k_surf,                                          & ! input: hydraulic conductivity at the surface (m s-1)
                       ! input (parameters)
                       theta_sat,                                       & ! input: soil porosity (-)
                       specificYield,                                   & ! input: specific yield (-)
                       kAnisotropic,                                    & ! input: anisotropy factor for lateral hydraulic conductivity (-)
                       zScale_TOPMODEL,                                 & ! input: scale factor for TOPMODEL-ish baseflow parameterization (m)
                       ! output (states)
                       scalarAquiferStorageNew,                         & ! output: aquifer storage at the end of the time step (m)
                       ! output (diagnostic variables and error control)
                       scalarWaterTableDepth,                           & ! output: water table depth at the start/end of the time step (m)
                       scalarAquiferBaseflow,                           & ! output: baseflow from the aquifer (m s-1)
                       err,message)                                       ! output: error control
 ! compute change in aquifer storage over the time step
 implicit none
 ! -----------------------------------------------------------------------------------------------------------------------------------------
 ! input (time step)
 real(dp),intent(in)                 :: dt                        ! time step (seconds)
 ! input (state variables)
 real(dp),intent(in)                 :: scalarAquiferStorage      ! aquifer storage at the start of the time step (m)
 real(dp),intent(in)                 :: scalarAquiferStorageTrial ! trial value of aquifer storage (m)
 ! input (fluxes)
 real(dp),intent(in)                 :: scalarAquiferRcharge      ! aquifer recharge averaged over the time step (m s-1)
 real(dp),intent(in)                 :: scalarAquiferTranspire    ! aquifer transpiration averaged over the time step (m s-1)
 ! input (available depth)
 real(dp),intent(in)                 :: mLayerVolFracLiqIter(:)   ! volumetric fraction of liquid water after itertations (-)
 real(dp),intent(in)                 :: mLayerVolFracIceIter(:)   ! volumetric fraction of ice after itertations (-)
 real(dp),intent(in)                 :: iLayerHeight(0:)          ! height at each layer interface (m)
 real(dp),intent(in)                 :: mLayerDepth(:)            ! depth of each soil layer (m)
 ! input (diagnostic variables)
 real(dp),intent(in)                 :: k_surf                    ! hydraulic conductivity at the surface (m s-1) 
 ! input (parameters)
 real(dp),intent(in)                 :: theta_sat                 ! total porosity (-)
 real(dp),intent(in)                 :: specificYield             ! fraction of water volume drained by gravity in an unconfined aquifer (-)
 real(dp),intent(in)                 :: kAnisotropic              ! anisotropy factor for lateral hydraulic conductivity (-)
 real(dp),intent(in)                 :: zScale_TOPMODEL           ! scale factor for TOPMODEL-ish baseflow parameterization (m)
 ! output (states)
 real(dp),intent(out)                :: scalarAquiferStorageNew   ! aquifer storage at the end of the time step (m)
 ! output (diagnostic variables and error control)
 real(dp),intent(out)                :: scalarWaterTableDepth     ! water table depth at the start/end of the time step (m)
 real(dp),intent(out)                :: scalarAquiferBaseflow     ! baseflow from the aquifer (m s-1)
 integer(i4b),intent(out)            :: err                       ! error code
 character(*),intent(out)            :: message                   ! error message
 ! -----------------------------------------------------------------------------------------------------------------------------------------
 ! algorithmic control parameters
 integer(i4b),parameter              :: maxiter=10                ! maximum number of iterations
 real(dp),parameter                  :: incTol= 1.e-6_dp          ! convergence tolerance for the iteration increment (m)
 real(dp),parameter                  :: resTol= 1.e-6_dp          ! convergence tolerance for the residual (m)
 ! define local variables
 character(len=256)                  :: cmessage                  ! error message for downwind routine
 integer(i4b)                        :: iter                      ! iteration index
 integer(i4b)                        :: iSoil                     ! index of each soil layer
 integer(i4b)                        :: nSoil                     ! number of soil layers
 integer(i4b)                        :: nUnsat                    ! number of unsaturated layers
 integer(i4b)                        :: nUnsat_init               ! number of unsaturated layers at the start of the iterations
 real(dp)                            :: baseflowMax               ! maximum baseflow rate from the aquifer (m s-1)
 real(dp),dimension(0:size(mLayerVolFracLiqIter)) :: waterReq2FillPore  ! water required to fill pore space (m)
 real(dp)                            :: drainablePorosity         ! drainable porosity (-)
 real(dp)                            :: aquiferStorageTrial       ! trial value of aquifer storage (m) -- local variable to preserve input
 real(dp)                            :: zWaterTrial               ! depth of the water table (m)
 real(dp)                            :: dBaseflow_dStorage        ! derivative in the baseflow term w.r.t. aquifer storage (s-1)
 real(dp)                            :: netFlux                   ! recharge minus baseflow (m s-1)
 real(dp)                            :: residual                  ! residual in aquifer storage (m)
 real(dp)                            :: dStorage                  ! iteration increment for aquifer storage (m)
 integer(i4b)                        :: idiff                     ! trial integers used to find the bracket
 integer(i4b),parameter              :: maxdiff=100               ! maximum trial values to find the bracket
 real(dp),parameter                  :: xdiff=0.00001_dp          ! increments used to find the bracket
 real(dp)                            :: xbeg,xend                 ! lower and upper points in the bracket
 real(dp)                            :: rbeg,rend                 ! residuals for the lower and upper points in the bracket
 integer(i4b)                        :: itest                     ! trial integers used to test along a line
 integer(i4b),parameter              :: ntest=100000              ! number of values to test along a line
 real(dp)                            :: xinc                      ! increment in storage to test along a line
 real(dp)                            :: xtry,rtry                 ! trial value and residual
 real(dp)                            :: zWaterTrial_init          ! initial value of the water table depth

 ! initialize error control
 err=0; message="groundwatr/"

 ! get the number of soil layers
 nSoil = size(mLayerDepth)

 ! initialize the auifer storage
 aquiferStorageTrial = scalarAquiferStorageTrial

 ! compute the maximum baseflow rate
 baseflowMax = kAnisotropic*k_surf

 ! compute the water required to fill pore space up to the level of each layer interface (m)
 call waterReq2FillLayer(theta_sat,                  & ! intent(in):  total porosity (-)
                         mLayerDepth,                & ! intent(in):  depth of each soil layer (m)
                         mLayerVolFracLiqIter,       & ! intent(in):  volumetric liquid water content in each soil layer (-)
                         mLayerVolFracIceIter,       & ! intent(in):  volumetric ice content in each soil layer (-)
                         waterReq2FillPore,          & ! intent(out): water required to fill pore space at each layer interface (m)
                         err,cmessage)                 ! intent(out): error control
 if(err/=0)then; message=trim(message)//trim(cmessage); return; endif

 ! compute the number of unsaturated layers, the drainable porosity, and the depth to the water table
 call waterTablePosition(mLayerVolFracLiqIter,       & ! intent(in):  volumetric liquid water content in each soil layer (-)
                         mLayerVolFracIceIter,       & ! intent(in):  volumetric ice content in each soil layer (-)
                         aquiferStorageTrial,        & ! intent(in):  aquifer storage (m)
                         waterReq2FillPore,          & ! intent(in):  water required to fill pore space up to the level of each layer interface (m(
                         iLayerHeight,               & ! intent(in):  height of each interface (m)
                         mLayerDepth,                & ! intent(in):  depth of each soil layer (m)
                         theta_sat,                  & ! intent(in):  soil porosity (-)
                         specificYield,              & ! intent(in):  fraction of water volume drained by gravity in an unconfined aquifer (-)
                         nUnsat_init,                & ! intent(out): number of unsaturated layers at the start of the iterations
                         drainablePorosity,          & ! intent(out): drainable porosity (-)
                         zWaterTrial,                & ! intent(out): water table depth (m)
                         err,cmessage)                 ! intent(out): error control
 if(err/=0)then; message=trim(message)//trim(cmessage); return; endif

 ! save the inital value of water table depth
 zWaterTrial_init = zWaterTrial

 ! *****
 ! * start iterations...
 ! *********************
 do iter=1,maxiter

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
  ! print progress
  if(iter == 1) write(*,'(a)')    'iter, aquiferStorageTrial, zWaterTrial, scalarAquiferBaseflow, scalarAquiferRcharge, scalarAquiferTranspire, dBaseflow_dStorage, residual, dStorage = '
  write(*,'(i4,1x,10(e20.10,1x))') iter, aquiferStorageTrial, zWaterTrial, scalarAquiferBaseflow, scalarAquiferRcharge, scalarAquiferTranspire, dBaseflow_dStorage, residual, dStorage

  ! ** update the aquifer storage (m) 
  aquiferStorageTrial = aquiferStorageTrial + dStorage

  ! ** update the number of unsaturated layers, the drainable porosity, and the depth to the water table
  call waterTablePosition(mLayerVolFracLiqIter,       & ! intent(in):  volumetric liquid water content in each soil layer (-)
                          mLayerVolFracIceIter,       & ! intent(in):  volumetric ice content in each soil layer (-)
                          aquiferStorageTrial,        & ! intent(in):  aquifer storage (m)
                          waterReq2FillPore,          & ! intent(in):  water required to fill pore space up to the level of each layer interface (m(
                          iLayerHeight,               & ! intent(in):  height of each interface (m)
                          mLayerDepth,                & ! intent(in):  depth of each soil layer (m)
                          theta_sat,                  & ! intent(in):  soil porosity (-)
                          specificYield,              & ! intent(in):  fraction of water volume drained by gravity in an unconfined aquifer (-)
                          nUnsat,                     & ! intent(out): number of unsaturated layers
                          drainablePorosity,          & ! intent(out): drainable porosity (-)
                          zWaterTrial,                & ! intent(out): water table depth (m)
                          err,cmessage)                 ! intent(out): error control
  if(err/=0)then; message=trim(message)//trim(cmessage); return; endif

  ! ** check for convergence and return desired quantities
  if(abs(dStorage) < incTol .and. abs(residual) < resTol)then
   ! (test that the rise/rall in the water table depth is tolerable)
   if(nUnsat_init < nSoil)then
    if(abs(zWaterTrial - zWaterTrial_init) > mLayerDepth(nUnsat_init+1)/2._dp)then
     write(message,'(a,2(f12.5,a))')trim(message)//'rise/fall in water table is intolerable [water table depth = ', zWaterTrial,&
                                    '; height of recharge point = ', iLayerHeight(nUnsat_init) - mLayerDepth(nUnsat_init)/2._dp, ']'
     err=-20; return  ! will cause a reduction in the time step and another trial
    endif  ! (if rise/fall in water table is intolerable)
   endif  ! (if water table is within the soil profile)
   ! update aquifer storage (m)
   scalarAquiferStorageNew = aquiferStorageTrial
   scalarWaterTableDepth   = zWaterTrial
   exit
  endif  ! convergence check
  ! check for non-convergence
  if(iter==maxiter)then
   ! return with a negative error
   ! NOTE: testing code below will not be activated
   message=trim(message)//'failed to converge'
   err=-20; return ! negative error forced time step reduction and another trial
   ! get brackets
   do idiff=1,maxdiff
    ! (values to test as brackets)
    xbeg = aquiferStorageTrial - xdiff*real(idiff,kind(xdiff))
    xend = aquiferStorageTrial + xdiff*real(idiff,kind(xdiff))
    ! (corresponding residuals)
    rbeg = resid(xbeg)
    rend = resid(xend)
    ! (check if we have a bracket)
    write(*,'(a,i4,1x,4(f13.5,1x))') 'checking for a bracket ', idiff, xbeg, xend, rbeg, rend
    if(rbeg*rend < 0._dp) exit
    if(idiff==maxdiff)then; err=20; message=trim(message)//'insufficient width given by maxdiff'; return; endif
   end do
   ! evaluate along a line
   xinc = (xend - xbeg)/real(ntest,kind(xinc))
   do itest=56600,56700
    xtry = xbeg + xinc*real(itest,kind(xinc))
    rtry = resid(xtry)
    write(*,'(a,i8,1x,f13.9,1x,e20.10)') 'checking for a root ', itest, xtry, rtry
   end do
   ! return convergence error
   message=trim(message)//'failed to converge'
   err=20; return
  endif
 end do  ! (loop through iterations)

 ! clean-up fluxes
 scalarAquiferBaseflow = scalarAquiferBaseflow + dBaseflow_dStorage*dStorage




 contains
 ! ---------------------------------------------------------------------------------------------------------------------------
 ! internal function: compute the residual given a value of aquifer storage (m)
 ! ---------------------------------------------------------------------------------------------------------------------------
 function resid(aquiferStorageTrial)
 implicit none
 ! input/output
 real(dp),intent(in)     :: aquiferStorageTrial  ! aquifer storage (m)
 real(dp)                :: resid                ! function value
 ! calculate the depth to the water table
 call waterTablePosition(mLayerVolFracLiqIter,  & ! intent(in):  volumetric liquid water content in each soil layer (-)
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
 if(err/=0)then; print*, 'in function resid... stopping '//trim(message)//trim(cmessage); stop; endif
 ! calculate the baseflow term (m s-1)
 scalarAquiferBaseflow = baseflowMax*exp(-zWaterTrial/zScale_TOPMODEL)
 ! compute the net flux (m s-1)
 netFlux = (scalarAquiferRcharge + scalarAquiferTranspire) - scalarAquiferBaseflow
 ! compute the residual in aquifer storage (m)
 resid   = aquiferStorageTrial - (scalarAquiferStorage + netFlux*dt)
 print*, 'nUnsat, zWaterTrial, scalarAquiferBaseflow, netFlux = ', nUnsat, zWaterTrial, scalarAquiferBaseflow, netFlux
 print*, 'aquiferStorageTrial = ', aquiferStorageTrial
 print*, 'waterReq2FillPore(:) = ', waterReq2FillPore(:)
 print*, 'iLayerHeight(nUnsat:nUnsat+1) = ', iLayerHeight(nUnsat:nUnsat+1)
 end function resid

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
  waterTableDepth = iLayerHeight(nUnsat+1) - fracFilled*mLayerDepth(nUnsat+1)
 endif
 end subroutine waterTablePosition


end module groundwatr_module

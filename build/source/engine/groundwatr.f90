module groundwatr_module
USE nrtype
implicit none
private
public::groundwatr
contains

 ! ************************************************************************************************
 ! new subroutine: compute change in snow density over the time step
 ! ************************************************************************************************
 subroutine groundwatr(&
                       ! input
                       dt,                                   &  ! input:  time step (seconds) 
                       scalarAquiferRcharge,                 &  ! input:  aquifer recharge (m s-1)
                       mLayerVolFracLiqIter,                 &  ! input:  volumetric fraction of liquid water after itertations (-)
                       mLayerVolFracIceIter,                 &  ! input:  volumetric fraction of ice after itertations (-)
                       ! input-output
                       scalarAquiferStorage,                 &  ! input-output: aquifer storage (m)
                       ! output
                       scalarWaterTableDepth,                &  ! output: water table depth at the end of the time step (m)
                       err,message)                             ! output: error control
 USE data_struc,only:mpar_meta,forc_meta,mvar_meta,indx_meta                    ! metadata
 USE data_struc,only:mpar_data,forc_data,mvar_data,indx_data,ix_soil,ix_snow    ! data structures
 USE var_lookup,only:iLookPARAM,iLookFORCE,iLookMVAR,iLookINDEX                 ! named variables for structure elements
 ! compute change in snow density over the time step
 implicit none
 ! input variables
 real(dp),intent(in)                 :: dt                       ! time step (seconds)
 real(dp),intent(in)                 :: scalarAquiferRcharge     ! aquifer recharge averaged over the time step (m s-1)
 real(dp),intent(in)                 :: mLayerVolFracLiqIter(:)  ! volumetric fraction of liquid water after itertations (-)
 real(dp),intent(in)                 :: mLayerVolFracIceIter(:)  ! volumetric fraction of ice after itertations (-)
 ! input-output variables
 real(dp),intent(inout)              :: scalarAquiferStorage     ! aquifer storage (m)
 ! output variables
 real(dp),intent(out)                :: scalarWaterTableDepth    ! water table depth at the end of the time step (m)
 integer(i4b),intent(out)            :: err                      ! error code
 character(*),intent(out)            :: message                  ! error message
 ! local pointers to soil parameters
 real(dp),pointer                    :: k_soil                   ! hydraulic conductivity (m s-1)
 real(dp),pointer                    :: kAnisotropic             ! anisotropy factor for lateral hydraulic conductivity (-) 
 real(dp),pointer                    :: zScale_TOPMODEL          ! scale factor for TOPMODEL-ish baseflow parameterization (m)
 real(dp),pointer                    :: theta_sat                ! soil porosity (-)
 real(dp),pointer                    :: specificYield            ! specific yield (-)
 ! algorithmic control parameters
 real(dp),pointer                    :: wimplicit                ! weight assigned to start-of-step fluxes (-)
 integer(i4b),parameter              :: maxiter=10               ! maximum number of iterations
 real(dp),parameter                  :: incTol= 1.e-6_dp         ! convergence tolerance for the iteration increment (m)
 real(dp),parameter                  :: resTol= 1.e-6_dp         ! convergence tolerance for the residual (m)
 ! coordinate variables
 integer(i4b)                        :: nSnow                    ! number of snow layers
 integer(i4b)                        :: nSoil                    ! number of soil layers
 integer(i4b),pointer                :: nLayers                  ! total number of layers
 real(dp),pointer                    :: mLayerDepth(:)           ! depth of each layer (m)
 real(dp),pointer                    :: iLayerHeight(:)          ! height of layer interfaces (m)
 ! define local variables
 character(len=256)                  :: cmessage                 ! error message for downwind routine
 integer(i4b)                        :: iter                     ! iteration index
 integer(i4b)                        :: iSoil                    ! index of each soil layer
 integer(i4b)                        :: nUnsat                   ! number of unsaturated layers
 real(dp)                            :: baseflowMax              ! maximum baseflow rate from the aquifer (m s-1)
 real(dp)                            :: fracFilled               ! fraction layer filled with water (-)
 real(dp),allocatable                :: waterReq2FillPore(:)     ! water required to fill pore space up to each layer interface (m)
 real(dp)                            :: drainablePorosity        ! drainable porosity (-)
 real(dp)                            :: aquiferStorageTrial      ! trial value for aquifer storage (m)
 real(dp)                            :: zWaterTrial              ! depth of the water table (m)
 real(dp)                            :: aquiferBaseflow          ! baseflow from the aquifer (m s-1)
 real(dp)                            :: dBaseflow_dStorage       ! derivative in the baseflow term w.r.t. aquifer storage (s-1)
 real(dp)                            :: netFlux                  ! recharge minus baseflow (m s-1)
 real(dp)                            :: residual                 ! residual in aquifer storage (m)
 real(dp)                            :: dStorage                 ! iteration increment for aquifer storage (m)

 ! initialize error control
 err=0; message="groundwatr/"

 ! assign pointers to soil parameters
 k_soil          => mpar_data%var(iLookPARAM%k_soil)             ! hydraulic conductivity (m s-1)
 kAnisotropic    => mpar_data%var(iLookPARAM%kAnisotropic)       ! anisotropy factor for lateral hydraulic conductivity (-)
 zScale_TOPMODEL => mpar_data%var(iLookPARAM%zScale_TOPMODEL)    ! scale factor for TOPMODEL-ish baseflow parameterization (m)
 theta_sat       => mpar_data%var(iLookPARAM%theta_sat)          ! soil porosity (-)
 specificYield   => mpar_data%var(iLookPARAM%specificYield)      ! specific yield (-)

 ! get the number of layers
 nLayers   => indx_data%var(iLookINDEX%nLayers)%dat(1)           ! total number of layers

 ! identify the number of snow and soil layers
 nSnow = count(indx_data%var(iLookINDEX%layerType)%dat == ix_snow)
 nSoil = count(indx_data%var(iLookINDEX%layerType)%dat == ix_soil)

 ! assign pointers to the height and depth vectors
 mLayerDepth  => mvar_data%var(iLookMVAR%mLayerDepth)%dat(nSnow+1:nLayers)
 iLayerHeight => mvar_data%var(iLookMVAR%iLayerHeight)%dat(nSnow+1:nLayers)

 ! compute the maximum baseflow rate
 baseflowMax = kAnisotropic*k_soil

 ! initialize the aquifer storage (m)
 aquiferStorageTrial = scalarAquiferStorage

 ! initialize the water required to fill pore space (m)
 allocate(waterReq2FillPore(0:nSoil),stat=err)
 if(err/=0)then; err=20; message=trim(message)//'problem allocating space for "waterReq2FillPore"'; return; endif
 waterReq2FillPore(nSoil) = 0._dp ! (no water equired to fill pore space at the bottom interface)

 ! compute the water required to fill pore space at each layer interface (m)
 do iSoil=nSoil-1,0,-1  ! (loop through soil layers, starting at the bottom)
  drainablePorosity = theta_sat - mLayerVolFracLiqIter(iSoil) - mLayerVolFracIceIter(iSoil)
  waterReq2FillPore(iSoil) = waterReq2FillPore(iSoil+1) + drainablePorosity*mLayerDepth(iSoil)
 end do

 ! *****
 ! * start iterations...
 ! *********************
 do iter=1,maxiter

  ! -----
  ! compute the depth to the water table and the drainable porosity...
  ! ---------------------------------------
  call waterTablePosition(aquiferStorageTrial,& ! intent(in):  aquifer storage (m)
                          nUnsat,             & ! intent(out): number of unsaturated layers
                          drainablePorosity,  & ! intent(out): drainable porosity (-)
                          zWaterTrial,        & ! intent(out): water table depth (m)
                          err,cmessage)         ! intent(out): error control
  if(err/=0)then; message=trim(message)//trim(cmessage); return; endif

  ! -----
  ! compute iteration increment...
  ! ------------------------------
  ! calculate the baseflow term (m s-1)
  aquiferBaseflow = baseflowMax*exp(-zWaterTrial/zScale_TOPMODEL)
  ! compute the net flux (m s-1)
  netFlux = scalarAquiferRcharge - aquiferBaseflow
  ! compute the residual in aquifer storage (m)
  residual = aquiferStorageTrial - (scalarAquiferStorage + netFlux*dt)
  ! compute the derivative in baseflow w.r.t. storage (s-1)
  dBaseflow_dStorage = aquiferBaseflow/(zScale_TOPMODEL*drainablePorosity)
  ! compute the iteration increment (m)
  dStorage = -residual/(1._dp + dBaseflow_dStorage*dt)
 
  ! check for convergence and update states
  if(abs(dStorage) < incTol .and. abs(residual) < resTol)then
   scalarAquiferStorage  = aquiferStorageTrial + dStorage
   call waterTablePosition(scalarAquiferStorage, & ! intent(in):  aquifer storage (m)
                           nUnsat,               & ! intent(out): number of unsaturated layers
                           drainablePorosity,    & ! intent(out): drainable porosity (-)
                           scalarWaterTableDepth,& ! intent(out): water table depth (m)
                           err,cmessage)           ! intent(out): error control
   if(err/=0)then; message=trim(message)//trim(cmessage); return; endif
   print*, 'baseflowMax, zScale_TOPMODEL = ', baseflowMax, zScale_TOPMODEL
   write(*,'(a)')            'scalarAquiferStorage, scalarWaterTableDepth, aquiferBaseflow, scalarAquiferRcharge, dBaseflow_dStorage, dStorage = '
   write(*,'(10(e20.10,1x))') scalarAquiferStorage, scalarWaterTableDepth, aquiferBaseflow, scalarAquiferRcharge, dBaseflow_dStorage, dStorage
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
 aquiferBaseflow = aquiferBaseflow + dBaseflow_dStorage*dStorage

 deallocate(waterReq2FillPore,stat=err)
 if(err/=0)then; err=20; message=trim(message)//'problem deallocating space for "waterReq2FillPore"'; return; endif


 contains

  ! ---------------------------------------------------------------------------------------------------------------------------
  ! internal subroutine: given aquifer storage,
  !                      compute the number of unsaturated layers, the drainable porosity, and the depth to the water table
  ! ---------------------------------------------------------------------------------------------------------------------------
  subroutine waterTablePosition(aquiferStorageTrial,& ! intent(in):  aquifer storage (m)
                                nUnsat,             & ! intent(out): number of unsaturated layers
                                drainablePorosity,  & ! intent(out): drainable porosity (-)
                                waterTableDepth,    & ! intent(out): water table depth (m)
                                err,message)          ! intent(out): error control
  implicit none
  ! input
  real(dp),intent(in)        :: aquiferStorageTrial   ! aquifer storage (m)
  ! output
  integer(i4b),intent(out)   :: nUnsat                ! number of unsaturated layers
  real(dp),intent(out)       :: drainablePorosity     ! drainable porosity (-)
  real(dp),intent(out)       :: waterTableDepth       ! water table depth (m)
  integer(i4b),intent(out)   :: err                   ! error code
  character(*),intent(out)   :: message               ! error message
  ! initialize error control
  err=0; message='waterTablePosition/'
  ! ***** water table below the soil column (note: waterReq2FillPore(nSoil) = zero)
  if(aquiferStorageTrial < waterReq2FillPore(nSoil))then
   nUnsat            = nSoil
   drainablePorosity = specificYield
   zWaterTrial       = iLayerHeight(nSoil) - aquiferStorageTrial/specificYield
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
   fracFilled  = (aquiferStorageTrial - waterReq2FillPore(nUnsat+1)) / (waterReq2FillPore(nUnsat) - waterReq2FillPore(nUnsat+1))  ! fraction layer filled with water
   zWaterTrial = iLayerHeight(nUnsat+1) + fracFilled*mLayerDepth(nUnsat+1)
   print*, 'zWaterTrial = ', zWaterTrial
   pause
  endif
  end subroutine waterTablePosition

 end subroutine groundwatr

end module groundwatr_module

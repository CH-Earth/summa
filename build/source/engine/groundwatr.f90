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
                       mLayerMatricHeadIter,                 &  ! input:  matric head (m)
                       mLayerVolFracLiqIter,                 &  ! input:  volumetric fraction of liquid water after itertations (-)
                       mLayerVolFracIceIter,                 &  ! input:  volumetric fraction of ice after itertations (-)
                       ! input-output
                       scalarAquiferStorage,                 &  ! input-output: aquifer storage (m)
                       ! output
                       mLayerMatricHeadNew,                  &  ! output: matric head (m)
                       mLayerVolFracLiqNew,                  &  ! output: volumetric fraction of liquid water after itertations (-)
                       mLayerVolFracIceNew,                  &  ! output: volumetric fraction of ice after itertations (-)
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
 real(dp),intent(in)                 :: mLayerVolFracLiqIter(:)  ! volumetric fraction of liquid water after iterations (-)
 real(dp),intent(in)                 :: mLayerVolFracIceIter(:)  ! volumetric fraction of ice after iterations (-)
 real(dp),intent(in)                 :: mLayerMatricHeadIter(:)  ! matric head after iterations (m)
 ! input-output variables
 real(dp),intent(inout)              :: scalarAquiferStorage     ! aquifer storage (m)
 ! output variables
 real(dp),intent(out)                :: mLayerVolFracLiqNew(:)   ! new volumetric fraction of liquid water (-)
 real(dp),intent(out)                :: mLayerVolFracIceNew(:)   ! new volumetric fraction of ice (-)
 real(dp),intent(out)                :: mLayerMatricHeadNew(:)   ! new matric head after (m)
 real(dp),intent(out)                :: scalarWaterTableDepth    ! water table depth at the end of the time step (m)
 integer(i4b),intent(out)            :: err                      ! error code
 character(*),intent(out)            :: message                  ! error message
 ! local pointers to soil parameters
 real(dp),pointer                    :: k_soil                   ! hydraulic conductivity (m s-1)
 real(dp),pointer                    :: kAnisotropic             ! anisotropy factor for lateral hydraulic conductivity (-) 
 real(dp),pointer                    :: zScale_TOPMODEL          ! scale factor for TOPMODEL-ish baseflow parameterization (m)
 real(dp),pointer                    :: theta_sat                ! soil porosity (-)
 real(dp),pointer                    :: specificYield            ! specific yield (-)
 ! local pointers to algorithmic control parameters
 real(dp),pointer                    :: wimplicit                ! weight assigned to start-of-step fluxes (-)
 ! local pointers to coordinate variables
 integer(i4b),pointer                :: nLayers                  ! total number of layers
 real(dp),pointer                    :: soilDepth                ! soil depth
 ! define local variables
 character(len=256)                  :: cmessage                 ! error message for downwind routine
 real(dp)                            :: baseflowMax              ! maximum baseflow rate from the aquifer (m s-1)
 real(dp)                            :: aquiferStorageTrial      ! trial value for aquifer storage (m)
 integer(i4b)                        :: iter                     ! iteration index
 integer(i4b),parameter              :: maxiter=10               ! maximum number of iterations
 real(dp)                            :: zWater                   ! depth of the water table (m)
 real(dp)                            :: aquiferBaseflow          ! baseflow from the aquifer (m s-1)
 real(dp)                            :: dBaseflow_dStorage       ! derivative in the baseflow term w.r.t. aquifer storage (s-1)
 real(dp)                            :: residual                 ! residual in aquifer storage (m)
 real(dp)                            :: dStorage                 ! iteration increment for aquifer storage (m)
 real(dp),parameter                  :: incTol= 1.e-6_dp         ! convergence tolerance for the iteration increment (m)
 real(dp),parameter                  :: resTol= 1.e-6_dp         ! convergence tolerance for the residual (m)

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

 ! assign pointers to soil depth
 soilDepth => mvar_data%var(iLookMVAR%iLayerHeight)%dat(nLayers) ! soil depth

 ! assign pointers to the height vector
 layerHeight => mvar_data%var(iLookMVAR%iLayerHeight)%dat(nSnow+1:nLayers)

 ! identify the maximum and minimum height
 zmin = layerHeight(nUnsat)
 if(nUnsat < nSoil)then; zmax = layerHeight(nUnsat+1)
 else;                   zmax = huge(soilDepth)
 endif

 ! compute the maximum baseflow rate
 baseflowMax = kAnisotropic*k_soil

 ! compute the drainable porosity
 if(nUnsat == nSoil)  ! (below the soil column)
  dPorosity = specificYield
 else                 ! (within the soil column)
  dPorosity = theta_sat - mLayerVolFracLiqIter(nUnsat+1) - mLayerVolFracIceIter(nUnsat+1)
 endif

 ! *****
 ! * start iterations...
 ! *********************
 do iter=1,maxiter

  ! -----
  ! compute iteration increment...
  ! ------------------------------
  ! calculate the baseflow term (m s-1)
  aquiferBaseflow = baseflowMax*exp(-zWaterTrial/zScale_TOPMODEL)
  ! compute the net flux (m s-1)
  netFlux = scalarAquiferRcharge - aquiferBaseflow
  ! compute the residual in the depth to the water table (-)
  residual = zWaterTrial - (zWaterInit - dt*netFlux/dPorosity)
  ! compute the derivative in the baseflow flux w.r.t. depth to the water table (s-1)
  dBaseflow_dWaterTable = -aquiferBaseflow/zScale_TOPMODEL
  ! compute the iteration increment in the depth to the water table (m)
  dWaterTable = residual/(-1._dp + dt*dBaseflow_dWaterTable/dPorosity)

  ! -----
  ! update water table depth...
  ! ---------------------------
  ! compute a temporary value for testing
  zTemp = zWaterTrial + dWaterTable

  ! identify movement of water table across layers
  if    (zTemp  > zMin .and. zTemp < zMax)then; iWaterCross = noCrossing ! water table does not cross a layer
  elseif(zTemp <= zMin)then;                    iWaterCross = iWaterRise ! water table rises into the layer above
  elseif(zTemp >= zMax)then;                    iWaterCross = iWaterFall ! water table falls into the layer below
  else; err=20; message=trim(message)//'cannot identify the type of water movement'; return; endif

  ! handle various cases associated with water table corssing soil layers
  select case(iWaterCross)

   ! ***** most common case, where water table does not cross a layer
   case(noCrossing)
    zWaterTrial = zTemp

   ! ***** case where water table rises into the next layer
   case(iWaterRise)   
    ! (update water table depth and liquid water content)
    zWaterTrial = layerHeight(nUnsat)  ! stop on the layer interface
    mLayerVolFracLiqNew(nUnsat) = theta_sat - mLayerVolFracIceIter(nUnsat) 
    ! (update number of unsaturated layers)
    nUnsat = nUnsat-1                  ! reduce the number of unsaturated layers
    if(nUnsat == 0)then; err=20; message=trim(message)//'soil column is completely saturated!'
    ! (update drainable porisity)
    dPorosity = theta_sat - mLayerVolFracLiqIter(nUnsat+1) - mLayerVolFracIceIter(nUnsat+1)


   ! ***** case where water table falls into the layer below
   case(iWaterFall)



 

    ! (get an initial estimate of the water table)
    zTemp = zWater + dStorage/fracAir

    ! (check that we did not cross a layer)
    do iAbsorb=1,maxAbsorb  ! (only absorb a certain number of layers [usually 2])
     if(zTemp < layerHeight(nUnsat)then
      err=20; message=trim(message)//'have not implemented specific yield within the soil column yet'; return; endif
     else
      zWater = zTemp
      exit
     endif
    end do

   endif  ! ** updating the water table depth


  ! check for convergence and update states
  if(abs(dStorage) < incTol .and. abs(residual) < resTol)then
   scalarAquiferStorage  = aquiferStorageTrial + dStorage
   scalarWaterTableDepth = soilDepth - scalarAquiferStorage/specificYield   ! water table depth at the end of the time step (m)
   print*, 'baseflowMax, zScale_TOPMODEL = ', baseflowMax, zScale_TOPMODEL
   write(*,'(a)')            'scalarAquiferStorage, scalarWaterTableDepth, aquiferBaseflow, scalarAquiferRcharge, dBaseflow_dStorage, dStorage = '
   write(*,'(10(e20.10,1x))') scalarAquiferStorage, scalarWaterTableDepth, aquiferBaseflow, scalarAquiferRcharge, dBaseflow_dStorage, dStorage
   exit
  ! ***** non-convergence -- update the aquifer storage and water table depth
  else
   ! ** update the aquifer storage (m) 
   aquiferStorageTrial = aquiferStorageTrial + dStorage
   ! ** update the water table depth (m)
   if(dStorage < 0._dp)then
    ! water table drops according to specific yield -- same specific yield in all layers
    zWater = zWater + dStorage/specificYield
   else
    ! water table rises according to the fraction of air in the soil -- different for all layers
    ! (compute fraction of air at the current water table position)
    if(nUnsat == nSoil)then
     fracAir = specificYield
    else
     fracAir = theta_sat - mLayerVolFracLiqIter(nUnsat+1) - mLayerVolFracIceIter(nUnsat+1)
    endif
    ! (get an initial estimate of the water table)
    zTemp = zWater + dStorage/fracAir
    ! (check that we did not cross a layer)
    do iAbsorb=1,maxAbsorb  ! (only absorb a certain number of layers [usually 2])
     if(zTemp < layerHeight(nUnsat)then
      err=20; message=trim(message)//'have not implemented specific yield within the soil column yet'; return; endif
     else
      zWater = zTemp
      exit
     endif
    end do
   endif  ! ** updating the water table depth
  endif  ! ***** non-convergence -- update
  ! check for non-convergence
  if(iter==maxiter)then; err=20; message=trim(message)//'failed to converge'; return; endif
 end do  ! (loop through iterations)

 ! clean-up fluxes
 aquiferBaseflow = aquiferBaseflow + dBaseflow_dStorage*dStorage

 ! assign states
 if(scalarWaterTableDepth < soilDepth)then
  err=20; message=trim(message)//'have not implemented specific yield within the soil column yet'; return
 else
  mLayerVolFracLiqNew(:) = mLayerVolFracLiqIter(:)  ! new volumetric fraction of liquid water (-)
  mLayerVolFracIceNew(:) = mLayerVolFracIceIter(:)  ! new volumetric fraction of ice (-)
  mLayerMatricHeadNew(:) = mLayerMatricHeadIter(:)  ! new matric head after (m)
 endif
 

 end subroutine groundwatr

end module groundwatr_module

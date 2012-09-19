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

 ! get the total number of model layers
 nLayers   => indx_data%var(iLookINDEX%nLayers)%dat(1)           ! total number of layers

 ! assign pointers to soil depth
 soilDepth => mvar_data%var(iLookMVAR%iLayerHeight)%dat(nLayers) ! soil depth

 ! compute the maximum baseflow rate
 baseflowMax = kAnisotropic*k_soil

 ! save the aquifer storage at the start of the step (m)
 aquiferStorageTrial = scalarAquiferStorage
 
 ! iterate
 do iter=1,maxiter
  ! calculate the depth to the water table (m)
  zWater = soilDepth - aquiferStorageTrial/specificYield
  if(zWater < soilDepth)then; err=20; message=trim(message)//'have not implemented specific yield within the soil column yet'; return; endif
  ! calculate the baseflow term (m s-1)
  aquiferBaseflow = baseflowMax*exp(-zWater/zScale_TOPMODEL)
  ! calculate the derivative in the baseflow term w.r.t. aquifer storage (s-1)
  dBaseflow_dStorage = aquiferBaseflow/(zScale_TOPMODEL*specificYield)
  ! compute the residual (m)
  residual = aquiferStorageTrial - (scalarAquiferStorage + (scalarAquiferRcharge - aquiferBaseflow)*dt)
  ! compute the iteration increment (m)
  dStorage = -residual/(1._dp + dt*dBaseflow_dStorage)
  ! update the aquifer storage
  aquiferStorageTrial = aquiferStorageTrial + dStorage
  ! check for convergence
  if(abs(dStorage) < incTol .and. abs(residual) < resTol)exit
  ! check for non-convergence
  if(iter==maxiter)then; err=20; message=trim(message)//'failed to converge'; return; endif
 end do  ! (loop through iterations)

 ! clean-up fluxes
 aquiferBaseflow = aquiferBaseflow + dBaseflow_dStorage*dStorage
 

 end subroutine groundwatr

end module groundwatr_module

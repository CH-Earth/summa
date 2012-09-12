module cmput_dvar_module
USE nrtype
implicit none
private
public::cmput_dvar
contains

 ! ************************************************************************************************
 ! new subroutine: compute diagnostic/ancillary variables
 ! ************************************************************************************************
 subroutine cmput_dvar(err,message)
 USE multiconst,only:iden_ice,iden_water                              ! intrinsic density of ice and water
 USE multiconst,only:lambda_air,lambda_ice,lambda_soil,lambda_water   ! thermal conductivity of individual constituents
 USE soil_utils_module,only:dTheta_dPsi                               ! compute derivative of the soil moisture characteristic
 USE snow_utils_module,only:tcond_snow                                ! compute thermal conductivity of snow
 USE data_struc,only:mpar_data,mvar_data,indx_data,ix_soil,ix_snow    ! data structures
 USE var_lookup,only:iLookPARAM,iLookMVAR,iLookINDEX                  ! named variables for structure elements
 implicit none
 ! dummy variables
 integer(i4b),intent(out)      :: err                    ! error code
 character(*),intent(out)      :: message                ! error message
 ! local pointers to model parameters
 real(dp),pointer              :: vGn_alpha              ! van Genutchen "alpha" parameter
 real(dp),pointer              :: vGn_n                  ! van Genutchen "n" parameter
 real(dp),pointer              :: theta_sat              ! soil porosity (-)
 real(dp),pointer              :: theta_res              ! soil residual volumetric water content (-)
 real(dp),pointer              :: soilAlbedo             ! soil albedo (-)
 ! local pointers to model state ariables
 real(dp),pointer              :: scalarAlbedo           ! surface albedo (-)
 real(dp),pointer              :: mLayerTemp(:)          ! temperature of each layer (K)
 real(dp),pointer              :: mLayerBulkDenIce(:)    ! bulk density of ice in each layer (kg m-3)
 real(dp),pointer              :: mLayerBulkDenLiq(:)    ! bulk density of liquid water in each layer (kg m-3)
 real(dp),pointer              :: mLayerMatricHead(:)    ! matric head of water in soil (m)
 real(dp),pointer              :: mLayerDepth(:)         ! depth of each layer
 ! local pointers to variable short-cuts
 real(dp),pointer              :: vGn_m                  ! van Genutchen "m" parameter
 real(dp),pointer              :: volHtCap_air           ! volumetric heat capacity of air (J m-3 K-1)
 real(dp),pointer              :: volHtCap_ice           ! volumetric heat capacity of ice (J m-3 K-1)
 real(dp),pointer              :: volHtCap_soil          ! volumetric heat capacity of dry soil (J m-3 K-1)
 real(dp),pointer              :: volHtCap_water         ! volumetric heat capacity of liquid water (J m-3 K-1)
 real(dp),pointer              :: volLatHt_fus           ! volumetric latent heat of fusion (J m-3 K-1)
 ! local pointers to model diagnostic variables
 real(dp),pointer              :: mLayerHeight(:)        ! height of the layer mid-point (top of soil = 0)
 real(dp),pointer              :: iLayerHeight(:)        ! height of the layer interface (top of soil = 0)
 real(dp),pointer              :: mLayerVolFracAir(:)    ! volumetric fraction of air in each layer (-)
 real(dp),pointer              :: mLayerVolFracIce(:)    ! volumetric fraction of ice in each layer (-)
 real(dp),pointer              :: mLayerVolFracLiq(:)    ! volumetric fraction of liquid water in each layer (-)
 real(dp),pointer              :: mLayerThermalC(:)      ! thermal conductivity at the mid-point of each layer (W m-1 K-1)
 real(dp),pointer              :: iLayerThermalC(:)      ! thermal conductivity at the interface of each layer (W m-1 K-1)
 real(dp),pointer              :: mLayerVolHtCapBulk(:)  ! volumetric heat capacity in each layer (J m-3 K-1)
 real(dp),pointer              :: mLayerdTheta_dPsi(:)   ! derivative in the soil water characteristic (m-1)
 ! local pointers to model index variables
 integer(i4b),pointer          :: nLayers                ! number of layers
 integer(i4b),pointer          :: layerType(:)           ! type of the layer
 ! local variables
 character(LEN=256)            :: cmessage               ! error message of downwind routine 
 integer(i4b)                  :: iLayer                 ! index of model layer
 real(dp),pointer              :: TCn                    ! thermal conductivity below the layer interface (W m-1 K-1)
 real(dp),pointer              :: TCp                    ! thermal conductivity above the layer interface (W m-1 K-1)
 real(dp)                      :: zdn                    ! height difference between interface and lower value (m)
 real(dp)                      :: zdp                    ! height difference between interface and upper value (m)
 ! initialize error control
 err=0; message="cmput_dvar/"
 ! assign pointers to model parameters
 vGn_alpha           => mpar_data%var(iLookPARAM%vGn_alpha)                  ! van Genutchen "alpha" parameter (m-1)
 vGn_n               => mpar_data%var(iLookPARAM%vGn_n)                      ! van Genutchen "n" parameter (-)
 theta_sat           => mpar_data%var(iLookPARAM%theta_sat)                  ! soil porosity (-)
 theta_res           => mpar_data%var(iLookPARAM%theta_res)                  ! soil residual volumetric water content (-)
 soilAlbedo          => mpar_data%var(iLookPARAM%soilAlbedo)                 ! soil albedo (-)
 ! assign pointers to model state variables
 scalarAlbedo        => mvar_data%var(iLookMVAR%scalarAlbedo)%dat(1)         ! surface albedo (-)
 mLayerTemp          => mvar_data%var(iLookMVAR%mLayerTemp)%dat              ! temperature of each layer (K)
 mLayerVolFracIce    => mvar_data%var(iLookMVAR%mLayerVolFracIce)%dat        ! volumetric fraction of ice in each layer (-)
 mLayerVolFracLiq    => mvar_data%var(iLookMVAR%mLayerVolFracLiq)%dat        ! volumetric fraction of liquid water in each layer (-)
 mLayerMatricHead    => mvar_data%var(iLookMVAR%mLayerMatricHead)%dat        ! matric head of water in soil (m)
 mLayerDepth         => mvar_data%var(iLookMVAR%mLayerDepth)%dat             ! depth of each layer (m)
 ! assign pointers to variable short-cuts
 vGn_m               => mvar_data%var(iLookMVAR%scalarVGn_m)%dat(1)          ! van Genutchen "m" parameter (-)
 volHtCap_air        => mvar_data%var(iLookMVAR%scalarVolHtCap_air)%dat(1)   ! volumetric heat capacity of air (J m-3 K-1)
 volHtCap_ice        => mvar_data%var(iLookMVAR%scalarVolHtCap_ice)%dat(1)   ! volumetric heat capacity of ice (J m-3 K-1)
 volHtCap_soil       => mvar_data%var(iLookMVAR%scalarVolHtCap_soil)%dat(1)  ! volumetric heat capacity of soil (J m-3 K-1)
 volHtCap_water      => mvar_data%var(iLookMVAR%scalarVolHtCap_water)%dat(1) ! volumetric heat capacity of water (J m-3 K-1)
 volLatHt_fus        => mvar_data%var(iLookMVAR%scalarvolLatHt_fus)%dat(1)   ! volumetric latent heat of fusion (J m-3)
 ! assign pointers to model diagnostic variables
 mLayerHeight        => mvar_data%var(iLookMVAR%mLayerHeight)%dat            ! height of the layer midpoint; top of soil = 0  (m)
 iLayerHeight        => mvar_data%var(iLookMVAR%iLayerHeight)%dat            ! height of the layer interface; top of soil = 0 (m)
 mLayerBulkDenIce    => mvar_data%var(iLookMVAR%mLayerBulkDenIce)%dat        ! bulk density of ice in each layer (kg m-3)
 mLayerBulkDenLiq    => mvar_data%var(iLookMVAR%mLayerBulkDenLiq)%dat        ! bulk density of liquid water in each layer (kg m-3)
 mLayerVolFracAir    => mvar_data%var(iLookMVAR%mLayerVolFracAir)%dat        ! volumetric fraction of air in each layer (-)
 mLayerThermalC      => mvar_data%var(iLookMVAR%mLayerThermalC)%dat          ! thermal conductivity at the mid-point of each layer (W m-1 K-1)
 iLayerThermalC      => mvar_data%var(iLookMVAR%iLayerThermalC)%dat          ! thermal conductivity at the interface of each layer (W m-1 K-1)
 mLayerVolHtCapBulk  => mvar_data%var(iLookMVAR%mLayerVolHtCapBulk)%dat      ! volumetric heat capacity in each layer (J m-3 K-1)
 mLayerdTheta_dPsi   => mvar_data%var(iLookMVAR%mLayerdTheta_dPsi)%dat       ! derivative in the soil water characteristic (m-1)
 ! assign pointers to index variables
 nLayers             => indx_data%var(iLookINDEX%nLayers)%dat(1)             ! number of layers
 layerType           => indx_data%var(iLookINDEX%layerType)%dat              ! layer type (ix_soil or ix_snow)
 ! loop through layers
 do iLayer=1,nLayers
  ! compute the derivative in the soil water characteristic (m-1)
  mLayerdTheta_dPsi(iLayer) = dTheta_dPsi(mLayerMatricHead(iLayer),vGn_alpha,theta_res,theta_sat,vGn_n,vGn_m) 
  ! compute the bulk density of ice and liquid water
  mLayerBulkDenIce(iLayer) = mLayerVolFracIce(iLayer) * iden_ice
  mLayerBulkDenLiq(iLayer) = mLayerVolFracLiq(iLayer) * iden_water
  ! compute the volumetric fraction of air
  mLayerVolFracAir(iLayer) = theta_sat - (mLayerVolFracIce(iLayer) + mLayerVolFracLiq(iLayer))
  ! compute the thermal conductivity of soil at the mid-point of each layer
  if(layerType(iLayer)==ix_soil)then
   mLayerThermalC(iLayer) = lambda_soil * (1._dp - theta_sat)      + & ! soil component 
                            lambda_ice  * mLayerVolFracIce(iLayer) + & ! ice component
                            lambda_water* mLayerVolFracLiq(iLayer) + & ! liquid water component
                            lambda_air  * mLayerVolFracAir(iLayer)     ! air component
  endif
  ! compute the thermal conductivity of snow at the mid-point of each layer
  if(layerType(iLayer)==ix_snow)then
   call tcond_snow(mLayerBulkDenIce(iLayer),mLayerThermalC(iLayer),err,cmessage)
   if(err/=0)then; message=trim(message)//trim(cmessage); return; endif
  endif
  ! compute volumetric heat capacity (J m-3 K-1)
  mLayerVolHtCapBulk(iLayer) = volHtCap_soil  * (1._dp - theta_sat)      + & ! soil component
                               volHtCap_ice   * mLayerVolFracIce(iLayer) + & ! ice component
                               volHtCap_water * mLayerVolFracLiq(iLayer) + & ! liquid water component
                               volHtCap_air   * mLayerVolFracAir(iLayer)     ! air component
 end do  ! looping through layers
 ! compute the thermal conductivity of snow at the interface of each layer
 iLayerThermalC(0)       = 0._dp ! prescribe the thermal conductivity at the domain boundaries (not used)
 iLayerThermalC(nLayers) = 0._dp ! prescribe the thermal conductivity at the domain boundaries (not used)
 ! loop through layers
 do iLayer=1,nLayers-1
  TCn => mLayerThermalC(iLayer)    ! thermal conductivity below the layer interface (W m-1 K-1)
  TCp => mLayerThermalC(iLayer+1)  ! thermal conductivity above the layer interface (W m-1 K-1)
  zdn =  iLayerHeight(iLayer)   - mLayerHeight(iLayer) ! height difference between interface and lower value (m)
  zdp =  mLayerHeight(iLayer+1) - iLayerHeight(iLayer) ! height difference between interface and upper value (m)
  iLayerThermalC(iLayer) = (TCn*TCp*(zdn + zdp)) / (TCn*zdp + TCp*zdn)
 end do
 ! identify the appropriate albedo
 if(layerType(nLayers)==ix_soil) scalarAlbedo = soilAlbedo
 end subroutine cmput_dvar

end module cmput_dvar_module

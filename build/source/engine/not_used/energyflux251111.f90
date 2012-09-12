module energyflux_module
USE nrtype
implicit none
private
public::diagn_evar
public::tempchange
public::energy_err
contains

 ! ************************************************************************************************
 ! new subroutine: compute change in temperature over the time step
 ! ************************************************************************************************
 subroutine tempchange(dt,mLayerTempIter,mLayerVolFracIceIter,mLayerTempNew,err,message)
 USE multiconst,only:gravity,  & ! acceleration of gravity     (m s-2)
                     Tfreeze,  & ! freezing point              (K)
                     LH_fus,   & ! latent heat of fusion       (J kg-1)
                     Em_Sno,   & ! emissivity of snow          (-)
                     sigma,    & ! Stefan Boltzman constant    (W m-2 K-4)
                     iden_ice, & ! intrinsic density of ice    (kg m-3)
                     iden_water  ! intrinsic density of water  (kg m-3)
 USE tridagSolv_module,only:tridag                                              ! solve tridiagonal system of equations
 USE data_struc,only:mpar_data,forc_data,mvar_data,indx_data,ix_soil,ix_snow    ! data structures
 USE var_lookup,only:iLookPARAM,iLookFORCE,iLookMVAR,iLookINDEX                 ! named variables for structure elements
 ! compute change in temperature over the time step
 implicit none
 ! dummy variables
 real(dp),intent(in)           :: dt                       ! time step (seconds)
 real(dp),intent(in)           :: mLayerTempIter(:)        ! trial temperature at the current iteration (K)
 real(dp),intent(in)           :: mLayerVolFracIceIter(:)  ! volumetric fraction of ice at the current iteration (-)
 real(dp),intent(out)          :: mLayerTempNew(:)         ! new temperature (K)
 integer(i4b),intent(out)      :: err                      ! error code
 character(*),intent(out)      :: message                  ! error message
 ! local pointers tp model parameters
 real(dp),pointer              :: Fabs_vis                 ! fraction of radiation in visible part of spectrum (-)
 real(dp),pointer              :: rad_ext                  ! extinction coefficient for penetrating sw radiation (m-1)
 real(dp),pointer              :: lowerBoundTemp           ! temperature of the lower boundary (K)
 ! local pointers to model forcing data
 real(dp),pointer              :: sw_down                  ! downward shortwave radiation (W m-2)
 real(dp),pointer              :: lw_down                  ! downward longwave radiation (W m-2)
 real(dp),pointer              :: airtemp                  ! air temperature at 2 meter height (K)
 real(dp),pointer              :: windspd                  ! wind speed at 10 meter height (m s-1)
 real(dp),pointer              :: airpres                  ! air pressure at 2 meter height (Pa)
 real(dp),pointer              :: spechum                  ! specific humidity at 2 meter height (g g-1)
 ! local pointers to model state variables
 real(dp),pointer              :: scalarAlbedo             ! surface albedo (-)
 real(dp),pointer              :: mLayerTemp(:)            ! temperature of each layer (K)
 real(dp),pointer              :: mLayerVolFracIce(:)      ! volumetric fraction of ice in each layer (-)
 real(dp),pointer              :: mLayerDepth(:)           ! depth of each layer (m)
 ! local pointers to variable short-cuts
 real(dp),pointer              :: volLatHt_fus             ! volumetric latent heat of fusion (J m-3 K-1)
 ! local pointers to model diagnostic variables
 real(dp),pointer              :: mLayerHeight(:)          ! height at the mid-point of each layer (m)
 real(dp),pointer              :: iLayerHeight(:)          ! height at the interface of each layer (m)
 real(dp),pointer              :: mLayerVolHtCapBulk(:)    ! bulk volumetric heat capacity (J m-3 K-1)
 real(dp),pointer              :: mLayerThermalC(:)        ! thermal conductivity at the mid-point of each layer (W m-1 K-1)
 real(dp),pointer              :: iLayerThermalC(:)        ! thermal conductivity at the interface of each layer (W m-1 K-1)
 real(dp),pointer              :: mLayerdTheta_dTk(:)      ! derivative in the freezing curve (K-1)
 real(dp),pointer              :: mLayerInfilFreeze(:)     ! increase in volumetric ice content caused by freezing infiltrating flux (-)
 ! local pointers to model index variables
 integer(i4b),pointer          :: nLayers                  ! number of layers
 integer(i4b),pointer          :: layerType(:)             ! type of the layer (ix_soil or ix_snow)
 ! define local variables
 character(LEN=256)            :: cmessage                 ! error message of downwind routine
 real(dp)                      :: radFlux                  ! radiative flux at the surface (J m-2 s-1)
 real(dp)                      :: lwderiv                  ! derivative in longwave radiation w.r.t. temperature (J m-2 s-1 K-1)
 real(dp),allocatable          :: dz_node(:)               ! distance between the mid-point of a layer and a neighbouring layer (m)
 real(dp),allocatable          :: dt_dudz(:)               ! the dt_dudz terms for the upper interface (s m-2)
 real(dp),allocatable          :: dt_dldz(:)               ! the dt_dudz terms for the lower interface (s m-2)
 real(dp),allocatable          :: hfusion(:)               ! the fusion term (J m-3 K-1)
 real(dp),allocatable          :: rv_temp(:)               ! common elements of the RHS vector (J m-3)
 real(dp),allocatable          :: d_m1(:)                  ! sub-diagonal elements of the tridiagonal system (J m-3 K-1)
 real(dp),allocatable          :: diag(:)                  ! diagonal elements of the tridiagonal system (J m-3 K-1)
 real(dp),allocatable          :: d_p1(:)                  ! super-diagonal elements of the tridiagonal system (J m-3 K-1)
 real(dp),allocatable          :: rvec(:)                  ! right-hand-side vector (J m-3)
 integer(i4b)                  :: iLayer                   ! index of model layer
 ! initialize error control
 err=0; message="tempchange/"
 ! assign pointers to model parameters
 Fabs_vis            => mpar_data%var(iLookPARAM%Fabs_vis)                   ! fraction of radiation in visible part of spectrum (-)
 rad_ext             => mpar_data%var(iLookPARAM%rad_ext)                    ! extinction coefficient for penetrating sw radiation (m-1)
 lowerBoundTemp      => mpar_data%var(iLookPARAM%lowerBoundTemp)             ! temperature of the lower boundary (K)
 ! assign pointers to model forcing variables
 sw_down             => forc_data%var(iLookFORCE%sw_down)                    ! downward shortwave radiation (W m-2)
 lw_down             => forc_data%var(iLookFORCE%lw_down)                    ! downward longwave radiation (W m-2)
 airtemp             => forc_data%var(iLookFORCE%airtemp)                    ! air temperature at 2 meter height (K)
 windspd             => forc_data%var(iLookFORCE%windspd)                    ! wind speed at 10 meter height (m s-1)
 airpres             => forc_data%var(iLookFORCE%airpres)                    ! air pressure at 2 meter height (Pa)
 spechum             => forc_data%var(iLookFORCE%spechum)                    ! specific humidity at 2 meter height (g g-1)
 ! assign pointers to model state variables
 scalarAlbedo        => mvar_data%var(iLookMVAR%scalarAlbedo)%dat(1)         ! surface albedo (-)
 mLayerTemp          => mvar_data%var(iLookMVAR%mLayerTemp)%dat              ! temperature of each layer (K)
 mLayerVolFracIce    => mvar_data%var(iLookMVAR%mLayerVolFracIce)%dat        ! volumetric fraction of ice in each layer (-)
 mLayerDepth         => mvar_data%var(iLookMVAR%mLayerDepth)%dat             ! depth of each layer (m)
 ! assin pointers to variable short-cuts
 volLatHt_fus        => mvar_data%var(iLookMVAR%scalarvolLatHt_fus)%dat(1)   ! volumetric latent heat of fusion (J m-3)
 ! assign pointers to model diagnostic variables
 mLayerHeight        => mvar_data%var(iLookMVAR%mLayerHeight)%dat            ! height at the mid-point of each layer (m)
 iLayerHeight        => mvar_data%var(iLookMVAR%iLayerHeight)%dat            ! height at the interface of each layer (m)
 mLayerThermalC      => mvar_data%var(iLookMVAR%mLayerThermalC)%dat          ! thermal conductivity at the mid-point of each layer (W m-1 K-1)
 iLayerThermalC      => mvar_data%var(iLookMVAR%iLayerThermalC)%dat          ! thermal conductivity at the interface of each layer (W m-1 K-1)
 mLayerVolHtCapBulk  => mvar_data%var(iLookMVAR%mLayerVolHtCapBulk)%dat      ! bulk volumetric heat capacity (J m-3 K-1)
 mLayerdTheta_dTk    => mvar_data%var(iLookMVAR%mLayerdTheta_dTk)%dat        ! derivative in the freezing curve (K-1)
 mLayerInfilFreeze   => mvar_data%var(iLookMVAR%mLayerInfilFreeze)%dat       ! increase in volumetric ice content caused by freezing infiltrating flux (-)
 ! assign pointers to index variables
 nLayers             => indx_data%var(iLookINDEX%nLayers)%dat(1)             ! number of layers
 layerType           => indx_data%var(iLookINDEX%layerType)%dat              ! layer type (ix_soil or ix_snow)

 ! allocate space for the dt/dzdz terms
 allocate(dz_node(nLayers-1),dt_dudz(nLayers),dt_dldz(nLayers),stat=err)
 if(err/=0)then; err=20; message=trim(message)//'problem allocating space for dt/dzdz vectors'; return; endif
 ! allocate space for temporary vectors
 allocate(hfusion(nLayers),rv_temp(nLayers),stat=err)
 if(err/=0)then; err=20; message=trim(message)//'problem allocating space for temporary vectors'; return; endif
 ! allocate space for the tri-diagonal vectors
 allocate(d_m1(nLayers-1),diag(nLayers),d_p1(nLayers-1),rvec(nLayers),stat=err)
 if(err/=0)then; err=20; message=trim(message)//'problem allocating space for tridag vectors'; return; endif

 ! compute the dt/dzdz terms for each layer(s m-2)
 dz_node(1:nLayers-1) = mLayerHeight(2:nLayers) - mLayerHeight(1:nLayers-1)
 dt_dudz(2:nLayers)   = dt/(mLayerDepth(2:nLayers)*dz_node)    ! upper distance, valid 2..nlayers
 dt_dldz(1:nLayers-1) = dt/(mLayerDepth(1:nLayers-1)*dz_node)  ! lower distance, valid 1..nlayers-1

 ! compute the dt/dzdz terms at the boundaries
 dt_dudz(1)       = dt/(mLayerDepth(1)       * mLayerDepth(1)*0.5_dp)
 dt_dldz(nLayers) = dt/(mLayerDepth(nLayers) * mLayerDepth(nLayers)*0.5_dp)

 ! compute the fusion term for each layer (J m-3 K-1)
 hfusion = LH_fus*iden_water*mLayerdTheta_dTk

 ! compute the radiative fluxes on the RHS -- note use temperature for current iteration (J m-2 s-1)
 radFlux = sw_down*(1._dp - scalarAlbedo) + lw_down - Em_Sno*sigma*mLayerTempIter(1)**4._dp

 ! compute the derivative in longwave radiation (J m-2 s-1 K-1)
 lwderiv = 4._dp*Em_Sno*sigma*mLayerTempIter(1)**3._dp

 ! compute the term on the RHS that is common to all layers (J m-3)
 rv_temp = mLayerVolHtCapBulk*mLayerTemp + hfusion*mLayerTempIter + &
            LH_fus*iden_ice*((mLayerVolFracIceIter+mLayerInfilFreeze) - mLayerVolFracIce)

 ! assemble the tri-diagonal system
 do iLayer=1,nLayers
  ! compute the off-diagonal elements
  if(iLayer<nLayers) d_p1(iLayer)   = -dt_dldz(iLayer)*iLayerThermalC(iLayer)
  if(iLayer>1)       d_m1(iLayer-1) = -dt_dudz(iLayer)*iLayerThermalC(iLayer-1)
  ! compute the diagonal elements (J m-3 K-1) and the right-hand-side vector (J m-3)
  if(iLayer==nLayers)then
   ! ** lower-most layer
   diag(iLayer) = mLayerVolHtCapBulk(iLayer) + hfusion(iLayer) + &
                   dt_dudz(iLayer)*iLayerThermalC(iLayer-1) + dt_dldz(iLayer)*iLayerThermalC(iLayer)
   rvec(iLayer) = rv_temp(iLayer) + dt_dldz(iLayer)*iLayerThermalC(iLayer)*lowerBoundTemp
  elseif(iLayer==1)then
   ! ** upper-most layer
   diag(iLayer) = mLayerVolHtCapBulk(iLayer) + hfusion(iLayer) + (dt/mLayerDepth(iLayer))*lwderiv + &
                   dt_dldz(iLayer)*iLayerThermalC(iLayer)
   rvec(iLayer) = rv_temp(iLayer) + (dt/mLayerDepth(iLayer))*(radFlux + lwderiv*mLayerTempIter(iLayer))
  else
   ! ** intermediate layers
   diag(iLayer) = mLayerVolHtCapBulk(iLayer) + hfusion(iLayer) + &
                   dt_dudz(iLayer)*iLayerThermalC(iLayer-1) + dt_dldz(iLayer)*iLayerThermalC(iLayer)
   rvec(iLayer) = rv_temp(iLayer)
  endif
 end do  ! (assembling the tri-diagonal system)

 ! solve the tridiagonal system of equations -- returns mLayerTempNew
 call tridag(d_m1,diag,d_p1,rvec,mLayerTempNew,err,cmessage)
 if(err/=0)then; message=trim(message)//trim(cmessage); return; endif
 print*, 'in energyflux, mLayerDelTemp(1:5) = ', mLayerTempNew(1:5) - mLayerTempIter(1:5)
 print*, 'in energyflux, mLayerTempNew(1:5) = ', mLayerTempNew(1:5)

 !print*, 'test tridag'
 ! check the tri-diagonal system
 !do iLayer=1,nLayers
 ! if(iLayer==1) then
 !  print*, iLayer, diag(iLayer)*mLayerTempNew(iLayer) + d_p1(iLayer)*mLayerTempNew(iLayer+1), rvec(iLayer)
 ! endif
 ! if(iLayer==nLayers) &
 !  print*, iLayer, d_m1(iLayer-1)*mLayerTempNew(iLayer-1) + diag(iLayer)*mLayerTempNew(iLayer), rvec(iLayer)
 ! if(iLayer>1 .and. iLayer<nLayers) &
 !  print*, iLayer, d_m1(iLayer-1)*mLayerTempNew(iLayer-1) + diag(iLayer)*mLayerTempNew(iLayer) + d_p1(iLayer)*mLayerTempNew(iLayer+1), rvec(iLayer)
 !end do

 ! deallocate vectors
 deallocate(dz_node,dt_dudz,dt_dldz,hfusion,rv_temp,d_m1,diag,d_p1,rvec,stat=err)
 if(err/=0)then; err=30; message=trim(message)//'problem deallocating vectors for tridag solution'; return; endif
 ! ====================================================================================================================

 end subroutine tempchange


 ! ************************************************************************************************
 ! new subroutine: compute error in the temperature equation
 ! ************************************************************************************************
 subroutine energy_err(dt,mLayerTempNew,mLayerVolFracIceNew,mLayerNrgError,err,message)
 USE multiconst,only:LH_fus,   & ! latent heat of fusion       (J kg-1)
                     Em_Sno,   & ! emissivity of snow          (-)
                     sigma,    & ! Stefan Boltzman constant    (W m-2 K-4)
                     iden_ice, & ! intrinsic density of ice    (kg m-3)
                     iden_water  ! intrinsic density of water  (kg m-3)
 USE data_struc,only:mpar_data,forc_data,mvar_data,indx_data,ix_soil,ix_snow    ! data structures
 USE var_lookup,only:iLookPARAM,iLookFORCE,iLookMVAR,iLookINDEX                 ! named variables for structure elements
 implicit none
 ! dummy variables
 real(dp),intent(in)           :: dt                       ! time step (seconds)
 real(dp),intent(in)           :: mLayerTempNew(:)         ! trial temperature at the current iteration (K)
 real(dp),intent(in)           :: mLayerVolFracIceNew(:)   ! volumetric fraction of ice at the current iteration (-)
 real(dp),intent(out)          :: mLayerNrgError(:)        ! energy error (J m-3)
 integer(i4b),intent(out)      :: err                      ! error code
 character(*),intent(out)      :: message                  ! error message
 ! internal variables
 integer(i4b)                              :: iLayer       ! loop through layers
 real(dp),dimension(0:size(mLayerTempNew)) :: mLayerFlux   ! fluxes at later interfaces
 ! local pointers to model forcing data
 real(dp),pointer              :: sw_down                  ! downward shortwave radiation (W m-2)
 real(dp),pointer              :: lw_down                  ! downward longwave radiation (W m-2)
 real(dp),pointer              :: airtemp                  ! air temperature at 2 meter height (K)
 real(dp),pointer              :: windspd                  ! wind speed at 10 meter height (m s-1)
 real(dp),pointer              :: airpres                  ! air pressure at 2 meter height (Pa)
 real(dp),pointer              :: spechum                  ! specific humidity at 2 meter height (g g-1)
 ! local pointers to model state variables
 real(dp),pointer              :: scalarAlbedo             ! surface albedo (-)
 real(dp),pointer              :: mLayerTemp(:)            ! temperature of each layer (K)
 real(dp),pointer              :: mLayerVolFracIce(:)      ! volumetric fraction of ice in each layer (-)
 real(dp),pointer              :: mLayerDepth(:)           ! depth of each layer (m)
 ! local pointers to model diagnostic variables
 real(dp),pointer              :: mLayerHeight(:)          ! height at the mid-point of each layer (m)
 real(dp),pointer              :: mLayerVolHtCapBulk(:)    ! bulk volumetric heat capacity (J m-3 K-1)
 real(dp),pointer              :: iLayerThermalC(:)        ! thermal conductivity at the interface of each layer (W m-1 K-1)
 real(dp),pointer              :: mLayerInfilFreeze(:)     ! increase in volumetric ice content caused by freezing infiltrating flux (-)
 ! local pointers to model index variables
 integer(i4b),pointer          :: nLayers                  ! number of layers
 ! initialize error control
 err=0; message="energy_err/"
 ! assign pointers to model forcing variables
 sw_down             => forc_data%var(iLookFORCE%sw_down)                    ! downward shortwave radiation (W m-2)
 lw_down             => forc_data%var(iLookFORCE%lw_down)                    ! downward longwave radiation (W m-2)
 airtemp             => forc_data%var(iLookFORCE%airtemp)                    ! air temperature at 2 meter height (K)
 windspd             => forc_data%var(iLookFORCE%windspd)                    ! wind speed at 10 meter height (m s-1)
 airpres             => forc_data%var(iLookFORCE%airpres)                    ! air pressure at 2 meter height (Pa)
 spechum             => forc_data%var(iLookFORCE%spechum)                    ! specific humidity at 2 meter height (g g-1)
 ! assign pointers to model state variables
 scalarAlbedo        => mvar_data%var(iLookMVAR%scalarAlbedo)%dat(1)         ! surface albedo (-)
 mLayerTemp          => mvar_data%var(iLookMVAR%mLayerTemp)%dat              ! temperature of each layer (K)
 mLayerVolFracIce    => mvar_data%var(iLookMVAR%mLayerVolFracIce)%dat        ! volumetric fraction of ice in each layer (-)
 mLayerDepth         => mvar_data%var(iLookMVAR%mLayerDepth)%dat             ! depth of each layer (m)
 ! assign pointers to model diagnostic variables
 mLayerHeight        => mvar_data%var(iLookMVAR%mLayerHeight)%dat            ! height at the mid-point of each layer (m)
 iLayerThermalC      => mvar_data%var(iLookMVAR%iLayerThermalC)%dat          ! thermal conductivity at the interface of each layer (W m-1 K-1)
 mLayerVolHtCapBulk  => mvar_data%var(iLookMVAR%mLayerVolHtCapBulk)%dat      ! bulk volumetric heat capacity (J m-3 K-1)
 mLayerInfilFreeze   => mvar_data%var(iLookMVAR%mLayerInfilFreeze)%dat       ! increase in volumetric ice content caused by freezing infiltrating flux (-)
 ! assign pointers to index variables
 nLayers             => indx_data%var(iLookINDEX%nLayers)%dat(1)             ! number of layers

 ! ***** compute fluxes at layer interfaces
 ! compute flux at the upper boundary -- positive downwards
 mLayerFlux(0) = sw_down*(1._dp - scalarAlbedo) + lw_down - Em_Sno*sigma*mLayerTempNew(1)**4._dp
 ! compute fluxes within the domain -- positive downwards
 do iLayer=1,nLayers-1
  mLayerFlux(iLayer) = -iLayerThermalC(iLayer)*(mLayerTempNew(iLayer+1) - mLayerTempNew(iLayer)) / &
                                               (mLayerHeight(iLayer+1) - mLayerHeight(iLayer))
 end do
 ! compute fluxes at the lower boundary (assume zero for now)
 mLayerFlux(nLayers) = 0._dp
 !print*, 'mLayerFlux = ', mLayerFlux
 !print*, 'nrg term = ', ((mLayerVolHtCapBulk*mLayerDepth)/dt)*(mLayerTempNew - mLayerTemp)
 !print*, 'phaseterm = ', (mLayerDepth/dt)*LH_fus*iden_ice*(mLayerVolFracIceNew - mLayerVolFracIce)
 !print *, mLayerVolFracIceNew - mLayerVolFracIce

 ! ***** compute the energy error
 mLayerNrgError = ((mLayerVolHtCapBulk*mLayerDepth)/dt)*(mLayerTempNew - mLayerTemp) - &
                   (mLayerDepth/dt)*LH_fus*iden_ice*((mLayerVolFracIceNew+mLayerInfilFreeze) - mLayerVolFracIce) - &
                   (-(mLayerFlux(1:nLayers) - mLayerFlux(0:nLayers-1)))
                    
 end subroutine energy_err


 ! ************************************************************************************************
 ! new subroutine: compute diagnostic energy variables (thermal conductivity and heat capacity) 
 ! ************************************************************************************************
 subroutine diagn_evar(err,message)
 USE multiconst,only:iden_ice,iden_water                              ! intrinsic density of ice and water
 USE multiconst,only:lambda_air,lambda_ice,lambda_soil,lambda_water   ! thermal conductivity of individual constituents
 USE snow_utils_module,only:tcond_snow                                ! compute thermal conductivity of snow
 USE data_struc,only:mpar_data,mvar_data,indx_data,ix_soil,ix_snow    ! data structures
 USE var_lookup,only:iLookPARAM,iLookMVAR,iLookINDEX                  ! named variables for structure elements
 implicit none
 ! dummy variables
 integer(i4b),intent(out)      :: err                    ! error code
 character(*),intent(out)      :: message                ! error message
 ! local pointers to model parameters
 real(dp),pointer              :: theta_sat              ! soil porosity (-)
 ! local pointers to model state variables
 real(dp),pointer              :: mLayerVolFracIce(:)    ! volumetric fraction of ice in each layer (-)
 real(dp),pointer              :: mLayerVolFracLiq(:)    ! volumetric fraction of liquid water in each layer (-)
 ! local pointers to variable short-cuts
 real(dp),pointer              :: volHtCap_air           ! volumetric heat capacity of air (J m-3 K-1)
 real(dp),pointer              :: volHtCap_ice           ! volumetric heat capacity of ice (J m-3 K-1)
 real(dp),pointer              :: volHtCap_soil          ! volumetric heat capacity of dry soil (J m-3 K-1)
 real(dp),pointer              :: volHtCap_water         ! volumetric heat capacity of liquid water (J m-3 K-1)
 real(dp),pointer              :: volLatHt_fus           ! volumetric latent heat of fusion (J m-3 K-1)
 ! local pointers to model diagnostic variables
 real(dp),pointer              :: mLayerHeight(:)        ! height of the layer mid-point (top of soil = 0)
 real(dp),pointer              :: iLayerHeight(:)        ! height of the layer interface (top of soil = 0)
 real(dp),pointer              :: mLayerVolFracAir(:)    ! volumetric fraction of air in each layer (-)
 real(dp),pointer              :: mLayerThermalC(:)      ! thermal conductivity at the mid-point of each layer (W m-1 K-1)
 real(dp),pointer              :: iLayerThermalC(:)      ! thermal conductivity at the interface of each layer (W m-1 K-1)
 real(dp),pointer              :: mLayerVolHtCapBulk(:)  ! volumetric heat capacity in each layer (J m-3 K-1)
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
 err=0; message="diagn_evar/"
 ! assign pointers to model parameters
 theta_sat           => mpar_data%var(iLookPARAM%theta_sat)                  ! soil porosity (-)
 ! assign pointers to model state variables
 mLayerVolFracIce    => mvar_data%var(iLookMVAR%mLayerVolFracIce)%dat        ! volumetric fraction of ice in each layer (-)
 mLayerVolFracLiq    => mvar_data%var(iLookMVAR%mLayerVolFracLiq)%dat        ! volumetric fraction of liquid water in each layer (-)
 ! assign pointers to variable short-cuts
 volHtCap_air        => mvar_data%var(iLookMVAR%scalarVolHtCap_air)%dat(1)   ! volumetric heat capacity of air (J m-3 K-1)
 volHtCap_ice        => mvar_data%var(iLookMVAR%scalarVolHtCap_ice)%dat(1)   ! volumetric heat capacity of ice (J m-3 K-1)
 volHtCap_soil       => mvar_data%var(iLookMVAR%scalarVolHtCap_soil)%dat(1)  ! volumetric heat capacity of soil (J m-3 K-1)
 volHtCap_water      => mvar_data%var(iLookMVAR%scalarVolHtCap_water)%dat(1) ! volumetric heat capacity of water (J m-3 K-1)
 volLatHt_fus        => mvar_data%var(iLookMVAR%scalarvolLatHt_fus)%dat(1)   ! volumetric latent heat of fusion (J m-3)
 ! assign pointers to model diagnostic variables
 mLayerHeight        => mvar_data%var(iLookMVAR%mLayerHeight)%dat            ! height of the layer midpoint; top of soil = 0  (m)
 iLayerHeight        => mvar_data%var(iLookMVAR%iLayerHeight)%dat            ! height of the layer interface; top of soil = 0 (m)
 mLayerVolFracAir    => mvar_data%var(iLookMVAR%mLayerVolFracAir)%dat        ! volumetric fraction of air in each layer (-)
 mLayerThermalC      => mvar_data%var(iLookMVAR%mLayerThermalC)%dat          ! thermal conductivity at the mid-point of each layer (W m-1 K-1)
 iLayerThermalC      => mvar_data%var(iLookMVAR%iLayerThermalC)%dat          ! thermal conductivity at the interface of each layer (W m-1 K-1)
 mLayerVolHtCapBulk  => mvar_data%var(iLookMVAR%mLayerVolHtCapBulk)%dat      ! volumetric heat capacity in each layer (J m-3 K-1)
 ! assign pointers to index variables
 nLayers             => indx_data%var(iLookINDEX%nLayers)%dat(1)             ! number of layers
 layerType           => indx_data%var(iLookINDEX%layerType)%dat              ! layer type (ix_soil or ix_snow)
 ! loop through layers
 do iLayer=1,nLayers
  ! compute the thermal conductivity of soil at the mid-point of each layer
  if(layerType(iLayer)==ix_soil)then
   mLayerThermalC(iLayer) = lambda_soil * (1._dp - theta_sat)      + & ! soil component
                            lambda_ice  * mLayerVolFracIce(iLayer) + & ! ice component
                            lambda_water* mLayerVolFracLiq(iLayer) + & ! liquid water component
                            lambda_air  * mLayerVolFracAir(iLayer)     ! air component
  endif
  ! compute the thermal conductivity of snow at the mid-point of each layer
  if(layerType(iLayer)==ix_snow)then
   call tcond_snow(mLayerVolFracIce(iLayer)*iden_ice,mLayerThermalC(iLayer),err,cmessage)
   if(err/=0)then; message=trim(message)//trim(cmessage); return; endif
  endif
  ! compute volumetric heat capacity (J m-3 K-1)
  mLayerVolHtCapBulk(iLayer) = volHtCap_soil  * (1._dp - theta_sat)      + & ! soil component
                               volHtCap_ice   * mLayerVolFracIce(iLayer) + & ! ice component
                               volHtCap_water * mLayerVolFracLiq(iLayer) + & ! liquid water component
                               volHtCap_air   * mLayerVolFracAir(iLayer)     ! air component
 end do  ! looping through layers
 ! compute the thermal conductivity of snow at the interface of each layer
 do iLayer=1,nLayers-1  ! (loop through layers)
  TCn => mLayerThermalC(iLayer)    ! thermal conductivity below the layer interface (W m-1 K-1)
  TCp => mLayerThermalC(iLayer+1)  ! thermal conductivity above the layer interface (W m-1 K-1)
  zdn =  iLayerHeight(iLayer)   - mLayerHeight(iLayer) ! height difference between interface and lower value (m)
  zdp =  mLayerHeight(iLayer+1) - iLayerHeight(iLayer) ! height difference between interface and upper value (m)
  iLayerThermalC(iLayer) = (TCn*TCp*(zdn + zdp)) / (TCn*zdp + TCp*zdn)
 end do
 ! assume the thermal conductivity at the domain boundaries is equal to the thermal conductivity of the layer
 iLayerThermalC(0)       = mLayerThermalC(1)
 iLayerThermalC(nLayers) = mLayerThermalC(nLayers)
 end subroutine diagn_evar


end module energyflux_module

module energyflux_module
USE nrtype
implicit none
private
public::tempchange
contains

 ! ************************************************************************************************
 ! new subroutine: compute change in temperature over the time step
 ! ************************************************************************************************
 subroutine tempchange(dt,mLayerTempIter,mLayerTotalNrg,mLayerDelTemp,err,message)
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
 real(dp),intent(in)           :: mLayerTempIter(:)        ! trial temperature of each layer (K)
 real(dp),intent(in)           :: mLayerTotalNrg(:)        ! total energy up to the current iteration (J m-3)
 real(dp),intent(out)          :: mLayerDelTemp(:)         ! temperature increment in each layer (K)
 integer(i4b),intent(out)      :: err                      ! error code
 character(*),intent(out)      :: message                  ! error message
 ! local pointers tp model parameters
 real(dp),pointer              :: Fabs_vis                 ! fraction of radiation in visible part of spectrum (-)
 real(dp),pointer              :: rad_ext                  ! extinction coefficient for penetrating sw radiation (m-1)
 ! local pointers to model forcing data
 real(dp),pointer              :: sw_down                  ! downward shortwave radiation (W m-2)
 ! local pointers to model state ariables
 real(dp),pointer              :: scalarAlbedo             ! surface albedo (-)
 real(dp),pointer              :: mLayerTemp(:)            ! temperature of each layer (K)
 real(dp),pointer              :: mLayerDepth(:)           ! depth of each layer (m)
 ! local pointers to variable short-cuts
 real(dp),pointer              :: volLatHt_fus             ! volumetric latent heat of fusion (J m-3 K-1)
 ! local pointers to model diagnostic variables
 real(dp),pointer              :: mLayerHeight(:)          ! height at the mid-point of each layer (m)
 real(dp),pointer              :: iLayerHeight(:)          ! height at the interface of each layer (m)
 real(dp),pointer              :: mLayerVolHtCapBulk(:)    ! bulk volumetric heat capacity (J m-3 K-1)
 real(dp),pointer              :: iLayerThermalC(:)        ! thermal conductivity at the interface of each layer (W m-1 K-1)
 real(dp),pointer              :: iLayerRadCondFlux(:)     ! energy flux at the interface of each layer (W m-2)
 real(dp),pointer              :: mLayerRadCondFlux(:)     ! change in energy in each layer (J m-3 s-1)
 real(dp),pointer              :: mLayerdTheta_dTk(:)      ! derivative in the freezing curve (K-1)
 ! local pointers to model index variables
 integer(i4b),pointer          :: nLayers                  ! number of layers
 integer(i4b),pointer          :: layerType(:)             ! type of the layer (ix_soil or ix_snow)
 ! define local variables
 character(LEN=256)            :: cmessage                 ! error message of downwind routine 
 real(dp)                      :: conv                     ! conversion factor -- depth (m) * volumetric heat capacity (J m-3 K) --> J m-2 K-1 
 real(dp),allocatable          :: mLayerVolHtCapApp(:)     ! apparent volumetric heat capacity (J m-3 K-1)
 real(dp),allocatable          :: mLayerRadCondFluxIter(:) ! net change in energy within each layer (J m-3 s-1)
 real(dp),allocatable          :: iLayerTempDeriv(:)       ! derivatives of fluxes at layer interfaces with respect to temperature (J m-2 K-1 s-1)
 real(dp),allocatable          :: a_tri(:)                 ! sub-diagonal elements of the tridiagonal system (-)
 real(dp),allocatable          :: b_tri(:)                 ! diagonal elements of the tridiagonal system (-)
 real(dp),allocatable          :: c_tri(:)                 ! super-diagonal elements of the tridiagonal system (-)
 real(dp),allocatable          :: r_vec(:)                 ! residual vector (K)
 real(dp),allocatable          :: jmat(:,:)                ! Jacobian matrix (-)
 integer(i4b)                  :: iLayer                   ! index of model layer
 ! initialize error control
 err=0; message="f-energyflux/"
 ! assign pointers to model parameters
 Fabs_vis            => mpar_data%var(iLookPARAM%Fabs_vis)                   ! fraction of radiation in visible part of spectrum (-)
 rad_ext             => mpar_data%var(iLookPARAM%rad_ext)                    ! extinction coefficient for penetrating sw radiation (m-1)
 ! assign pointers to model forcing variables
 sw_down             => forc_data%var(iLookFORCE%sw_down)                    ! downward shortwave radiation (W m-2)
 ! assign pointers to model state variables
 scalarAlbedo        => mvar_data%var(iLookMVAR%scalarAlbedo)%dat(1)         ! surface albedo (-)
 mLayerTemp          => mvar_data%var(iLookMVAR%mLayerTemp)%dat              ! temperature of each layer (K)
 mLayerDepth         => mvar_data%var(iLookMVAR%mLayerDepth)%dat             ! depth of each layer (m)
 ! assin pointers to variable short-cuts
 volLatHt_fus        => mvar_data%var(iLookMVAR%scalarvolLatHt_fus)%dat(1)   ! volumetric latent heat of fusion (J m-3)
 ! assign pointers to model diagnostic variables
 mLayerHeight        => mvar_data%var(iLookMVAR%mLayerHeight)%dat            ! height at the mid-point of each layer (m)
 iLayerHeight        => mvar_data%var(iLookMVAR%iLayerHeight)%dat            ! height at the interface of each layer (m)
 iLayerThermalC      => mvar_data%var(iLookMVAR%iLayerThermalC)%dat          ! thermal conductivity at the interface of each layer (W m-1 K-1)
 iLayerRadCondFlux   => mvar_data%var(iLookMVAR%iLayerRadCondFlux)%dat       ! energy flux at the interface of each layer (W m-2)
 mLayerRadCondFlux   => mvar_data%var(iLookMVAR%mLayerRadCondFlux)%dat       ! change in energy in each layer (J m-3 s-1)
 mLayerVolHtCapBulk  => mvar_data%var(iLookMVAR%mLayerVolHtCapBulk)%dat      ! bulk volumetric heat capacity (J m-3 K-1)
 mLayerdTheta_dTk    => mvar_data%var(iLookMVAR%mLayerdTheta_dTk)%dat        ! derivative in the freezing curve (K-1)
 ! assign pointers to index variables
 nLayers             => indx_data%var(iLookINDEX%nLayers)%dat(1)             ! number of layers
 layerType           => indx_data%var(iLookINDEX%layerType)%dat              ! layer type (ix_soil or ix_snow)


 ! allocate space for the apparent heat capacity (J m-3 K-1)
 allocate(mLayerVolHtCapApp(nLayers),stat=err)
 if(err/=0)then; err=20; message=trim(message)//'problem allocating space for the apparent volumetric heat capacity'; return; endif
 ! allocate space for the derivatives of fluxes at layer interfaces with respect to temperature
 allocate(mLayerRadCondFluxIter(nLayers),iLayerTempDeriv(0:nLayers),stat=err)
 if(err/=0)then; err=20; message=trim(message)//'problem allocating space for flux derivative vectors'; return; endif
 allocate(jmat(nLayers,nLayers),stat=err)
 if(err/=0)then; err=20; message=trim(message)//'problem allocating space for the Jacobian matrix'; return; endif
 ! allocate space for the tri-diagonal vectors, and the tridiagonal solution
 allocate(a_tri(nLayers-1),b_tri(nLayers),c_tri(nLayers-1),r_vec(nLayers),stat=err)
 if(err/=0)then; err=20; message=trim(message)//'problem allocating space for tridag vectors'; return; endif
 ! compute internal fluxes
 call intrnlFlux(mLayerTempIter,mLayerRadCondFlux,err,cmessage) 
 if(err/=0)then; err=10; message=trim(message)//trim(cmessage); return; endif
 ! compute the apparent heat capacity (J m-3 K-1)
 mlayerVolHtCapApp = mLayerVolHtCapBulk + iden_water*LH_fus*mLayerdTheta_dTk
 ! compute the residual vector in terms of temperature (K)
 r_vec = mLayerTempIter - (mLayerTemp + (dt*mLayerRadCondFlux + mLayerTotalNrg)/mlayerVolHtCapApp)
 print*, mlayerVolHtCapApp*(mLayerTempIter-mLayerTemp)
 print*, dt*mLayerRadCondFlux
 print*, mLayerTotalNrg
 print*, 'r_vec = ', r_vec
 ! compute finite-difference Jacobian
 call fDiffJacbn(dt,mLayerTempIter,jmat,err,cmessage)
 if(err/=0)then; err=10; message=trim(message)//trim(cmessage); return; endif
 print*, 'Jacobian = '
 do iLayer=1,nLayers
  write(*,'(10(e20.10))'), jmat(:,iLayer)
 end do
 print*, ' after Jacobian'

 ! compute derivative of the flux at the bottom of the domain, and the top of the surface layer
 iLayerTempDeriv(0)       = 0._dp 
 iLayerTempDeriv(nLayers) = 4._dp * Em_Sno * sigma * mLayerTempIter(nLayers)**3._dp   ! longwave radiation derivative (W m-2 K-1)
 ! compute derivative of the fluxes at all layer interfaces, and build the tri-diagonal system
 do iLayer=1,nLayers
  ! compute derivative of the fluxes at all layer interfaces
  if(iLayer<nLayers) & ! (flux at the top of the surface layer computed earlier)
   iLayerTempDeriv(iLayer) = iLayerThermalC(iLayer) / (mLayerHeight(iLayer+1) - mLayerHeight(iLayer)) ! J m-2 K-1 s-1
  ! compute the conversion factor (m-1 s J-1 m3 K = J-1 m2 K s)
  conv = (1._dp/mLayerDepth(iLayer))*(dt/mlayerVolHtCapApp(iLayer))   ! J-1
  ! assemble the elements of the tri-diagonal matrix (dimensionless)
  if(iLayer>1)       a_tri(iLayer-1) = 0._dp   - (                             iLayerTempDeriv(iLayer-1))*conv
                     b_tri(iLayer+0) = 1._dp   - (-iLayerTempDeriv(iLayer)   - iLayerTempDeriv(iLayer-1))*conv 
  if(iLayer<nLayers) c_tri(iLayer+0) = 0._dp   - ( iLayerTempDeriv(iLayer)                              )*conv
 end do  ! (looping through layers)
 write(*,'(a8,1x,10(e20.10,1x))'), 'a_tri = ', a_tri
 write(*,'(a8,1x,10(e20.10,1x))'), 'b_tri = ', b_tri
 write(*,'(a8,1x,10(e20.10,1x))'), 'c_tri = ', c_tri
 ! solve the tridiagonal system of equations -- returns mLayerDelTemp
 call tridag(a_tri,b_tri,c_tri,-r_vec,mLayerDelTemp,err,cmessage)
 if(err/=0)then; message=trim(message)//trim(cmessage); return; endif
 print*, 'in energyflux, mLayerDelTemp = ', mLayerDelTemp
 ! deallocate vectors
 deallocate(mLayerRadCondFluxIter,iLayerTempDeriv,jmat,a_tri,b_tri,c_tri,r_vec,stat=err)
 if(err/=0)then; err=30; message=trim(message)//'problem deallocating vectors for tridag solution'; return; endif
 ! ====================================================================================================================

 end subroutine tempchange


 ! ************************************************************************************************
 ! new subroutine: compute finite-difference approximation of the Jacobian matrix
 ! ************************************************************************************************
 SUBROUTINE fDiffJacbn(dt,x_in,jmat,err,message)
 IMPLICIT NONE
 REAL(DP), INTENT(IN)    :: dt
 REAL(DP), DIMENSION(:), INTENT(IN)    :: x_in
 REAL(DP), DIMENSION(:,:), INTENT(OUT) :: jmat
 REAL(DP), PARAMETER :: EPS=-1.0e-8_dp
 INTEGER(I4B) :: j,n
 REAL(DP), DIMENSION(size(x_in)) :: fvec,ftest,xsav,xtry,xph,h
 integer(i4b),intent(out)        :: err                    ! error code
 character(*),intent(out)        :: message                ! error message
 character(LEN=256)              :: cmessage               ! error message of downwind routine 
 err=0; message="fDiffJacbn/"
 ! check that the size of the vectors match
 if ( all((/size(fvec),size(jmat,1),size(jmat,2)/) == size(x_in)) ) then
  n=size(x_in)
 else
  err=20; message=trim(message)//"sizeMismatch"; return
 endif
 ! compute initial function vector
 call res_discrete(dt,x_in,fvec,err,cmessage)
 if(err/=0)then; err=10; message=trim(message)//trim(cmessage); return; endif
 ! compute Jacobian
 xtry=x_in
 xsav=xtry
 h=EPS*abs(xsav)
 where (h == 0.0) h=EPS
 xph=xsav+h
 h=xph-xsav
 do j=1,n
  xtry(j)=xph(j)
  call res_discrete(dt,xtry,ftest,err,cmessage)
  if(err/=0)then; err=10; message=trim(message)//trim(cmessage); return; endif
  jmat(:,j)=(ftest(:)-fvec(:))/h(j)
  xtry(j)=xsav(j)
 end do
 end subroutine fDiffJacbn


 ! ************************************************************************************************
 ! new subroutine: compute residual vector
 ! ************************************************************************************************
 subroutine res_discrete(dt,mLayerTempTrial,rvec,err,message)
 ! compute the residual temperature vector
 USE multiconst,only:&
                     Tfreeze,    &                         ! freezing point              (K)
                     gravity,    &                         ! acceleration of gravity     (m s-2)
                     iden_water, &                         ! intrinsic density of water  (kg m-3)
                     LH_fus                                ! latent heat of fusion       (J kg-1)
 USE data_struc,only:mvar_data,indx_data                   ! data structures
 USE var_lookup,only:iLookMVAR                             ! named variables for structure elements
 implicit none
 ! dummy variables
 real(dp),intent(in)            :: dt                      ! time step (seconds)
 real(dp),intent(in)            :: mLayerTempTrial(:)      ! temperature vector (K)
 real(dp),intent(out)           :: rvec(:)                 ! residual vector (K)
 integer(i4b),intent(out)       :: err                     ! error code
 character(*),intent(out)       :: message                 ! error message
 ! local pointers to model diagnostic variables
 real(dp),pointer              :: volLatHt_fus             ! volumetric latent heat of fusion (J m-3 K-1)
 real(dp),pointer              :: mLayerTemp(:)            ! temperature of each layer (K)
 real(dp),pointer              :: mLayerVolHtCapBulk(:)    ! bulk volumetric heat capacity (J m-3 K-1)
 real(dp),pointer              :: mLayerdTheta_dTk(:)      ! derivative in the freezing curve (K-1)
 ! local variables
 real(dp),dimension(size(mLayerTempTrial)) :: dEnergy_dt   ! change in energy over layer (J m-3 s-1)
 real(dp),dimension(size(mLayerTempTrial)) :: AppHtCap     ! apparent heat capacity (J m-3 K-1)
 character(LEN=256)                        :: cmessage     ! error message of downwind routine
 err=0; message="res_discrete/"

 ! assign pointers to model diagnostic variables
 volLatHt_fus        => mvar_data%var(iLookMVAR%scalarvolLatHt_fus)%dat(1)   ! volumetric latent heat of fusion (J m-3)
 mLayerTemp          => mvar_data%var(iLookMVAR%mLayerTemp)%dat              ! temperature of each layer (K)
 mLayerVolHtCapBulk  => mvar_data%var(iLookMVAR%mLayerVolHtCapBulk)%dat      ! bulk volumetric heat capacity (J m-3 K-1)
 mLayerdTheta_dTk    => mvar_data%var(iLookMVAR%mLayerdTheta_dTk)%dat        ! derivative in the freezing curve (K-1)
 ! compute change in energy over the layer (J m-3 s-1)
 call intrnlFlux(mLayerTempTrial,dEnergy_dt,err,cmessage)
 if(err/=0)then; err=10; message=trim(message)//trim(cmessage); return; endif
 ! compute apparent heat capacity (J m-3 K-1)
 AppHtCap = mLayerVolHtCapBulk + iden_water*LH_fus*mLayerdTheta_dTk
 ! compute the residual temperature (K)
 rvec = mLayerTempTrial - (mLayerTemp + (dt*dEnergy_dt)/AppHtCap)
 end subroutine res_discrete


 ! ************************************************************************************************
 ! new subroutine: compute fluxes within the snow-soil system, and the residual energy
 ! ************************************************************************************************
 subroutine intrnlFlux(mLayerTempIter,mLayerRadCondFluxIter,err,message)
 USE multiconst,only:Tfreeze,    & ! temperature at the freezing point (K)
                     gravity,    & ! gravitational acceleration (m s-2)
                     LH_fus        ! latent heat of fusion (J kg-1)
 USE data_struc,only:mpar_data,forc_data,mvar_data,indx_data,ix_soil,ix_snow    ! data structures
 USE var_lookup,only:iLookPARAM,iLookFORCE,iLookMVAR,iLookINDEX                 ! named variables for structure elements
 ! compute energy fluxes within the domain
 implicit none
 ! dummy variables
 real(dp),intent(in)           :: mLayerTempIter(:)        ! trial temperature of each layer (K)
 real(dp),intent(out)          :: mLayerRadCondFluxIter(:) ! radiation and conductive flux (J m-3 s-1)
 integer(i4b),intent(out)      :: err                      ! error code
 character(*),intent(out)      :: message                  ! error message
 ! local pointers tp model parameters
 real(dp),pointer              :: Fabs_vis                 ! fraction of radiation in visible part of spectrum (-)
 real(dp),pointer              :: rad_ext                  ! extinction coefficient for penetrating sw radiation (m-1)
 ! local pointers to model forcing data
 real(dp),pointer              :: sw_down                  ! downward shortwave radiation (W m-2)
 ! local pointers to model state ariables  
 real(dp),pointer              :: scalarAlbedo             ! surface albedo (-)
 real(dp),pointer              :: mLayerTemp(:)            ! temperature of each layer (K)
 ! local pointers to variable short-cuts
 real(dp),pointer              :: volLatHt_fus             ! volumetric latent heat of fusion (J m-3 K-1)
 ! local pointers to model diagnostic variables
 real(dp),pointer              :: mLayerHeight(:)          ! height at the mid-point of each layer (m)
 real(dp),pointer              :: iLayerHeight(:)          ! height at the interface of each layer (m)
 real(dp),pointer              :: mLayerVolHtCapBulk(:)    ! bulk volumetric heat capacity (J m-3 K-1)
 real(dp),pointer              :: iLayerThermalC(:)        ! thermal conductivity at the interface of each layer (W m-1 K-1)
 real(dp),pointer              :: iLayerRadCondFlux(:)     ! energy flux at the interface of each layer (W m-2)
 real(dp),pointer              :: mLayerdTheta_dPsi(:)     ! derivative in the soil water characteristic (m-1)
 ! local pointers to model index variables
 integer(i4b),pointer          :: nLayers                  ! number of layers
 integer(i4b),pointer          :: layerType(:)             ! type of the layer (ix_soil or ix_snow)
 ! define local variables
 character(LEN=256)            :: cmessage                 ! error message of downwind routine 
 real(dp)                      :: sw__net                  ! net shortwave radiation (W m-2)
 real(dp)                      :: sw__vis                  ! shortwave radiation in the visible part of the spectrum
 real(dp)                      :: depth_layr               ! depth at the layer interface (m)
 real(dp)                      :: condv_flux               ! conductive heat flux (W m-2)
 real(dp)                      :: swrad_flux               ! radiative heat flux (W m-2)
 integer(i4b)                  :: iLayer                   ! index of model layer
 ! initialize error control
 err=0; message="f-energyflux/"
 ! assign pointers to model parameters
 Fabs_vis            => mpar_data%var(iLookPARAM%Fabs_vis)                   ! fraction of radiation in visible part of spectrum (-)
 rad_ext             => mpar_data%var(iLookPARAM%rad_ext)                    ! extinction coefficient for penetrating sw radiation (m-1)
 ! assign pointers to model forcing variables
 sw_down             => forc_data%var(iLookFORCE%sw_down)                    ! downward shortwave radiation (W m-2)
 ! assign pointers to model state variables
 scalarAlbedo        => mvar_data%var(iLookMVAR%scalarAlbedo)%dat(1)         ! surface albedo (-)
 mLayerTemp          => mvar_data%var(iLookMVAR%mLayerTemp)%dat              ! temperature of each layer (K)
 ! assign pointers to variable short-cuts
 volLatHt_fus        => mvar_data%var(iLookMVAR%scalarvolLatHt_fus)%dat(1)   ! volumetric latent heat of fusion (J m-3)
 ! assign pointers to model diagnostic variables
 mLayerHeight        => mvar_data%var(iLookMVAR%mLayerHeight)%dat            ! height at the mid-point of each layer (m)
 iLayerHeight        => mvar_data%var(iLookMVAR%iLayerHeight)%dat            ! height at the interface of each layer (m)
 iLayerThermalC      => mvar_data%var(iLookMVAR%iLayerThermalC)%dat          ! thermal conductivity at the interface of each layer (W m-1 K-1)
 iLayerRadCondFlux   => mvar_data%var(iLookMVAR%iLayerRadCondFlux)%dat       ! energy flux at the interface of each layer (W m-2)
 mLayerVolHtCapBulk  => mvar_data%var(iLookMVAR%mLayerVolHtCapBulk)%dat      ! bulk volumetric heat capacity (J m-3 K-1)
 mLayerdTheta_dPsi   => mvar_data%var(iLookMVAR%mLayerdTheta_dPsi)%dat       ! derivative in the soil water characteristic (m-1)
 ! assign pointers to index variables
 nLayers             => indx_data%var(iLookINDEX%nLayers)%dat(1)             ! number of layers
 layerType           => indx_data%var(iLookINDEX%layerType)%dat              ! layer type (ix_soil or ix_snow)

 ! prescribe zero flux at the bottom of the domain
 iLayerRadCondFlux(0) = 0._dp
 ! compute diagnostic energy variables (thermal conductivvity and volumetric heat capacity)
 call diagn_evar(err,cmessage)
 if(err/=0)then; err=10; message=trim(message)//trim(cmessage); return; endif
 ! compute the flux at the top of the domain
 call surfceFlux(mLayerTempIter(nLayers),iLayerRadCondFlux(nLayers),err,cmessage)
 if(err/=0)then; err=10; message=trim(message)//trim(cmessage); return; endif
 ! separate absorbed shortwave radiation into visible and near infra-red components
 sw__net = sw_down * (1._dp - scalarAlbedo)
 sw__vis = Fabs_vis * sw__net  ! sw radiation in the visible part of the spectrum
 ! compute the energy fluxes at layer interfaces (W m-2) -- the **top** of each layer
 !   NOTE: fluxes are positive upward
 do iLayer=1,nLayers
  !  compute fluxes at **top** of each layer (flux at top of the domain is computed earlier)
  if(ilayer<nLayers)then
   ! compute conductive flux (W m-2)
   condv_flux = iLayerThermalC(iLayer) * ( (mLayerTempIter(iLayer+1)  - mLayerTempIter(iLayer)) /  &
                                           (mLayerHeight(iLayer+1)    - mLayerHeight(iLayer))   )
   ! compute shortwave radiative flux (W m-2)
   if(mLayerHeight(iLayer+1)>0._dp)then ! snow is above iLayerHeight(iLayer)
    depth_layr = iLayerHeight(nLayers) - iLayerHeight(iLayer) ! depth of the top of the layer
    swrad_flux = sw__vis * exp(-rad_ext*depth_layr)           ! flux at the top of the layer
   else
    swrad_flux = 0._dp  ! soil is above iLayerHeight(iLayer)
   endif
   ! compute total energy flux at the top of each layer (W m-2)
   iLayerRadCondFlux(iLayer) = condv_flux + swrad_flux
  endif
  !print*, 'flux = ', iLayer, iLayerRadCondFlux(iLayer)
  ! compute temporal change in layer energy (J m-3 s-1)
  mLayerRadCondFluxIter(iLayer)  = (iLayerRadCondFlux(iLayer) - iLayerRadCondFlux(iLayer-1)) &
                                 / (iLayerHeight(iLayer) - iLayerHeight(iLayer-1))
 end do  ! looping through layers
 end subroutine intrnlFlux





 ! ************************************************************************************************
 ! new subroutine: compute surface energy balance
 ! ************************************************************************************************
 subroutine surfceFlux(surfTemp,surfFlux,err,message)
 USE multiconst,only:ave_slp,  & ! mean sea level pressure     (Pa) 
                     R_da,     & ! gas constant for dry air    (Pa m3 kg-1 K-1; [J = Pa m3])
                     Cp_air,   & ! specific heat of air        (J kg-1 K-1)
                     iden_air, & ! intrinsic density of air    (kg m-3)
                     gravity,  & ! acceleration of gravity     (m s-2)
                     LH_sub,   & ! latent heat of sublimation  (J kg-1)
                     LH_vap,   & ! latent heat of vaporization (J kg-1)
                     Em_Sno,   & ! emissivity of snow          (-)
                     sigma       ! Stefan Boltzman constant    (W m-2 K-4)
 USE conv_funcs_module,only:relhm2sphm                                          ! compute specific humidity
 USE snow_utils_module,only:astability                                          ! compute atmospheric stability
 USE soil_utils_module,only:RH_soilair                                          ! compute relative humidity of air in soil pores
 USE tridagSolv_module,only:tridag                                              ! solve tridiagonal system of equations
 USE data_struc,only:mpar_data,forc_data,mvar_data,indx_data,ix_soil,ix_snow    ! data structures
 USE var_lookup,only:iLookPARAM,iLookFORCE,iLookMVAR,iLookINDEX                 ! named variables for structure elements
 ! compute energy fluxes within the domain
 implicit none
 ! dummy variables
 real(dp),intent(in)           :: surfTemp               ! trial surface temperature (K)
 real(dp),intent(out)          :: surfFlux               ! computed surface flux (W m-2)
 integer(i4b),intent(out)      :: err                    ! error code
 character(*),intent(out)      :: message                ! error message
 ! local pointers tp model parameters
 real(dp),pointer              :: mheight                ! measurement height (m)
 real(dp),pointer              :: minwind                ! minimum windspeed (m s-1)
 ! local pointers to model forcing data
 real(dp),pointer              :: sw_down                ! downward shortwave radiation (W m-2)
 real(dp),pointer              :: lw_down                ! downward longwave radiation (W m-2)
 real(dp),pointer              :: airtemp                ! air temperature at 2 meter height (K)
 real(dp),pointer              :: windspd                ! wind speed at 10 meter height (m s-1)
 real(dp),pointer              :: airpres                ! air pressure at 2 meter height (Pa) 
 real(dp),pointer              :: spechum                ! specific humidity at 2 meter height (g g-1)
 ! local pointers to model state ariables
 real(dp),pointer              :: scalarAlbedo           ! surface albedo (-)
 real(dp),pointer              :: mLayerMatricHead(:)    ! matric head of each layer (m)
 ! local pointers to variable short-cuts
 real(dp),pointer              :: volHtCap_air           ! volumetric heat capacity of air (J m-3 K-1)
 ! local pointers to model diagnostic variables
 real(dp),pointer              :: scalarSenHeat          ! sensible heat flux at the surface (W m-2)
 real(dp),pointer              :: scalarLatHeat          ! latent heat flux at the surface (W m-2)
 real(dp),pointer              :: scalarMassLiquid       ! evaporation/dew (kg m-2 s-1)
 real(dp),pointer              :: scalarMassSolid        ! sublimation/frost (kg m-2 s-1)
 real(dp),pointer              :: mLayerVolFracLiq(:)    ! volumetric fraction of liquid water (-)
 real(dp),pointer              :: scalarExCoef           ! turbulent exchange coefficient (-)
 ! local pointers to model index variables
 integer(i4b),pointer          :: nLayers
 integer(i4b),pointer          :: layerType(:)           ! type of the layer (ix_soil or ix_snow)
 ! define local variables
 character(LEN=256)            :: cmessage               ! error message of downwind routine 
 real(dp)                      :: ptcorr                 ! potential temperature correction
 real(dp)                      :: PTsurf                 ! surface *potential* temperaure
 real(dp)                      :: PT_air                 ! air *potential* temperature
 real(dp)                      :: RHsurf                 ! relative humidity of the top layer (-)
 real(dp)                      :: Q_surf                 ! specific humidity at the surface (g/g)
 real(dp)                      :: T_grad                 ! potential temperature gradient between air and surface
 real(dp)                      :: T_mean                 ! mean temperature between air and surface
 real(dp)                      :: RiBulk                 ! bulk Richardson number
 real(dp)                      :: LW__up                 ! upward longwave radiation, positive downward (W m-2)
 ! initialize error control
 err=0; message="f-energyflux/"
 ! assign pointers to model parameters
 mheight             => mpar_data%var(iLookPARAM%mheight)                    ! measurement height (m)
 minwind             => mpar_data%var(iLookPARAM%minwind)                    ! minimum windspeed (m s-1)
 ! assign pointers to model forcing variables
 sw_down             => forc_data%var(iLookFORCE%sw_down)                    ! downward shortwave radiation (W m-2)
 lw_down             => forc_data%var(iLookFORCE%lw_down)                    ! downward longwave radiation (W m-2)
 airtemp             => forc_data%var(iLookFORCE%airtemp)                    ! air temperature at 2 meter height (K)
 windspd             => forc_data%var(iLookFORCE%windspd)                    ! wind speed at 10 meter height (m s-1)
 airpres             => forc_data%var(iLookFORCE%airpres)                    ! air pressure at 2 meter height (Pa)
 spechum             => forc_data%var(iLookFORCE%spechum)                    ! specific humidity at 2 meter height (g g-1)
 ! assign pointers to model state variables
 scalarAlbedo        => mvar_data%var(iLookMVAR%scalarAlbedo)%dat(1)         ! surface albedo (-)
 mLayerMatricHead    => mvar_data%var(iLookMVAR%mLayerMatricHead)%dat        ! matric head of each layer (m)
 ! assin pointers to variable short-cuts
 volHtCap_air        => mvar_data%var(iLookMVAR%scalarVolHtCap_air)%dat(1)   ! volumetric heat capacity of air (J m-3 K-1)
 ! assign pointers to model diagnostic variables
 scalarSenHeat       => mvar_data%var(iLookMVAR%scalarSenHeat)%dat(1)        ! sensible heat flux at the surface (W m-2)
 scalarLatHeat       => mvar_data%var(iLookMVAR%scalarLatHeat)%dat(1)        ! latent heat flux at the surface (W m-2)
 scalarMassLiquid    => mvar_data%var(iLookMVAR%scalarMassLiquid)%dat(1)     ! evaporation/dew (kg m-2 s-1)
 scalarMassSolid     => mvar_data%var(iLookMVAR%scalarMassSolid)%dat(1)      ! sublimation/frost (kg m-2 s-1)
 mLayerVolFracLiq    => mvar_data%var(iLookMVAR%mLayerVolFracLiq)%dat        ! volumetric fraction of liquid water (-)
 scalarExcoef        => mvar_data%var(iLookMVAR%scalarExCoef)%dat(1)         ! turbulent exchange coefficient (-)
 ! assign pointers to index variables
 nLayers             => indx_data%var(iLookINDEX%nLayers)%dat(1)             ! number of layers
 layerType           => indx_data%var(iLookINDEX%layerType)%dat              ! layer type (ix_soil or ix_snow)

 ! compute potential temperature in the atmospheric layer
 ptcorr = (ave_slp/airpres)**(R_da/Cp_air)    ! potential temperature correction
 PTsurf = surfTemp * ptcorr
 PT_air = airtemp  * ptcorr
 ! compute relative humidity for the top layer (assume saturated if snow)
 select case(layerType(nLayers))
  case(ix_soil)
   ! compute the relative humity of air in the soil pores in the top layer (-)
   RHsurf     = RH_soilair(mLayerMatricHead(nLayers),surfTemp)
  case(ix_snow); RHsurf = 1._dp
  case default; err=10; message=trim(message)//"unknownOption"; return
 end select
 ! compute specific humity of the top layer
 Q_surf = relhm2sphm(RHsurf,airpres,surfTemp)
 ! compute stability corrections
 T_grad = PT_air - PTsurf
 T_mean = 0.5_dp * (PT_air + PTsurf)
 RiBulk = 0._dp
 if(windspd>minwind) RiBulk = (T_grad/T_mean) * ((gravity*mheight)/(windspd*windspd))
 call astability(RiBulk,scalarExCoef,err,cmessage)
 if(err/=0)then; err=10; message=trim(message)//trim(cmessage); return; endif
 ! compute sublimation/frost (kg m-2 s-1)
 if(mLayerVolFracLiq(nLayers)==0._dp) then
  scalarMassLiquid= 0._dp
  scalarMassSolid = iden_air * windspd * scalarExCoef * (spechum - Q_surf)
 ! compute evaporation/dew (kg m-2 s-1)
 else
  scalarMassLiquid= iden_air * windspd * scalarExCoef * (spechum - Q_surf)
  scalarMassSolid = 0._dp
 endif
 ! compute sensible and latent heat (W m-2)
 scalarSenHeat = volHtCap_air * windspd * scalarExCoef * (PT_air - PTsurf)      ! W m-2
 scalarLatHeat = LH_vap*scalarMassLiquid + LH_sub*scalarMassSolid         ! W m-2
 ! compute upward longwave radiation -- note sign convention is positive downward
 LW__up = -Em_Sno * sigma * surfTemp**4._dp                               ! W m-2
 ! compute flux at the top of the domain
 surfFlux      = scalarSenHeat + scalarLatHeat + sw_down * (1._dp - scalarAlbedo) + lw_down + LW__up  ! W m-2
 end subroutine surfceFlux 

 
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
 end subroutine diagn_evar


end module energyflux_module

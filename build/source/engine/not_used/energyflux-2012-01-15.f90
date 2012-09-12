module energyflux_module
USE nrtype
implicit none
private
public::diagn_evar
public::tempchange
contains

 ! ************************************************************************************************
 ! new subroutine: compute change in temperature over the time step
 ! ************************************************************************************************
 subroutine tempchange(dt,                   & ! time step (seconds)
                       mLayerTempIter,       & ! trial temperature at the current iteration (K)
                       mLayerVolFracIceIter, & ! volumetric fraction of ice at the current iteration (-)
                       mLayerVolFracLiqIter, & ! volumetric fraction of liquid water at the current iteration (-)
                       mLayerMatricHeadIter, & ! matric head at the current iteration (m)
                       mLayerTempDiff,       & ! iteration increment for temperature (K)
                       mLayerTempNew,        & ! new temperature (K)
                       mLayerVolFracIceNew,  & ! new volumetric fraction of ice (-)
                       mLayerVolFracLiqNew,  & ! new volumetric fraction of liquid water (-)
                       mLayerMatricHeadNew,  & ! new matric head (m)
                       err,message)            ! error control
 USE multiconst,only:gravity,  & ! acceleration of gravity     (m s-2)
                     Tfreeze,  & ! freezing point              (K)
                     LH_fus,   & ! latent heat of fusion       (J kg-1)
                     LH_vap,   & ! latent heat of vaporization (J kg-1)
                     LH_sub,   & ! latent heat of sublimation  (J kg-1)
                     Em_Sno,   & ! emissivity of snow          (-)
                     sigma,    & ! Stefan Boltzman constant    (W m-2 K-4)
                     iden_air, & ! intrinsic density of air    (kg m-3)
                     iden_ice, & ! intrinsic density of ice    (kg m-3)
                     iden_water  ! intrinsic density of water  (kg m-3)
 USE conv_funcs_module,only:relhm2sphm    ! compute specific humidity
 USE snow_utils_module,only:astability    ! compute atmospheric stability
 USE snow_utils_module,only:fracliquid    ! compute volumetric fraction of liquid water
 USE snow_utils_module,only:dFracLiq_dTk  ! differentiate the freezing curve w.r.t. temperature (snow)
 USE soil_utils_module,only:RH_soilair    ! compute relative humidity of air in soil pores
 USE soil_utils_module,only:volFracLiq    ! compute volumetric fraction of liquid water based on matric head
 USE soil_utils_module,only:matricHead    ! compute the matric head based on volumetric liquid water content
 USE soil_utils_module,only:crit_soilT    ! compute the critical soil temperature above which all water is unfrozen
 USE soil_utils_module,only:dTheta_dTk    ! differentiate the freezing curve w.r.t. temperature (soil)
 USE tridagSolv_module,only:tridag                                              ! solve tridiagonal system of equations
 USE data_struc,only:mpar_data,forc_data,mvar_data,indx_data,ix_soil,ix_snow    ! data structures
 USE var_lookup,only:iLookPARAM,iLookFORCE,iLookMVAR,iLookINDEX                 ! named variables for structure elements
 ! compute change in temperature over the time step
 implicit none
 ! dummy variables
 real(dp),intent(in)           :: dt                       ! time step (seconds)
 real(dp),intent(in)           :: mLayerTempIter(:)        ! trial temperature at the current iteration (K)
 real(dp),intent(in)           :: mLayerVolFracIceIter(:)  ! volumetric fraction of ice at the current iteration (-)
 real(dp),intent(in)           :: mLayerVolFracLiqIter(:)  ! volumetric fraction of liquid water at the current iteration (-)
 real(dp),intent(in)           :: mLayerMatricHeadIter(:)  ! matric head at the current iteration (m)
 real(dp),intent(out)          :: mLayerTempDiff(:)        ! iteration increment for temperature (K) 
 real(dp),intent(out)          :: mLayerTempNew(:)         ! new temperature (K)
 real(dp),intent(out)          :: mLayerVolFracIceNew(:)   ! new volumetric fraction of ice (-)
 real(dp),intent(out)          :: mLayerVolFracLiqNew(:)   ! new volumetric fraction of liquid water (-)
 real(dp),intent(out)          :: mLayerMatricHeadNew(:)   ! new matric head (m)
 integer(i4b),intent(out)      :: err                      ! error code
 character(*),intent(out)      :: message                  ! error message
 ! local pointers to model parameters
 real(dp),pointer              :: mheight                  ! measurement height (m)
 real(dp),pointer              :: minwind                  ! minimum windspeed (m s-1)
 real(dp),pointer              :: Fabs_vis                 ! fraction of radiation in visible part of spectrum (-)
 real(dp),pointer              :: rad_ext                  ! extinction coefficient for penetrating sw radiation (m-1)
 real(dp),pointer              :: snowfrz_scale            ! scaling parameter for the snow freezing curve (K-1)
 real(dp),pointer              :: lowerBoundTemp           ! temperature of the lower boundary (K)
 ! local pointers to soil parameters
 real(dp),pointer              :: vGn_alpha                ! van Genutchen "alpha" parameter
 real(dp),pointer              :: vGn_n                    ! van Genutchen "n" parameter
 real(dp),pointer              :: theta_sat                ! soil porosity (-)
 real(dp),pointer              :: theta_res                ! soil residual volumetric water content (-)
 real(dp),pointer              :: k_soil                   ! hydraulic conductivity (m s-1)
 real(dp),pointer              :: spec_storage             ! specific storage coefficient (m-1)
 real(dp),pointer              :: f_impede                 ! ice impedence factor (-)
 ! local pointers to vegetation parameters
 real(dp),pointer              :: LAI                      ! leaf area index (m2 m-2)
 real(dp),pointer              :: maxStomatalConduct       ! maximum stomatal conmductance (m s-1)
 real(dp),pointer              :: psi_closed               ! matric head when the stomata are fully closed (m)
 ! local pointers to derived model variables that are constant over the simulation period
 real(dp),pointer              :: vGn_m                    ! van Genutchen "m" parameter (-)
 real(dp),pointer              :: kappa                    ! constant in the freezing curve function (m K-1)
 real(dp),pointer              :: mLayerRootDensity(:)     ! fraction of roots in each soil layer (-)
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
 real(dp),pointer              :: iLayerHeight(:)          ! height at the interface of each layer (m)
 real(dp),pointer              :: mLayerVolHtCapBulk(:)    ! bulk volumetric heat capacity (J m-3 K-1)
 real(dp),pointer              :: mLayerThermalC(:)        ! thermal conductivity at the mid-point of each layer (W m-1 K-1)
 real(dp),pointer              :: iLayerThermalC(:)        ! thermal conductivity at the interface of each layer (W m-1 K-1)
 real(dp),pointer              :: mLayerdTheta_dTk(:)      ! derivative in the freezing curve (K-1)
 real(dp),pointer              :: mLayerEvapLimitFac(:)    ! moisture avail factor limiting transpiration in each layer (-)
 real(dp),pointer              :: mLayerTranspire(:)       ! transpiration loss from each soil layer (kg m-2 s-1)
 ! local pointers to diagnostic scalar variables
 real(dp),pointer              :: scalarPotentialET        ! potential ET (kg m-2 s-1)
 real(dp),pointer              :: scalarMassLiquid         ! evaporation/dew (kg m-2 s-1)
 real(dp),pointer              :: scalarMassSolid          ! sublimation/frost (kg m-2 s-1)
 real(dp),pointer              :: scalarSenHeat            ! sensible heat flux at the surface (W m-2)
 real(dp),pointer              :: scalarLatHeat            ! latent heat flux at the surface (W m-2)
 real(dp),pointer              :: scalarExCoef             ! turbulent exchange coefficient (-)
 ! local pointers to model variables that are constant over the simulation period
 real(dp),pointer              :: volHtCap_air             ! volumetric heat capacity of air (J m-3 K-1)
 real(dp),pointer              :: volLatHt_fus             ! volumetric latent heat of fusion (J m-3 K-1)
 ! local pointers to model index variables
 integer(i4b),pointer          :: nLayers                  ! number of layers
 integer(i4b),pointer          :: layerType(:)             ! type of the layer (ix_soil or ix_snow)
 ! define local variables
 character(LEN=256)            :: cmessage                 ! error message of downwind routine
 real(dp)                      :: surfTemp                 ! short-cut to surface temperature -- assume temperature of the top layer (K)
 real(dp)                      :: T_grad                   ! gradient in temperature between the atmosphere and surface (K)
 real(dp)                      :: T_mean                   ! mean of the atmosphere and surface temperature (K)
 real(dp)                      :: RiBulk                   ! bulk Richardson number (-)
 real(dp),parameter            :: RHsurf=1._dp             ! relative humidity for the top layer (assume saturated)
 real(dp)                      :: aerodynConduct           ! aerodynamic conductance (m s-1)
 real(dp)                      :: stomatalConduct          ! stomatal conductance (m s-1)
 real(dp)                      :: plantWiltFactor          ! plant wilting factor (-)
 real(dp)                      :: evapRatio                ! ratio of actual:potential evapotranspiration (-)
 real(dp)                      :: Qderiv                   ! derivative in specific humidity w.r.t. temperature (kg kg-1 K-1)
 real(dp)                      :: ExSen                    ! exchange of sensible heat (J m-2 s-1 K-1
 real(dp)                      :: ExLat                    ! exchange of latent heat (J m-2 s-1)
 real(dp)                      :: LW_deriv                 ! derivative in longwave radiation w.r.t. temperature (J m-2 s-1 K-1)
 real(dp)                      :: Qh_deriv                 ! derivative in sensible heat w.r.t. temperature (J m-2 s-1 K-1)
 real(dp)                      :: Qe_deriv                 ! derivative in latent heat w.r.t. temperature (J m-2 s-1 K-1)
 real(dp),parameter            :: epsT=1.d-10              ! offset from Tcrit when re-setting iterations at the critical temperature (K)
 real(dp)                      :: theta                    ! liquid water equivalent of total water [liquid water + ice] (-)
 real(dp)                      :: Tcrit                    ! critical temperature above which all water is unfrozen (K)
 real(dp)                      :: del_liq                  ! change in volumetric fraction of liquid water over the iteration (-)
 real(dp),dimension(size(mLayerTempIter)-1) :: dz_node     ! distance between the mid-point of a layer and a neighbouring layer (m)
 real(dp),dimension(size(mLayerTempIter))   :: dt_dudz     ! the dt_dudz terms for the upper interface (s m-2)
 real(dp),dimension(size(mLayerTempIter))   :: dt_dldz     ! the dt_dudz terms for the lower interface (s m-2)
 real(dp),dimension(size(mLayerTempIter))   :: hfusion     ! the fusion term (J m-3 K-1)
 real(dp),dimension(0:size(mLayerTempIter)) :: mLayerFlux  ! fluxes at later interfaces (J m-2 s-1)
 real(dp),dimension(size(mLayerTempIter)-1) :: d_m1        ! sub-diagonal elements of the tridiagonal system (J m-3 K-1)
 real(dp),dimension(size(mLayerTempIter))   :: diag        ! diagonal elements of the tridiagonal system (J m-3 K-1)
 real(dp),dimension(size(mLayerTempIter)-1) :: d_p1        ! super-diagonal elements of the tridiagonal system (J m-3 K-1)
 real(dp),dimension(size(mLayerTempIter))   :: rvec        ! right-hand-side vector (J m-3)
 integer(i4b)                  :: nSoil,nSnow              ! number of soil and snow layers
 integer(i4b)                  :: iLayer                   ! index of model layer
 ! initialize error control
 err=0; message="tempchange/"
 ! assign pointers to model parameters
 mheight             => mpar_data%var(iLookPARAM%mheight)                    ! measurement height (m)
 minwind             => mpar_data%var(iLookPARAM%minwind)                    ! minimum windspeed (m s-1)
 Fabs_vis            => mpar_data%var(iLookPARAM%Fabs_vis)                   ! fraction of radiation in visible part of spectrum (-)
 rad_ext             => mpar_data%var(iLookPARAM%rad_ext)                    ! extinction coefficient for penetrating sw radiation (m-1)
 snowfrz_scale       => mpar_data%var(iLookPARAM%snowfrz_scale)              ! scaling parameter for the snow freezing curve (K-1)
 lowerBoundTemp      => mpar_data%var(iLookPARAM%lowerBoundTemp)             ! temperature of the lower boundary (K)
 ! assign pointers to soil parameters
 vGn_alpha           => mpar_data%var(iLookPARAM%vGn_alpha)                  ! van Genutchen "alpha" parameter (m-1)
 vGn_n               => mpar_data%var(iLookPARAM%vGn_n)                      ! van Genutchen "n" parameter (-)
 theta_sat           => mpar_data%var(iLookPARAM%theta_sat)                  ! soil porosity (-)
 theta_res           => mpar_data%var(iLookPARAM%theta_res)                  ! soil residual volumetric water content (-)
 k_soil              => mpar_data%var(iLookPARAM%k_soil)                     ! hydraulic conductivity (m s-1)
 spec_storage        => mpar_data%var(iLookPARAM%spec_storage)               ! specific storage coefficient (m-1)
 f_impede            => mpar_data%var(iLookPARAM%f_impede)                   ! ice impedence factor (-)
 ! assign pointers to vegetation parameters
 LAI                 => mpar_data%var(iLookPARAM%LAI)                        ! leaf area index (m2 m-2)
 maxStomatalConduct  => mpar_data%var(iLookPARAM%maxStomatalConduct)         ! maximum stomatal conmductance (m s-1)
 psi_closed          => mpar_data%var(iLookPARAM%psi_closed)                 ! matric head when the stomata are fully closed (m)
 ! assign pointers to model variables that are constant over the simulation period
 vGn_m               => mvar_data%var(iLookMVAR%scalarVGn_m)%dat(1)          ! van Genutchen "m" parameter (-)
 kappa               => mvar_data%var(iLookMVAR%scalarKappa)%dat(1)          ! constant in the freezing curve function (m K-1)
 mLayerRootDensity   => mvar_data%var(iLookMVAR%mLayerRootDensity)%dat       ! fraction of roots in each soil layer (-)
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
 iLayerHeight        => mvar_data%var(iLookMVAR%iLayerHeight)%dat            ! height at the interface of each layer (m)
 mLayerThermalC      => mvar_data%var(iLookMVAR%mLayerThermalC)%dat          ! thermal conductivity at the mid-point of each layer (W m-1 K-1)
 iLayerThermalC      => mvar_data%var(iLookMVAR%iLayerThermalC)%dat          ! thermal conductivity at the interface of each layer (W m-1 K-1)
 mLayerVolHtCapBulk  => mvar_data%var(iLookMVAR%mLayerVolHtCapBulk)%dat      ! bulk volumetric heat capacity (J m-3 K-1)
 mLayerdTheta_dTk    => mvar_data%var(iLookMVAR%mLayerdTheta_dTk)%dat        ! derivative in the freezing curve (K-1)
 mLayerEvapLimitFac  => mvar_data%var(iLookMVAR%mLayerEvapLimitFac)%dat      ! moisture avail factor limiting transpiration in each layer (-) 
 mLayerTranspire     => mvar_data%var(iLookMVAR%mLayerTranspire)%dat         ! transpiration loss from each soil layer (kg m-2 s-1)
 ! assign pointers to diagnostic scalar variables
 scalarPotentialET   => mvar_data%var(iLookMVAR%scalarPotentialET)%dat(1)    ! potential ET (kg m-2 s-1)
 scalarMassLiquid    => mvar_data%var(iLookMVAR%scalarMassLiquid)%dat(1)     ! transpiration (kg m-2 s-1)
 scalarMassSolid     => mvar_data%var(iLookMVAR%scalarMassSolid)%dat(1)      ! sublimation/frost (kg m-2 s-1)
 scalarSenHeat       => mvar_data%var(iLookMVAR%scalarSenHeat)%dat(1)        ! sensible heat flux at the surface (W m-2)
 scalarLatHeat       => mvar_data%var(iLookMVAR%scalarLatHeat)%dat(1)        ! latent heat flux at the surface (W m-2)
 scalarExcoef        => mvar_data%var(iLookMVAR%scalarExCoef)%dat(1)         ! turbulent exchange coefficient (-)
 ! assign pointers to model variables that are constant over the simulation period
 volHtCap_air        => mvar_data%var(iLookMVAR%scalarVolHtCap_air)%dat(1)   ! volumetric heat capacity of air (J m-3 K-1)
 volLatHt_fus        => mvar_data%var(iLookMVAR%scalarvolLatHt_fus)%dat(1)   ! volumetric latent heat of fusion (J m-3)
 ! assign pointers to index variables
 nLayers             => indx_data%var(iLookINDEX%nLayers)%dat(1)             ! number of layers
 layerType           => indx_data%var(iLookINDEX%layerType)%dat              ! layer type (ix_soil or ix_snow)
 ! identify the number of snow and soil layers
 nSnow = count(layerType==ix_snow)
 nSoil = count(layerType==ix_soil)

 !print*, 'in tempchange'
 ! print forcing
 !print*, 'sw_down, lw_down, airtemp, windspd, airpres, spechum = ', &
 !         sw_down, lw_down, airtemp, windspd, airpres, spechum

 ! ***** compute the derivative in the freezing curve w.r.t. temperature (K-1)
 do iLayer=1,nLayers
  select case(layerType(iLayer))
   ! ***** process snow layers *****
   case(ix_snow)
    mLayerdTheta_dTk(iLayer) = dFracLiq_dTk(mLayerTempIter(iLayer),snowfrz_scale)
   ! ***** process soil layers *****
   case(ix_soil)
    if(mLayerVolFracIceIter(iLayer)>0._dp)then
     mLayerdTheta_dTk(iLayer) = dTheta_dTk(mLayerTempIter(iLayer),theta_res,theta_sat,vGn_alpha,vGn_n,vGn_m)
    else
     mLayerdTheta_dTk(iLayer) = 0._dp
    endif
   ! **** check that the case was identified correctly
   case default
    err=40; message=trim(message)//"cannot identify the layer as snow or soil"; return
  endselect
 end do

 ! compute the dt/dzdz terms for each layer (s m-2)
 dz_node(1:nLayers-1) = mLayerHeight(2:nLayers) - mLayerHeight(1:nLayers-1)
 dt_dudz(2:nLayers)   = dt/(mLayerDepth(2:nLayers)*dz_node)    ! upper distance, valid 2..nlayers
 dt_dldz(1:nLayers-1) = dt/(mLayerDepth(1:nLayers-1)*dz_node)  ! lower distance, valid 1..nlayers-1

 ! compute the dt/dzdz terms at the boundaries
 dt_dudz(1)       = dt/(mLayerDepth(1)       * mLayerDepth(1)*0.5_dp)
 dt_dldz(nLayers) = dt/(mLayerDepth(nLayers) * mLayerDepth(nLayers)*0.5_dp)

 ! define short-cut to surface temperature -- assume temperature of the top layer (K)
 surfTemp= mLayerTempIter(1)

 ! compute stability corrections
 T_grad = airtemp - surfTemp
 T_mean = 0.5_dp*(airtemp + surfTemp)
 if(windspd<minwind) windspd=minwind
 RiBulk = (T_grad/T_mean) * ((gravity*mheight)/(windspd*windspd))
 call astability(RiBulk,scalarExCoef,err,cmessage)
 if(err/=0)then; err=10; message=trim(message)//trim(cmessage); return; endif

 ! compute aerodynamic resistance
 aerodynConduct = scalarExCoef*windspd

 ! compute the ratio of actual:potential evapotranspiration
 if(nSnow>0)then
  evapRatio = 1._dp
 else
  ! ** compute the factor limiting evaporation for each soil layer (-)
  do iLayer=nSnow+1,nLayers
   ! compute the stomatal conductance of the canopy (m s-1)
   plantWiltFactor = min(psi_closed - mLayerMatricHeadIter(iLayer-nSnow), 0._dp)/psi_closed
   stomatalConduct = LAI*maxStomatalConduct*plantWiltFactor
   ! compute the factor limiting evaporation for a given soil layer (-)
   mLayerEvapLimitFac(iLayer-nSnow) = stomatalConduct / (stomatalConduct + aerodynConduct)
  end do ! (looping through soil layers)
  ! ** compute the ratio of actual:potential evapotranspiration
  evapRatio = sum(mLayerEvapLimitFac*mLayerRootDensity)
 endif  ! (if surface is snow-free)

 ! compute the derivative in specific humidity w.r.t. temperature (kg kg-1 K-1)
 Qderiv = dsphum_dTk(RHsurf,airpres,surfTemp) 

 ! compute the exchange factor for sensible and latent heat
 ExSen = volHtCap_air * windspd * scalarExCoef       ! J m-2 s-1 K-1
 if(nSnow >0) ExLat = LH_sub * iden_air * windspd * scalarExCoef              ! J m-2 s-1
 if(nSnow==0) ExLat = LH_vap * iden_air * windspd * scalarExCoef * evapRatio  ! J m-2 s-1

 ! compute sensible and latent heat at the current iteration (W m-2) --> positive downwards
 scalarSenHeat = ExSen*(airtemp - mLayerTempIter(1))
 scalarLatHeat = ExLat*(spechum - relhm2sphm(RHsurf,airpres,mLayerTempIter(1)))

 ! compute the fusion term for each layer (J m-3 K-1)
 hfusion = LH_fus*iden_water*mLayerdTheta_dTk

 ! compute the derivative in longwave radiation (J m-2 s-1 K-1)
 LW_deriv = 4._dp*Em_Sno*sigma*surfTemp**3._dp

 ! compute the derivative in sensible and latent heat (J m-2 s-1 K-1)
 Qh_deriv = ExSen
 Qe_deriv = ExLat*Qderiv

 ! ***** compute fluxes at layer interfaces
 ! compute flux at the upper boundary -- positive downwards
 mLayerFlux(0) = sw_down*(1._dp - scalarAlbedo) + lw_down - Em_Sno*sigma*mLayerTempIter(1)**4._dp + scalarSenHeat + scalarLatHeat
 ! compute fluxes within the domain -- positive downwards
 do iLayer=1,nLayers-1
  mLayerFlux(iLayer) = -iLayerThermalC(iLayer)*(mLayerTempIter(iLayer+1) - mLayerTempIter(iLayer)) / &
                                               (mLayerHeight(iLayer+1) - mLayerHeight(iLayer))
 end do
 ! compute fluxes at the lower boundary [NOTE: use mLayerThermalC(iLayer)]
 mLayerFlux(nLayers) = -mLayerThermalC(iLayer)*(lowerBoundTemp - mLayerTempIter(iLayer))/(mLayerDepth(iLayer)*0.5_dp)

 ! assemble the tri-diagonal system
 do iLayer=1,nLayers
  ! compute the off-diagonal elements
  if(iLayer<nLayers) d_p1(iLayer)   = -dt_dldz(iLayer)*iLayerThermalC(iLayer)
  if(iLayer>1)       d_m1(iLayer-1) = -dt_dudz(iLayer)*iLayerThermalC(iLayer-1)
  ! compute the diagonal elements (J m-3 K-1)
  ! ** lower-most layer
  if(iLayer==nLayers)then
   diag(iLayer) = mLayerVolHtCapBulk(iLayer) + hfusion(iLayer) + &
                  dt_dldz(iLayer)*mLayerThermalC(iLayer) +       &    ! NOTE: use mLayerThermalC
                  dt_dudz(iLayer)*iLayerThermalC(iLayer-1)
  ! ** upper-most layer
  elseif(iLayer==1)then
   diag(iLayer) = mLayerVolHtCapBulk(iLayer) + hfusion(iLayer) + &
                  dt_dldz(iLayer)*iLayerThermalC(iLayer) +       &
                  LW_deriv + Qh_deriv + Qe_deriv
  ! ** intermediate layers
  else
   diag(iLayer) = mLayerVolHtCapBulk(iLayer) + hfusion(iLayer) + &
                  dt_dldz(iLayer)*iLayerThermalC(iLayer) +       &
                  dt_dudz(iLayer)*iLayerThermalC(iLayer-1)
  endif  ! computing diagonal elements
  ! compute the residual vector (J m-3)
  rvec(iLayer) = -(mLayerFlux(iLayer) - mLayerFlux(iLayer-1))*(dt/mLayerDepth(iLayer)) -    &
                 mLayerVolHtCapBulk(iLayer)*(mLayerTempIter(iLayer) - mLayerTemp(iLayer)) + &
                 LH_fus*iden_ice*(mLayerVolFracIceIter(iLayer) - mLayerVolFracIce(iLayer))
 end do  ! (assembling the tri-diagonal system)

 print*, 'temperature residual (J m-3) = ', rvec(1:5)

 ! solve the tridiagonal system of equations -- returns mLayerTempNew
 call tridag(d_m1,diag,d_p1,rvec,mLayerTempDiff,err,cmessage)
 if(err/=0)then; message=trim(message)//trim(cmessage); return; endif

 ! update temperature (K)
 mLayerTempNew = mLayerTempIter + mLayerTempDiff

 ! update volumetric liquid and ice content (-)
 do iLayer=1,nLayers  ! (process snow and soil separately)
  select case(layerType(iLayer))
   ! ** snow
   case(ix_snow)
    ! for crossing cases: re-set iterations at the critical temperature
    if(mLayerTempNew(iLayer)>Tfreeze) mLayerTempNew(iLayer)=Tfreeze
    ! compute the volumetric fraction of liquid water and ice (-)
    mLayerVolFracLiqNew(iLayer) = fracliquid(mLayerTempNew(iLayer),snowfrz_scale)
    mLayerVolFracIceNew(iLayer) = mLayerVolFracIceIter(iLayer) - & 
                                  ( mLayerVolFracLiqNew(iLayer) - mLayerVolFracLiqIter(iLayer) )*(iden_water/iden_ice)
   ! ** soil
   case(ix_soil)
    ! calculate the critical soil temperature above which all water is unfrozen (K)
    theta = mLayerVolFracIceIter(iLayer)*(iden_ice/iden_water) + mLayerVolFracLiqIter(iLayer)
    Tcrit = crit_soilT(theta,theta_res,theta_sat,vGn_alpha,vGn_n,vGn_m)
    ! for crossing cases: re-set iterations at the critical temperature
    if(mLayerTempIter(iLayer)>Tcrit .and. mLayerTempNew(iLayer)<Tcrit) print*,'-ve crossing ', iLayer
    if(mLayerTempIter(iLayer)<Tcrit .and. mLayerTempNew(iLayer)>Tcrit) print*,'+ve crossing ', iLayer
    if(mLayerTempIter(iLayer)>Tcrit .and. mLayerTempNew(iLayer)<Tcrit) mLayerTempNew(iLayer)=Tcrit-epsT
    if(mLayerTempIter(iLayer)<Tcrit .and. mLayerTempNew(iLayer)>Tcrit) mLayerTempNew(iLayer)=Tcrit+epsT
    ! compute the matric head (m) volumetric fraction of liquid water and ice (-)
    if(mLayerTempNew(iLayer)<Tcrit)then
     mLayerMatricHeadNew(iLayer-nSnow) = kappa*(mLayerTempNew(iLayer) - Tfreeze)
     mLayerVolFracLiqNew(iLayer)       = volFracLiq(mLayerMatricHeadNew(iLayer-nSnow),&
                                                    vGn_alpha,theta_res,theta_sat,vGn_n,vGn_m)
     mLayerVolFracIceNew(iLayer)       = mLayerVolFracIceIter(iLayer) - &
                                         ( mLayerVolFracLiqNew(iLayer) - mLayerVolFracLiqIter(iLayer) )*(iden_water/iden_ice)
    else
     mLayerMatricHeadNew(iLayer-nSnow) = matricHead(theta,vGn_alpha,theta_res,theta_sat,vGn_n,vGn_m)
     mLayerVolFracLiqNew(iLayer)       = theta
     mLayerVolFracIceNew(iLayer)       = 0._dp
    endif
   case default; err=10; message=trim(message)//'unknown case for model layer'; return
  endselect
  ! sanity check
  if(mLayerVolFracIceNew(iLayer) < 0._dp)then; err=10; message=trim(message)//'volumetric ice content < 0'; return; endif
 end do ! (looping through layers)

 ! compute un-stressed ET (kg m-2 s-1)
 scalarPotentialET = iden_air * windspd * scalarExCoef * (spechum - relhm2sphm(RHsurf,airpres,mLayerTempNew(1)))

 ! compute transpiration (kg m-2 s-1)
 if(nSnow==0)then
  ! ** compute actual transpiration from each soil layer (kg m-2 s-1)
  mLayerTranspire(1:nSoil) =  mLayerEvapLimitFac(1:nSoil)*mLayerRootDensity(1:nSoil)*scalarPotentialET
  ! ** compute total transpiration (kg m-2 s-1)
  scalarMassLiquid = sum(mLayerTranspire)
 else
  mLayerTranspire(1:nSoil) = 0._dp
  scalarMassLiquid = 0._dp
 endif

 ! compute sublimation/frost (kg m-2 s-1)
 if(nSnow==0) scalarMassSolid  = 0._dp
 if(nsnow >0) scalarMassSolid  = scalarPotentialET

 ! update sensible and latent heat (W m-2) --> positive downwards
 scalarSenHeat = ExSen*(airtemp - mLayerTempNew(1))
 scalarLatHeat = ExLat*(spechum - relhm2sphm(RHsurf,airpres,mLayerTempNew(1)))

 ! check that total transpiration matches the latent heat flux
 !print*, 'check ET ', scalarMassLiquid*LH_vap, scalarLatHeat

 ! ====================================================================================================================

 end subroutine tempchange


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
 ! compute the volumetric fraction of air
 where(layerType==ix_snow)
  mLayerVolFracAir = 1._dp - (mLayerVolFracIce + mLayerVolFracLiq)
 elsewhere
  mLayerVolFracAir = theta_sat - (mLayerVolFracIce + mLayerVolFracLiq)
 endwhere
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
  if(layerType(iLayer)==ix_snow)then
   mLayerVolHtCapBulk(iLayer) = volHtCap_ice   * mLayerVolFracIce(iLayer) + & ! ice component
                                volHtCap_water * mLayerVolFracLiq(iLayer) + & ! liquid water component
                                volHtCap_air   * mLayerVolFracAir(iLayer)     ! air component
  else
   mLayerVolHtCapBulk(iLayer) = volHtCap_soil  * (1._dp - theta_sat)      + & ! soil component
                                volHtCap_ice   * mLayerVolFracIce(iLayer) + & ! ice component
                                volHtCap_water * mLayerVolFracLiq(iLayer) + & ! liquid water component
                                volHtCap_air   * mLayerVolFracAir(iLayer)     ! air component
  endif
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


 ! ************************************************************************************************
 ! new function: compute derivative in specific humidity w.r.t. temperature (kg kg-1 K-1)
 ! ************************************************************************************************
 function dsphum_dTk(r_hum, apres, Tk)
 ! compute derivative in specific humidity w.r.t. temperature (kg kg-1 K-1)
 ! (based on Teten's formula))
 USE multiconst,only:TFreeze, &            ! temperature at freezing              (K)
                     satvpfrz,&            ! sat vapour pressure at 273.16K       (Pa)
                     w_ratio               ! molecular ratio water to dry air     (-)
 implicit none
 ! dummies
 real(dp),intent(in)         :: r_hum       ! relative humidity (fraction)
 real(dp),intent(in)         :: apres       ! atmospheric pressure (Pa)
 real(dp),intent(in)         :: Tk          ! temperature (K)
 real(dp)                    :: dsphum_dTk  ! derivative in specific humidity w.r.t. temperature (kg kg-1 K-1)
 ! locals
 real(dp)                    :: Tc          ! temperature (oC)
 real(dp)                    :: pvp         ! partial vapour pressure (Pa)
 real(dp)                    :: pvp_p       ! derivative in partial vapor pressure w.r.t. temperature (Pa K-1)
 real(dp),parameter          :: a= 17.27_dp ! 1st parameter in Teten's formula for partial vapor pressure (-)
 real(dp),parameter          :: b=237.30_dp ! 2nd parameter in Teten's formula for partial vapor pressure (K))
 ! convert temperature to deg C
 Tc  = Tk - Tfreeze
 ! compute the partial vapour pressure (Pa)
 pvp = r_hum * satVpFrz * exp( a*Tc/(b + Tc) )
 ! compute the derivative in partial vapor pressure w.r.t. temperature (Pa/K)
 pvp_p = pvp * (a/(b + Tc) - (a*Tc)/(b + Tc)**2._dp)
 ! compute the derivative in specific humidity w.r.t. temperature (kg kg-1 K-1)
 dsphum_dTk = w_ratio*pvp_p/(apres - (1._dp - w_ratio)*pvp) + &
              w_ratio*(1._dp - w_ratio)*pvp*pvp_p/(apres - (1._dp - w_ratio)*pvp)**2._dp
 end function dsphum_dTk


end module energyflux_module

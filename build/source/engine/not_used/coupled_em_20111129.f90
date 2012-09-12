module coupled_em_module
USE nrtype
implicit none
private
public::coupled_em
contains

 ! ************************************************************************************************
 ! new subroutine: run the coupled energy-mass model for one timestep
 ! ************************************************************************************************
 subroutine coupled_em(dt_init,err,message)
 ! provide access to subroutines
 USE soil_utils_module,only:hydCond         ! compute hydraulic conductivity
 USE soil_utils_module,only:volFracLiq      ! compute volumetric fraction of liquid water
 USE soil_utils_module,only:dTheta_dPsi     ! compute derivative of the soil moisture characteristic
 USE soil_utils_module,only:dTheta_dTk      ! differentiate the freezing curve w.r.t. temperature (soil)
 USE soil_utils_module,only:crit_soilT      ! compute the critical temperature above which all water is unfrozen
 USE snow_utils_module,only:fracliquid      ! fractional liquid water content based on temperature (snow)
 USE snow_utils_module,only:templiquid      ! temperature based on fractional liquid water content (snow)
 USE snow_utils_module,only:dFracLiq_dTk    ! differentiate the freezing curve w.r.t. temperature (snow)
 USE energyflux_module,only:diagn_evar      ! compute diagnostic energy variables -- thermal conductivity and heat capacity
 USE energyflux_module,only:tempchange      ! compute change in temperature over the time step
 USE energyflux_module,only:energy_err      ! error in energy equation (J m-3)
 USE energyflux_module,only:sfcMassFlx      ! compute the mass flux at the surface (kg m-2 s-1)
 USE liquidflux_module,only:liquidflow      ! compute liquid water flow through the snowpack
 USE liquidflux_module,only:masschange      ! compute change in mass over the time step for the soil
 ! provide access to data
 USE multiconst,only:Tfreeze,      & ! temperature at freezing              (K)
                     gravity,      & ! acceleration of gravity              (m s-2)
                     LH_fus,       & ! latent heat of fusion                (J kg-1)
                     lambda_air,   & ! thermal conductivity of air          (W m-1 K-1) 
                     lambda_ice,   & ! thermal conductivity of ice          (W m-1 K-1)
                     lambda_soil,  & ! thermal conductivity of soil         (W m-1 K-1)
                     lambda_water, & ! thermal conductivity of liquid water (W m-1 K-1)
                     Cp_air,       & ! specific heat of air                 (J kg-1 K-1)
                     Cp_ice,       & ! specific heat of ice                 (J kg-1 K-1)
                     Cp_soil,      & ! specific heat of soil                (J kg-1 K-1)
                     Cp_water,     & ! specific heat of liquid water        (J kg-1 K-1)
                     iden_air,     & ! intrinsic density of air             (kg m-3)
                     iden_ice,     & ! intrinsic density of ice             (kg m-3)
                     iden_water      ! intrinsic density of liquid water    (kg m-3)
 USE data_struc,only:forcFileInfo                                     ! extract time step of forcing data
 USE data_struc,only:mpar_data,mvar_data,indx_data,ix_soil,ix_snow    ! data structures
 USE var_lookup,only:iLookPARAM,iLookMVAR,iLookINDEX                  ! named variables for structure elements
 implicit none
 ! define output
 real(dp),intent(inout)               :: dt_init                  ! used to initialize the size of the sub-step
 integer(i4b),intent(out)             :: err                      ! error code
 character(*),intent(out)             :: message                  ! error message
 ! local pointers to snow parameters
 real(dp),pointer                     :: Fcapil                   ! capillary retention as a fraction of the total pore volume (-)
 real(dp),pointer                     :: snowfrz_scale            ! scaling parameter for the snow freezing curve (K-1)
 ! local pointers to soil parameters
 real(dp),pointer                     :: vGn_alpha                ! van Genutchen "alpha" parameter
 real(dp),pointer                     :: vGn_n                    ! van Genutchen "n" parameter
 real(dp),pointer                     :: theta_sat                ! soil porosity (-)
 real(dp),pointer                     :: theta_res                ! soil residual volumetric water content (-)
 real(dp),pointer                     :: k_soil                   ! hydraulic conductivity (m s-1)
 real(dp),pointer                     :: spec_storage             ! specific storage coefficient (m-1)
 real(dp),pointer                     :: f_impede                 ! ice impedence factor (-)
 ! local pointers to derived model variables that are constant over the simulation period
 real(dp),pointer                     :: vGn_m                    ! van Genutchen "m" parameter (-)
 real(dp),pointer                     :: volHtCap_air             ! volumetric heat capacity of air (J m-3 K-1)
 real(dp),pointer                     :: volHtCap_ice             ! volumetric heat capacity of ice (J m-3 K-1)
 real(dp),pointer                     :: volHtCap_soil            ! volumetric heat capacity of dry soil (J m-3 K-1)
 real(dp),pointer                     :: volHtCap_water           ! volumetric heat capacity of liquid water (J m-3 K-1)
 real(dp),pointer                     :: volLatHt_fus             ! volumetric latent heat of fusion (J m-3 K-1)
 ! local pointers to model variables that are constant over the time step
 real(dp),pointer                     :: mLayerDepth(:)           ! depth of the layer (m)
 real(dp),pointer                     :: mLayerHeight(:)          ! height of the layer mid-point (m)
 real(dp),pointer                     :: iLayerHeight(:)          ! height of the layer interface (m)
 ! local pointers to model state variables -- all layers
 real(dp),pointer                     :: mLayerTemp(:)            ! temperature of each layer (K)
 real(dp),pointer                     :: mLayerVolFracIce(:)      ! volumetric fraction of ice in each layer (-)
 real(dp),pointer                     :: mLayerVolFracLiq(:)      ! volumetric fraction of liquid water in each layer (-)
 ! local pointers to model state variables -- soil layers
 real(dp),pointer                     :: mLayerMatricHead(:)      ! matric head in each ***soil*** layer (m)
 ! local pointers to model diagnostic variables -- surface scalars
 real(dp),pointer                     :: scalarFracLiqTop         ! fraction of total water in the top snow/soil layer that is liquid (-)
 real(dp),pointer                     :: scalarMassLiquid         ! evaporation/dew (kg m-2 s-1)
 real(dp),pointer                     :: scalarMassSolid          ! sublimation/frost (kg m-2 s-1)
 real(dp),pointer                     :: scalarSenHeat            ! sensible heat flux at the surface (W m-2)
 real(dp),pointer                     :: scalarLatHeat            ! latent heat flux at the surface (W m-2)
 ! local pointers to model diagnostic variables -- all layers
 real(dp),pointer                     :: mLayerVolFracAir(:)      ! volumetric fraction of air in each layer (-)
 real(dp),pointer                     :: iLayerThermalC(:)        ! thermal conductivity at the interface of each layer (W m-1 K-1)
 real(dp),pointer                     :: mLayerRadCondFlux(:)     ! change in layer energy (J m-3 s-1)
 real(dp),pointer                     :: mLayerVolHtCapBulk(:)    ! bulk volumetric heat capacity (J m-3 K-1)
 real(dp),pointer                     :: mLayerdTheta_dTk(:)      ! derivative in the freezing curve (K-1)
 real(dp),pointer                     :: mLayerMeltFreeze(:)      ! melt/freeze in each layer (kg m-3)
 real(dp),pointer                     :: mLayerInfilFreeze(:)     ! increase in volumetric ice content caused by freezing infiltrating flux (-)
 ! local pointers to model diagnostic variables -- soil layers
 real(dp),pointer                     :: mLayerdTheta_dPsi(:)     ! derivative in the soil water characteristic (m-1)
 ! local pointers to model diagnostic variables -- snow layers
 real(dp),pointer                     :: mLayerPoreSpace(:)       ! pore space in each snow layer (-)
 real(dp),pointer                     :: mLayerRelSat(:)          ! relative saturation in each snow layer (-)
 ! local pointers to model index variables
 integer(i4b),pointer                 :: nLayers                  ! number of layers
 integer(i4b),pointer                 :: layerType(:)             ! type of the layer (ix_soil or ix_snow)



 ! define local model state variables
 real(dp),allocatable                 :: mLayerTempIter(:)        ! temperature vector at the current iteration (K)
 real(dp),allocatable                 :: mLayerTempNew(:)         ! temperature vector at the next iteration (K)
 real(dp),allocatable                 :: mLayerVolFracIceIter(:)  ! volumetric fraction of ice at the current iteration (-)
 real(dp),allocatable                 :: mLayerVolFracIceNew(:)   ! volumetric fraction of ice at the next iteration (-)
 real(dp),allocatable                 :: mLayerVolFracLiqIter(:)  ! volumetric fraction of liquid water at the current iteration (-)
 real(dp),allocatable                 :: mLayerVolFracLiqNew(:)   ! volumetric fraction of liquid water at the next iteration (-)
 real(dp),allocatable                 :: mLayerMatricHeadIter(:)  ! matric head at the current iteration (-)
 real(dp),allocatable                 :: mLayerMatricHeadNew(:)   ! matric head at the next iteration (-)
 
 ! define local model diagnostic variables
 real(dp),allocatable                 :: mLayerInfilFreezeNew(:)  ! increase in volumetric ice content caused by freezing infiltrating flux (-)

 ! define local error monitoring variables
 real(dp),allocatable                 :: mLayerNrgIncrement(:)    ! change in energy from one iteration to the next (J m-3)
 real(dp),allocatable                 :: mLayerMatricIncrement(:) ! relative change in matric head from one iteration to the next (m)
 real(dp),allocatable                 :: mLayerNrgError(:)        ! energy error in each layer (J m-3)

 ! control the length of the sub-step
 real(dp)                             :: dt                   ! length of time step (seconds)
 real(dp)                             :: dt_sub               ! length of the sub-step (seconds)
 real(dp)                             :: dt_done              ! length of time step completed (seconds)
 integer(i4b),parameter               :: n_inc=2              ! minimum number of iterations to increase time step
 integer(i4b),parameter               :: n_dec=3              ! maximum number of iterations to increase time step
 real(dp),parameter                   :: F_inc = 1.02_dp      ! factor used to increase time step
 real(dp),parameter                   :: F_dec = 0.90_dp      ! factor used to decrease time step
 real(dp),parameter                   :: eps   = 1.d-10       ! small increment used at the freezing point
 
 ! define local variables
 character(len=256)                   :: cmessage             ! error message
 logical(lgt)                         :: printflag            ! .true. if want to print
 integer(i4b)                         :: nSnow                ! number of snow layers
 integer(i4b)                         :: nSoil                ! number of soil layers
 integer(i4b)                         :: nsub                 ! number of sub-steps
 integer(i4b)                         :: iter                 ! iteration index
 integer(i4b)                         :: niter                ! number of iterations
 integer(i4b),parameter               :: maxiter=100           ! maximum number of iterations
 real(dp),dimension(1)                :: nrg_max,matric_max   ! maximum change in energy/matric head for a given iteration
 real(dp),parameter                   :: nrg_tol=1.d-0        ! iteration tolerance for energy (J m-3) -- convergence when nrg_max < nrg_tol
 real(dp),parameter                   :: matric_tol=1.d-4     ! iteration tolerance for matric head -- convergence when matric_max < matric_tol
 integer(i4b)                         :: iLayer               ! loop through model layers
 real(dp)                             :: dIce_dT              ! derivative in volumetric ice content w.r.t. temperature (K-1)
 real(dp)                             :: del_ice              ! change in volumetric ice content over the iteration (-)
 real(dp)                             :: del_h2o              ! change in volumetric liquid water content over the iteration (-)
 real(dp)                             :: kappa                ! constant in the freezing curve function (m K-1)
 real(dp)                             :: liqMx                ! maximum possible volumetric fraction of liquid water at a given temperature (-)
 real(dp)                             :: theta                ! volumetric fraction of total water, liquid plus ice (-)
 real(dp)                             :: Tcrit                ! critical soil temperature above which all water is unfrozen (K)
 real(dp)                             :: scalarMassSolidIter  ! sublimation/frost at current iteration (kg m-2 s-1)
 real(dp)                             :: scalarMassSolidNew   ! sublimation/frost at new iteration (kg m-2 s-1) 
 real(dp)                             :: scalarMassLiquidIter ! evaporation/dew at current iteration (kg m-2 s-1)
 real(dp)                             :: scalarMassLiquidNew  ! evaporation/dew at new iteration (kg m-2 s-1) 


 ! initialize error control
 err=0; message="coupled_em/"

 ! assign local pointers to the model index structures
 nLayers           => indx_data%var(iLookINDEX%nLayers)%dat(1)             ! number of layers
 layerType         => indx_data%var(iLookINDEX%layerType)%dat              ! layer type (ix_soil or ix_snow)

 ! identify the number of snow and soil layers
 nSnow = count(layerType==ix_snow)
 nSoil = count(layerType==ix_soil)

 ! assign pointers to snow parameters
 snowfrz_scale     => mpar_data%var(iLookPARAM%snowfrz_scale)  ! scaling parameter for the snow freezing curve (K-1)
 Fcapil            => mpar_data%var(iLookPARAM%Fcapil)         ! capillary retention as a fraction of the total pore volume (-)

 ! assign pointers to soil parameters
 vGn_alpha         => mpar_data%var(iLookPARAM%vGn_alpha)      ! van Genutchen "alpha" parameter (m-1)
 vGn_n             => mpar_data%var(iLookPARAM%vGn_n)          ! van Genutchen "n" parameter (-)
 theta_sat         => mpar_data%var(iLookPARAM%theta_sat)      ! soil porosity (-)
 theta_res         => mpar_data%var(iLookPARAM%theta_res)      ! soil residual volumetric water content (-)
 k_soil            => mpar_data%var(iLookPARAM%k_soil)         ! hydraulic conductivity (m s-1)
 spec_storage      => mpar_data%var(iLookPARAM%spec_storage)   ! specific storage coefficient (m-1)
 f_impede          => mpar_data%var(iLookPARAM%f_impede)       ! ice impedence factor (-)

 ! assign pointers to model variables that are constant over the simulation period
 vGn_m             => mvar_data%var(iLookMVAR%scalarVGn_m)%dat(1)           ! van Genutchen "m" parameter (-)
 volHtCap_air      => mvar_data%var(iLookMVAR%scalarVolHtCap_air)%dat(1)    ! volumetric heat capacity of air (J m-3 K-1)
 volHtCap_ice      => mvar_data%var(iLookMVAR%scalarVolHtCap_ice)%dat(1)    ! volumetric heat capacity of ice (J m-3 K-1)
 volHtCap_soil     => mvar_data%var(iLookMVAR%scalarVolHtCap_soil)%dat(1)   ! volumetric heat capacity of soil (J m-3 K-1)
 volHtCap_water    => mvar_data%var(iLookMVAR%scalarVolHtCap_water)%dat(1)  ! volumetric heat capacity of water (J m-3 K-1)
 volLatHt_fus      => mvar_data%var(iLookMVAR%scalarvolLatHt_fus)%dat(1)    ! volumetric latent heat of fusion (J m-3)

 ! assign pointers to model variables that are constant over the time step
 mLayerDepth       => mvar_data%var(iLookMVAR%mLayerDepth)%dat              ! depth of the layer (m)
 mLayerHeight      => mvar_data%var(iLookMVAR%mLayerHeight)%dat             ! height of the layer mid-point (m)
 iLayerHeight      => mvar_data%var(iLookMVAR%iLayerHeight)%dat             ! height of the layer interface (m)

 ! assign pointers to model state variables -- all layers
 mLayerTemp        => mvar_data%var(iLookMVAR%mLayerTemp)%dat               ! temperature of each layer (K)
 mLayerVolFracIce  => mvar_data%var(iLookMVAR%mLayerVolFracIce)%dat         ! volumetric fraction of ice in each layer (-)
 mLayerVolFracLiq  => mvar_data%var(iLookMVAR%mLayerVolFracLiq)%dat         ! volumetric fraction of liquid water in each layer (-)

 ! assign pointers to model state variables -- soil layers
 mLayerMatricHead  => mvar_data%var(iLookMVAR%mLayerMatricHead)%dat         ! matric head in each **soil** layer (m)

 ! assign pointers to model diagnostic variables -- surface scalars
 scalarFracLiqTop  => mvar_data%var(iLookMVAR%scalarFracLiqTop)%dat(1)      ! fraction of total water in the top snow/soil layer that is liquid (-)
 scalarMassLiquid  => mvar_data%var(iLookMVAR%scalarMassLiquid)%dat(1)      ! evaporation/dew (kg m-2 s-1)
 scalarMassSolid   => mvar_data%var(iLookMVAR%scalarMassSolid)%dat(1)       ! sublimation/frost (kg m-2 s-1)
 scalarSenHeat     => mvar_data%var(iLookMVAR%scalarSenHeat)%dat(1)         ! sensible heat flux at the surface (W m-2)
 scalarLatHeat     => mvar_data%var(iLookMVAR%scalarLatHeat)%dat(1)         ! latent heat flux at the surface (W m-2)

 ! assign pointers to model diagnostic variables -- all layers
 mLayerVolFracAir  => mvar_data%var(iLookMVAR%mLayerVolFracAir)%dat         ! volumetric fraction of air in each layer (-)
 iLayerThermalC    => mvar_data%var(iLookMVAR%iLayerThermalC)%dat           ! thermal conductivity at the interface of each layer (W m-1 K-1)
 mLayerVolHtCapBulk=> mvar_data%var(iLookMVAR%mLayerVolHtCapBulk)%dat       ! bulk volumetric heat capacity (J m-3 K-1)
 mLayerRadCondFlux => mvar_data%var(iLookMVAR%mLayerRadCondFlux)%dat        ! change in layer energy (J m-3 s-1)
 mLayerdTheta_dTk  => mvar_data%var(iLookMVAR%mLayerdTheta_dTk)%dat         ! derivative in the freezing curve (K-1)
 mLayerMeltFreeze  => mvar_data%var(iLookMVAR%mLayerMeltFreeze)%dat         ! melt/freeze in each layer (kg m-3)
 mLayerInfilFreeze => mvar_data%var(iLookMVAR%mLayerInfilFreeze)%dat        ! increase in volumetric ice content caused by freezing infiltrating flux (-)

 ! assign pointers to model diagnostic variables -- soil only
 mLayerdTheta_dPsi => mvar_data%var(iLookMVAR%mLayerdTheta_dPsi)%dat        ! derivative in the soil water characteristic (m-1)

 ! assign pointers to model diagnostic variables -- snow only
 mLayerPoreSpace   => mvar_data%var(iLookMVAR%mLayerPoreSpace)%dat          ! pore space in each **snow** layer
 mLayerRelsat      => mvar_data%var(iLookMVAR%mLayerRelsat)%dat             ! relative saturation in each **snow** layer

 ! initialize print flag
 printflag=.false.

 ! get the length of the time step (seconds)
 dt = forcFileInfo%data_step

 ! define the constant in the freezing curve function (m K-1)
 kappa  = (iden_ice/iden_water)*(LH_fus/(gravity*Tfreeze))  ! NOTE: J = kg m2 s-2

 ! allocate space for state variables at the start and end of the iteration
 allocate(mLayerTempIter(nLayers),      mLayerTempNew(nLayers),       &  ! all layers
          mLayerVolFracIceIter(nLayers),mLayerVolFracIceNew(nLayers), &  ! all layers
          mLayerVolFracLiqIter(nLayers),mLayerVolFracLiqNew(nLayers), &  ! all layers
          mLayerMatricHeadIter(nSoil),  mLayerMatricHeadNew(nSoil),   &  ! **soil layers only**
          stat=err)
 if(err/=0)then; err=20; message='problem allocating space for state variable vectors - 1'; return; endif

 ! allocate space for diagnostic variables
 allocate(mLayerInfilFreezeNew(nLayers),stat=err)
 if(err/=0)then; err=20; message='problem allocating space for diagnostic variable vectors - 1'; return; endif

 ! allocate space for error monitoring
 allocate(mLayerNrgIncrement(nLayers), &  ! (all layers)
          mLayerNrgError(nLayers),     &  ! (all layers)
          mLayerMatricIncrement(nSoil),&  ! NOTE: nSoil
          stat=err)
 if(err/=0)then; err=20; message='problem allocating space for error monitoring'; return; endif

 ! initialize the length of the sub-step
 dt_sub  = min(dt_init,dt)
 dt_done = 0._dp

 ! initialize the number of sub-steps
 nsub=0

 ! initialize state variables at the next iteration
 mLayerTempIter       = mLayerTemp
 mLayerVolFracIceIter = mLayerVolFracIce
 mLayerMatricHeadIter = mLayerMatricHead
 mLayerVolFracLiqIter = mLayerVolFracLiq

 ! loop through sub-steps
 do  ! continuous do statement with exit clause (alternative to "while")

  ! increment the number of sub-steps
  nsub = nsub+1

  ! initialize number of iterations
  niter=0

  ! initialize melt/freeze (kg m-3)
  mLayerMeltFreeze = 0._dp

  ! initialize the change in volumetric ice content associated with freezing the infiltrating liquid water flux (-)
  mLayerInfilFreeze = 0._dp

  ! compute the fraction of liquid water in the top layer
  scalarFracLiqTop = mLayerVolFracLiqIter(1)/(mLayerVolFracLiqIter(1)+mLayerVolFracIceIter(1))
  
  ! compute an initial estimate of vapor transport, for use in the mass equations
  call sfcMassFlx(mLayerTempIter(1),scalarMassSolidNew,scalarMassLiquidNew,err,cmessage)
  if(err/=0)then; err=10; message=trim(message)//trim(cmessage); return; endif
  print*, 'scalarMassLiquidNew, scalarMassSolidNew = ', scalarMassLiquidNew, scalarMassSolidNew
  
  ! update volumetric liquid water and volumetric ice content
  mLayerVolFracIceIter(1) = mLayerVolFracIceIter(1) + (scalarMassSolidNew/iden_ice) * dt/mLayerDepth(1)
  mLayerVolFracLiqIter(1) = mLayerVolFracLiqIter(1) + (scalarMassLiquidNew/iden_water) * dt/mLayerDepth(1)

  ! initialize iterations
  scalarMassSolidIter  = scalarMassSolidNew
  scalarMassLiquidIter = scalarMassLiquidNew
  pause


  !print*, '**************************************************************************'
  !print*, '**************************************************************************'
  !print*, '**************************************************************************'
  !print*, 'new set of iterations....'

  printflag=.false.

  ! iterate
  do iter=1,maxiter

   ! increment number of iterations
   niter=niter+1

   if(niter>10) printflag=.true.
   if(printflag) print*, '****************************'
   if(printflag) print*, 'new iteration, iter = ', iter
   
   !print*, '***** new iteration, iter = ', iter

   ! compute the volumetric fraction of air
   where(layerType==ix_snow)
    mLayerVolFracAir = 1._dp - (mLayerVolFracIceIter + mLayerVolFracLiqIter)
   elsewhere
    mLayerVolFracAir = theta_sat - (mLayerVolFracIceIter + mLayerVolFracLiqIter)
   endwhere
   !print*, 'mLayerVolFracIceIter(1:5) = ', mLayerVolFracIceIter(1:5)
   !print*, 'mLayerVolFracLiqIter(1:5) = ', mLayerVolFracLiqIter(1:5)

   ! compute diagnostic energy variables (thermal conductivity and volumetric heat capacity)
   call diagn_evar(err,cmessage)
   if(err/=0)then; err=10; message=trim(message)//trim(cmessage); return; endif

   ! compute pore space
   mLayerPoreSpace = mLayerVolFracAir + mLayerVolFracLiq

   ! compute the relative saturation *** in each snow layer *** at the start of the sub-step
   if(iter==1)then
    do iLayer=1,nSnow
     if(mLayerVolFracLiq(iLayer)/mLayerPoreSpace(iLayer) > Fcapil) then
      mLayerRelsat(iLayer)    = (mLayerVolFracLiq(iLayer)/mLayerPoreSpace(iLayer) - Fcapil) / (1._dp - Fcapil)
     else
      mLayerRelsat(iLayer) = 0._dp
     endif
    end do
   endif

   ! print temperatures
   !print*, 'mLayerTempIter(1:5) = ', mLayerTempIter(1:5)

   ! compute liquid water flow through the snowpack (also energy gained by freezing the infiltrating flux)
   call liquidflow(dt_sub,&                            ! time step (seconds)
                   mLayerTempIter(1:nSnow),          & ! layer temperature (K)
                   mLayerVolFracLiqIter(1:nSnow),    & ! volumetric fraction of liquid water at the current iteration (-)
                   mLayerVolFracIceIter(1:nSnow),    & ! volumetric fraction of ice at the current iteration (-)
                   mLayerInfilFreezeNew(1:nSnow),    & ! volumetric fraction of infiltrating liquid water flux that freezes (-)
                   mLayerVolFracLiqNew(1:nSnow),     & ! volumetric liquid water content at the next iteration (-)
                   err,cmessage)
   ! check if we need to adjust the time step
   if(err<0)then
    cycle
   endif
   if(err>0)then; err=10; message=trim(message)//trim(cmessage); return; endif
   ! adjust volumetric liquid water and ice content
   mLayerVolFracLiqIter(1:nSnow) = mLayerVolFracLiqNew(1:nSnow)
   mLayerVolFracIceIter(1:nSnow) = mLayerVolFracIceIter(1:nSnow) + (mLayerInfilFreezeNew(1:nSnow) - mLayerInfilFreeze(1:nSnow))

   ! compute melt/freeze in each layer (kg m-3)
   mLayerMeltFreeze(1:nSnow) = mLayerMeltFreeze(1:nSnow) + iden_ice*(mLayerInfilFreezeNew(1:nSnow) - mLayerInfilFreeze(1:nSnow))
   !print*, 're-freeze volume = ', mLayerInfilFreezeNew(1:5) - mLayerInfilFreeze(1:5)
   !pause

   ! re-set mLayerInfilFreeze
   mLayerInfilFreeze = mLayerInfilFreezeNew

   ! compute the derivative in the soil water characteristic (m-1)
   do iLayer=1,nSoil ! NOTE: nSoil
    mLayerdTheta_dPsi(iLayer) = dTheta_dPsi(mLayerMatricHeadIter(iLayer),vGn_alpha,theta_res,theta_sat,vGn_n,vGn_m)
   end do

   ! compute the matric head at the next iteration (note liquid water and ice vectors are defined for all layers)
   call masschange(dt_sub,&                                ! time step (seconds)
                   mLayerMatricHeadIter,                 & ! matric head in each layer at the current iteration (m)
                   mLayerVolFracIceIter(nSnow+1:nLayers),& ! volumetric fraction of ice at the current iteration (-)
                   mLayerVolFracLiqIter(nSnow+1:nLayers),& ! volumetric fraction of liquid water at the current iteration (-)
                   mLayerMatricHeadNew,                  & ! matric head in each layer at the next iteration (m)
                   mLayerVolFracLiqNew(nSnow+1:nLayers), & ! volumetric fraction of liquid water at the next iteration (-)
                   err,cmessage)
   if(err/=0)then; err=10; message=trim(message)//trim(cmessage); return; endif

   !print*, 'mLayerVolFracLiqIter = ', mLayerVolFracLiqIter(1:10)
   !print*, 'mLayerVolFracLiqNew = ',  mLayerVolFracLiqNew(1:10)
   !print*, 'mLayerMatricHeadIter = ', mLayerMatricHeadIter(1:10)
   !print*, 'mLayerMatricHeadNew = ',  mLayerMatricHeadNew(1:10)

   ! compute volumetric liquid water content for the soil layers
   !do iLayer=1,nSoil
   ! mLayerVolFracLiqNew(iLayer+nSnow) = volFracLiq(mLayerMatricHeadNew(iLayer),vGn_alpha,theta_res,theta_sat,vGn_n,vGn_m)
   !end do
   !mLayerVolFracLiqNew(nSnow+1:nLayers) = mLayerVolFracLiqIter(nSnow+1:nLayers) + mLayerdTheta_dPsi*(mLayerMatricHeadNew - mLayerMatricHeadIter)

   ! compute the relative change in matric head
   mLayerMatricIncrement = abs(mLayerMatricHeadNew - mLayerMatricHeadIter)

   if(printflag)then
    print*, 'mLayerVolFracLiqNew, before phase change = '
    do iLayer=51,70
     write(*,'(i4,1x,6(f20.5))'), iLayer, mLayerVolFracLiqIter(iLayer), mLayerVolFracLiqNew(iLayer), mLayerVolFracLiqNew(iLayer) - mLayerVolFracLiqIter(iLayer), &
                                    mLayerMatricHeadNew(iLayer-nSnow),mLayerMatricHeadIter(iLayer-nSnow), mLayerMatricHeadNew(iLayer-nSnow) - mLayerMatricHeadIter(iLayer-nSnow)
    end do
   endif

   ! update mass vectors
   mLayerMatricHeadIter = mLayerMatricHeadNew
   mLayerVolFracLiqIter = mLayerVolFracLiqNew

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

   ! compute the fraction of liquid water in the top layer
   scalarFracLiqTop = mLayerVolFracLiqIter(1)/(mLayerVolFracLiqIter(1)+mLayerVolFracIceIter(1))
  
   ! compute the change in temperature over the iteration
   call tempchange(dt_sub,                 & ! time step (seconds)
                   mLayerTempIter,         & ! temperature at the current iteration (K)
                   mLayerVolFracIceIter,   & ! volumetric ice content at the current iteration (-)
                   mLayerTempNew,          & ! temperature at the new iteration (K)
                   err,cmessage)
   if(err/=0)then; err=10; message=trim(message)//trim(cmessage); return; endif

   ! add/subtract evaporation/dew [liquid water] and sublimation/frost [ice] from the top layer (kg m-2 s-1 --> volume over time step)
   mLayerVolFracIceIter(1) = mLayerVolFracIceIter(1) + ((scalarMassSolidNew - scalarMassSolidIter)/iden_ice) * dt/mLayerDepth(1)
   mLayerVolFracLiqIter(1) = mLayerVolFracLiqIter(1) + ((scalarMassLiquidNew - scalarMassLiquidIter)/iden_water) * dt/mLayerDepth(1)

   ! re-set iterations
   scalarMassSolidIter  = scalarMassSolidNew
   scalarMassLiquidIter = scalarMassLiquidNew


   if(printflag)then
    print*, 'temp_new = '
    do iLayer=50,70
     print*, iLayer, mLayerTempNew(iLayer)
    end do
   endif

   ! print temperatures
   !print*, 'mLayerTempNew(1:5) = ', mLayerTempNew(1:5)

   ! compute new liquid water and ice
   do iLayer=1,nLayers
    ! check ice content exists (always true for snow)
    if(mLayerVolFracIceIter(iLayer)>0._dp)then
     ! compute change in volumetric ice content and volumetric liquid water content
     del_ice = -mLayerdTheta_dTk(iLayer)*(iden_water/iden_ice)*(mLayerTempNew(iLayer) - mLayerTempIter(iLayer))
     del_h2o = -del_ice*(iden_ice/iden_water)
     ! ensure that the change in ice content does not exceed available water or ice
     if(-del_ice>mLayerVolFracIceIter(iLayer)) del_ice = -mLayerVolFracIceIter(iLayer)
     if(-del_h2o>mLayerVolFracLiqIter(iLayer)) del_h2o = -mLayerVolFracLiqIter(iLayer)
     ! compute the new ice content and liquid water content
     mLayerVolFracIceNew(iLayer) = mLayerVolFracIceIter(iLayer) + del_ice
     mLayerVolFracLiqNew(iLayer) = mLayerVolFracLiqIter(iLayer) + del_h2o
     if(mLayerVolFracLiqNew(iLayer)>1._dp) then
      print*, mLayerTempNew(iLayer) - mLayerTempIter(iLayer), mLayerdTheta_dTk
      pause ' excessive water'
     endif
    endif
    ! process snow and soil
    select case(layerType(iLayer))
     ! **************************************************
     ! ***** process snow *******************************
     ! **************************************************
     case(ix_snow)
      ! check snow
      if(mLayerVolFracIceNew(iLayer) < 0.1_dp) stop ' small volumetric ice content'
      ! if temperature is above freezing, then estimate temperature based on the freezing curve
      if(mLayerTempNew(iLayer)>Tfreeze)then
       ! case where no liquid water exists
       if(mLayerVolFracLiqNew(iLayer)<epsilon(mLayerVolFracLiqNew))then
        mLayerTempNew(iLayer) = 0.5_dp*(Tfreeze + mLayerTempIter(iLayer))
        del_h2o = fracliquid(mLayerTempNew(iLayer),snowfrz_scale)*mLayerVolFracIceNew(iLayer) - mLayerVolFracLiqNew(iLayer)
        del_ice = -del_h2o*(iden_water/iden_ice)
        mLayerVolFracIceNew(iLayer) = mLayerVolFracIceIter(iLayer) + del_ice
        mLayerVolFracLiqNew(iLayer) = mLayerVolFracLiqIter(iLayer) + del_h2o
        printflag=.true.
        pause 'test correction'
       ! case where liquid water exists
       else
        print*, mLayerTempNew(iLayer), mLayerVolFracLiqNew(iLayer), mLayerVolFracIceNew(iLayer)
        mLayerTempNew(iLayer) = templiquid(mLayerVolFracLiqNew(iLayer)/mLayerVolFracIceNew(iLayer),snowfrz_scale)
        print*, mLayerTempNew(iLayer)
        pause ' test templiquid '
       endif  ! (if no liquid water exists)
      endif  ! (if snow temperature is greater than Tfreeze)
      ! test for active flow
      if(mLayerVolFracLiqNew(iLayer)>0.1_dp) then
       printflag=.true.
       pause 'active flow'
      endif
     ! **************************************************
     ! ***** process soil *******************************
     ! **************************************************
     case(ix_soil)
      ! ***** check if temperatures drop below the "freezing point" (NOTE: Tcrit<Tfreeze)
      if(mLayerVolFracIceIter(iLayer)<epsilon(mLayerVolFracIceIter))then  ! (no ice content at start of iteration)
       ! compute the "liquid equivalent" volumetric fraction of total water (liquid plus ice)
       theta = mLayerVolFracLiqIter(iLayer) + mLayerVolFracIceIter(iLayer)*(iden_ice/iden_water)
       ! compute the critical soil temperature above which all water is unfrozen (K)
       Tcrit = crit_soilT(theta,theta_res,theta_sat,vGn_alpha,vGn_n,vGn_m)
       ! process cases where temperatures drops below the freezing point
       if(mLayerTempIter(iLayer)>Tcrit .and. mLayerTempNew(iLayer)<=Tcrit)then
        mLayerTempNew(iLayer)       = Tcrit
        mLayerVolFracIceNew(iLayer) = eps
        mLayerVolFracLiqNew(iLayer) = mLayerVolFracLiqIter(iLayer) - eps*(iden_ice/iden_water)
       else
        mLayerVolFracIceNew(iLayer) = mLayerVolFracIceIter(iLayer)
        mLayerVolFracLiqNew(iLayer) = mLayerVolFracLiqIter(iLayer)
       endif
      endif  ! (no ice content at start of iteration)
     ! check case is found
     case default; err=35; message="f-fuse/def_output/varTypeNotFound"; return
    endselect
   end do ! (looping through layers)

   ! compute melt/freeze in each layer (kg m-3)
   mLayerMeltFreeze = mLayerMeltFreeze + iden_ice*(mLayerVolFracIceNew - mLayerVolFracIceIter)

   if(printflag)then
    print*, 'mLayerVolFracLiqNew, after phase change = '
    do iLayer=50,70
     write(*,'(i4,1x,3(f20.10))'), iLayer, mLayerVolFracLiqIter(iLayer), mLayerVolFracLiqNew(iLayer), mLayerVolFracLiqNew(iLayer) - mLayerVolFracLiqIter(iLayer)
    end do
   endif

   ! compute residual error
   call energy_err(dt_sub,mLayerTempNew,mLayerVolFracIceNew,mLayerNrgError,err,cmessage)
   if(err/=0)then; err=10; message=trim(message)//trim(cmessage); return; endif
 
   if(printflag)then
    do iLayer=50,70
     print*, 'mLayerNrgError = ', iLayer, mLayerNrgError(iLayer)
    end do
   endif

   ! compute the increment in energy (J m-3)
   mLayerNrgIncrement = mLayerVolHtCapBulk*(mLayerTempNew - mLayerTempIter) - LH_fus*iden_ice*(mLayerVolFracIceNew - mLayerVolFracIceIter)
   if(printflag)then
    print*, 'mLayerNrgIncrement = '
    do iLayer=50,70
     write(*,'(i4,1x,3(e20.10,1x))') iLayer, mLayerNrgIncrement(iLayer), mLayerVolHtCapBulk(iLayer)*(mLayerTempNew(iLayer) - mLayerTempIter(iLayer)), &
               LH_fus*iden_ice*(mLayerVolFracIceNew(iLayer) - mLayerVolFracIceIter(iLayer))
    end do
   endif

   ! update temperature vector
   mLayerTempIter = mLayerTempNew

   ! update volumetric liquid and ice content
   mLayerVolFracLiqIter = mLayerVolFracLiqNew
   mLayerVolFracIceIter = mLayerVolFracIceNew

   ! compute maximum derivatives
   nrg_max    = maxval(abs(mLayerNrgIncrement))
   matric_max = maxval(abs(mLayerMatricIncrement) - matric_tol*abs(mLayerMatricHeadNew))

   !print*, 'abs(mLayerMatricIncrement)', abs(mLayerMatricIncrement)
   !print*, 'matric_tol*abs(mLayerMatricHeadNew)', matric_tol*abs(mLayerMatricHeadNew)
   !print*, 'abs(mLayerMatricIncrement) - matric_tol*abs(mLayerMatricHeadNew)', abs(mLayerMatricIncrement) - matric_tol*abs(mLayerMatricHeadNew)
   !pause

   print*, 'nrg_max, matric_max', nrg_max, matric_max , maxloc(abs(mLayerNrgIncrement))
   pause
   !if(printflag) pause

   ! update mass vectors
   mLayerMatricHeadIter = mLayerMatricHeadNew
   mLayerVolFracLiqIter = mLayerVolFracLiqNew

   ! check for convergence
   if (nrg_max(1) < nrg_tol .and. matric_max(1) < 0._dp) exit

   ! check for lack of convergence
   if(niter==maxiter)then; err=20; message=trim(message)//'failed to converge'; return; endif

  end do  ! (iterating)

  ! update the state vectors
  mLayerTemp       = mLayerTempNew
  mLayerVolFracIce = mLayerVolFracIceNew
  mLayerMatricHead = mLayerMatricHeadNew
  mLayerVolFracLiq = mLayerVolFracLiqNew

  !print*, 'after  iteration loop'
  !print*, 'mLayerTempIter(1:5) = ', mLayerTempIter(1:5)


  ! test matric head
  !print*,'mLayerMatricHeadNew  = ', mLayerMatricHeadNew

  ! increment the time step increment
  dt_done = dt_done + dt_sub
  print*, dt_done, dt_sub, niter

  ! modify the length of the time step
  if(niter<n_inc) dt_sub = dt_sub*F_inc
  if(niter>n_dec) dt_sub = dt_sub*F_dec

  ! save the time step to initialize the subsequent step
  if(dt_done<dt .or. nsub==1) dt_init = dt_sub

  ! exit do-loop if finished
  if(dt_done>=dt)exit

  ! make sure that we don't exceed the step
  dt_sub = min(dt-dt_done, dt_sub)
  !pause

 end do  ! (sub-step loop)

 !print*, 'nsub = ', nsub
 !pause

 

 ! deallocate space for state variables at the current iteration
 deallocate(mLayerTempIter,mLayerVolFracIceIter,mLayerMatricHeadIter,mLayerVolFracLiqIter,stat=err)
 if(err/=0)then; err=20; message='problem deallocating space for state variable vectors - 1'; return; endif

 ! deallocate space for state variables at the current iteration
 deallocate(mLayerTempNew,mLayerVolFracIceNew,mLayerMatricHeadNew,mLayerVolFracLiqNew,stat=err)
 if(err/=0)then; err=20; message='problem deallocating space for state variable vectors - 2'; return; endif

 ! deallocate space for diagnostic variables
 deallocate(mLayerInfilFreezeNew,stat=err)
 if(err/=0)then; err=20; message='problem deallocating space for diagnostic vectors'; return; endif

 ! deallocate space for error monitoring
 deallocate(mLayerNrgIncrement,mLayerNrgError,mLayerMatricIncrement,stat=err)
 if(err/=0)then; err=20; message='problem deallocating space for error monitoring vectors'; return; endif

 ! estimate the mass balance error
 !print*, 'mass = ', sum((mLayerVolFracLiqNew*iden_water + mLayerVolFracIceNew*iden_ice)*mLayerDepth), &
 !                   sum((mLayerVolFracLiq*iden_water + mLayerVolFracIce*iden_ice)*mLayerDepth)


 end subroutine coupled_em

end module coupled_em_module

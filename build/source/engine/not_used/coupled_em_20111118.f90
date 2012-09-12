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
 USE soil_utils_module,only:dTheta_dTk      ! differentiate the freezing curve w.r.t. temperature
 USE soil_utils_module,only:crit_soilT      ! compute the critical temperature above which all water is unfrozen
 USE energyflux_module,only:tempchange      ! compute change in temperature over the time step
 USE energyflux_module,only:energy_err      ! error in energy equation (J m-3)
 USE liquidflux_module,only:masschange      ! compute change in mass over the time step
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
 real(dp),intent(inout)               :: dt_init                 ! used to initialize the size of the sub-step
 integer(i4b),intent(out)             :: err                     ! error code
 character(*),intent(out)             :: message                 ! error message
 ! local pointers to model parameters
 real(dp),pointer                     :: vGn_alpha               ! van Genutchen "alpha" parameter
 real(dp),pointer                     :: vGn_n                   ! van Genutchen "n" parameter
 real(dp),pointer                     :: theta_sat               ! soil porosity (-)
 real(dp),pointer                     :: theta_res               ! soil residual volumetric water content (-)
 real(dp),pointer                     :: k_soil                  ! hydraulic conductivity (m s-1)
 real(dp),pointer                     :: spec_storage            ! specific storage coefficient (m-1)
 real(dp),pointer                     :: f_impede                ! ice impedence factor (-)
 ! local pointers to derived model variables that are constant over the simulation period
 real(dp),pointer                     :: vGn_m                   ! van Genutchen "m" parameter (-)
 real(dp),pointer                     :: volHtCap_air            ! volumetric heat capacity of air (J m-3 K-1)
 real(dp),pointer                     :: volHtCap_ice            ! volumetric heat capacity of ice (J m-3 K-1)
 real(dp),pointer                     :: volHtCap_soil           ! volumetric heat capacity of dry soil (J m-3 K-1)
 real(dp),pointer                     :: volHtCap_water          ! volumetric heat capacity of liquid water (J m-3 K-1)
 real(dp),pointer                     :: volLatHt_fus            ! volumetric latent heat of fusion (J m-3 K-1)
 ! local pointers to model variables that are constant over the time step
 real(dp),pointer                     :: mLayerDepth(:)          ! depth of the layer (m)
 real(dp),pointer                     :: mLayerHeight(:)         ! height of the layer mid-point (m)
 real(dp),pointer                     :: iLayerHeight(:)         ! height of the layer interface (m)
 ! local pointers to model state variables
 real(dp),pointer                     :: mLayerTemp(:)           ! temperature of each layer (K)
 real(dp),pointer                     :: mLayerVolFracIce(:)     ! volumetric fraction of ice in each layer (-)
 real(dp),pointer                     :: mLayerVolFracLiq(:)     ! volumetric fraction of liquid water in each layer (-)
 real(dp),pointer                     :: mLayerMatricHead(:)     ! matric head in each layer (m)
 ! local pointers to model diagnostic variables
 real(dp),pointer                     :: iLayerThermalC(:)       ! thermal conductivity at the interface of each layer (W m-1 K-1)
 real(dp),pointer                     :: mLayerRadCondFlux(:)    ! change in layer energy (J m-3 s-1)
 real(dp),pointer                     :: mLayerVolHtCapBulk(:)   ! bulk volumetric heat capacity (J m-3 K-1)
 real(dp),pointer                     :: mLayerdTheta_dPsi(:)    ! derivative in the soil water characteristic (m-1)
 real(dp),pointer                     :: mLayerdTheta_dTk(:)     ! derivative in the freezing curve (K-1)


 ! local pointers to model index variables
 integer(i4b),pointer                 :: nLayers                 ! number of layers
 integer(i4b),pointer                 :: layerType(:)            ! type of the layer (ix_soil or ix_snow)



 ! define local model state variables
 real(dp),allocatable                 :: mLayerTempIter(:)       ! temperature vector at the current iteration (K)
 real(dp),allocatable                 :: mLayerTempNew(:)        ! temperature vector at the next iteration (K)
 real(dp),allocatable                 :: mLayerVolFracIceIter(:) ! volumetric fraction of ice at the current iteration (-)
 real(dp),allocatable                 :: mLayerVolFracIceNew(:)  ! volumetric fraction of ice at the next iteration (-)
 real(dp),allocatable                 :: mLayerVolFracLiqIter(:) ! volumetric fraction of liquid water at the current iteration (-)
 real(dp),allocatable                 :: mLayerVolFracLiqNew(:)  ! volumetric fraction of liquid water at the next iteration (-)
 real(dp),allocatable                 :: mLayerMatricHeadIter(:) ! matric head at the current iteration (-)
 real(dp),allocatable                 :: mLayerMatricHeadNew(:)  ! matric head at the next iteration (-)
 



 real(dp),allocatable                 :: mLayerDelTemp(:)     ! change in temperatre over the time step (K)
 real(dp),allocatable                 :: mLayerDelIce(:)      ! change in volumetric ice content over the time step (-)
 real(dp),allocatable                 :: nrg_diff(:)          ! change in energy from one iteration to the next (J m-3)
 real(dp),allocatable                 :: matric_diff(:)       ! relative change in matric head from one iteration to the next (m)

 real(dp),allocatable                 :: mLayerNrgMelt(:)     ! energy available for melting (J m-3)

 logical(lgt),allocatable             :: mLayerCrossFlag(:)   ! flag when temperature crosses zero 
 real(dp)                             :: dIce_dT              ! derivative in volumetric ice content w.r.t. temperature (K-1)
 real(dp)                             :: del_ice              ! change in volumetric ice content over the iteration (-)
 real(dp)                             :: kappa                ! constant in the freezing curve function (m K-1)
 real(dp)                             :: liqMx                ! maximum possible volumetric fraction of liquid water at a given temperature (-)
 real(dp),allocatable                 :: theta(:)             ! volumetric fraction of total water, liquid plus ice (-)
 real(dp),allocatable                 :: Tcrit(:)             ! critical soil temperature above which all water is unfrozen (K)
 real(dp),allocatable                 :: mLayerNrgError(:)    ! energy error in each layer (J m-3)

 ! length of the time step
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
 integer(i4b)                         :: nsub                 ! number of sub-steps
 integer(i4b)                         :: iter                 ! iteration index
 integer(i4b)                         :: niter                ! number of iterations
 integer(i4b),parameter               :: maxiter=20           ! maximum number of iterations
 real(dp),dimension(1)                :: nrg_max,matric_max   ! maximum change in energy/matric head for a given iteration
 real(dp),parameter                   :: nrg_tol=1.d-0        ! iteration tolerance for energy (J m-3) -- convergence when nrg_max < nrg_tol
 real(dp),parameter                   :: matric_tol=1.d-2     ! iteration tolerance for matric head -- convergence when matric_max < matric_tol
 integer(i4b)                         :: iLayer               ! loop through model layers

 ! initialize error control
 err=0; message="coupled_em/"

 ! assign pointers to model parameters
 vGn_alpha         => mpar_data%var(iLookPARAM%vGn_alpha)      ! van Genutchen "alpha" parameter (m-1)
 vGn_n             => mpar_data%var(iLookPARAM%vGn_n)          ! van Genutchen "n" parameter (-)
 theta_sat         => mpar_data%var(iLookPARAM%theta_sat)      ! soil porosity (-)
 theta_res         => mpar_data%var(iLookPARAM%theta_res)      ! soil residual volumetric water content (-)
 k_soil            => mpar_data%var(iLookPARAM%k_soil)         ! hydraulic conductivity (m s-1)
 spec_storage      => mpar_data%var(iLookPARAM%spec_storage)   ! specific storage coefficient (m-1)
 f_impede          => mpar_data%var(iLookPARAM%f_impede)       ! ice impedence factor (-)

 ! assign pointers to model variables that are constant over the simulation period
 vGn_m             => mvar_data%var(iLookMVAR%scalarVGn_m)%dat(1)          ! van Genutchen "m" parameter (-)
 volHtCap_air      => mvar_data%var(iLookMVAR%scalarVolHtCap_air)%dat(1)   ! volumetric heat capacity of air (J m-3 K-1)
 volHtCap_ice      => mvar_data%var(iLookMVAR%scalarVolHtCap_ice)%dat(1)   ! volumetric heat capacity of ice (J m-3 K-1)
 volHtCap_soil     => mvar_data%var(iLookMVAR%scalarVolHtCap_soil)%dat(1)  ! volumetric heat capacity of soil (J m-3 K-1)
 volHtCap_water    => mvar_data%var(iLookMVAR%scalarVolHtCap_water)%dat(1) ! volumetric heat capacity of water (J m-3 K-1)
 volLatHt_fus      => mvar_data%var(iLookMVAR%scalarvolLatHt_fus)%dat(1)   ! volumetric latent heat of fusion (J m-3)

 ! assign pointers to model variables that are constant over the time step
 mLayerDepth       => mvar_data%var(iLookMVAR%mLayerDepth)%dat             ! depth of the layer (m)
 mLayerHeight      => mvar_data%var(iLookMVAR%mLayerHeight)%dat            ! height of the layer mid-point (m)
 iLayerHeight      => mvar_data%var(iLookMVAR%iLayerHeight)%dat            ! height of the layer interface (m)

 ! assign pointers to model state variables
 mLayerTemp        => mvar_data%var(iLookMVAR%mLayerTemp)%dat              ! temperature of each layer (K)
 mLayerVolFracIce  => mvar_data%var(iLookMVAR%mLayerVolFracIce)%dat        ! volumetric fraction of ice in each layer (-)
 mLayerVolFracLiq  => mvar_data%var(iLookMVAR%mLayerVolFracLiq)%dat        ! volumetric fraction of liquid water in each layer (-)
 mLayerMatricHead  => mvar_data%var(iLookMVAR%mLayerMatricHead)%dat        ! matric head in each layer (m)

 ! assign pointers to model diagnostic variables
 iLayerThermalC    => mvar_data%var(iLookMVAR%iLayerThermalC)%dat          ! thermal conductivity at the interface of each layer (W m-1 K-1)
 mLayerVolHtCapBulk=> mvar_data%var(iLookMVAR%mLayerVolHtCapBulk)%dat      ! bulk volumetric heat capacity (J m-3 K-1)
 mLayerRadCondFlux => mvar_data%var(iLookMVAR%mLayerRadCondFlux)%dat       ! change in layer energy (J m-3 s-1)
 mLayerdTheta_dPsi => mvar_data%var(iLookMVAR%mLayerdTheta_dPsi)%dat       ! derivative in the soil water characteristic (m-1)
 mLayerdTheta_dTk  => mvar_data%var(iLookMVAR%mLayerdTheta_dTk)%dat        ! derivative in the freezing curve (K-1)

 ! assign local pointers to the model index structures
 nLayers           => indx_data%var(iLookINDEX%nLayers)%dat(1)             ! number of layers
 layerType         => indx_data%var(iLookINDEX%layerType)%dat              ! layer type (ix_soil or ix_snow)

 !print*, 'geometry'
 !do iLayer=1,nLayers
 ! write(*,'(i4,1x,10(f12.6,1x))'), iLayer, iLayerHeight(iLayer-1:iLayer), mLayerHeight(iLayer), mLayerDepth(iLayer)
 !end do

 ! get the length of the time step (seconds)
 dt = forcFileInfo%data_step
 !print*, 'dt = ', dt

 ! define the constant in the freezing curve function (m K-1)
 kappa  = (iden_ice/iden_water)*(LH_fus/(gravity*Tfreeze))  ! NOTE: J = kg m2 s-2

 ! allocate space for state variables at the current iteration
 allocate(mLayerTempIter(nLayers),mLayerVolFracIceIter(nLayers),mLayerMatricHeadIter(nLayers),&
          mLayerVolFracLiqIter(nLayers),stat=err)
 if(err/=0)then; err=20; message='problem allocating space for state variable vectors - 1'; return; endif

 ! allocate space for state variables at the current iteration
 allocate(mLayerTempNew(nLayers),mLayerVolFracIceNew(nLayers),mLayerMatricHeadNew(nLayers),&
          mLayerVolFracLiqNew(nLayers),stat=err)
 if(err/=0)then; err=20; message='problem allocating space for state variable vectors - 2'; return; endif

 ! allocate space for other crap
 allocate(mLayerDelTemp(nLayers),mLayerCrossFlag(nLayers),&
          mLayerNrgMelt(nLayers),mLayerDelIce(nLayers),&
          nrg_diff(nLayers),&
          matric_diff(nLayers),&
          Tcrit(nLayers),&
          theta(nLayers),&
          mLayerNrgError(nLayers),&
          stat=err)
 if(err/=0)then; err=20; message='problem allocating space for crap vectors'; return; endif

 ! initialize state variables at the next iteration
 mLayerTempIter       = mLayerTemp
 mLayerVolFracIceIter = mLayerVolFracIce
 mLayerMatricHeadIter = mLayerMatricHead
 mLayerVolFracLiqIter = mLayerVolFracLiq
 
 !print*, 'initial mass = ', sum(mLayerVolFracLiq*iden_water*mLayerDepth)
 !pause

 ! compute the derivative in the soil water characteristic (m-1)
 do iLayer=1,nLayers
  mLayerdTheta_dPsi(iLayer) = dTheta_dPsi(mLayerMatricHeadIter(iLayer),vGn_alpha,theta_res,theta_sat,vGn_n,vGn_m)
 end do

 ! initialize the length of the sub-step
 dt_sub  = min(dt_init,dt)
 dt_done = 0._dp

 ! initialize the number of sub-steps
 nsub=0

 ! loop through sub-steps
 do  ! continuous do statement with exit clause (alternative to "while")

  ! increment the number of sub-steps
  nsub = nsub+1

  ! initialize number of iterations
  niter=0

  print*, '*********************************************************'
  print*, '****************** start iterating **********************'
  print*, '*********************************************************'

  ! iterate
  do iter=1,maxiter

   ! increment number of iterations
   niter=niter+1
   if(niter>5) pause 'number of iterations is greater than 5'
   print*, '***** new iteration *****', niter
   print*, 'mLayerVolFracLiqIter = ', mLayerVolFracLiqIter
   print*, 'mLayerVolFracIceIter = ', mLayerVolFracIceIter


   ! ***** compute the derivative in the freezing curve w.r.t. temperature(K-1)
   do iLayer=1,nLayers
    if(mLayerVolFracIceIter(iLayer)>0._dp)then
     mLayerdTheta_dTk(iLayer) = dTheta_dTk(mLayerTempIter(iLayer),theta_res,theta_sat,vGn_alpha,vGn_n,vGn_m)
    else
     mLayerdTheta_dTk(iLayer) = 0._dp
    endif
   end do
   print*,'mLayerdTheta_dTk = ', mLayerdTheta_dTk
  
   ! compute the change in temperature over the iteration
   call tempchange(dt_sub,mLayerTempIter,mLayerVolFracIceIter,mLayerTempNew,err,cmessage)
   if(err/=0)then; err=10; message=trim(message)//trim(cmessage); return; endif
   ! compute residual error
   call energy_err(dt_sub,mLayerTempNew,mLayerVolFracIceNew,mLayerNrgError,err,cmessage)
   if(err/=0)then; err=10; message=trim(message)//trim(cmessage); return; endif
   print*, 'mLayerNrgError = ', mLayerNrgError

   ! process cases where temperature crosses the critical temperature
   do iLayer=1,nLayers
    ! compute the "liquid equivalent" volumetric fraction of total water (liquid plus ice)
    theta(iLayer) = mLayerVolFracLiqIter(iLayer) + mLayerVolFracIceIter(iLayer)*(iden_ice/iden_water)
    ! compute the critical soil temperature above which all water is unfrozen (K)
    Tcrit(iLayer) = crit_soilT(theta(iLayer),theta_res,theta_sat,vGn_alpha,vGn_n,vGn_m)
    ! process cases where temperatures drops below the freezing point
    if(mLayerTempIter(iLayer)>Tcrit(iLayer) .and. mLayerTempNew(iLayer)<=Tcrit(iLayer))then
     mLayerTempNew(iLayer)=Tcrit(iLayer) - eps
     print*, 'epsilon = ', eps
     write(*,'(a16,1x,f30.25)') 'Tcrit         = ', Tcrit(iLayer)
     write(*,'(a16,1x,f30.25)') 'mLayerTempNew = ', mLayerTempNew(iLayer)
    endif
   end do ! (looping through layers)  
   print*, 'mLayerTempIter = ', mLayerTempIter
   print*, 'mLayerTempNew = ', mLayerTempNew

   ! compute new liquid water and ice
   do iLayer=1,nLayers
    if(mLayerTempNew(iLayer)<Tfreeze)then
     write(*,'(a16,1x,f30.25)') 'Tcrit         = ', Tcrit(iLayer)
     write(*,'(a16,1x,f30.25)') 'mLayerTempNew = ', mLayerTempNew(iLayer)
     LiqMx = theta_res + (theta_sat - theta_res) * &
                           (1._dp + (vGn_alpha*kappa*(Tcrit(iLayer) - Tfreeze))**vGn_n)**(-vGn_m)
     print*, 'LiqMx | Tcrit = ', LiqMx
     ! compute the maximum volumetric liquid water content possible at a given temperature
     LiqMx = theta_res + (theta_sat - theta_res) * &
                           (1._dp + (vGn_alpha*kappa*(mLayerTempNew(iLayer) - Tfreeze))**vGn_n)**(-vGn_m)
     print*, 'LiqMx | temp  = ', LiqMx

     ! compute the actual liquid water content
     mLayerVolFracLiqNew(iLayer) = min(theta(iLayer),LiqMx)
     print*, iLayer, theta(iLayer),LiqMx
     ! compute the ice content
     mLayerVolFracIceNew(iLayer) = (theta(iLayer) - mLayerVolFracLiqNew(iLayer))!*(iden_water/iden_ice)
     ! see what we get by using the derivative
     del_ice = -mLayerdTheta_dTk(iLayer)*(iden_water/iden_ice)*(mLayerTempNew(iLayer) - mLayerTempIter(iLayer))
     print*, 'del_ice = ', del_ice

    else
     mLayerVolFracLiqNew(iLayer) = mLayerVolFracLiqIter(iLayer)
     mLayerVolFracIceNew(iLayer) = mLayerVolFracIceIter(iLayer)
    endif
   end do
   print*, 'mLayerVolFracLiqNew = ', mLayerVolFracLiqNew
   print*, 'mLayerVolFracIceNew = ', mLayerVolFracIceNew
   print*, 'change in ice content = ', mLayerVolFracIceNew - mLayerVolFracIceIter

   ! compute the increment in energy (J m-3)
   nrg_diff = mLayerVolHtCapBulk*(mLayerTempNew - mLayerTempIter) - LH_fus*iden_ice*(mLayerVolFracIceNew - mLayerVolFracIceIter)
   print*, 'nrg_diff (temp) = ', mLayerVolHtCapBulk*(mLayerTempNew - mLayerTempIter)
   print*, 'nrg_diff (ice) = ', LH_fus*iden_ice*(mLayerVolFracIceNew - mLayerVolFracIceIter)
   print*, 'nrg_diff (all) = ', mLayerVolHtCapBulk*(mLayerTempNew - mLayerTempIter) - LH_fus*iden_ice*(mLayerVolFracIceNew - mLayerVolFracIceIter)

   ! update temperature vector
   mLayerTempIter = mLayerTempNew

   !print*,'mLayerTempIter = '
   !print*,mLayerTempIter
   !print*,'mLayerDelTemp = '
   !print*,mLayerDelTemp

   !print*, 'mLayerDelTemp*mLayerVolHtCapBulk = ', mLayerDelTemp*mLayerVolHtCapBulk
   !print*, 'LH_fus*iden_ice*(mLayerVolFracIceNew - mLayerVolFracIceIter) = ', &
   !         LH_fus*iden_ice*(mLayerVolFracIceNew - mLayerVolFracIceIter)
   !print*, 'nrg_diff = ', nrg_diff

   ! update volumetric liquid and ice content
   mLayerVolFracLiqIter = mLayerVolFracLiqNew
   mLayerVolFracIceIter = mLayerVolFracIceNew

   ! compute the derivative in the soil water characteristic (m-1)
   do iLayer=1,nLayers
    mLayerdTheta_dPsi(iLayer) = dTheta_dPsi(mLayerMatricHeadIter(iLayer),vGn_alpha,theta_res,theta_sat,vGn_n,vGn_m)
   end do

   ! compute the matric head at the next iteration
   call masschange(dt_sub,&               ! time step (seconds)
                   mLayerMatricHeadIter,& ! matric head in each layer at the current iteration (m)
                   mLayerVolFracIceIter,& ! volumetric fraction of ice at the current iteration (-)
                   mLayerVolFracLiqIter,& ! volumetric fraction of liquid water at the current iteration (-)
                   mLayerMatricHeadNew, & ! matric head in each layer at the next iteration (m)
                   err,cmessage)
   if(err/=0)then; err=10; message=trim(message)//trim(cmessage); return; endif

   ! compute the relative change in matric head
   matric_diff = abs(mLayerMatricHeadNew - mLayerMatricHeadIter) / abs(mLayerMatricHeadNew)

   ! compute maximum derivatives
   nrg_max    = maxval(abs(nrg_diff))
   matric_max = maxval(abs(matric_diff)/abs(mLayerMatricHeadIter))
   !print*, 'dt_sub = ', dt_sub, 'iter = ', iter, 'nrg_max = ', nrg_max, 'matric_max = ', matric_max

   !print*,'mLayerMatricHead     = ', mLayerMatricHead
   !print*,'mLayerMatricHeadIter = ', mLayerMatricHeadIter
   !print*,'mLayerMatricHeadNew  = ', mLayerMatricHeadNew
   !print*,'matric_diff          = ', matric_diff
   
   ! compute volumetric liquid water content
   do iLayer=1,nLayers
    mLayerVolFracLiqNew(iLayer) = volFracLiq(mLayerMatricHeadNew(iLayer),vGn_alpha,theta_res,theta_sat,vGn_n,vGn_m)
   end do
   !mLayerVolFracLiqNew = mLayerVolFracLiqIter + mLayerdTheta_dPsi*(mLayerMatricHeadNew - mLayerMatricHeadIter)

   !print*,'mLayerVolFracLiq    = ',mLayerVolFracLiq
   !print*,'mLayerVolFracLiqNew = ',mLayerVolFracLiqNew
   !print*,'mLayerVolFracLiqDiff= ',mLayerVolFracLiqNew - mLayerVolFracLiqIter
   !print*, '***** mass check *****'
   !print*, sum(mLayerdTheta_dPsi*(mLayerMatricHeadNew - mLayerMatricHeadIter)*iden_water*mLayerDepth)
   !print*, 'mass = ', sum((mLayerVolFracLiqNew*iden_water + mLayerVolFracIceNew*iden_ice)*mLayerDepth) 

   ! update mass vectors
   mLayerMatricHeadIter = mLayerMatricHeadNew
   mLayerVolFracLiqIter = mLayerVolFracLiqNew

   ! compute energy available for melting and freezing
  


   ! identify cases where the temperature crosses zero
   where((mLayerTempIter>Tfreeze .and. mLayerTempNew<Tfreeze) .or. &
         (mLayerTempIter<Tfreeze .and. mLayerTempNew>Tfreeze) )
    ! compute energy available for melting and freezing (J m-3 s-1)
    mLayerNrgMelt  = mLayerVolHtCapBulk*(mLayerTempNew-Tfreeze)/dt
    mLayerCrossFlag=.true.
   elsewhere
    mLayerCrossFlag=.false.
   endwhere
   !print*, mLayerCrossFlag


   ! check for convergence
   if (nrg_max(1) < nrg_tol .and. matric_max(1) < matric_tol) exit

  end do  ! (iterating)


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

  ! update the state vectors
  mLayerTemp       = mLayerTempNew
  mLayerVolFracIce = mLayerVolFracIceNew
  mLayerMatricHead = mLayerMatricHeadNew
  mLayerVolFracLiq = mLayerVolFracLiqNew
 
 end do  ! (sub-step loop)

 print*, 'nsub = ', nsub

 ! deallocate space for state variables at the current iteration
 deallocate(mLayerTempIter,mLayerVolFracIceIter,mLayerMatricHeadIter,mLayerVolFracLiqIter,stat=err)
 if(err/=0)then; err=20; message='problem deallocating space for state variable vectors - 1'; return; endif

 ! deallocate space for state variables at the current iteration
 deallocate(mLayerTempNew,mLayerVolFracIceNew,mLayerMatricHeadNew,mLayerVolFracLiqNew,stat=err)
 if(err/=0)then; err=20; message='problem deallocating space for state variable vectors - 2'; return; endif

 ! deallocate space for other crap
 deallocate(mLayerDelTemp,mLayerCrossFlag,mLayerNrgMelt,mLayerDelIce,&
            nrg_diff,matric_diff,Tcrit,theta,stat=err)
 if(err/=0)then; err=20; message='problem deallocating space for crap vectors'; return; endif

 ! estimate the mass balance error
 !print*, 'mass = ', sum((mLayerVolFracLiqNew*iden_water + mLayerVolFracIceNew*iden_ice)*mLayerDepth), &
 !                   sum((mLayerVolFracLiq*iden_water + mLayerVolFracIce*iden_ice)*mLayerDepth)


 end subroutine coupled_em

end module coupled_em_module

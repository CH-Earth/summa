module energyflux_module
USE nrtype
implicit none
private
public::diagn_evar
public::tempchange
public::phsechange
! local pointers to model parameters
real(dp),pointer              :: mheight                  ! measurement height (m)
real(dp),pointer              :: wimplicit                ! weight assigned to start-of-step fluxes (-)
real(dp),pointer              :: snowfrz_scale            ! scaling parameter for the snow freezing curve (K-1)
real(dp),pointer              :: lowerBoundTemp           ! temperature of the lower boundary (K)
! local pointers to soil parameters
real(dp),pointer              :: vGn_alpha                ! van Genutchen "alpha" parameter
real(dp),pointer              :: vGn_n                    ! van Genutchen "n" parameter
real(dp),pointer              :: theta_sat                ! soil porosity (-)
real(dp),pointer              :: theta_res                ! soil residual volumetric water content (-)
! local pointers to vegetation parameters
real(dp),pointer              :: LAI                      ! leaf area index (m2 m-2)
real(dp),pointer              :: minStomatalResist        ! minimum stomatal resistance (s m-1)
real(dp),pointer              :: plantWiltPsi             ! critical matric head when stomatal resitance 2 x min (m)
real(dp),pointer              :: plantWiltExp             ! empirical exponent in plant wilting factor expression (-)
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
real(dp),pointer              :: mLayerMatricHead(:)      ! matric head in each layer (m)
real(dp),pointer              :: mLayerVolFracIce(:)      ! volumetric fraction of ice in each layer (-)
real(dp),pointer              :: mLayerDepth(:)           ! depth of each layer (m)
! local pointers to model diagnostic variables
real(dp),pointer              :: mLayerHeight(:)          ! height at the mid-point of each layer (m)
real(dp),pointer              :: iLayerHeight(:)          ! height at the interface of each layer (m)
real(dp),pointer              :: mLayerVolHtCapBulk(:)    ! bulk volumetric heat capacity (J m-3 K-1)
real(dp),pointer              :: mLayerThermalC(:)        ! thermal conductivity at the mid-point of each layer (W m-1 K-1)
real(dp),pointer              :: iLayerThermalC(:)        ! thermal conductivity at the interface of each layer (W m-1 K-1)
real(dp),pointer              :: mLayerTcrit(:)           ! critical soil temperature above which all water is unfrozen (K)
real(dp),pointer              :: mLayerdTheta_dTk(:)      ! derivative in the freezing curve (K-1)
real(dp),pointer              :: mLayerTranspireLim(:)    ! moisture avail factor limiting transpiration in each layer (-)
real(dp),pointer              :: mLayerInitTranspire(:)   ! transpiration loss from each soil layer at the start of the step (m s-1)
real(dp),pointer              :: mLayerTranspire(:)       ! transpiration loss from each soil layer (m s-1)
real(dp),pointer              :: iLayerInitNrgFlux(:)     ! energy flux at layer interfaces at the start of the time step (W m-2)
! local pointers to diagnostic scalar variables
real(dp),pointer              :: scalarTranspireLim       ! aggregate soil moist & veg limit on transpiration, weighted by root density (-)
real(dp),pointer              :: scalarPotentialET        ! potential ET (kg m-2 s-1)
real(dp),pointer              :: scalarMassLiquid         ! evaporation/dew (kg m-2 s-1)
real(dp),pointer              :: scalarMassSolid          ! sublimation/frost (kg m-2 s-1)
real(dp),pointer              :: scalarSenHeat            ! sensible heat flux at the surface (W m-2)
real(dp),pointer              :: scalarLatHeat            ! latent heat flux at the surface (W m-2)
real(dp),pointer              :: scalarExCoef             ! turbulent exchange coefficient (-)
real(dp),pointer              :: scalarExSen              ! exchange factor for sensible heat (J m-2 s-1 K-1)
real(dp),pointer              :: scalarExLat              ! exchange factor for latent heat (J m-2 s-1)
! local pointers to model index variables
integer(i4b),pointer          :: nLayers                  ! number of layers
integer(i4b),pointer          :: layerType(:)             ! type of the layer (ix_soil or ix_snow)
integer(i4b)                  :: nSoil,nSnow              ! number of snow and soil layers
! local parameters
real(dp),parameter            :: RHsurf=1._dp             ! relative humidity of the surface (-)
! look-up values for method used to compute derivative
integer(i4b)                  :: fDerivMeth               ! index for the method used to calculate flux derivatives
integer(i4b),parameter        :: numerical=1001           ! look-up value for numerical solution
integer(i4b),parameter        :: analytical=1002          ! look-up value for analytical solution
contains


 ! ************************************************************************************************
 ! new subroutine: compute change in temperature over the time step
 ! ************************************************************************************************
 subroutine tempchange(dt,                   & ! time step (seconds)
                       iter,                 & ! current iteration count
                       mLayerTempIter,       & ! trial temperature at the current iteration (K)
                       mLayerVolFracIceIter, & ! volumetric fraction of ice at the current iteration (-)
                       mLayerVolFracLiqIter, & ! volumetric fraction of liquid water at the current iteration (-)
                       mLayerMatricHeadIter, & ! matric head at the current iteration (m)
                       mLayerTempDiffOld,    & ! iteration increment for temperature at the last iteration (K)
                       mLayerTempDiff,       & ! iteration increment for temperature (K)
                       mLayerTempNew,        & ! new temperature (K)
                       mLayerVolFracIceNew,  & ! new volumetric fraction of ice (-)
                       mLayerVolFracLiqNew,  & ! new volumetric fraction of liquid water (-)
                       mLayerMatricHeadNew,  & ! new matric head (m)
                       err,message)            ! error control
 ! physical constants
 USE multiconst,only:&
                     LH_fus,&                  ! latent heat of fusion       (J kg-1)
                     LH_vap,&                  ! latent heat of vaporization (J kg-1)
                     LH_sub,&                  ! latent heat of sublimation (J kg-1)
                     Tfreeze,&                 ! freezing point of pure water (K)
                     iden_air,&                ! intrinsic density of air (kg m-3)
                     iden_ice,&                ! intrinsic density of ice (kg m-3)
                     iden_water                ! intrinsic density of water (kg m-3)
 ! model decisions
 USE data_struc,only:model_decisions           ! model decision structure
 USE data_struc,only:ix_soil,ix_snow           ! names variables for snow and soil
 USE var_lookup,only:iLookDECISIONS            ! named variables for elements of the decision structure
 ! utility modules
 USE snow_utils_module,only:dFracLiq_dTk       ! differentiate the freezing curve w.r.t. temperature (snow)
 USE soil_utils_module,only:dTheta_dTk         ! differentiate the freezing curve w.r.t. temperature (soil)
 USE conv_funcs_module,only:relhm2sphm         ! compute specific humidity 
 USE tridagSolv_module,only:tridag             ! solve tridiagonal system of equations
 ! compute change in temperature over the time step
 implicit none
 ! dummy variables
 real(dp),intent(in)           :: dt                       ! time step (seconds)
 integer(i4b),intent(in)       :: iter                     ! iteration count
 real(dp),intent(in)           :: mLayerTempIter(:)        ! trial temperature at the current iteration (K)
 real(dp),intent(in)           :: mLayerVolFracIceIter(:)  ! volumetric fraction of ice at the current iteration (-)
 real(dp),intent(in)           :: mLayerVolFracLiqIter(:)  ! volumetric fraction of liquid water at the current iteration (-)
 real(dp),intent(in)           :: mLayerMatricHeadIter(:)  ! matric head at the current iteration (m)
 real(dp),intent(in)           :: mLayerTempDiffOld(:)     ! iteration increment for temperature at the last iteration (K) 
 real(dp),intent(out)          :: mLayerTempDiff(:)        ! iteration increment for temperature (K) 
 real(dp),intent(out)          :: mLayerTempNew(:)         ! new temperature (K)
 real(dp),intent(out)          :: mLayerVolFracIceNew(:)   ! new volumetric fraction of ice (-)
 real(dp),intent(out)          :: mLayerVolFracLiqNew(:)   ! new volumetric fraction of liquid water (-)
 real(dp),intent(out)          :: mLayerMatricHeadNew(:)   ! new matric head (m)
 integer(i4b),intent(out)      :: err                      ! error code
 character(*),intent(out)      :: message                  ! error message
 ! define general local variables
 character(LEN=256)            :: cmessage                 ! error message of downwind routine
 logical(lgt)                  :: printflag                ! .true. if print progress to the screen
 logical(lgt)                  :: fTranspire               ! .true. if computing transpiration
 real(dp)                      :: theta                    ! total volumetric water content (liquid plus ice)
 ! define the local variables for the solution
 real(dp),dimension(0:size(mLayerTempIter)) :: dFlux_dTempAbove ! derivative in flux w.r.t. temperature in the layer above (J m-2 s-1 K-1)
 real(dp),dimension(0:size(mLayerTempIter)) :: dFlux_dTempBelow ! derivative in flux w.r.t. temperature in the layer below (J m-2 s-1 K-1)
 real(dp)                                   :: wtim        ! weighted time (s-1)
 real(dp),dimension(size(mLayerTempIter)-1) :: d_m1        ! sub-diagonal elements of the tridiagonal system (J m-3 K-1)
 real(dp),dimension(size(mLayerTempIter))   :: diag        ! diagonal elements of the tridiagonal system (J m-3 K-1)
 real(dp),dimension(size(mLayerTempIter)-1) :: d_p1        ! super-diagonal elements of the tridiagonal system (J m-3 K-1)
 real(dp),dimension(size(mLayerTempIter))   :: rvec        ! residual vector (J m-3)
 ! define the local variables for the line search
 real(dp),dimension(size(mLayerTempIter))   :: g           ! gradient of the function vector (J m-3 J m-3 K-1)
 real(dp)                                   :: fold,fnew   ! function values (J m-3 J m-3)
 real(dp),parameter            :: STPMX=5._dp              ! maximum step size in line search (K)
 real(dp)                      :: stpmax                   ! scaled maximum step size
 integer(i4b)                  :: iLayer                   ! index of model layers
 real(dp)                      :: critDiff                 ! temperature difference from critical temperature (K)
 real(dp),parameter            :: epsT=1.d-10              ! offset from Tcrit when re-setting iterations at the critical temperature (K)
 logical(lgt)                  :: crossFlag                ! .true. if temperature crosses the critical temperature
 ! define the local variables for the Jacobian
 real(dp),dimension(size(mLayerTempIter),size(mLayerTempIter)) :: jmat  ! jacobian matrix
 real(dp),dimension(size(mLayerTempIter))                      :: ftest,xsav,xtry ! compute jacobian: function vector and dependent variables
 real(dp),dimension(size(mLayerTempIter))                      :: xph,h           ! compute jacobian: perturbed vector and finite difference increment
 real(dp),parameter                                            :: eps=-1.0e-8_dp  ! compute jacobian: finite difference increment
 integer(i4b)                                                  :: ijac            ! cpmpute jacobian: index of columns
 logical(lgt)                                                  :: test_jac        ! .true. if desire to test the Jacobian
 ! define the numerical method
 integer(i4b)                         :: num_method               ! numerical method
 integer(i4b),parameter               :: itertive=2001            ! named index for the iterative method
 integer(i4b),parameter               :: non_iter=2002            ! named index for the non-iterative methof
 integer(i4b),parameter               :: itersurf=2003            ! named index for the case where iterate only on the surface energy balance

 ! initialize error control
 err=0; message="tempchange/"

 ! initialize print flag
 printflag=.false.

 ! initialize the Jacobian test
 test_jac=.false.

 ! identify the numerical method
 select case(trim(model_decisions(iLookDECISIONS%num_method)%decision))
  case('itertive'); num_method=itertive  ! iterative
  case('non_iter'); num_method=non_iter  ! non-iterative
  case('itersurf'); num_method=itersurf  ! iterate only on the surface energy balance
  case default
   err=10; message=trim(message)//"unknown option for the numerical method"; return
 end select

 ! identify if there is a need to transpire
 fTranspire=.true.
 select case(trim(model_decisions(iLookDECISIONS%bound_cond)%decision))
  case('headflux'); fTranspire=.false. 
  case('headhead'); fTranspire=.false.
 end select
 if(nSnow>0) fTranspire=.false.

 ! identify the method used to calculate flux derivatives
 select case(trim(model_decisions(iLookDECISIONS%fDerivMeth)%decision))
  case('numericl'); fDerivMeth=numerical
  case('analytic'); fDerivMeth=analytical
  case default
   err=10; message=trim(message)//"unknown method used to calculate flux derivatives [option="//trim(model_decisions(iLookDECISIONS%fDerivMeth)%decision)//"]"; return
 end select

 ! assign variables to data structures
 call assignVars(err,cmessage)
 if(err/=0)then; message=trim(message)//trim(cmessage); return; endif

 ! compute initial fluxes 
 if(iter==1)then


  ! compute the surface flux (J m-2 s-1) and the derivative in the surface flux w.r.t. surface temperature (J m-2 s-1 K-1)
  call 




aerodynamic and canopy conductance (lagged, based on values at iter=m)
  call flxConduct(mLayerTemp(1),      & ! intent(in): surface temperature (K)
                  mLayerMatricHead,   & ! intent(in): matric head in each layer (m)
                  err,cmessage)         ! intent(out): error control
  if(err/=0)then; message=trim(message)//trim(cmessage); return; endif
  ! compute the energy flux at layer interfaces at the start of the time step
  !  --> populates vector iLayerInitNrgFlux
  call iLayer_nrg(mLayerTemp,         & ! intent(in): trial temperature (K)
                  iLayerInitNrgFlux,  & ! intent(out): initial energy flux at layer interfaces (W m-2)
                  err,cmessage)         ! intent(out): error control
  if(err/=0)then; message=trim(message)//trim(cmessage); return; endif
  ! compute transpiration for each layer at the start of the time step (m s-1)
  scalarPotentialET = iden_air * windspd * scalarExCoef * (spechum - relhm2sphm(RHsurf,airpres,mLayerTemp(1))) ! kg m-2 s-1
  if(fTranspire)then; mLayerInitTranspire(1:nSoil) = (mLayerTranspireLim(1:nSoil)*mLayerRootDensity(1:nSoil)*scalarPotentialET)/iden_water
  else; mLayerInitTranspire(1:nSoil) = 0._dp; endif
 endif

 ! compute the aerodynamic and canopy conductance (lagged, based on values at iter=m)
 call flxConduct(mLayerTempIter(1),   & ! intent(in): surface temperature (K)
                 mLayerMatricHeadIter,& ! intent(in): matric head in each layer (m)
                 err,cmessage)          ! intent(out): error control
 if(err/=0)then; message=trim(message)//trim(cmessage); return; endif

 ! compute the residual vector (J m-3)
 print*, 'main residual calculation'
 call e_residual(dt,                  & ! intent(in): time step (seconds)
                 mLayerTempIter,      & ! intent(in): trial temperature (K)
                 mLayerVolFracIceIter,& ! intent(in): trial volumetric fraction of ice (-)
                 rvec,                & ! intent(out): residual vector (J m-3)
                 err,cmessage)          ! intent(out): error control
 if(err/=0)then; message=trim(message)//trim(cmessage); return; endif

 ! ***** compute the derivative in the freezing curve w.r.t. temperature (K-1)
 do iLayer=1,nLayers
  select case(layerType(iLayer))
   ! ***** process snow layers *****
   case(ix_snow)
    theta = mLayerVolFracIceIter(iLayer)*(iden_ice/iden_water) + mLayerVolFracLiqIter(iLayer)
    mLayerdTheta_dTk(iLayer) = dFracLiq_dTk(mLayerTempIter(iLayer),snowfrz_scale)*theta
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

 ! compute the derivative in flux at layer interfaces w.r.t. temperature in the layer above and in the layer below
 call dFlux_dTemp(fDerivMeth,            & ! intent(in): method used to compute derivatives (analytical or numerical)
                  mLayerTempIter,        & ! intent(in): trial temperature (K)
                  mLayerVolFracIceIter,  & ! intent(in): trial volumetric fraction of ice (-)
                  mLayerVolFracLiqIter,  & ! intent(in): trial volumetric fraction of liquid water (-)
                  dFlux_dTempAbove,      & ! intent(out): derivatives in the flux w.r.t. temperature in the layer above (W m-2 K-1)
                  dFlux_dTempBelow,      & ! intent(out): derivatives in the flux w.r.t. temperature in the layer below (W m-2 K-1)
                  err,cmessage)            ! intent(out): error control
 if(err/=0)then; message=trim(message)//trim(cmessage); return; endif
 ! check against the numerical solution
 !call dFlux_dTemp(numerical,             & ! intent(in): method used to compute derivatives (analytical or numerical)
 !                 mLayerTempIter,        & ! intent(in): trial temperature (K)
 !                 mLayerVolFracIceIter,  & ! intent(in): trial volumetric fraction of ice (-)
 !                 mLayerVolFracLiqIter,  & ! intent(in): trial volumetric fraction of liquid water (-)
 !                 dFlux_dTempAbove,      & ! intent(out): derivatives in the flux w.r.t. temperature in the layer above (W m-2 K-1)
 !                 dFlux_dTempBelow,      & ! intent(out): derivatives in the flux w.r.t. temperature in the layer below (W m-2 K-1)
 !                 err,cmessage)            ! intent(out): error control
 !if(err/=0)then; message=trim(message)//trim(cmessage); return; endif

 ! assemble the tri-diagonal matrix
 wtim = (1._dp - wimplicit)*dt  ! weighted time
 diag = (wtim/mLayerDepth)*(-dFlux_dTempBelow(0:nLayers-1) + dFlux_dTempAbove(1:nLayers)) + mLayerdTheta_dTk*LH_fus*iden_water + mLayerVolHtCapBulk
 d_m1 = (wtim/mLayerDepth(1:nLayers-1))*(-dFlux_dTempAbove(1:nLayers-1) )
 d_p1 = (wtim/mLayerDepth(1:nLayers-1))*( dFlux_dTempBelow(1:nLayers-1) )
 !write(*,'(a,1x,5(e20.10,1x))') 'diag = ', diag(46:50)
 !write(*,'(a,1x,5(e20.10,1x))') 'd_m1 = ', d_m1(46:50)
 !write(*,'(a,1x,5(e20.10,1x))') 'd_p1 = ', d_p1(46:50)

 ! assemble the tri-diagonal matrix
 !call getTridiag(dt,                  & ! intent(in): time step (s)
 !                mLayerTempIter,      & ! intent(in): trial temperature (K)
 !                mLayerVolFracIceIter,& ! intent(in): trial volumetric fraction of ice (-)
 !                mLayerVolFracLiqIter,& ! intent(in): trial volumetric fraction of liquid water (-)
 !                d_m1,                & ! intent(out): sub-diagonal elements of the tridiagonal system (J m-3 K-1)
 !                diag,                & ! intent(out): diagonal elements of the tridiagonal system (J m-3 K-1)
 !                d_p1,                & ! intent(out): super-diagonal elements of the tridiagonal system (J m-3 K-1)
 !                err,cmessage)          ! intent(out): error control
 !if(err/=0)then; message=trim(message)//trim(cmessage); return; endif
 !write(*,'(a,1x,5(e20.10,1x))') 'diag = ', diag(46:50)
 !write(*,'(a,1x,5(e20.10,1x))') 'd_m1 = ', d_m1(46:50)
 !write(*,'(a,1x,5(e20.10,1x))') 'd_p1 = ', d_p1(46:50)

 ! ***** compute the Jacobian using one-sided finite differences
 if(test_jac)then
  print*, 'starting to compute Jacobian'
  xtry=mLayerTempIter
  xsav=xtry
  h=EPS*abs(xsav)
  where (h == 0.0) h=EPS
  xph=xsav+h
  h=xph-xsav
  do ijac=1,5
   xtry(ijac)=xph(ijac)
   call e_residual(dt,                  & ! intent(in): time step (seconds)
                   xtry,                & ! intent(in): trial temperature (K)
                   mLayerVolFracIceIter,& ! intent(in): trial volumetric fraction of ice (-)
                   ftest,               & ! intent(out): residual vector (J m-3)
                   err,cmessage)          ! intent(out): error control
   if(err/=0)then; err=10; message=trim(message)//trim(cmessage); return; endif
   jmat(:,ijac)=(ftest(:)-rvec(:))/h(ijac)
   xtry(ijac)=xsav(ijac)
  end do
  print*, 'Jacobian = '
  do ijac=1,5
   write(*,'(i4,1x,5(e20.10,1x))') ijac, jmat(1:5,ijac)
  end do
  write(*,'(a,1x,5(e20.10,1x))') 'rvec = ', rvec(1:5) 
  write(*,'(a,1x,5(e20.10,1x))') 'diag = ', diag(1:5)
  write(*,'(a,1x,5(e20.10,1x))') 'd_m1 = ', d_m1(1:5)
  write(*,'(a,1x,5(e20.10,1x))') 'd_p1 = ', d_p1(1:5)
  print*, 'mLayerTempIter = ', mLayerTempIter(1:5)
  stop 'after computing Jacobian '
 endif
 
 ! solve the tridiagonal system of equations -- returns mLayerTempDiff
 call tridag(d_m1,                    & ! intent(in): sub-diagonal elements of the tridiagonal system (J m-3 K-1)
             diag,                    & ! intent(in): diagonal elements of the tridiagonal system (J m-3 K-1)
             d_p1,                    & ! intent(in): super-diagonal elements of the tridiagonal system (J m-3 K-1)
             -rvec,                   & ! intent(in): residual vector (J m-3)
             mLayerTempDiff,          & ! intent(out): temperature increment (K)
             err,cmessage)              ! intent(out): error control
 if(err/=0)then; message=trim(message)//trim(cmessage); return; endif
 print*, 'mLayerTempDiff = ', mLayerTempDiff(1:5)
 !print*, 'mLayerTempIter = ', mLayerTempIter(1:50)
 !print*, 'lowerBoundTemp = ', lowerBoundTemp
 !print*, 'temp gradient = ', mLayerTempIter(2:50) - mLayerTempIter(1:49)

 ! check for oscillations (borrowed from the SHAW model)
 if(iter > 10)then
  ! if delta is jumping between positive and negative, then cut the delta in half
  do iLayer=1,nLayers
   ! (case of oscillatory behavior)
   if(mLayerTempDiff(iLayer)*mLayerTempDiffOld(iLayer) < 0._dp)then
    mLayerTempDiff(iLayer) = mLayerTempDiff(iLayer)/2._dp
    !message=trim(message)//'check fix to oscillatory behavior'
    !err=20; return
   endif  ! (if there is an oscillation)
  end do
 endif

 ! check
 !if(printflag)then
 ! do iLayer=45,50
 !  write(*,'(a,2(i4,1x),10(f15.8,1x))') 'temp start, temp iter, ice iter, temp increment, nrg increment',&
 !   iLayer, layerType(iLayer), mLayerTemp(iLayer), mLayerTempIter(iLayer), mLayerVolFracIceIter(iLayer), mLayerTempDiff(iLayer), mLayerTempDiff(iLayer)*mLayerVolHtCapBulk(iLayer)
 ! end do
 !endif

 ! initialize flag for excessive temperature extrapolation
 crossFlag=.false.

 ! adjust del temperature in cases where snow temperature exceeds Tfreeze -- use bi-section
 if(nSnow>0)then
  do iLayer=1,nSnow
   if(mLayerTempIter(iLayer)+mLayerTempDiff(iLayer) > Tfreeze)then
    mLayerTempDiff(iLayer) = (Tfreeze-mLayerTempIter(iLayer))*0.5_dp
    crossFlag=.true.
   endif
  end do
 endif

 ! adjust del temperature in cases where soil temperature crosses the critical temperature
 do iLayer=nSnow+1,nLayers
  !if(iLayer > 45)then
  ! write(*,'(a,i4,2(f20.15,1x),e20.10)') 'iLayer, mLayerTcrit(iLayer-nSnow), mLayerTempIter(iLayer), mLayerTempDiff(iLayer) = ',&
  !                                        iLayer, mLayerTcrit(iLayer-nSnow), mLayerTempIter(iLayer), mLayerTempDiff(iLayer)
  !endif
  ! get the difference from the critical temperature (K)
  critDiff = mLayerTcrit(iLayer-nSnow) - mLayerTempIter(iLayer)
  ! set temperature to Tcrit in cases where temperatures cross Tcrit
  if(critDiff > 0._dp)then  ! (mLayerTempIter < Tcrit)
   if(mLayerTempDiff(iLayer) > critDiff)then; mLayerTempDiff(iLayer) = critDiff + epsT; crossFlag=.true.; endif
  else                      ! (mLayerTempIter > Tcrit)
   if(mLayerTempDiff(iLayer) < critDiff)then; mLayerTempDiff(iLayer) = critDiff - epsT; crossFlag=.true.; endif
  endif
 end do

 ! update temperature and compute phase change if the soil temperature crosses the critical temperature
 if(crossFlag .or. num_method==non_iter)then

  ! update temperature
  mLayerTempNew = mLayerTempIter + mLayerTempDiff
  ! compute phase change
  call phsechange(mLayerTempNew,       & ! intent(in): new temperature vector (K)
                  mLayerMatricHeadIter,& ! intent(in): matric head at the current iteration (m)
                  mLayerVolFracLiqIter,& ! intent(in): volumetric fraction of liquid water at the current iteration (-)
                  mLayerVolFracIceIter,& ! intent(in): volumetric fraction of ice at the current iteration (-)
                  mLayerMatricHeadNew, & ! intent(out): new matric head (m)
                  mLayerVolFracLiqNew, & ! intent(out): new volumetric fraction of liquid water (-)
                  mLayerVolFracIceNew, & ! intent(out): new volumetric fraction of ice (-)
                  err,cmessage)          ! intent(out): error control
  if(err/=0)then; message=trim(message)//trim(cmessage); return; endif

 ! otherwise, conduct a line search
 else

  ! compute the gradient of the function vector (J m-3 J m-3 K-1)
  do iLayer=1,nLayers
   if(iLayer==1)then;           g(iLayer) =                                 diag(iLayer)*rvec(iLayer) + d_p1(iLayer)*rvec(iLayer+1)
   elseif(iLayer==nLayers)then; g(iLayer) = d_m1(iLayer-1)*rvec(iLayer-1) + diag(iLayer)*rvec(iLayer)
   else;                        g(iLayer) = d_m1(iLayer-1)*rvec(iLayer-1) + diag(iLayer)*rvec(iLayer) + d_p1(iLayer)*rvec(iLayer+1)
   endif
  end do

  ! compute the function value (J m-3 J m-3)
  fold = 0.5_dp*dot_product(rvec,rvec)
 
  ! compute maximum step size (K)
  stpmax=STPMX*real(nLayers,dp)
  !stpmax=STPMX*max(sqrt(dot_product(mLayerTempIter,mLayerTempIter)),real(nLayers,dp))

  ! conduct the line search
  call lnsrch(dt,                      & ! intent(in): time step (seconds)
              mLayerTempIter,          & ! intent(in): trial temperature vector (K)
              stpmax,                  & ! intent(in): maximum step size (K)
              fold,                    & ! intent(in): function value for trial temperature vector (J m-3 J m-3)
              g,                       & ! intent(in): gradient of the function vector (J m-3 J m-3 K-1)
              mLayerTempDiff,          & ! intent(in): iteration increment (K)
              mLayerMatricHeadIter,    & ! intent(in): matric head at the current iteration (m)
              mLayerVolFracLiqIter,    & ! intent(in): volumetric fraction of liquid water at the current iteration (-)
              mLayerVolFracIceIter,    & ! intent(in): volumetric fraction of ice at the current iteration (-)
              mLayerTempNew,           & ! intent(out): new temperature vector (K)
              rvec,                    & ! intent(out): new residual vector (J m-3)
              fnew,                    & ! intent(out): new function value (J m-3 J m-3)
              mLayerMatricHeadNew,     & ! intent(out): new matric head (m)
              mLayerVolFracLiqNew,     & ! intent(out): new volumetric fraction of liquid water (-)
              mLayerVolFracIceNew,     & ! intent(out): new volumetric fraction of ice (-)
              err,cmessage)              ! intent(out): error control
  ! negative error code requires convergence check (done in upwind routine), so just check positive errors
  if(err>0)then; message=trim(message)//trim(cmessage); return; endif

 endif

 ! compute un-stressed ET (kg m-2 s-1)
 scalarPotentialET = iden_air * windspd * scalarExCoef * (spechum - relhm2sphm(RHsurf,airpres,mLayerTempNew(1)))

 ! compute transpiration (m s-1)
 if(fTranspire)then
  ! ** compute actual transpiration from each soil layer (m s-1)
  mLayerTranspire(1:nSoil) =  (mLayerTranspireLim(1:nSoil)*mLayerRootDensity(1:nSoil)*scalarPotentialET)/iden_water
  ! ** compute total transpiration (kg m-2 s-1)
  scalarMassLiquid = ( wimplicit*(sum(mLayerInitTranspire)) + (1._dp - wimplicit)*sum(mLayerTranspire) )*iden_water
 else
  mLayerTranspire(1:nSoil) = 0._dp
  scalarMassLiquid = 0._dp
 endif

 ! compute sublimation/frost (kg m-2 s-1)
 if(nSnow==0) scalarMassSolid  = 0._dp
 if(nsnow >0) scalarMassSolid  = scalarPotentialET

 ! update sensible and latent heat (W m-2) --> positive downwards
 scalarSenHeat = scalarExSen*(airtemp - mLayerTempNew(1))
 if(nSnow==0) scalarLatHeat = scalarMassLiquid*LH_vap
 if(nsnow >0) scalarLatHeat = scalarMassSolid*LH_sub

 ! check that total transpiration matches the latent heat flux
 !print*, 'check ET ', scalarMassLiquid*LH_vap, scalarLatHeat

 ! ====================================================================================================================

 end subroutine tempchange


 ! ************************************************************************************************
 ! new subroutine: compute phase change impacts on matric head and volumetric liquid water and ice
 ! ************************************************************************************************
 subroutine phsechange(mLayerTempNew,       & ! intent(in): new temperature vector (K)
                       mLayerMatricHeadIter,& ! intent(in): matric head at the current iteration (m)
                       mLayerVolFracLiqIter,& ! intent(in): volumetric fraction of liquid water at the current iteration (-)
                       mLayerVolFracIceIter,& ! intent(in): volumetric fraction of ice at the current iteration (-)
                       mLayerMatricHeadNew, & ! intent(out): new matric head (m)
                       mLayerVolFracLiqNew, & ! intent(out): new volumetric fraction of liquid water (-)
                       mLayerVolFracIceNew, & ! intent(out): new volumetric fraction of ice (-)
                       err,message)           ! intent(out): error control
 ! physical constants
 USE multiconst,only:&
                     Tfreeze,  & ! freezing point              (K)
                     iden_ice, & ! intrinsic density of ice    (kg m-3)
                     iden_water  ! intrinsic density of water  (kg m-3)
 ! utility routines
 USE snow_utils_module,only:fracliquid    ! compute volumetric fraction of liquid water
 USE soil_utils_module,only:volFracLiq    ! compute volumetric fraction of liquid water based on matric head
 USE soil_utils_module,only:matricHead    ! compute the matric head based on volumetric liquid water content
 ! data structures
 USE data_struc,only:mpar_data,mvar_data,indx_data,ix_soil,ix_snow    ! data structures
 USE var_lookup,only:iLookPARAM,iLookMVAR,iLookINDEX                  ! named variables for structure elements
 implicit none
 ! input variables
 real(dp),intent(in)           :: mLayerTempNew(:)         ! new estimate of temperature (K)
 real(dp),intent(in)           :: mLayerMatricHeadIter(:)  ! before phase change: matric head (m)
 real(dp),intent(in)           :: mLayerVolFracLiqIter(:)  ! before phase change: volumetric fraction of liquid water (-)
 real(dp),intent(in)           :: mLayerVolFracIceIter(:)  ! before phase change: volumetric fraction of ice (-)
 ! output variables
 real(dp),intent(out)          :: mLayerMatricHeadNew(:)   ! after phase change: matric head (m)
 real(dp),intent(out)          :: mLayerVolFracLiqNew(:)   ! after phase change: volumetric fraction of liquid water (-)
 real(dp),intent(out)          :: mLayerVolFracIceNew(:)   ! after phase change: volumetric fraction of ice (-)
 integer(i4b),intent(out)      :: err                      ! error code
 character(*),intent(out)      :: message                  ! error message
 ! local pointers to model parameters
 real(dp),pointer              :: snowfrz_scale            ! scaling parameter for the snow freezing curve (K-1)
 real(dp),pointer              :: vGn_alpha                ! van Genutchen "alpha" parameter
 real(dp),pointer              :: vGn_n                    ! van Genutchen "n" parameter
 real(dp),pointer              :: theta_sat                ! soil porosity (-)
 real(dp),pointer              :: theta_res                ! soil residual volumetric water content (-)
 real(dp),pointer              :: vGn_m                    ! van Genutchen "m" parameter (-)
 real(dp),pointer              :: kappa                    ! constant in the freezing curve function (m K-1)
 ! local pointers to model variables
 real(dp),pointer              :: mLayerTcrit(:)           ! critical soil temperature above which all water is unfrozen (K)
 integer(i4b),pointer          :: layerType(:)             ! type of the layer (ix_soil or ix_snow)
 ! define local variables
 real(dp)                      :: theta                    ! liquid water equivalent of total water (-)
 integer(i4b)                  :: nSnow                    ! number of snow layers
 integer(i4b)                  :: iLayer                   ! index of model layer
 logical(lgt)                  :: printflag                ! flag to print debug information
 ! initialize error control
 err=0; message="phsechange/"

 ! initialize print flag
 printflag=.false.

 ! assign pointers to model parameters
 snowfrz_scale       => mpar_data%var(iLookPARAM%snowfrz_scale)              ! scaling parameter for the snow freezing curve (K-1)
 vGn_alpha           => mpar_data%var(iLookPARAM%vGn_alpha)                  ! van Genutchen "alpha" parameter (m-1)
 vGn_n               => mpar_data%var(iLookPARAM%vGn_n)                      ! van Genutchen "n" parameter (-)
 theta_sat           => mpar_data%var(iLookPARAM%theta_sat)                  ! soil porosity (-)
 theta_res           => mpar_data%var(iLookPARAM%theta_res)                  ! soil residual volumetric water content (-)
 vGn_m               => mvar_data%var(iLookMVAR%scalarVGn_m)%dat(1)          ! van Genutchen "m" parameter (-)
 kappa               => mvar_data%var(iLookMVAR%scalarKappa)%dat(1)          ! constant in the freezing curve function (m K-1)

 ! assign pointers to index variables
 mLayerTcrit         => mvar_data%var(iLookMVAR%mLayerTcrit)%dat             ! critical soil temperature above which all water is unfrozen (K) 
 layerType           => indx_data%var(iLookINDEX%layerType)%dat              ! layer type (ix_soil or ix_snow)

 ! identify the number of snow layers
 nSnow = count(layerType==ix_snow)

 ! update volumetric liquid and ice content (-)
 do iLayer=1,size(layerType)  ! (process snow and soil separately)
  ! compute liquid water equivalent of total water (liquid plus ice)
  theta = mLayerVolFracIceIter(iLayer)*(iden_ice/iden_water) + mLayerVolFracLiqIter(iLayer)
  select case(layerType(iLayer))
   ! ** snow
   case(ix_snow)
    ! compute the volumetric fraction of liquid water and ice (-)
    mLayerVolFracLiqNew(iLayer) = fracliquid(mLayerTempNew(iLayer),snowfrz_scale)*theta
    mLayerVolFracIceNew(iLayer) = (theta - mLayerVolFracLiqNew(iLayer))*(iden_water/iden_ice)
   ! ** soil
   case(ix_soil)
    ! compute the matric head (m) volumetric fraction of liquid water and ice (-)
    if(mLayerTempNew(iLayer)<mLayerTcrit(iLayer-nSnow))then
     mLayerMatricHeadNew(iLayer-nSnow) = kappa*(mLayerTempNew(iLayer) - Tfreeze)
     mLayerVolFracLiqNew(iLayer)       = volFracLiq(mLayerMatricHeadNew(iLayer-nSnow),&
                                                    vGn_alpha,theta_res,theta_sat,vGn_n,vGn_m)
     mLayerVolFracIceNew(iLayer)       = (theta - mLayerVolFracLiqNew(iLayer))*(iden_water/iden_ice)
    else
     ! update matric head when all water is **unfrozen** -- if matric head > 0 at iter=m then no change in matric head
     if(mLayerMatricHeadIter(iLayer-nSnow) > 0._dp)then ! saturated at the start of the iteration
      mLayerMatricHeadNew(iLayer-nSnow) = mLayerMatricHeadIter(iLayer-nSnow)
     else
      ! some water is frozen at the start of the iteration      
      mLayerMatricHeadNew(iLayer-nSnow) = matricHead(min(theta,theta_sat),vGn_alpha,theta_res,theta_sat,vGn_n,vGn_m)
     endif
     ! update liquid water and ice content 
     mLayerVolFracLiqNew(iLayer)       = min(theta,theta_sat)
     mLayerVolFracIceNew(iLayer)       = 0._dp
    endif
   case default; err=10; message=trim(message)//'unknown case for model layer'; return
  endselect
  ! sanity check
  if(iLayer > 40 .and. printflag) &
   write(*,'(a,i4,1x,10(f20.10,1x))') 'in phase change, liquid (iter, new), ice (iter, new, diff)', &
    iLayer, mLayerVolFracLiqIter(iLayer), mLayerVolFracLiqNew(iLayer), mLayerVolFracIceIter(iLayer), mLayerVolFracIceNew(iLayer), &
    mLayerVolFracIceNew(iLayer) - mLayerVolFracIceIter(iLayer)
  if(mLayerVolFracIceNew(iLayer) < 0._dp)then
   write(message,'(a,i0,a,e20.10,a)')trim(message)//"volumetric ice content < 0 [iLayer=",iLayer,&
                                     &"; mLayerVolFracIceNew(iLayer)=",mLayerVolFracIceNew(iLayer),"]"
   err=10; return
  endif
 end do ! (looping through layers)
 endsubroutine phsechange


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
 real(dp),pointer              :: lambda_drysoil    ! thermal conductivity of dry soil (W m-1)
 real(dp),pointer              :: lambda_wetsoil    ! thermal conductivity of wet soil (W m-1)
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
 real(dp)                      :: lambda_wet             ! thermal conductivity of the wet material
 real(dp)                      :: kerstenNum             ! the Kersten number (-), defining weight applied to conductivity of the wet medium

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
 lambda_drysoil      => mvar_data%var(iLookMVAR%scalarLambda_drysoil)%dat(1) ! thermal conductivity of dry soil (W m-1)
 lambda_wetsoil      => mvar_data%var(iLookMVAR%scalarLambda_wetsoil)%dat(1) ! thermal conductivity of wet soil (W m-1)
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
   !mLayerThermalC(iLayer) = lambda_soil * (1._dp - theta_sat)      + & ! soil component
   !                         lambda_ice  * mLayerVolFracIce(iLayer) + & ! ice component
   !                         lambda_water* mLayerVolFracLiq(iLayer) + & ! liquid water component
   !                         lambda_air  * mLayerVolFracAir(iLayer)     ! air component
   ! compute the thermal conductivity of the wet material (W m-1)
   lambda_wet = lambda_wetsoil**(1._dp - theta_sat) * lambda_water**theta_sat * lambda_ice**(theta_sat - mLayerVolFracLiq(iLayer))
   ! compute the Kersten number (-)
   kerstenNum = log10( (mLayerVolFracIce(iLayer) + mLayerVolFracLiq(iLayer))/theta_sat ) + 1._dp
   ! ...and, compute the thermal conductivity
   mLayerThermalC(iLayer) = kerstenNum*lambda_wet + (1._dp - kerstenNum)*lambda_drysoil
   !print*, 'iLayer, mLayerVolFracLiq(iLayer), lambda_wet, lambda_drysoil, kerstenNum, mLayerThermalC(iLayer) = ', &
   !         iLayer, mLayerVolFracLiq(iLayer), lambda_wet, lambda_drysoil, kerstenNum, mLayerThermalC(iLayer)
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
 ! ************************************************************************************************
 ! ************************************************************************************************
 ! ***** PRIVATE SUBROUTINES **********************************************************************
 ! ************************************************************************************************
 ! ************************************************************************************************
 ! ************************************************************************************************

 ! ************************************************************************************************
 ! private subroutine: compute the aerodynamic & canopy conductance (lagged, based on values at iter=m) 
 ! ************************************************************************************************
 subroutine flxConduct(surfTemp,&             ! intent(in): surface temperature (K)
                       mLayerMatricHeadIter,& ! intent(in): matric head in each layer (m)
                       err,message)           ! intent(out): error control
 ! physical constants
 USE multiconst,only:&
                     Cp_air,   &            ! specific heat of air        (J kg-1 K-1)
                     LH_fus,   &            ! latent heat of fusion       (J kg-1)
                     LH_vap,   &            ! latent heat of vaporization (J kg-1)
                     LH_sub,   &            ! latent heat of sublimation  (J kg-1)
                     iden_air               ! intrinsic density of air    (kg m-3)
 ! utility modules
 USE snow_utils_module,only:bulkRichardson  ! compute the bulk Richardson number (-) and its derivative w.r.t. temperature (K-1)
 USE snow_utils_module,only:astability      ! compute atmospheric stability
 ! compute the aerodynamic and canopy conductance (lagged, based on values at iter=m)
 implicit none
 ! input
 real(dp),intent(in)           :: surfTemp                 ! surface temperature (K)
 real(dp),intent(in)           :: mLayerMatricHeadIter(:)  ! matric head in each layer (m)





 ! output
 ! (populates shared variables: scalarExSen, scalarExLat, scalarTranspireLim)
 integer(i4b),intent(out)      :: err                      ! error code
 character(*),intent(out)      :: message                  ! error message
 ! local variables
 character(LEN=256)            :: cmessage                 ! error message of downwind routine 
 real(dp)                      :: RiBulk                   ! bulk Richardson number (-)
 real(dp)                      :: dRiBulk_dTemp            ! derivative in the bulk Richardson number w.r.t. temperature (K-1) 
 real(dp)                      :: aerodynResist            ! aerodynamic resistance (s m-1)
 real(dp)                      :: stomatalResist           ! stomatal resistance (s m-1)
 real(dp)                      :: plantWiltFactor          ! plant wilting factor (-) 
 integer(i4b)                  :: iLayer                   ! index of model layers
 ! initialize error control
 err=0; message="flxConduct/"

 ! compute the bulk Richardson number and its derivative
 call bulkRichardson(airtemp,surfTemp,windspd,mheight, & ! (input)
                     RiBulk,dRiBulk_dTemp,err,cmessage)  ! (output)
 if(err/=0)then; err=10; message=trim(message)//trim(cmessage); return; endif

 ! compute stability corrections
 call astability(RiBulk,scalarExCoef,err,cmessage)
 if(err/=0)then; err=10; message=trim(message)//trim(cmessage); return; endif

 ! compute aerodynamic resistance (s m-1)
 aerodynResist = 1._dp / (scalarExCoef*windspd)

 ! compute the ratio of actual:potential evapotranspiration
 call evapResistance(aerodynResist,        &  ! intent(in): aerodynamic resistance (s m-1)
                     mLayerMatricHeadIter, &  ! intent(in): matric head in each layer (m)
                     scalarTranspireLim,   &  ! intent(out): resistance to evaporation at the surface (-)
                     mLayerTranspireLim,   &  ! intent(out): resistance to evaporation in each layer (-)
                     err,cmessage)            ! intent(out): error control
 if(err/=0)then; err=10; message=trim(message)//trim(cmessage); return; endif
 
 ! compute the exchange factor for sensible and latent heat
 scalarExSen = Cp_air * iden_air * windspd * scalarExCoef       ! J m-2 s-1 K-1
 if(nSnow >0) scalarExLat = LH_sub * iden_air * windspd * scalarExCoef                       ! J m-2 s-1
 if(nSnow==0) scalarExLat = LH_vap * iden_air * windspd * scalarExCoef * scalarTranspireLim  ! J m-2 s-1

 end subroutine flxConduct

 ! ************************************************************************************************
 ! private subroutine: compute the resistance to evaporation (-)
 ! ************************************************************************************************
 subroutine evapResistance(aerodynResist,        &  ! intent(in): aerodynamic resistance (s m-1)
                           mLayerMatricHeadIter, &  ! intent(in): matric head in each layer (m)
                           evapResist,           &  ! intent(out): resistance to evaporation at the surface (-)
                           mLayerEvapResist,     &  ! intent(out): resistance to evaporation in each layer (-)
                           err,message)             ! intent(out): error control
 implicit none
 ! input
 real(dp),intent(in)           :: aerodynResist            ! aerodynamic resistance
 real(dp),intent(in)           :: mLayerMatricHeadIter(:)  ! matric head in each layer (m)
 ! output
 real(dp),intent(out)          :: evapResist               ! resistance to evaporation (-) 
 real(dp),intent(out)          :: mLayerEvapResist(:)
 integer(i4b),intent(out)      :: err                      ! error code
 character(*),intent(out)      :: message                  ! error message
 ! local variables
 real(dp)                      :: aerodynResist            ! aerodynamic resistance (s m-1)
 real(dp)                      :: stomatalResist           ! stomatal resistance (s m-1)
 real(dp)                      :: plantWiltFactor          ! plant wilting factor (-) 
 integer(i4b)                  :: iLayer                   ! index of model layers
 ! initialize error control
 err=0; message='evapResistance/'
 ! compute the ratio of actual:potential evapotranspiration
 if(nSnow>0)then
  evapResist = 1._dp
 else
  ! ** compute the factor limiting evaporation for each soil layer (-)
  evapResist = 0._dp  ! (initialize the weighted average)
  do iLayer=1,nSoil
   ! compute the stomatal resistance of the canopy (m s-1)
   plantWiltFactor = 1._dp + ( min(mLayerMatricHeadIter(iLayer),0._dp) / plantWiltPsi )**plantWiltExp
   stomatalResist  = (minStomatalResist/LAI)*plantWiltFactor
   ! compute the factor limiting evaporation for a given soil layer (-)
   mLayerEvapResist(iLayer) = aerodynResist / (stomatalResist + aerodynResist)
   ! commpute the weighted average (weighted by root density)
   evapResist = evapResist + mLayerEvapResist(iLayer)*mLayerRootDensity(iLayer)
  end do ! (looping through soil layers)
 endif  ! (if surface is snow-free)
 end function evapResistance


 ! ************************************************************************************************
 ! private function: compute surface energy flux
 ! ************************************************************************************************
 function surfaceFlx(surfaceTempIter)  ! intent(in): trial surface temperature (K)
 ! physical constants
 USE multiconst,only:&
                     Em_Sno,   & ! emissivity of snow          (-)
                     sigma       ! Stefan Boltzman constant    (W m-2 K-4)
 ! utility modules
 USE conv_funcs_module,only:relhm2sphm    ! compute specific humidity
 implicit none
 real(dp),intent(in)           :: surfaceTempIter      ! trial temperature at the current iteration (K)
 real(dp)                      :: surfaceFlx           ! energy flux at the surface (W m-2)

 ! compute the bulk Richardson number and its derivative
 call bulkRichardson(airtemp,surfTemp,windspd,mheight, & ! (input)
                     RiBulk,dRiBulk_dTemp,err,cmessage)  ! (output)
 if(err/=0)then; err=10; message=trim(message)//trim(cmessage); return; endif

 ! compute stability corrections
 call astability(RiBulk,scalarExCoef,err,cmessage)
 if(err/=0)then; err=10; message=trim(message)//trim(cmessage); return; endif

 ! compute aerodynamic resistance (s m-1)
 aerodynResist = 1._dp / (scalarExCoef*windspd)

 ! compute the ratio of actual:potential evapotranspiration
 call evapResistance(aerodynResist,        &  ! intent(in): aerodynamic resistance (s m-1)
                     mLayerMatricHeadIter, &  ! intent(in): matric head in each layer (m)
                     scalarTranspireLim,   &  ! intent(out): resistance to evaporation at the surface (-)
                     mLayerTranspireLim,   &  ! intent(out): resistance to evaporation in each layer (-)
                     err,cmessage)            ! intent(out): error control
 if(err/=0)then; err=10; message=trim(message)//trim(cmessage); return; endif

 ! compute the exchange factor for sensible and latent heat
 scalarExSen = Cp_air * iden_air * windspd * scalarExCoef       ! J m-2 s-1 K-1
 if(nSnow >0) scalarExLat = LH_sub * iden_air * windspd * scalarExCoef                       ! J m-2 s-1
 if(nSnow==0) scalarExLat = LH_vap * iden_air * windspd * scalarExCoef * scalarTranspireLim  ! J m-2 s-1





 ! compute sensible and latent heat at the current iteration (W m-2) --> positive downwards
 scalarSenHeat = scalarExSen*(airtemp - surfaceTempIter)
 scalarLatHeat = scalarExLat*(spechum - relhm2sphm(RHsurf,airpres,surfaceTempIter))

 ! compute the surface energy flux -- positive downwards
 surfaceFlx    = sw_down*(1._dp - scalarAlbedo) + lw_down - Em_Sno*sigma*surfaceTempIter**4._dp + scalarSenHeat + scalarLatHeat

 end function surfaceFlx


 ! ************************************************************************************************
 ! private subroutine: compute energy fluxes at layer interfaces
 ! ************************************************************************************************
 subroutine iLayer_nrg(mLayerTempIter,      & ! intent(in): trial temperature (K)
                       iLayerNrgFlux,       & ! intent(out): energy flux at the layer interfaces (W m-2)
                       err,message)           ! intent(out): error control
 ! compute residual vector in the energy equation
 implicit none
 ! input
 real(dp),intent(in)           :: mLayerTempIter(:)        ! trial temperature at the current iteration (K)
 ! output
 real(dp),intent(out)          :: iLayerNrgFlux(0:)         ! energy flux at the layer interfaces (W m-2)
 integer(i4b),intent(out)      :: err                      ! error code
 character(*),intent(out)      :: message                  ! error message
 ! local variables
 integer(i4b)                  :: iLayer                   ! index of model layers
 ! initialize error control
 err=0; message='iLayer_nrg/'

 ! ***** compute fluxes at layer interfaces
 ! compute flux at the upper boundary -- positive downwards
 iLayerNrgFlux(0) = surfaceFlx(mLayerTempIter(1))
 ! compute fluxes within the domain -- positive downwards
 do iLayer=1,nLayers-1
  iLayerNrgFlux(iLayer) = -iLayerThermalC(iLayer)*(mLayerTempIter(iLayer+1) - mLayerTempIter(iLayer)) / &
                                                  (mLayerHeight(iLayer+1) - mLayerHeight(iLayer))
 end do
 ! compute fluxes at the lower boundary [NOTE: use mLayerThermalC(iLayer)]
 iLayerNrgFlux(nLayers) = -mLayerThermalC(iLayer)*(lowerBoundTemp - mLayerTempIter(iLayer))/(mLayerDepth(iLayer)*0.5_dp)

 end subroutine iLayer_nrg


 ! ************************************************************************************************
 ! private subroutine: compute residual vector in the energy equation 
 ! ************************************************************************************************
 subroutine e_residual(dt,                  & ! intent(in): time step (seconds)
                       mLayerTempIter,      & ! intent(in): trial temperature (K)
                       mLayerVolFracIceIter,& ! intent(in): trial volumetric fraction of ice (-)
                       rvec,                & ! intent(out): residual vector (J m-3)
                       err,message)           ! intent(out): error control
 ! physical constants
 USE multiconst,only:&
                     Cp_air,   & ! specific heat of air        (J kg-1 K-1)
                     LH_fus,   & ! latent heat of fusion       (J kg-1)
                     Em_Sno,   & ! emissivity of snow          (-)
                     sigma,    & ! Stefan Boltzman constant    (W m-2 K-4)
                     iden_air, & ! intrinsic density of air    (kg m-3)
                     iden_ice    ! intrinsic density of ice    (kg m-3)
 ! utility modules
 USE conv_funcs_module,only:relhm2sphm    ! compute specific humidity
 ! compute residual vector in the energy equation
 implicit none
 ! input
 real(dp),intent(in)           :: dt                       ! time step (seconds)
 real(dp),intent(in)           :: mLayerTempIter(:)        ! trial temperature at the current iteration (K)
 real(dp),intent(in)           :: mLayerVolFracIceIter(:)  ! volumetric fraction of ice at the current iteration (-)
 ! output
 real(dp),intent(out)          :: rvec(:)                  ! residual vector (J m-3)
 integer(i4b),intent(out)      :: err                      ! error code
 character(*),intent(out)      :: message                  ! error message
 ! local variables
 character(LEN=256)            :: cmessage                 ! error message of downwind routine
 real(dp),parameter            :: RHsurf=1._dp             ! relative humidity of the surface (-)
 real(dp)                      :: iLayerNrgFlux(0:size(mLayerTempIter)) ! energy flux at layer interfaces (J m-2 s-1) 
 integer(i4b)                  :: iLayer                   ! index of model layers
 real(dp)                      :: nrg0,nrg1                ! energy content at the start of the time step / current iteration (J m-3)
 real(dp)                      :: flx0,flx1                ! fluxes at the start of the time step / current iteration (J m-3)
 real(dp)                      :: phse                     ! phase change term (J m-3)
 ! initialize error control
 err=0; message='e_residual/'

 ! compute energy fluxes at layer interfaces
 call iLayer_nrg(mLayerTempIter,      & ! intent(in): trial temperature (K)
                 iLayerNrgFlux,       & ! intent(out): energy flux at the layer interfaces (W m-2)
                 err,cmessage)          ! intent(out): error control
 if(err/=0)then; message=trim(message)//trim(cmessage); return; endif


 ! print header
 print*, 'iLayer, rvec(iLayer), nrg0, nrg1, (nrg1 - nrg0), flx0, flx1, phse, temp'

 ! compute the residual vector (J m-3)
 do iLayer=1,nLayers
  ! (compute individual terms)
  nrg0 = mLayerVolHtCapBulk(iLayer)*mLayerTemp(iLayer)      ! energy content at the start of the time step (J m-3)
  nrg1 = mLayerVolHtCapBulk(iLayer)*mLayerTempIter(iLayer)  ! energy content at the current iteration (J m-3)
  flx0 = -(iLayerInitNrgFlux(iLayer) - iLayerInitNrgFlux(iLayer-1))*(dt/mLayerDepth(iLayer))  ! flux at the start of the time step (J m-3)
  flx1 = -(iLayerNrgFlux(iLayer) - iLayerNrgFlux(iLayer-1))*(dt/mLayerDepth(iLayer))          ! flux at the current iteration (J m-3)
  phse = LH_fus*iden_ice*(mLayerVolFracIceIter(iLayer) - mLayerVolFracIce(iLayer))   ! phase change term (J m-3)
  ! (compute residuals)
  rvec(iLayer) = nrg1 - (nrg0 + flx0*wimplicit + flx1*(1._dp - wimplicit) + phse)
  ! (print progress)
  if(iLayer < 5) write(*,'(a,1x,i4,1x,f9.3,1x,10(e20.10,1x))') 'residuals = ', iLayer, mLayerTempIter(iLayer), rvec(iLayer), nrg0, nrg1, (nrg1 - nrg0), flx0, flx1, phse
  !rvec(iLayer) = mLayerVolHtCapBulk(iLayer)*mLayerTempIter(iLayer) - &
  !               (mLayerVolHtCapBulk(iLayer)*mLayerTemp(iLayer) + &
  !               -(iLayerInitNrgFlux(iLayer) - iLayerInitNrgFlux(iLayer-1))*(dt/mLayerDepth(iLayer))*wimplicit + &
  !               -(iLayerNrgFlux(iLayer) - iLayerNrgFlux(iLayer-1))*(dt/mLayerDepth(iLayer))*(1._dp - wimplicit) + &
  !               LH_fus*iden_ice*(mLayerVolFracIceIter(iLayer) - mLayerVolFracIce(iLayer)) )
  !if(iLayer > 45) write(*,'(a,1x,i4,1x,10(e20.10,1x))') 'residuals = ', iLayer, rvec(iLayer)
 end do
 endsubroutine e_residual


 ! ************************************************************************************************
 ! private subroutine: compute derivative in fluxes at layer interfaces w.r.t.
 !                      temperature in the layer above and the layer below
 ! ************************************************************************************************
 subroutine dFlux_dTemp(dMethod,               & ! intent(in): method used to compute derivatives (analytical or numerical)
                        mLayerTempIter,        & ! intent(in): trial temperature (K)
                        mLayerVolFracIceIter,  & ! intent(in): trial volumetric fraction of ice (-)
                        mLayerVolFracLiqIter,  & ! intent(in): trial volumetric fraction of liquid water (-)
                        dFlux_dTempAbove,      & ! intent(out): derivatives in the flux w.r.t. temperature in the layer above (W m-2 K-1)
                        dFlux_dTempBelow,      & ! intent(out): derivatives in the flux w.r.t. temperature in the layer below (W m-2 K-1)
                        err,message)            ! intent(out): error control
 ! compute derivative in fluxes at layer interfaces w.r.t. matric head in the layer above and the layer below
 USE multiconst,only:          & ! (physical constants)
                     Em_Sno,   & ! emissivity of snow          (-)
                     sigma       ! Stefan Boltzman constant    (W m-2 K-4)
 implicit none
 ! input
 integer(i4b),intent(in)       :: dMethod                  ! method used to compute derivatives (analytical or numerical)
 real(dp),intent(in)           :: mLayerTempIter(:)        ! trial temperature (K)
 real(dp),intent(in)           :: mLayerVolFracIceIter(:)  ! trial volumetric fraction of ice (-)
 real(dp),intent(in)           :: mLayerVolFracLiqIter(:)  ! trial volumetric fraction of liquid water (-)
 ! output
 real(dp),intent(out)          :: dFlux_dTempAbove(0:)     ! derivatives in the flux w.r.t. matric head in the layer above (s-1)
 real(dp),intent(out)          :: dFlux_dTempBelow(0:)     ! derivatives in the flux w.r.t. matric head in the layer below (s-1)
 integer(i4b),intent(out)      :: err                      ! error code
 character(*),intent(out)      :: message                  ! error message
 ! local variables
 character(LEN=256)            :: cmessage                 ! error message of downwind routine 
 real(dp),parameter            :: dx=1.e-8_dp              ! finite difference increment
 integer(i4b)                  :: iLayer                   ! layer index
 real(dp)                      :: Qderiv                   ! derivative in specific humidity w.r.t. temperature (kg kg-1 K-1)
 real(dp)                      :: LW_deriv                 ! derivative in longwave radiation w.r.t. temperature (J m-2 s-1 K-1)
 real(dp)                      :: Qh_deriv                 ! derivative in sensible heat w.r.t. temperature (J m-2 s-1 K-1)
 real(dp)                      :: Qe_deriv                 ! derivative in latent heat w.r.t. temperature (J m-2 s-1 K-1)
 real(dp)                      :: dz                       ! spatial difference in height
 real(dp)                      :: flux0,flux1,flux2        ! used to test the numerical derivatives
 ! initialize error control
 err=0; message="dFlux_dTemp/"

 ! ***** compute derivatives in surface fluxes...
 ! **********************************************

 ! compute the derivative in specific humidity w.r.t. temperature (kg kg-1 K-1)
 Qderiv = dsphum_dTk(RHsurf,airpres,mLayerTempIter(1))

 ! compute the derivative in sensible and latent heat (J m-2 s-1 K-1)
 Qh_deriv = scalarExSen
 Qe_deriv = scalarExLat*Qderiv

 ! compute the derivative in longwave radiation (J m-2 s-1 K-1)
 LW_deriv = 4._dp*Em_Sno*sigma*mLayerTempIter(1)**3._dp

 ! ***** loop through layers...
 ! ****************************

 ! compute the derivative in flux at layer interfaces w.r.t. matric head
 do iLayer=0,nLayers  ! loop through interfaces

  ! ***** the upper boundary
  if(iLayer==0)then  ! (upper boundary)

   dFlux_dTempAbove(iLayer) = -huge(Qderiv)  ! don't expect this to be used, so deliberately set to a ridiculous value to cause problems
   ! derivatives in the flux w.r.t. temperature
   if(dMethod==analytical)then
    dFlux_dTempBelow(iLayer) = -(LW_deriv + Qh_deriv + Qe_deriv)
   ! compute numerical derivatives
   else
    flux0 = surfaceFlx(mLayerTempIter(1)) 
    flux1 = surfaceFlx(mLayerTempIter(1)+dx) 
    dFlux_dTempBelow(iLayer) = (flux1 - flux0)/dx
   endif

  ! ***** the lower boundary
  elseif(iLayer==nLayers)then  ! (lower boundary)

   dFlux_dTempBelow(iLayer) = -huge(Qderiv)  ! don't expect this to be used, so deliberately set to a ridiculous value to cause problems
   ! derivatives in the flux w.r.t. temperature
   if(dMethod==analytical)then
    dFlux_dTempAbove(iLayer) = mLayerThermalC(iLayer)/(mLayerDepth(iLayer)*0.5_dp)
   else
    flux0 = -mLayerThermalC(iLayer)*(lowerBoundTemp - (mLayerTempIter(iLayer)   ))/(mLayerDepth(iLayer)*0.5_dp)
    flux1 = -mLayerThermalC(iLayer)*(lowerBoundTemp - (mLayerTempIter(iLayer)+dx))/(mLayerDepth(iLayer)*0.5_dp)
    dFlux_dTempAbove(iLayer) = (flux1 - flux0)/dx
   endif

  ! ***** internal layers
  else

   if(dMethod==analytical)then
    dFlux_dTempAbove(iLayer) =  iLayerThermalC(iLayer)/(mLayerHeight(iLayer+1) - mLayerHeight(iLayer))
    dFlux_dTempBelow(iLayer) = -iLayerThermalC(iLayer)/(mLayerHeight(iLayer+1) - mLayerHeight(iLayer))
   else
    dz    = (mLayerHeight(iLayer+1) - mLayerHeight(iLayer))
    flux0 = -iLayerThermalC(iLayer)*( mLayerTempIter(iLayer+1)     -  mLayerTempIter(iLayer)    ) / dz
    flux1 = -iLayerThermalC(iLayer)*( mLayerTempIter(iLayer+1)     - (mLayerTempIter(iLayer)+dx)) / dz
    flux2 = -iLayerThermalC(iLayer)*((mLayerTempIter(iLayer+1)+dx) -  mLayerTempIter(iLayer)    ) / dz
    dFlux_dTempAbove(iLayer) = (flux1 - flux0)/dx
    dFlux_dTempBelow(iLayer) = (flux2 - flux0)/dx
   endif  ! case for the numerical method

  endif  ! type of layer (upper, internal, or lower)

  ! print results
  !if(iLayer>45) &
  ! write(*,'(2(i4,1x),2(e20.10,1x))') dMethod, iLayer, dFlux_dTempAbove(iLayer), dFlux_dTempBelow(iLayer)

 end do  ! looping through layers

 end subroutine dFlux_dTemp


 ! ************************************************************************************************
 ! private subroutine: assemble the tri-diagonal matrix
 ! ************************************************************************************************
 subroutine getTridiag(dt,                  & ! intent(in): time step (s)
                       mLayerTempIter,      & ! intent(in): trial temperature (K)
                       mLayerVolFracIceIter,& ! intent(in): trial volumetric fraction of ice (-)
                       mLayerVolFracLiqIter,& ! intent(in): trial volumetric fraction of liquid water (-)
                       d_m1,                & ! intent(out): sub-diagonal elements of the tridiagonal system (J m-3 K-1)
                       diag,                & ! intent(out): diagonal elements of the tridiagonal system (J m-3 K-1)
                       d_p1,                & ! intent(out): super-diagonal elements of the tridiagonal system (J m-3 K-1)
                       err,message)           ! intent(out): error control
 ! physical constants
 USE multiconst,only:&
                     Cp_air,   & ! specific heat of air        (J kg-1 K-1)
                     LH_fus,   & ! latent heat of fusion       (J kg-1)
                     Em_Sno,   & ! emissivity of snow          (-)
                     sigma,    & ! Stefan Boltzman constant    (W m-2 K-4)
                     iden_air, & ! intrinsic density of air    (kg m-3)
                     iden_ice, & ! intrinsic density of ice    (kg m-3)
                     iden_water  ! intrinsic density of water  (kg m-3)
 ! utility modules
 USE snow_utils_module,only:dFracLiq_dTk  ! differentiate the freezing curve w.r.t. temperature (snow)
 USE soil_utils_module,only:dTheta_dTk    ! differentiate the freezing curve w.r.t. temperature (soil)
 ! names variables to define the layer type
 USE data_struc,only:ix_soil,ix_snow
 ! assemble the tridiagonal matrix
 implicit none
 ! input variables
 real(dp),intent(in)           :: dt                       ! time step (s)
 real(dp),intent(in)           :: mLayerTempIter(:)        ! trial temperature (K)
 real(dp),intent(in)           :: mLayerVolFracIceIter(:)  ! trial volumetric fraction of ice (-)
 real(dp),intent(in)           :: mLayerVolFracLiqIter(:)  ! trial volumetric fraction of liquid water (-)
 ! output variables
 real(dp),intent(out)          :: d_m1(:)                  ! sub-diagonal elements of the tridiagonal system (J m-3 K-1)
 real(dp),intent(out)          :: diag(:)                  ! diagonal elements of the tridiagonal system (J m-3 K-1)
 real(dp),intent(out)          :: d_p1(:)                  ! super-diagonal elements of the tridiagonal system (J m-3 K-1)
 integer(i4b),intent(out)      :: err                      ! error code
 character(*),intent(out)      :: message                  ! error message
 ! define local variables
 real(dp)                      :: theta                    ! volumetric liquid wareter equivalent of total water [liquid + ice] (-)
 real(dp),parameter            :: RHsurf=1._dp             ! relative humidity of the surface (-)
 real(dp)                      :: Qderiv                   ! derivative in specific humidity w.r.t. temperature (kg kg-1 K-1)
 real(dp)                      :: LW_deriv                 ! derivative in longwave radiation w.r.t. temperature (J m-2 s-1 K-1)
 real(dp)                      :: Qh_deriv                 ! derivative in sensible heat w.r.t. temperature (J m-2 s-1 K-1)
 real(dp)                      :: Qe_deriv                 ! derivative in latent heat w.r.t. temperature (J m-2 s-1 K-1)
 real(dp)                      :: dz_node(nLayers)         ! distance between the mid-point of a layer and a neighbouring layer (m)
 real(dp)                      :: dt_dudz(nLayers)         ! the dt_dudz terms for the upper interface (s m-2)
 real(dp)                      :: dt_dldz(nLayers)         ! the dt_dudz terms for the lower interface (s m-2)
 real(dp)                      :: hfusion(nLayers)         ! the fusion term (J m-3 K-1)
 integer(i4b)                  :: iLayer                   ! index of model layer
 ! initialize error control
 err=0; message="getTridiag/"

 ! compute the dt/dzdz terms for each layer (s m-2)
 dz_node(1:nLayers-1) = mLayerHeight(2:nLayers) - mLayerHeight(1:nLayers-1)
 dt_dudz(2:nLayers)   = (1._dp - wimplicit)*dt/(mLayerDepth(2:nLayers)*dz_node)    ! upper distance, valid 2..nlayers
 dt_dldz(1:nLayers-1) = (1._dp - wimplicit)*dt/(mLayerDepth(1:nLayers-1)*dz_node)  ! lower distance, valid 1..nlayers-1

 ! compute the dt/dzdz terms at the boundaries
 dt_dudz(1)       = (1._dp - wimplicit)*dt/(mLayerDepth(1)       * mLayerDepth(1)*0.5_dp)
 dt_dldz(nLayers) = (1._dp - wimplicit)*dt/(mLayerDepth(nLayers) * mLayerDepth(nLayers)*0.5_dp)

 ! compute the derivative in specific humidity w.r.t. temperature (kg kg-1 K-1)
 Qderiv = dsphum_dTk(RHsurf,airpres,mLayerTempIter(1)) 

 ! compute the derivative in sensible and latent heat (J m-2 s-1 K-1)
 Qh_deriv = scalarExSen
 Qe_deriv = scalarExLat*Qderiv

 ! compute the derivative in longwave radiation (J m-2 s-1 K-1)
 LW_deriv = 4._dp*Em_Sno*sigma*mLayerTempIter(1)**3._dp

 ! compute the fusion term for each layer (J m-3 K-1)
 hfusion = LH_fus*iden_water*mLayerdTheta_dTk

 ! assemble the tri-diagonal system
 do iLayer=1,nLayers
  ! compute the off-diagonal elements (J m-3 K-1)
  if(iLayer<nLayers) d_p1(iLayer)   = -dt_dldz(iLayer)*iLayerThermalC(iLayer)
  if(iLayer>1)       d_m1(iLayer-1) = -dt_dudz(iLayer)*iLayerThermalC(iLayer-1)
  ! compute the diagonal elements (J m-3 K-1)
  ! ** lower-most layer
  if(iLayer==nLayers)then
   diag(iLayer) = mLayerVolHtCapBulk(iLayer) + hfusion(iLayer) + &
                   (dt_dldz(iLayer)*mLayerThermalC(iLayer) +       &    ! NOTE: use mLayerThermalC
                    dt_dudz(iLayer)*iLayerThermalC(iLayer-1))
  ! ** upper-most layer
  elseif(iLayer==1)then
   diag(iLayer) = mLayerVolHtCapBulk(iLayer) + hfusion(iLayer) + &
                   dt_dldz(iLayer)*iLayerThermalC(iLayer) +       &
                  (1._dp - wimplicit)*dt*(LW_deriv + Qh_deriv + Qe_deriv)/mLayerDepth(iLayer)
  ! ** intermediate layers
  else
   diag(iLayer) = mLayerVolHtCapBulk(iLayer) + hfusion(iLayer) + &
                  dt_dldz(iLayer)*iLayerThermalC(iLayer) +       &
                  dt_dudz(iLayer)*iLayerThermalC(iLayer-1)
  endif  ! computing diagonal elements
 end do  ! (assembling the tri-diagonal system)

 end subroutine getTridiag


 ! ************************************************************************************************
 ! private subroutine: line search
 ! ************************************************************************************************
 SUBROUTINE lnsrch(dt,               & ! intent(in): time step (seconds)
            xold,                    & ! intent(in): trial temperature vector (K)
            stpmax,                  & ! intent(in): maximum step size (K)
            fold,                    & ! intent(in): function value for trial temperature vector (J m-3 J m-3)
            g,                       & ! intent(in): gradient of the function vector (J m-3 J m-3 K-1)
            p,                       & ! intent(in): iteration increment (K)
            mLayerMatricHeadIter,    & ! intent(in): matric head at the current iteration (m)
            mLayerVolFracLiqIter,    & ! intent(in): volumetric fraction of liquid water at the current iteration (-)
            mLayerVolFracIceIter,    & ! intent(in): volumetric fraction of ice at the current iteration (-)
            x,                       & ! intent(out): new temperature vector (K)
            rvec,                    & ! intent(out): new residual vector (J m-3)
            f,                       & ! intent(out): new function value (J m-3 J m-3)
            mLayerMatricHeadNew,     & ! intent(out): new matric head (m)
            mLayerVolFracLiqNew,     & ! intent(out): new volumetric fraction of liquid water (-)
            mLayerVolFracIceNew,     & ! intent(out): new volumetric fraction of ice (-)
            err,message)               ! intent(out): error control
 IMPLICIT NONE
 ! input variables
 real(dp) :: dt
 REAL(DP), DIMENSION(:), INTENT(IN) :: xold,g
 REAL(DP), DIMENSION(:), INTENT(INOUT) :: p
 REAL(DP), INTENT(IN) :: stpmax,fold
 real(dp),intent(in)           :: mLayerMatricHeadIter(:)  ! before phase change: matric head (m)
 real(dp),intent(in)           :: mLayerVolFracLiqIter(:)  ! before phase change: volumetric fraction of liquid water (-)
 real(dp),intent(in)           :: mLayerVolFracIceIter(:)  ! before phase change: volumetric fraction of ice (-)
 ! output variables
 REAL(DP), DIMENSION(:), INTENT(OUT) :: x,rvec
 REAL(DP), INTENT(OUT) :: f
 real(dp),intent(out)          :: mLayerMatricHeadNew(:)   ! after phase change: matric head (m)
 real(dp),intent(out)          :: mLayerVolFracLiqNew(:)   ! after phase change: volumetric fraction of liquid water (-)
 real(dp),intent(out)          :: mLayerVolFracIceNew(:)   ! after phase change: volumetric fraction of ice (-)
 integer(i4b),intent(out)  :: err         ! error code
 character(*),intent(out)  :: message     ! error message
 ! local variables
 character(LEN=256)            :: cmessage                 ! error message of downwind routine
 REAL(DP), PARAMETER :: ALF=1.0e-4_dp,TOLX=epsilon(x)
 INTEGER(I4B) :: ndum
 REAL(DP) :: a,alam,alam2,alamin,b,disc,f2,fold2,pabs,rhs1,rhs2,slope,&
     tmplam
 ! initialize error control
 err=0; message="lnsrch/"
 ! check arguments
 if ( all((/size(g),size(p),size(x)/) == size(xold)) ) then
  ndum=size(xold)
 else
  err=20; message=trim(message)//"sizeMismatch"; return
 endif
 ! start procedure
 pabs=sqrt(dot_product(p,p))
 if (pabs > stpmax) p(:)=p(:)*stpmax/pabs
 slope=dot_product(g,p)
 alamin=TOLX/maxval(abs(p(:))/max(abs(xold(:)),1.0_dp))
 alam=1.0_dp
 do
  ! update the temperature vector (K)
  x(:)=xold(:)+alam*p(:)
  ! compute the aerodynamic and canopy conductance
  call flxConduct(x(1),                & ! intent(in): surface temperature (K)
                  mLayerMatricHeadIter,& ! intent(in): matric head in each layer (m)
                  err,cmessage)          ! intent(out): error control
  if(err/=0)then; message=trim(message)//trim(cmessage); return; endif
  ! compute phase change
  call phsechange(x,                   & ! intent(in): new temperature vector (K)
                  mLayerMatricHeadIter,& ! intent(in): matric head at the current iteration (m)
                  mLayerVolFracLiqIter,& ! intent(in): volumetric fraction of liquid water at the current iteration (-)
                  mLayerVolFracIceIter,& ! intent(in): volumetric fraction of ice at the current iteration (-)
                  mLayerMatricHeadNew, & ! intent(out): new matric head (m)
                  mLayerVolFracLiqNew, & ! intent(out): new volumetric fraction of liquid water (-)
                  mLayerVolFracIceNew, & ! intent(out): new volumetric fraction of ice (-)
                  err,cmessage)          ! intent(out): error control
  if(err/=0)then; message=trim(message)//trim(cmessage); return; endif
  ! compute the residual vector (J m-3)
  print*, 'call from lnsrch'
  call e_residual(dt,                  & ! intent(in): time step (seconds)
                  x,                   & ! intent(in): trial temperature (K)
                  mLayerVolFracIceNew, & ! intent(in): trial volumetric fraction of ice (-)
                  rvec,                & ! intent(out): residual vector (J m-3)
                  err,cmessage)          ! intent(out): error control
  if(err/=0)then; message=trim(message)//trim(cmessage); return; endif
  ! compute the function evaluation
  f=0.5_dp*dot_product(rvec,rvec)
  ! additional exit criteria
  if(f<0.5_dp)return
  ! check if backtracked all the way to the original value
  if (alam < alamin) then
   x(:)=xold(:)
   err=-10; message=trim(message)//'warning: check convergence'
   RETURN
  ! check if improved the solution sufficiently
  else if (f <= fold+ALF*alam*slope) then
   RETURN
  ! build another trial vector
  else
   if (alam == 1.0_dp) then
    tmplam=-slope/(2.0_dp*(f-fold-slope))
   else
    rhs1=f-fold-alam*slope
    rhs2=f2-fold2-alam2*slope
    a=(rhs1/alam**2-rhs2/alam2**2)/(alam-alam2)
    b=(-alam2*rhs1/alam**2+alam*rhs2/alam2**2)/&
        (alam-alam2)
    if (a == 0.0_dp) then
     tmplam=-slope/(2.0_dp*b)
    else
     disc=b*b-3.0_dp*a*slope
     if (disc < 0.0_dp)then; err=-10; message=trim(message)//'warning: roundoff problem in lnsrch'; return; endif
     tmplam=(-b+sqrt(disc))/(3.0_dp*a)
    end if
    if (tmplam > 0.5_dp*alam) tmplam=0.5_dp*alam
   end if
  end if
  alam2=alam
  f2=f
  fold2=fold
  alam=max(tmplam,0.1_dp*alam)
 end do
 END SUBROUTINE lnsrch



 ! ************************************************************************************************
 ! private subroutine: assign variables to elements in the data structures
 ! ************************************************************************************************
 subroutine assignVars(err,message)
 ! data structures
 USE data_struc,only:mpar_data,forc_data,mvar_data,indx_data,ix_soil,ix_snow    ! data structures
 USE var_lookup,only:iLookPARAM,iLookFORCE,iLookMVAR,iLookINDEX                 ! named variables for structure elements
 ! assign variables to elements in the data structures
 implicit none
 ! dummy variables
 integer(i4b),intent(out)      :: err                      ! error code
 character(*),intent(out)      :: message                  ! error message
 ! initialize error control
 err=0; message="assignVars/"
 ! assign pointers to model parameters
 mheight             => mpar_data%var(iLookPARAM%mheight)                    ! measurement height (m)
 wimplicit           => mpar_data%var(iLookPARAM%wimplicit)                  ! weight assigned to start-of-step fluxes (-)
 snowfrz_scale       => mpar_data%var(iLookPARAM%snowfrz_scale)              ! scaling parameter for the snow freezing curve (K-1)
 lowerBoundTemp      => mpar_data%var(iLookPARAM%lowerBoundTemp)             ! temperature of the lower boundary (K)
 ! assign pointers to soil parameters
 vGn_alpha           => mpar_data%var(iLookPARAM%vGn_alpha)                  ! van Genutchen "alpha" parameter (m-1)
 vGn_n               => mpar_data%var(iLookPARAM%vGn_n)                      ! van Genutchen "n" parameter (-)
 theta_sat           => mpar_data%var(iLookPARAM%theta_sat)                  ! soil porosity (-)
 theta_res           => mpar_data%var(iLookPARAM%theta_res)                  ! soil residual volumetric water content (-)
 ! assign pointers to vegetation parameters
 LAI                 => mpar_data%var(iLookPARAM%LAI)                        ! leaf area index (m2 m-2)
 minStomatalResist   => mpar_data%var(iLookPARAM%minStomatalResist)          ! minimum stomatal resistance (s m-1)
 plantWiltPsi        => mpar_data%var(iLookPARAM%plantWiltPsi)               ! critical matric head when stomatal resitance 2 x min (m)
 plantWiltExp        => mpar_data%var(iLookPARAM%plantWiltExp)               ! empirical exponent in plant wilting factor expression (-)
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
 mLayerMatricHead    => mvar_data%var(iLookMVAR%mLayerMatricHead)%dat        ! matric head in each layer (m)
 mLayerVolFracIce    => mvar_data%var(iLookMVAR%mLayerVolFracIce)%dat        ! volumetric fraction of ice in each layer (-)
 mLayerDepth         => mvar_data%var(iLookMVAR%mLayerDepth)%dat             ! depth of each layer (m)
 ! assign pointers to model diagnostic variables
 mLayerHeight        => mvar_data%var(iLookMVAR%mLayerHeight)%dat            ! height at the mid-point of each layer (m)
 iLayerHeight        => mvar_data%var(iLookMVAR%iLayerHeight)%dat            ! height at the interface of each layer (m)
 mLayerThermalC      => mvar_data%var(iLookMVAR%mLayerThermalC)%dat          ! thermal conductivity at the mid-point of each layer (W m-1 K-1)
 iLayerThermalC      => mvar_data%var(iLookMVAR%iLayerThermalC)%dat          ! thermal conductivity at the interface of each layer (W m-1 K-1)
 mLayerVolHtCapBulk  => mvar_data%var(iLookMVAR%mLayerVolHtCapBulk)%dat      ! bulk volumetric heat capacity (J m-3 K-1)
 mLayerTcrit         => mvar_data%var(iLookMVAR%mLayerTcrit)%dat             ! critical soil temperature above which all water is unfrozen (K)
 mLayerdTheta_dTk    => mvar_data%var(iLookMVAR%mLayerdTheta_dTk)%dat        ! derivative in the freezing curve (K-1)
 mLayerTranspireLim  => mvar_data%var(iLookMVAR%mLayerTranspireLim)%dat      ! soil moist & veg limit on transpiration for each layer (-) 
 mLayerInitTranspire => mvar_data%var(iLookMVAR%mLayerInitTranspire)%dat      ! transpiration loss from each soil layer at the start of the step (m s-1)
 mLayerTranspire     => mvar_data%var(iLookMVAR%mLayerTranspire)%dat         ! transpiration loss from each soil layer (m s-1)
 iLayerInitNrgFlux   => mvar_data%var(iLookMVAR%iLayerInitNrgFlux)%dat       ! energy flux at layer interfaces at the start of the time step (W m-2)
 ! assign pointers to diagnostic scalar variables
 scalarTranspireLim  => mvar_data%var(iLookMVAR%scalarTranspireLim)%dat(1)   ! aggregate soil moist & veg limit on transpiration, weighted by root density (-)
 scalarPotentialET   => mvar_data%var(iLookMVAR%scalarPotentialET)%dat(1)    ! potential ET (kg m-2 s-1)
 scalarMassLiquid    => mvar_data%var(iLookMVAR%scalarMassLiquid)%dat(1)     ! transpiration (kg m-2 s-1)
 scalarMassSolid     => mvar_data%var(iLookMVAR%scalarMassSolid)%dat(1)      ! sublimation/frost (kg m-2 s-1)
 scalarSenHeat       => mvar_data%var(iLookMVAR%scalarSenHeat)%dat(1)        ! sensible heat flux at the surface (W m-2)
 scalarLatHeat       => mvar_data%var(iLookMVAR%scalarLatHeat)%dat(1)        ! latent heat flux at the surface (W m-2)
 scalarExcoef        => mvar_data%var(iLookMVAR%scalarExCoef)%dat(1)         ! turbulent exchange coefficient (-)
 scalarExSen         => mvar_data%var(iLookMVAR%scalarExSen)%dat(1)          ! exchange factor for sensible heat (J m-2 s-1 K-1)
 scalarExLat         => mvar_data%var(iLookMVAR%scalarExLat)%dat(1)          ! exchange factor for latent heat (J m-2 s-1)
 ! assign pointers to index variables
 nLayers             => indx_data%var(iLookINDEX%nLayers)%dat(1)             ! number of layers
 layerType           => indx_data%var(iLookINDEX%layerType)%dat              ! layer type (ix_soil or ix_snow)
 ! identify the number of snow and soil layers
 nSnow = count(layerType==ix_snow)
 nSoil = count(layerType==ix_soil)
 end subroutine assignVars


 ! ************************************************************************************************
 ! private function: compute derivative in specific humidity w.r.t. temperature (kg kg-1 K-1)
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

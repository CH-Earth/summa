module liquidflux_module
USE nrtype
implicit none
private
public::liquidflow
public::masschange
contains

 ! ************************************************************************************************
 ! new subroutine: compute liquid water flux through the snowpack
 ! ************************************************************************************************
 subroutine liquidflow(dt,&                       ! time step (seconds)
                       mLayerVolFracLiqIter,    & ! volumetric fraction of liquid water at the current iteration (-)
                       mLayerVolFracIceIter,    & ! volumetric fraction of ice at the current iteration (-)
                       mLayerVolFracLiqNew,     & ! volumetric fraction of liquid water at the next iteration (-)
                       err,message)
 USE multiconst,only:iden_ice,iden_water                                        ! intrinsic density of ice and water (kg m-3)
 USE data_struc,only:mpar_data,forc_data,mvar_data,indx_data,ix_soil,ix_snow    ! data structures
 USE var_lookup,only:iLookPARAM,iLookFORCE,iLookMVAR,iLookINDEX,iLookDECISIONS  ! named variables for structure elements
 ! compute energy fluxes within the domain
 implicit none
 ! dummy variables
 real(dp),intent(in)           :: dt                         ! time step (seconds)
 real(dp),intent(in)           :: mLayerVolFracLiqIter(:)    ! volumetric fraction of liquid water at the current iteration (-)
 real(dp),intent(in)           :: mLayerVolFracIceIter(:)    ! volumetric fraction of ice at the current iteration (-)
 real(dp),intent(out)          :: mLayerVolFracLiqNew(:)     ! volumetric fraction of liquid water at the next iteration (-)
 integer(i4b),intent(out)      :: err                        ! error code
 character(*),intent(out)      :: message                    ! error message
 ! local pointers to model parameters
 real(dp),pointer              :: Fcapil                     ! capillary retention as a fraction of the total pore volume (-)
 real(dp),pointer              :: k_snow                     ! hydraulic conductivity of snow (m s-1), 0.0055 = approx. 20 m/hr, from UEB
 real(dp),pointer              :: mw_exp                     ! exponent for meltwater flow (-)
 ! local pointers to model forcing data
 real(dp),pointer              :: rainfall                   ! rainfall (kg m-2 s-1) 
 ! local pointers to model state variables 
 real(dp),pointer              :: mLayerDepth(:)             ! depth of the layer (m)
 real(dp),pointer              :: mLayerVolFracLiq(:)        ! volumetric fraction of liquid water in each snow layer (-)
 real(dp),pointer              :: mLayerVolFracIce(:)        ! volumetric fraction of ice in each snow layer (-)
 ! local pointers to model diagnostic variables
 real(dp),pointer              :: iLayerLiquidWaterFlux(:)   ! vertical liquid water flux at layer interfaces (m s-1)
 ! local pointers to model index variables
 integer(i4b),pointer          :: layerType(:)               ! type of the layer (ix_soil or ix_snow)
 ! define local variables
 integer(i4b)                  :: iLayer                     ! layer index
 integer(i4b)                  :: nSnow                      ! number of snow layers
 real(dp)                      :: xmin,xmax                  ! bounds for bi-section (-)
 real(dp)                      :: volFracLiqTrial            ! trial value of volumetric fraction of liquid water (-)
 real(dp)                      :: volFracLiqNew              ! new value of volumetric fraction of liquid water (-)
 real(dp)                      :: poreSpace                  ! pore space (-)
 real(dp)                      :: volLiqRes                  ! residual volumetric liquid water content (-)
 real(dp)                      :: relPSpace                  ! relative pore space; min=res, max=porespace (-)
 real(dp)                      :: phseChnge                  ! volumetric liquid water equivalent associated with phase change (-)
 real(dp)                      :: relSaturn                  ! relative saturation [0,1] (-)
 real(dp)                      :: vDrainage                  ! vertical drainage (m s-1)
 real(dp)                      :: dflw_dliq                  ! derivative in vertical drainage (m s-1)
 real(dp)                      :: dt_dz                      ! dt/dz (s m-1)
 real(dp)                      :: lin_error                  ! linearization error [-residual] (-)
 real(dp)                      :: increment                  ! iteration increment (-)
 real(dp),parameter            :: atol=1.d-4                 ! absolute iteration tolerance (-)
 integer(i4b),parameter        :: maxiter=10                 ! maximum number of iterations
 integer(i4b)                  :: iter                       ! iteration index
 logical(lgt)                  :: printflag                  ! flag to print crap to the screen
 ! initialize error control
 err=0; message="liquidflow/"

 ! initialize printflag
 printflag=.false.
 
 ! assign local pointers to the model index structures
 layerType => indx_data%var(iLookINDEX%layerType)%dat              ! layer type (ix_soil or ix_snow)

 ! get the number of snow layers
 nSnow = count(layerType==ix_snow)

 ! check that the input vectors match nSnow
 if(size(mLayerVolFracLiqIter)/=nSnow .or. size(mLayerVolFracLiqNew)/=nSnow) then
  err=20; message=trim(message)//'size mismatch of input vectors'; return
 endif

 ! assign pointers to model parameters
 Fcapil        => mpar_data%var(iLookPARAM%Fcapil)           ! capillary retention as a fraction of the total pore volume (-)
 k_snow        => mpar_data%var(iLookPARAM%k_snow)           ! hydraulic conductivity of snow (m s-1), 0.0055 = approx. 20 m/hr, from UEB
 mw_exp        => mpar_data%var(iLookPARAM%mw_exp)           ! exponent for meltwater flow (-)

 ! assign pointers to model forcing data
 rainfall => mvar_data%var(iLookMVAR%scalarRainfall)%dat(1)  ! computed rainfall rate (kg m-2 s-1)

 ! assign pointers to model state variables
 mLayerDepth           => mvar_data%var(iLookMVAR%mLayerDepth)%dat(1:nSnow)         ! depth of the layer (m)
 mLayerVolFracLiq      => mvar_data%var(iLookMVAR%mLayerVolFracLiq)%dat(1:nSnow)    ! volumetric fraction of liquid water in each snow layer (-)
 mLayerVolFracIce      => mvar_data%var(iLookMVAR%mLayerVolFracIce)%dat(1:nSnow)    ! volumetric fraction of ice in each snow layer (-)

 ! assign pointers to model diagnostic variables
 iLayerLiquidWaterFlux => mvar_data%var(iLookMVAR%iLayerLiquidWaterFlux)%dat        ! vertical liquid water flux at layer interfaces (m s-1)

 ! define the liquid flux at the upper boundary -- include evaporation/dew (m s-1)
 iLayerLiquidWaterFlux(0) = rainfall/iden_water

 ! check the meltwater exponent is >=1
 if(mw_exp<1._dp)then; err=20; message=trim(message)//'meltwater exponent < 1'; return; endif

 do iLayer=1,nSnow
  ! initialize the volumetric fraction of liquid water (-)
  volFracLiqTrial = mLayerVolFracLiqIter(iLayer)
  ! compute dt/dz
  dt_dz     = dt/mLayerDepth(iLayer)
  ! compute the pore space (-)
  poreSpace = 1._dp - mLayerVolFracIceIter(iLayer)
  ! compute the residual volumetric liquid water content
  volLiqRes = Fcapil*poreSpace
  ! compute the relative pore space
  relPSpace = poreSpace - volLiqRes
  ! compute the liquid water equivalent associated with phase change
  phseChnge = (iden_ice/iden_water)*(mLayerVolFracIceIter(iLayer) - mLayerVolFracIce(iLayer))
  ! initialize bounds
  xmin = 0._dp
  xmax = poreSpace
  ! ****************************
  ! ***** begin iterations *****
  ! ****************************
  do iter=1,maxiter
   ! ***** check that flow occurs
   if(volFracLiqTrial > volLiqRes)then
    ! compute the relative saturation (-)
    relSaturn = (volFracLiqTrial - volLiqRes)/relPSpace
    ! compute the initial flow estimate (m s-1)
    vDrainage = k_snow*relSaturn**mw_exp
    ! compute the change in the vertical liquid water flux with volumetric liquid water content (m s-1)
    dflw_dliq = ( (k_snow*mw_exp)/relPSpace ) * relSaturn**(mw_exp - 1._dp)
   ! **** default = flow does not occur, and the derivative is zero
   else
    vDrainage = 0._dp
    dflw_dliq = 0._dp
   endif
   ! compute the residual (-)
   lin_error = dt_dz*(iLayerLiquidWaterFlux(iLayer-1) - vDrainage) - phseChnge - &
                (volFracLiqTrial - mLayerVolFracLiq(iLayer))
   ! compute the iteration increment (-) and new value
   increment = lin_error/(1._dp + dt_dz*dflw_dliq)
   volFracLiqNew = volFracLiqTrial + increment
   ! update bounds
   if(increment> 0._dp) xmin = volFracLiqTrial
   if(increment<=0._dp) xmax = volFracLiqTrial
   ! use bi-section if outside bounds
   if(volFracLiqNew<xmin .or. volFracLiqNew>xmax) then
    volFracLiqNew = 0.5_dp*(xmin+xmax)
    increment = volFracLiqNew - volFracLiqTrial
   endif
   if(printflag)&
    write(*,'(a,2(i4,1x),10(f20.10,1x))') 'in watersnow, iter, iLayer, volFracLiqNew, volLiqRes = ', iter, iLayer, volFracLiqNew, volLiqRes
   ! check for convergence
   if(abs(increment) < atol) exit
   ! get ready for next iteration
   volFracLiqTrial = volFracLiqNew
   ! check for non-convergence
   if(iter==maxiter)then; err=20; message=trim(message)//'failed to converge'; return; endif
  end do  ! (end iterations)
  ! save state
  mLayerVolFracLiqNew(iLayer) = volFracLiqNew
  ! get ready to process the next snow layer
  iLayerLiquidWaterFlux(iLayer) = vDrainage + dflw_dliq*increment  ! second term will be zero if converge completely
  ! *** now process the next layer
 end do  ! (looping through snow layers)

 end subroutine liquidflow

 ! ************************************************************************************************
 ! new subroutine: compute liquid water flux through the soil
 ! ************************************************************************************************
 subroutine masschange(dt,&                   ! time step (seconds)
                       iter,&
                       mLayerMatricHeadIter,& ! matric head in each layer at the current iteration (m)
                       mLayerVolFracIceIter,& ! volumetric fraction of ice at the current iteration (-)
                       mLayerVolFracLiqIter,& ! volumetric fraction of liquid water at the current iteration (-)
                       mLayerMatricHeadDiff,& ! iteration increment for matric head (m)
                       mLayerMatricHeadNew, & ! matric head in each layer at the next iteration (m)
                       mLayerVolFracLiqNew,& ! volumetric fraction of liquid water at the next iteration (-)
                       err,message)
 USE multiconst,only:iden_ice,iden_water                                        ! intrinsic density of ice and water (kg m-3)
 USE data_struc,only:model_decisions                                            ! model decision structure
 USE data_struc,only:mpar_data,forc_data,mvar_data,indx_data,ix_soil,ix_snow    ! data structures
 USE var_lookup,only:iLookPARAM,iLookFORCE,iLookMVAR,iLookINDEX,iLookDECISIONS  ! named variables for structure elements
 USE soil_utils_module,only:volFracLiq      ! compute volumetric fraction of liquid water
 USE soil_utils_module,only:matricHead      ! compute matric head (m)
 USE soil_utils_module,only:dTheta_dPsi     ! compute derivative of the soil moisture characteristic
 USE tridagSolv_module,only:tridag          ! solve tridiagonal system of equations
 ! compute energy fluxes within the domain
 implicit none
 ! dummy variables
 real(dp),intent(in)           :: dt                       ! time step (seconds)
 integer(i4b),intent(in)       :: iter
 real(dp),intent(in)           :: mLayerMatricHeadIter(:)  ! matric head in each layer at the current iteration (m)
 real(dp),intent(in)           :: mLayerVolFracIceIter(:)  ! volumetric fraction of ice at the current iteration (-)
 real(dp),intent(in)           :: mLayerVolFracLiqIter(:)  ! volumetric fraction of liquid water at the current iteration (-)
 real(dp),intent(out)          :: mLayerMatricHeadDiff(:)  ! iteration increment for matric head (m)
 real(dp),intent(out)          :: mLayerMatricHeadNew(:)   ! matric head in each layer at the next iteration (m)
 real(dp),intent(out)          :: mLayerVolFracLiqNew(:)   ! volumetric fraction of liquid water at the next iteration (m)
 integer(i4b),intent(out)      :: err                      ! error code
 character(*),intent(out)      :: message                  ! error message
 ! local pointers to model parameters
 real(dp),pointer              :: vGn_alpha                ! van Genutchen "alpha" parameter
 real(dp),pointer              :: vGn_n                    ! van Genutchen "n" parameter
 real(dp),pointer              :: theta_sat                ! soil porosity (-)
 real(dp),pointer              :: theta_res                ! soil residual volumetric water content (-)
 real(dp),pointer              :: k_soil                   ! hydraulic conductivity (m s-1)
 real(dp),pointer              :: spec_storage             ! specific storage coefficient (m-1)
 real(dp),pointer              :: f_impede                 ! ice impedence factor (-)
 real(dp),pointer              :: upperBoundHead           ! upper boundary condition for matric head (m)
 real(dp),pointer              :: lowerBoundHead           ! lower boundary condition for matric head (m)
 ! local pointers to model forcing data
 real(dp),pointer              :: rainfall                 ! rainfall (kg m-2 s-1) 
 ! local pointers to model variables that are constant over the simulation period 
 real(dp),pointer              :: vGn_m                    ! van Genutchen "m" parameter (-)
 ! local pointers to model state variables
 real(dp),pointer              :: mLayerVolFracIce(:)      ! volumetric fraction of ice at the start of the time step (-)
 real(dp),pointer              :: mLayerVolFracLiq(:)      ! volumetric fraction of liquid water at the start of the time step (-)
 real(dp),pointer              :: mLayerMatricHead(:)      ! matric head in each layer at the start of the time step (m)
 ! local pointers to variable short-cuts
 real(dp),pointer              :: volLatHt_fus             ! volumetric latent heat of fusion (J m-3 K-1)
 ! local pointers to diagnostic scalar variables
 real(dp),pointer              :: scalarSfcMeltPond        ! ponded water caused by melt of the "snow without a layer" (kg m-2)
 ! local pointers to model diagnostic variables
 real(dp),pointer              :: mLayerDepth(:)           ! depth of the layer (m)
 real(dp),pointer              :: mLayerHeight(:)          ! height of the layer mid-point (m)
 real(dp),pointer              :: mLayerdTheta_dPsi(:)     ! derivative in the soil water characteristic (m-1)
 real(dp),pointer              :: mLayerTranspire(:)       ! transpiration loss from each soil layer (kg m-2 s-1)
 real(dp),pointer              :: surfaceRunoff            ! surface runoff (m s-1)
 ! local pointers to model index variables
 integer(i4b),pointer          :: layerType(:)             ! type of the layer (ix_soil or ix_snow)
 ! define local variables
 character(LEN=256)            :: cmessage                 ! error message of downwind routine
 integer(i4b)                  :: nSnow                    ! number of snow layers
 integer(i4b)                  :: nSoil                    ! number of soil layers
 integer(i4b)                  :: ibeg,iend                ! start and end indices of the soil layers in concatanated vector
 integer(i4b)                  :: iLayer                   ! layer index
 real(dp)                      :: q_top,q_bot              ! boundary fluxes (m s-1)
 real(dp),dimension(size(mLayerMatricHeadIter))    :: mLayerHydCond  ! hydraulic conductivity at layer mid-point (m s-1)
 real(dp),dimension(0:size(mLayerMatricHeadIter))  :: iLayerHydCond  ! hydraulic conductivity at layer interface (m s-1)
 real(dp),dimension(size(mLayerMatricHeadIter)-1)  :: dz_node        ! distance between the mid-point of a layer and a neighbouring layer (m)
 real(dp),dimension(size(mLayerMatricHeadIter))    :: dt_dudz        ! the dt_dudz terms for the upper interface (s m-2)
 real(dp),dimension(size(mLayerMatricHeadIter))    :: dt_dldz        ! the dt_dudz terms for the lower interface (s m-2)
 real(dp),dimension(size(mLayerMatricHeadIter)-1)  :: d_m1           ! sub-diagonal elements of the tridiagonal system (m-1)
 real(dp),dimension(size(mLayerMatricHeadIter))    :: diag           ! diagonal elements of the tridiagonal system (m-1)
 real(dp),dimension(size(mLayerMatricHeadIter)-1)  :: d_p1           ! super-diagonal elements of the tridiagonal system (m-1)
 real(dp),dimension(size(mLayerMatricHeadIter))    :: rvec           ! right-hand-side vector (-)
 real(dp),dimension(size(mLayerMatricHeadIter))    :: sink           ! sink term associated with transpiration (-)
 integer(i4b),parameter        :: diriclet=1001            ! look-up value for diriclet boundary conditions
 integer(i4b),parameter        :: neumann=1002             ! look-up value for neumann boundary conditions
 integer(i4b)                  :: bc_upper,bc_lower        ! boundary condition index (diriclet or neumann)
 real(dp)                      :: fracFrozen               ! fraction of frozen ground (-)
 real(dp),dimension(size(mLayerMatricHeadIter))   :: g     ! gradient of the function vector (m-1)
 real(dp)                      :: fold,fnew                ! function values (-)
 real(dp),parameter            :: STPMX=5._dp              ! maximum step size in line search (m)
 real(dp)                      :: stpmax                   ! scaled maximum step size
 logical(lgt)                  :: printflag                ! flag to print crap to the screen

 ! initialize error control
 err=0; message="masschange/"

 ! initilaize printflag
 printflag=.false.

 ! assign local pointers to the model index structures
 layerType => indx_data%var(iLookINDEX%layerType)%dat ! layer type (ix_soil or ix_snow)

 ! get the number of snow and soil layers
 nSnow = count(layerType==ix_snow)
 nSoil = count(layerType==ix_soil)

 ! check the size of the input arguments
 if(any((/size(mLayerMatricHeadIter),size(mLayerVolFracIceIter),size(mLayerVolFracLiqIter),size(mLayerMatricHeadNew)/) /= nSoil)) then
  err=20; message=trim(message)//'size mis-match for the input arguments'; return
 endif

 ! get indices for the data structures
 ibeg = count(layerType==ix_snow)+1
 iend = ibeg + (nSoil-1)

 ! assign pointers to model parameters
 vGn_alpha         => mpar_data%var(iLookPARAM%vGn_alpha)                      ! van Genutchen "alpha" parameter (m-1)
 vGn_n             => mpar_data%var(iLookPARAM%vGn_n)                          ! van Genutchen "n" parameter (-)
 theta_sat         => mpar_data%var(iLookPARAM%theta_sat)                      ! soil porosity (-)
 theta_res         => mpar_data%var(iLookPARAM%theta_res)                      ! soil residual volumetric water content (-)
 k_soil            => mpar_data%var(iLookPARAM%k_soil)                         ! hydraulic conductivity (m s-1)
 spec_storage      => mpar_data%var(iLookPARAM%spec_storage)                   ! specific storage coefficient (m-1)
 f_impede          => mpar_data%var(iLookPARAM%f_impede)                       ! ice impedence factor (-)

 ! assign pointers to head boundary conditions (included as parameters for convenience)
 lowerBoundHead    => mpar_data%var(iLookPARAM%lowerBoundHead)
 upperBoundHead    => mpar_data%var(iLookPARAM%upperBoundHead)

 ! assign pointers to model forcing data
 rainfall          => mvar_data%var(iLookMVAR%scalarRainfall)%dat(1)             ! computed rainfall rate (kg m-2 s-1)

 ! assign pointers to model variables that are constant over the simulation period
 vGn_m             => mvar_data%var(iLookMVAR%scalarVGn_m)%dat(1)              ! van Genutchen "m" parameter (-)
 volLatHt_fus      => mvar_data%var(iLookMVAR%scalarvolLatHt_fus)%dat(1)       ! volumetric latent heat of fusion (J m-3)

 ! assign pointers to model state variables
 mLayerVolFracIce  => mvar_data%var(iLookMVAR%mLayerVolFracIce)%dat(ibeg:iend) ! volumetric fraction of ice in each layer (-)
 mLayerVolFracLiq  => mvar_data%var(iLookMVAR%mLayerVolFracLiq)%dat(ibeg:iend) ! volumetric fraction of liquid water in each layer (-)
 mLayerMatricHead  => mvar_data%var(iLookMVAR%mLayerMatricHead)%dat(ibeg:iend) ! matric head in each layer (m)

 ! assign local pointers to diagnostic scalar variables
 scalarSfcMeltPond => mvar_data%var(iLookMVAR%scalarSfcMeltPond)%dat(1)        ! ponded water caused by melt of the "snow without a layer" (kg m-2)

 ! assign pointers to model diagnostic variables
 mLayerDepth       => mvar_data%var(iLookMVAR%mLayerDepth)%dat(ibeg:iend)      ! depth of the layer (m)
 mLayerHeight      => mvar_data%var(iLookMVAR%mLayerHeight)%dat(ibeg:iend)     ! height of the layer mid-point (m)
 mLayerdTheta_dPsi => mvar_data%var(iLookMVAR%mLayerdTheta_dPsi)%dat  ! (soil only) ! derivative in the soil water characteristic (m-1)
 mLayerTranspire   => mvar_data%var(iLookMVAR%mLayerTranspire)%dat    ! (soil only) ! transpiration loss from each soil layer (kg m-2 s-1)
 surfaceRunoff     => mvar_data%var(iLookMVAR%scalarSurfaceRunoff)%dat(1)      ! surface runoff (m s-1)

 ! identify the boundary conditions
 select case(trim(model_decisions(iLookDECISIONS%bound_cond)%decision))
  case('fluxflux'); bc_upper=neumann;  bc_lower=neumann
  case('fluxhead'); bc_upper=neumann;  bc_lower=diriclet
  case('headflux'); bc_upper=diriclet; bc_lower=neumann
  case('headhead'); bc_upper=diriclet; bc_lower=diriclet
  case default
   err=10; message=trim(message)//"unknownBoundaryConditionOption[option="//trim(model_decisions(iLookDECISIONS%bound_cond)%decision)//"]"; return
 end select

 ! check the b/c make sense
 if(nSnow>0 .and. bc_upper==diriclet)then
  err=20; message=trim(message)//'using diriclet bc for the top of the soil zone when snow is present'; return
 endif

 ! define boundary fluxes (m s-1)
 if(bc_upper==neumann) then
  if(nSnow==0) q_top = rainfall/iden_water !+ (scalarSfcMeltPond/dt)/iden_water
  if(nSnow>0)  q_top = mvar_data%var(iLookMVAR%iLayerLiquidWaterFlux)%dat(nSnow)  ! liquid water flux from the base of the snowpack (m s-1)
  !if(nSnow>0)  q_top = 0._dp 
  ! compute fraction of frozen ground -- note element 1 is top soil element, because just soil layers passed
  !fracFrozen = min(mLayerVolFracIceIter(1),theta_sat)/theta_sat
  ! compute surface runoff
  !surfaceRunoff = fracFrozen*q_top
  ! update boundary condition
  !q_top = surfaceRunoff - q_top
 endif
 if(bc_lower==neumann) q_bot = 0._dp


 ! compute the dt/dzdz terms for each layer(s m-2)
 dz_node(1:nSoil-1) = mLayerHeight(2:nSoil) - mLayerHeight(1:nSoil-1)
 dt_dudz(2:nSoil)   = dt/(mLayerDepth(2:nSoil)*dz_node)    ! upper distance, valid 2..nSoil
 dt_dldz(1:nSoil-1) = dt/(mLayerDepth(1:nSoil-1)*dz_node)  ! lower distance, valid 1..nSoil-1

 ! compute the dt/dzdz terms at the boundaries
 dt_dudz(1)     = dt/(mLayerDepth(1)     * mLayerDepth(1)*0.5_dp)
 dt_dldz(nSoil) = dt/(mLayerDepth(nSoil) * mLayerDepth(nSoil)*0.5_dp)

 ! compute the sink term
 if(nSnow==0)then
  sink(1:nSoil) = (mLayerTranspire(1:nSoil)/iden_water) *dt/mLayerDepth(1:nSoil)
 else
  sink(1:nSoil) = 0._dp
 endif
 !print*, 'sink = ', sink

 ! *****
 ! compute the hydraulic conductivity for all layers
 call hydCond_all(mLayerMatricHeadIter, & ! intent(in): matric head (m)
                  mLayerVolFracIceIter, & ! intent(in): volumetric fraction of ice (-)
                  mLayerVolFracLiqIter, & ! intent(in): volumetric fraction of liquid water (-)
                  err,cmessage)           ! intent(out): error control
 if(err/=0)then; message=trim(message)//trim(cmessage); return; endif

 ! *****
 ! compute the residual vector
 call vlFrcLiqRes(mLayerMatricHeadIter, & ! intent(in): matric head (m)
                  mLayerVolFracIceIter, & ! intent(in): volumetric fraction of ice (-)
                  mLayerVolFracLiqIter, & ! intent(in): volumetric fraction of liquid water (-)
                  rvec,                 & ! intent(out): residual vector (-)
                  err,cmessage)           ! intent(out): error control
 if(err/=0)then; message=trim(message)//trim(cmessage); return; endif
 if(printflag) print*, 'residuals -- take1 ', rvec
 pause

 ! *****
 ! compute the tri-diagonal matrix
 call get_tridiag(mLayerMatricHeadIter,  & ! intent(in): matric head (m)
                  d_m1,                  & ! intent(out): sub-diagonal vector
                  diag,                  & ! intent(out): diagonal vector
                  d_p1,                  & ! intent(out): super-diagonal vector
                  err,cmessage)            ! intent(out): error control
 if(err/=0)then; message=trim(message)//trim(cmessage); return; endif

 ! *****
 ! solve the tridiagonal system of equations -- returns mLayerMatricHeadNew
 call tridag(d_m1,diag,d_p1,rvec,mLayerMatricHeadDiff,err,cmessage)
 if(err/=0)then; message=trim(message)//trim(cmessage); return; endif

 ! test tridag
 if(printflag)then
  write(*,'(a)') 'test tridag: iLayer, d_m1(iLayer-1), diag(iLayer), d_p1(iLayer), mLayerMatricHeadDiff(iLayer), rvec(iLayer)'
  do iLayer=nSoil-5,nSoil
   if(iLayer==1) then
    write(*,'(i4,1x,10(f20.10,1x))') iLayer, diag(iLayer), d_p1(iLayer), mLayerMatricHeadDiff(iLayer), rvec(iLayer)
    print*, iLayer, diag(iLayer)*mLayerMatricHeadDiff(iLayer) + d_p1(iLayer)*mLayerMatricHeadDiff(iLayer+1), rvec(iLayer)
   endif
   if(iLayer==nSoil) then
    write(*,'(i4,1x,10(f20.10,1x))') iLayer, d_m1(iLayer-1), diag(iLayer), mLayerMatricHeadDiff(iLayer), rvec(iLayer)
    print*, iLayer, d_m1(iLayer-1)*mLayerMatricHeadDiff(iLayer-1) + diag(iLayer)*mLayerMatricHeadDiff(iLayer), rvec(iLayer)
    print*, iLayer, d_m1(iLayer-1)*mLayerMatricHeadDiff(iLayer-1), diag(iLayer)*mLayerMatricHeadDiff(iLayer)
   endif
   if(iLayer>1 .and. iLayer<nSoil) then
    write(*,'(i4,1x,10(f20.10,1x))') iLayer, d_m1(iLayer-1), diag(iLayer), d_p1(iLayer), mLayerMatricHeadDiff(iLayer), rvec(iLayer)
    print*, iLayer, d_m1(iLayer-1)*mLayerMatricHeadDiff(iLayer-1) + diag(iLayer)*mLayerMatricHeadDiff(iLayer) + d_p1(iLayer)*mLayerMatricHeadDiff(iLayer+1), rvec(iLayer)
   endif
  end do
 endif ! (if printflag)

 ! compute the gradient of the function vector (m-1)
 do iLayer=1,nSoil
  if(iLayer==1)then;         g(iLayer) =                                 diag(iLayer)*rvec(iLayer) + d_p1(iLayer)*rvec(iLayer+1)
  elseif(iLayer==nSoil)then; g(iLayer) = d_m1(iLayer-1)*rvec(iLayer-1) + diag(iLayer)*rvec(iLayer)
  else;                      g(iLayer) = d_m1(iLayer-1)*rvec(iLayer-1) + diag(iLayer)*rvec(iLayer) + d_p1(iLayer)*rvec(iLayer+1)
  endif
 end do

 ! compute the function value (-)
 fold = 0.5_dp*dot_product(rvec,rvec)

 ! compute maximum step size (K)
 stpmax=STPMX*real(nSoil,dp)

 ! *****
 ! compute line search
 call lnsrch(mLayerMatricHeadIter,    & ! intent(in): matric head at the current iteration (m)
             stpmax,                  & ! intent(in): maximum step size (m)
             fold,                    & ! intent(in): function value for trial matric head vector (-)
             g,                       & ! intent(in): gradient of the function vector (m-1)
             mLayerMatricHeadDiff,    & ! intent(in): iteration increment (m)
             mLayerVolFracLiqIter,    & ! intent(in): volumetric fraction of liquid water at the current iteration (-)
             mLayerVolFracIceIter,    & ! intent(in): volumetric fraction of ice at the current iteration (-)
             mLayerMatricHeadNew,     & ! intent(out): new matric head vector (K)
             rvec,                    & ! intent(out): new residual vector (-)
             fnew,                    & ! intent(out): new function value (-)
             mLayerVolFracLiqNew,     & ! intent(out): new volumetric fraction of liquid water (-)
             err,cmessage)               ! intent(out): error control
 if(err/=0)then; message=trim(message)//trim(cmessage); return; endif

 ! update the matric head increment
 mLayerMatricHeadDiff = mLayerMatricHeadNew - mLayerMatricHeadIter

 ! compute the boundary fluxes
 if(bc_upper==diriclet)&
  q_top = iLayerHydCond(0)*(upperBoundHead - mLayerMatricHeadNew(1))/(0.5_dp*mLayerDepth(1)) + iLayerHydCond(0)
 if(bc_lower==diriclet)&
  q_bot = iLayerHydCond(nSoil)*(mLayerMatricHeadNew(nSoil) - lowerBoundHead)/(0.5_dp*mLayerDepth(nSoil)) + iLayerHydCond(nSoil)
 !print*, 'q_top, q_bot = ', q_top, q_bot
 !print*, 'change in mass = ', q_top - q_bot


 contains

  ! ************************************************************************************************
  ! internal subroutine: compute hydraulic conductivity
  ! ************************************************************************************************
  subroutine hydCond_all(mLayerMatricHeadTrial, & ! intent(in): matric head (m)
                         mLayerVolFracIceTrial, & ! intent(in): volumetric fraction of ice (-)
                         mLayerVolFracLiqTrial, & ! intent(in): volumetric fraction of liquid water (-)
                         err,message)             ! intent(out): error control
  USE multiconst,only:iden_ice,iden_water    ! intrinsic density of ice and water
  USE soil_utils_module,only:hydCond         ! compute hydraulic conductivity
  ! compute hydraulic conductivity for all layers
  implicit none
  ! input
  real(dp),intent(in)           :: mLayerMatricHeadTrial(:)  ! matric head in each layer (m)
  real(dp),intent(in)           :: mLayerVolFracIceTrial(:)  ! volumetric fraction of ice in each layer (-)
  real(dp),intent(in)           :: mLayerVolFracLiqTrial(:)  ! volumetric fraction of liquid water in each layer (-)
  ! output
  integer(i4b),intent(out)      :: err                       ! error code
  character(*),intent(out)      :: message                   ! error message
  ! local variables
  character(LEN=256)            :: cmessage                  ! error message of downwind routine 
  real(dp),dimension(size(mLayerMatricHeadTrial))  :: relIce ! relative ice content in each layer soil layer (-)
  real(dp)                      :: iceImpede                 ! ice impedence factor (-)
  real(dp)                      :: test                      ! scalar for testing
  ! initialize error control
  err=0; message="hydCond_all/"
 
  ! loop through layers
  do iLayer=1,nSoil
   ! compute the relative volumetric fraction of ice
   relIce(iLayer) = (mLayerVolFracIceTrial(iLayer)*iden_ice) / ( mLayerVolFracIceTrial(iLayer)*iden_ice + mLayerVolFracLiqTrial(iLayer)*iden_water )
   ! compute the ice impedence factor (-)
   iceImpede = ( exp(-f_impede*(1._dp - relIce(iLayer))) - exp(-f_impede) ) / (1._dp - exp(-f_impede))
   ! compute the hydraulic conductivity for a given layer (m s-1)
   mLayerHydCond(iLayer) = hydCond(mLayerMatricHeadTrial(iLayer),k_soil,vGn_alpha,vGn_n,vGn_m) * (1._dp - iceImpede)
   ! compute the hydraulic conductivity for layer interfaces (m s-1)
   if(iLayer>1) iLayerHydCond(iLayer-1) = (mLayerHydCond(iLayer-1) + mLayerHydCond(iLayer))/2._dp
  end do

  ! if diriclet, compute hydraulic conductivity at the boundaries (m s-1)
  if(bc_upper==diriclet) iLayerHydCond(0)     = hydCond(upperBoundHead,k_soil,vGn_alpha,vGn_n,vGn_m) * 10._dp**(f_impede*relIce(1))
  if(bc_lower==diriclet) iLayerHydCond(nSoil) = hydCond(lowerBoundHead,k_soil,vGn_alpha,vGn_n,vGn_m) * 10._dp**(f_impede*relIce(nSoil))

  ! if newmann, set hydraulic conductivity at boundaries to boundary fluxes
  if(bc_upper==neumann) iLayerHydCond(0)     = q_top
  if(bc_lower==neumann) iLayerHydCond(nSoil) = q_bot

  end subroutine hydCond_all


  ! ************************************************************************************************
  ! internal subroutine: compute the residual vector
  ! ************************************************************************************************
  subroutine vlFrcLiqRes(mLayerMatricHeadTrial, & ! intent(in): matric head (m)
                         mLayerVolFracIceTrial, & ! intent(in): volumetric fraction of ice (-)
                         mLayerVolFracLiqTrial, & ! intent(in): volumetric fraction of liquid water (-)
                         rvec,                  & ! intent(out): residual vector (-)
                         err,message)             ! intent(out): error control
  ! compute residual vector for the volumetric fraction of liquid water
  implicit none
  ! input
  real(dp),intent(in)           :: mLayerMatricHeadTrial(:)  ! matric head in each layer (m)
  real(dp),intent(in)           :: mLayerVolFracIceTrial(:)  ! volumetric fraction of ice in each layer (-)
  real(dp),intent(in)           :: mLayerVolFracLiqTrial(:)  ! volumetric fraction of liquid water in each layer (-)
  ! output
  real(dp),intent(out)          :: rvec(:)                   ! residual vector (-)
  integer(i4b),intent(out)      :: err                       ! error code
  character(*),intent(out)      :: message                   ! error message
  ! local variables
  character(LEN=256)            :: cmessage                  ! error message of downwind routine 
  real(dp)                      :: matricDiffUpper           ! spatial difference in matric head (m)
  real(dp)                      :: matricDiffLower           ! spatial difference in matric head (m)
  ! initialize error control
  err=0; message="vlFrcLiqRes/"
  ! loop through soil layers
  do iLayer=1,nSoil
   ! compute the spatial differences in matric head at the bottom of the layer
   if(iLayer<nSoil)then
    matricDiffLower = mLayerMatricHeadTrial(iLayer+1) - mLayerMatricHeadTrial(iLayer)
   else  ! (special case of the bottom layer)
    select case(bc_lower)
     case(neumann);  matricDiffLower = 0._dp
     case(diriclet); matricDiffLower = lowerBoundHead - mLayerMatricHeadTrial(iLayer)
    endselect
   endif
   ! compute the spatial differences in matric head at the top of the layer
   if(iLayer>1)then
    matricDiffUpper = mLayerMatricHeadTrial(iLayer) - mLayerMatricHeadTrial(iLayer-1)
   else  ! (special case of the top layer)
    select case(bc_upper)
     case(neumann);  matricDiffUpper = 0._dp
     case(diriclet); matricDiffUpper = mLayerMatricHeadTrial(iLayer) - mLayerMatricHeadTrial(iLayer-1)
    endselect
   endif
   ! compute the residuals
   rvec(iLayer) = dt_dldz(iLayer)*(iLayerHydCond(iLayer)  *matricDiffLower) - &
                  dt_dudz(iLayer)*(iLayerHydCond(iLayer-1)*matricDiffUpper) - &
                  (dt/mLayerDepth(iLayer))*(iLayerHydCond(iLayer) - iLayerHydCond(iLayer-1)) -    &
                  (iden_ice/iden_water)*(mLayerVolFracIceTrial(iLayer)-mLayerVolFracIce(iLayer)) - &
                  (mLayerVolFracLiqTrial(iLayer)-mLayerVolFracLiq(iLayer)) + sink(iLayer)
   write(*,'(a,3(f15.12,1x))') 'mlayerMatricHeadTrial(iLayer), matricDiffLower, rvec(iLayer) = ', &
    mLayerMatricHeadTrial(iLayer), matricDiffLower, rvec(iLayer)
   if(printflag)then
    if(iLayer==nSoil)then
     print*, 'matricDiffLower, matricDiffUpper = ', matricDiffLower, matricDiffUpper
     print*, 'capil lower, capil upper = ', dt_dldz(iLayer)*(iLayerHydCond(iLayer)  *matricDiffLower), &
                                            dt_dudz(iLayer)*(iLayerHydCond(iLayer-1)*matricDiffUpper)
     print*, 'hydcon lower,upper = ', iLayerHydCond(iLayer), iLayerHydCond(iLayer-1)
     print*, 'gravity = ', (dt/mLayerDepth(iLayer))*(iLayerHydCond(iLayer) - iLayerHydCond(iLayer-1))
     print*, 'phase change = ', (iden_ice/iden_water)*(mLayerVolFracIceTrial(iLayer)-mLayerVolFracIce(iLayer))
     print*, 'sink = ', sink(iLayer)
     print*, 'rvec = ', rvec(iLayer)
    endif
   endif
  end do  ! (looping through layers)
  pause 'old routine: after rvec'
  end subroutine 


  ! ************************************************************************************************
  ! internal subroutine: assemble the tri-diagonal matrix
  ! ************************************************************************************************
  subroutine get_tridiag(mLayerMatricHeadTrial, & ! intent(in): matric head (m)
                         d_m1,                  & ! intent(out): sub-diagonal vector
                         diag,                  & ! intent(out): diagonal vector
                         d_p1,                  & ! intent(out): super-diagonal vector
                         err,message)             ! intent(out): error control
  ! compute residual vector for the volumetric fraction of liquid water
  implicit none
  ! input
  real(dp),intent(in)           :: mLayerMatricHeadTrial(:)  ! matric head in each layer (m)
  ! output
  real(dp),intent(out)          :: d_m1(:)                   ! sub-diagonal vector (m-1) 
  real(dp),intent(out)          :: diag(:)                   ! diagonal vector (m-1) 
  real(dp),intent(out)          :: d_p1(:)                   ! super-diagonal vector (m-1) 
  integer(i4b),intent(out)      :: err                       ! error code
  character(*),intent(out)      :: message                   ! error message
  ! local variables
  character(LEN=256)            :: cmessage                  ! error message of downwind routine 
  ! initialize error control
  err=0; message="get_tridiag/"

  ! loop through soil layers
  do iLayer=1,nSoil

   ! compute the derivative in the soil water characteristic (m-1)
   mLayerdTheta_dPsi(iLayer) = dTheta_dPsi(mLayerMatricHeadTrial(iLayer),vGn_alpha,theta_res,theta_sat,vGn_n,vGn_m)

   ! compute the off-diagonal elements
   if(iLayer<nSoil) d_p1(iLayer)   = -dt_dldz(iLayer)*iLayerHydCond(iLayer)
   if(iLayer>1)     d_m1(iLayer-1) = -dt_dudz(iLayer)*iLayerHydCond(iLayer-1)

   ! compute the diagonal elements
   if(bc_lower==neumann .and. iLayer==nSoil)then  ! Neumann boundary conditions for the lower-most layer
    diag(iLayer) = mLayerdTheta_dPsi(iLayer) + dt_dudz(iLayer)*iLayerHydCond(iLayer-1)
   elseif(bc_lower==neumann .and. iLayer==1)then  ! Neumann boundary conditions for the upper-most layer
    diag(iLayer) = mLayerdTheta_dPsi(iLayer) + dt_dldz(iLayer)*iLayerHydCond(iLayer)
   else ! interior layers, or diriclet boundary conditions 
    diag(iLayer) = mLayerdTheta_dPsi(iLayer) + dt_dldz(iLayer)*iLayerHydCond(iLayer) + dt_dudz(iLayer)*iLayerHydCond(iLayer-1)
   endif

  end do  ! (looping through layers)
  
  end subroutine get_tridiag


  ! ************************************************************************************************
  ! internal subroutine: line search
  ! ************************************************************************************************
  SUBROUTINE lnsrch(xold,                    & ! intent(in): matric head at the current iteration (m)
                    stpmax,                  & ! intent(in): maximum step size (m)
                    fold,                    & ! intent(in): function value for trial matric head vector (J m-3 J m-3)
                    g,                       & ! intent(in): gradient of the function vector (J m-3 J m-3 K-1)
                    p,                       & ! intent(in): iteration increment (m)
                    mLayerVolFracLiqIter,    & ! intent(in): volumetric fraction of liquid water at the current iteration (-)
                    mLayerVolFracIceIter,    & ! intent(in): volumetric fraction of ice at the current iteration (-)
                    x,                       & ! intent(out): new matric head vector (K)
                    rvec,                    & ! intent(out): new residual vector (J m-3)
                    f,                       & ! intent(out): new function value (J m-3 J m-3)
                    mLayerVolFracLiqNew,     & ! intent(out): new volumetric fraction of liquid water (-)
                    err,message)               ! intent(out): error control
  IMPLICIT NONE
  ! input variables
  real(dp) :: dt
  REAL(DP), DIMENSION(:), INTENT(IN) :: xold,g
  REAL(DP), DIMENSION(:), INTENT(INOUT) :: p
  REAL(DP), INTENT(IN) :: stpmax,fold
  real(dp),intent(in)           :: mLayerVolFracLiqIter(:)  ! volumetric fraction of liquid water (-)
  real(dp),intent(in)           :: mLayerVolFracIceIter(:)  ! volumetric fraction of ice (-)
  ! output variables
  REAL(DP), DIMENSION(:), INTENT(OUT) :: x,rvec
  REAL(DP), INTENT(OUT) :: f
  real(dp),intent(out)          :: mLayerVolFracLiqNew(:)   ! volumetric fraction of liquid water (-)
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
   ! update the matric head vector (m)
   x(:)=xold(:)+alam*p(:)
   ! compute matric head and volumetric fraction of liquid water
   do iLayer=1,nSoil
    mLayerVolFracLiqNew(iLayer) = volFracLiq(mLayerMatricHeadNew(iLayer),vGn_alpha,theta_res,theta_sat,vGn_n,vGn_m)
   end do
   ! compute the hydraulic conductivity for all layers
   call hydCond_all(x,                    & ! intent(in): matric head (m)
                    mLayerVolFracIceIter, & ! intent(in): volumetric fraction of ice (-)
                    mLayerVolFracLiqNew,  & ! intent(in): volumetric fraction of liquid water (-)
                    err,cmessage)           ! intent(out): error control
   if(err/=0)then; message=trim(message)//trim(cmessage); return; endif
   ! compute the residuals again
   call vlFrcLiqRes(x,                     & ! intent(in): matric head (m)
                    mLayerVolFracIceIter,  & ! intent(in): volumetric fraction of ice (-)
                    mLayerVolFracLiqNew,   & ! intent(in): volumetric fraction of liquid water (-)
                    rvec,                  & ! intent(out): residual vector (-)
                    err,cmessage)            ! intent(out): error control
   if(err/=0)then; message=trim(message)//trim(cmessage); return; endif
   ! compute the function evaluation
   f=0.5_dp*dot_product(rvec,rvec)
   ! additional exit criteria
   if(f<1.d-4)return
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
     if (tmplam > 0.5_dp*alam) tmplam=0.5_dp*alam
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


 end subroutine masschange

end module liquidflux_module

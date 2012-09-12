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
 subroutine liquidflow(dt,                      & ! time step (seconds)
                       iter,                    & ! current iteration count 
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
 integer(i4b),intent(in)       :: iter                       ! current iteration count
 real(dp),intent(in)           :: mLayerVolFracLiqIter(:)    ! volumetric fraction of liquid water at the current iteration (-)
 real(dp),intent(in)           :: mLayerVolFracIceIter(:)    ! volumetric fraction of ice at the current iteration (-)
 real(dp),intent(out)          :: mLayerVolFracLiqNew(:)     ! volumetric fraction of liquid water at the next iteration (-)
 integer(i4b),intent(out)      :: err                        ! error code
 character(*),intent(out)      :: message                    ! error message
 ! local pointers to model parameters
 real(dp),pointer              :: Fcapil                     ! capillary retention as a fraction of the total pore volume (-)
 real(dp),pointer              :: k_snow                     ! hydraulic conductivity of snow (m s-1), 0.0055 = approx. 20 m/hr, from UEB
 real(dp),pointer              :: mw_exp                     ! exponent for meltwater flow (-)
 ! local pointers to algorithmic control parameters
 real(dp),pointer              :: wimplicit                  ! weight assigned to start-of-step fluxes (-)
 ! local pointers to model forcing data
 real(dp),pointer              :: rainfall                   ! rainfall (kg m-2 s-1) 
 ! local pointers to model state variables 
 real(dp),pointer              :: mLayerDepth(:)             ! depth of the layer (m)
 real(dp),pointer              :: mLayerVolFracLiq(:)        ! volumetric fraction of liquid water in each snow layer (-)
 real(dp),pointer              :: mLayerVolFracIce(:)        ! volumetric fraction of ice in each snow layer (-)
 ! local pointers to model diagnostic variables
 real(dp),pointer              :: mLayerPoreSpace(:)         ! pore space in each snow layer (-)
 real(dp),pointer              :: mLayerThetaResid(:)        ! residual volumetric liquid water content in each snow layer (-)
 real(dp),pointer              :: iLayerInitLiqFluxSnow(:)   ! vertical liquid water flux at layer interfaces at the start of the time step (m s-1)
 real(dp),pointer              :: iLayerLiqFluxSnow(:)       ! vertical liquid water flux at layer interfaces at the end of the time step (m s-1)
 ! local pointers to model index variables
 integer(i4b),pointer          :: layerType(:)               ! type of the layer (ix_soil or ix_snow)
 ! define local variables
 integer(i4b)                  :: iLayer                     ! layer index
 integer(i4b)                  :: nSnow                      ! number of snow layers
 real(dp)                      :: xmin,xmax                  ! bounds for bi-section (-)
 real(dp)                      :: volFracLiqTrial            ! trial value of volumetric fraction of liquid water (-)
 real(dp)                      :: volFracLiqNew              ! new value of volumetric fraction of liquid water (-)
 real(dp)                      :: multResid                  ! multiplier for the residual water content (-)
 real(dp),parameter            :: residThrs=550._dp          ! ice density threshold to reduce residual liquid water content (kg m-3)
 real(dp),parameter            :: residScal=10._dp           ! scaling factor for residual liquid water content reduction factor (kg m-3)
 real(dp)                      :: volLiqRes                  ! residual volumetric liquid water content (-)
 real(dp)                      :: relPSpace                  ! relative pore space; min=res, max=porespace (-)
 real(dp)                      :: phseChnge                  ! volumetric liquid water equivalent associated with phase change (-)
 real(dp)                      :: relSaturn                  ! relative saturation [0,1] (-)
 real(dp)                      :: vDrainage                  ! vertical drainage (m s-1)
 real(dp)                      :: dflw_dliq                  ! derivative in vertical drainage (m s-1)
 real(dp)                      :: dt_dz                      ! dt/dz (s m-1)
 real(dp)                      :: lin_error                  ! linearization error [-residual] (-)
 real(dp)                      :: increment                  ! iteration increment (-)
 real(dp),parameter            :: atol=1.d-6                 ! absolute iteration tolerance (-)
 integer(i4b),parameter        :: maxiter=10                 ! maximum number of iterations
 integer(i4b)                  :: jiter                      ! internal iteration index
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

 ! assign pointers to algorithmic control parameters
 wimplicit     => mpar_data%var(iLookPARAM%wimplicit)        ! weight assigned to start-of-step fluxes (-)

 ! assign pointers to model forcing data
 rainfall => mvar_data%var(iLookMVAR%scalarRainfall)%dat(1)  ! computed rainfall rate (kg m-2 s-1)

 ! assign pointers to model state variables
 mLayerDepth       => mvar_data%var(iLookMVAR%mLayerDepth)%dat(1:nSnow)         ! depth of the layer (m)
 mLayerVolFracLiq  => mvar_data%var(iLookMVAR%mLayerVolFracLiq)%dat(1:nSnow)    ! volumetric fraction of liquid water in each snow layer (-)
 mLayerVolFracIce  => mvar_data%var(iLookMVAR%mLayerVolFracIce)%dat(1:nSnow)    ! volumetric fraction of ice in each snow layer (-)

 ! assign pointers to model diagnostic variables
 mLayerPoreSpace       => mvar_data%var(iLookMVAR%mLayerPoreSpace)%dat          ! pore space in each snowlayer (-)
 mLayerThetaResid      => mvar_data%var(iLookMVAR%mLayerThetaResid)%dat         ! residual volumetric liquid water content in each snow layer (-)
 iLayerInitLiqFluxSnow => mvar_data%var(iLookMVAR%iLayerInitLiqFluxSnow)%dat    ! liquid flux at layer interfaces at the start of the time step (m s-1)
 iLayerLiqFluxSnow     => mvar_data%var(iLookMVAR%iLayerLiqFluxSnow)%dat        ! liquid flux at layer interfaces at the end of the time step (m s-1)

 ! define the liquid flux at the upper boundary -- include evaporation/dew (m s-1)
 iLayerInitLiqFluxSnow(0) = rainfall/iden_water
 iLayerLiqFluxSnow(0)     = rainfall/iden_water

 ! check the meltwater exponent is >=1
 if(mw_exp<1._dp)then; err=20; message=trim(message)//'meltwater exponent < 1'; return; endif

 ! compute properties fixed over the time step, and initial fluxes
 if(iter==1)then
  ! loop through snow layers
  do iLayer=1,nSnow
   ! **(1)** compute properties fixed over the time step
   ! compute the reduction in liquid water holding capacity at high snow density (-)
   multResid = 1._dp / ( 1._dp + exp( (mLayerVolFracIce(iLayer)*iden_ice - residThrs) / residScal) )
   ! compute the pore space (-)
   mLayerPoreSpace(iLayer)  = 1._dp - mLayerVolFracIce(iLayer)
   ! compute the residual volumetric liquid water content (-)
   mLayerThetaResid(iLayer) = Fcapil*mLayerPoreSpace(iLayer) * multResid
   ! **(2)** compute the initial drainage flux
   call mw_func( mLayerVolFracLiq(iLayer),mLayerThetaResid(iLayer),mLayerPoreSpace(iLayer), & ! input
                 iLayerInitLiqFluxSnow(iLayer) )                                              ! output
  end do  ! (looping through snow layers)
 endif  ! (if the first iteration)

 ! compute fluxes at the end of the step
 do iLayer=1,nSnow
  ! compute dt/dz
  dt_dz     = dt/mLayerDepth(iLayer)
  ! compute the liquid water equivalent associated with phase change
  phseChnge = (iden_ice/iden_water)*(mLayerVolFracIceIter(iLayer) - mLayerVolFracIce(iLayer))
  ! initialize bounds
  xmin = 0._dp
  xmax = mLayerPoreSpace(iLayer)
  ! ****************************
  ! ***** begin iterations *****
  ! ****************************
  ! initialize the volumetric fraction of liquid water (-)
  volFracLiqTrial = mLayerVolFracLiqIter(iLayer)
  do jiter=1,maxiter
   ! compute the drainage flux and its derivative
   call mw_func(volFracLiqTrial,mLayerThetaResid(iLayer),mLayerPoreSpace(iLayer), & ! input
                vDrainage, dflw_dliq)                                               ! output
   ! compute the residual (-)
   lin_error = dt_dz*(          wimplicit *(iLayerInitLiqFluxSnow(iLayer-1) - iLayerInitLiqFluxSnow(iLayer)) +   &
                       (1._dp - wimplicit)*(iLayerLiqFluxSnow(iLayer-1)     - vDrainage                    ) ) - &
               phseChnge - (volFracLiqTrial - mLayerVolFracLiq(iLayer))
   ! compute the iteration increment (-) and new value
   increment = lin_error/(1._dp + dt_dz*(1._dp - wimplicit)*dflw_dliq)
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
    write(*,'(a,2(i4,1x),10(f20.10,1x))') 'in watersnow, iter, iLayer, volFracLiqNew, volLiqRes = ', iter, iLayer, increment, volFracLiqNew, volLiqRes
   ! check for convergence
   if(abs(increment) < atol) exit
   ! get ready for next iteration
   volFracLiqTrial = volFracLiqNew
  end do  ! (end iterations)
  ! save state
  mLayerVolFracLiqNew(iLayer) = volFracLiqNew
  ! get ready to process the next snow layer
  iLayerLiqFluxSnow(iLayer) = vDrainage + dflw_dliq*increment ! second term will be zero if converge completely
  ! *** now process the next layer
 end do  ! (looping through snow layers)

 contains

  ! calculate vertical drainage and its derivative
  subroutine mw_func(volFracLiqTrial, &   ! intent(in): volumetric fraction of liquid water (-)
                     volFracLiqRes,   &   ! intent(in): residual volumetric fraction of liquid water (-)
                     poreSpace,       &   ! intent(in): pore space (-)
                     vDrainage,       &   ! intent(out): vertical drainage (m s-1)
                     dflw_dliq)           ! intent(out): optional: function derivative (m s-1)
  implicit none
  ! input variables
  real(dp),intent(in)           :: volFracLiqTrial   ! volumetric fraction of liquid water (-)
  real(dp),intent(in)           :: volFracLiqRes     ! residual volumetric fraction of liquid water (-)
  real(dp),intent(in)           :: poreSpace         ! pore space (-)
  ! output variables
  real(dp),intent(out)          :: vDrainage         ! vertical drainage (m s-1)
  real(dp),intent(out),optional :: dflw_dliq         ! function derivative (m s-1)
  ! ***** check that flow occurs
  if(volFracLiqTrial > volFracLiqRes)then
   ! compute the relative saturation (-)
   relSaturn = (volFracLiqTrial - volFracLiqRes)/(poreSpace - volFracLiqRes)
   ! compute the initial flow estimate (m s-1)
   vDrainage = k_snow*relSaturn**mw_exp
   ! compute the change in the vertical liquid water flux with volumetric liquid water content (m s-1)
   if(present(dflw_dliq)) dflw_dliq = ( (k_snow*mw_exp)/(poreSpace - volFracLiqRes) ) * relSaturn**(mw_exp - 1._dp)
  ! **** default = flow does not occur, and the derivative is zero
  else
   vDrainage = 0._dp
   if(present(dflw_dliq)) dflw_dliq = 0._dp
  endif
  end subroutine mw_func

 end subroutine liquidflow

 ! ************************************************************************************************
 ! new subroutine: compute liquid water flux through the soil
 ! ************************************************************************************************
 subroutine masschange(dt,&                   ! time step (seconds)
                       iter,&                 ! iteration index
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
 integer(i4b),intent(in)       :: iter                       ! iteration index
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
 real(dp),pointer              :: specficStorage           ! specific storage coefficient (m-1)
 real(dp),pointer              :: f_impede                 ! ice impedence factor (-)
 real(dp),pointer              :: upperBoundHead           ! upper boundary condition for matric head (m)
 real(dp),pointer              :: lowerBoundHead           ! lower boundary condition for matric head (m)
 real(dp),pointer              :: wimplicit                ! weight assigned to start-of-step fluxes (-)
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
 real(dp),pointer              :: scalarRainPlusMelt       ! rain plus melt, used as input to the soil zone before computing surface runoff (m s-1)
 real(dp),pointer              :: scalarSurfaceRunoff      ! surface runoff (m s-1)
 ! local pointers to model diagnostic variables
 real(dp),pointer              :: mLayerDepth(:)           ! depth of the layer (m)
 real(dp),pointer              :: mLayerHeight(:)          ! height of the layer mid-point (m)
 real(dp),pointer              :: mLayerdTheta_dPsi(:)     ! derivative in the soil water characteristic (m-1)
 real(dp),pointer              :: mLayerInitTranspire(:)   ! transpiration loss from each soil layer at the start of the time step (m s-1)
 real(dp),pointer              :: mLayerTranspire(:)       ! transpiration loss from each soil layer (m s-1)
 real(dp),pointer              :: iLayerInitLiqFluxSoil(:) ! liquid flux at soil layer interfaces at the start of the time step (m s-1)
 real(dp),pointer              :: iLayerLiqFluxSoil(:)     ! liquid flux at soil layer interfaces at the end of the time step (m s-1)
 ! local pointers to model index variables
 integer(i4b),pointer          :: layerType(:)             ! type of the layer (ix_soil or ix_snow)
 ! define local variables
 character(LEN=256)            :: cmessage                 ! error message of downwind routine
 integer(i4b)                  :: nSnow                    ! number of snow layers
 integer(i4b)                  :: nSoil                    ! number of soil layers
 integer(i4b)                  :: ibeg,iend                ! start and end indices of the soil layers in concatanated vector
 integer(i4b)                  :: iLayer                   ! layer index
 real(dp)                      :: satArea                  ! saturated area (-)
 real(dp)                      :: iceImpedeSurface         ! ice impedence factor at the surface (-)
 real(dp),parameter            :: vic_bpar=0.1_dp          ! the VIC "b" parameter, defining the fraction of saturated area
 real(dp),dimension(size(mLayerMatricHeadIter))    :: theta          ! volumetric fraction of total water, if all is melted (-)
 real(dp),dimension(size(mLayerMatricHeadIter))    :: mLayerHydCond  ! hydraulic conductivity at layer mid-point (m s-1)
 real(dp),dimension(0:size(mLayerMatricHeadIter))  :: iLayerHydCond  ! hydraulic conductivity at layer interface (m s-1)
 real(dp),dimension(size(mLayerMatricHeadIter)-1)  :: d_m1           ! sub-diagonal elements of the tridiagonal system (m-1)
 real(dp),dimension(size(mLayerMatricHeadIter))    :: diag           ! diagonal elements of the tridiagonal system (m-1)
 real(dp),dimension(size(mLayerMatricHeadIter)-1)  :: d_p1           ! super-diagonal elements of the tridiagonal system (m-1)
 real(dp),dimension(size(mLayerMatricHeadIter))    :: rvec           ! right-hand-side vector (-)
 integer(i4b),parameter        :: diriclet=1001            ! look-up value for diriclet boundary conditions
 integer(i4b),parameter        :: neumann=1002             ! look-up value for neumann boundary conditions
 integer(i4b)                  :: bc_upper,bc_lower        ! boundary condition index (diriclet or neumann)
 real(dp)                      :: fracFrozen               ! fraction of frozen ground (-)
 real(dp),dimension(size(mLayerMatricHeadIter))   :: g     ! gradient of the function vector (m-1)
 real(dp)                      :: fold,fnew                ! function values (-)
 real(dp),parameter            :: STPMX=5._dp              ! maximum step size in line search (m)
 real(dp)                      :: stpmax                   ! scaled maximum step size
 logical(lgt)                  :: printflag                ! flag to print crap to the screen
 real(dp),dimension(size(mLayerMatricHeadIter),size(mLayerMatricHeadIter)) :: jmat  ! jacobian matrix
 real(dp),dimension(size(mLayerMatricHeadIter))                            :: ftest,xsav,xtry ! compute jacobian: function vector and dependent variables
 real(dp),dimension(size(mLayerMatricHeadIter))                            :: xph,h           ! compute jacobian: perturbed vector and finite difference increment
 real(dp),parameter                                                        :: eps=-1.0e-8_dp  ! compute jacobian: finite difference increment
 integer(i4b)                                                              :: ijac            ! compute jacobian: index of columns
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
 specficStorage    => mpar_data%var(iLookPARAM%specficStorage)                 ! specific storage coefficient (m-1)
 f_impede          => mpar_data%var(iLookPARAM%f_impede)                       ! ice impedence factor (-)

 ! assign pointers to head boundary conditions (included as parameters for convenience)
 lowerBoundHead    => mpar_data%var(iLookPARAM%lowerBoundHead)
 upperBoundHead    => mpar_data%var(iLookPARAM%upperBoundHead)

 ! assign pointers to algorithmic control parameters
 wimplicit         => mpar_data%var(iLookPARAM%wimplicit)                      ! weight assigned to start-of-step fluxes (-)

 ! assign pointers to model forcing data
 rainfall          => mvar_data%var(iLookMVAR%scalarRainfall)%dat(1)           ! computed rainfall rate (kg m-2 s-1)

 ! assign pointers to model variables that are constant over the simulation period
 vGn_m             => mvar_data%var(iLookMVAR%scalarVGn_m)%dat(1)              ! van Genutchen "m" parameter (-)
 volLatHt_fus      => mvar_data%var(iLookMVAR%scalarvolLatHt_fus)%dat(1)       ! volumetric latent heat of fusion (J m-3)

 ! assign pointers to model state variables
 mLayerVolFracIce  => mvar_data%var(iLookMVAR%mLayerVolFracIce)%dat(ibeg:iend) ! volumetric fraction of ice in each layer (-)
 mLayerVolFracLiq  => mvar_data%var(iLookMVAR%mLayerVolFracLiq)%dat(ibeg:iend) ! volumetric fraction of liquid water in each layer (-)
 mLayerMatricHead  => mvar_data%var(iLookMVAR%mLayerMatricHead)%dat      ! (soil only) ! matric head in each layer (m)

 ! assign local pointers to diagnostic scalar variables
 scalarSfcMeltPond     => mvar_data%var(iLookMVAR%scalarSfcMeltPond)%dat(1)    ! ponded water caused by melt of the "snow without a layer" (kg m-2)
 scalarRainPlusMelt    => mvar_data%var(iLookMVAR%scalarRainPlusMelt)%dat(1)   ! rain plus melt (m s-1)
 scalarSurfaceRunoff   => mvar_data%var(iLookMVAR%scalarSurfaceRunoff)%dat(1)  ! surface runoff (m s-1)

 ! assign pointers to model diagnostic variables
 mLayerDepth           => mvar_data%var(iLookMVAR%mLayerDepth)%dat(ibeg:iend)  ! depth of the layer (m)
 mLayerHeight          => mvar_data%var(iLookMVAR%mLayerHeight)%dat(ibeg:iend) ! height of the layer mid-point (m)
 mLayerdTheta_dPsi     => mvar_data%var(iLookMVAR%mLayerdTheta_dPsi)%dat       ! (soil only) ! derivative in the soil water characteristic (m-1)
 mLayerInitTranspire   => mvar_data%var(iLookMVAR%mLayerInitTranspire)%dat     ! (soil only) ! transpiration loss from each soil layer at start-of-step (m s-1)
 mLayerTranspire       => mvar_data%var(iLookMVAR%mLayerTranspire)%dat         ! (soil only) ! transpiration loss from each soil layer (m s-1)
 iLayerInitLiqFluxSoil => mvar_data%var(iLookMVAR%iLayerInitLiqFluxSoil)%dat   ! (soil only) ! liquid flux at layer interfaces at the start of the time step (m s-1)
 iLayerLiqFluxSoil     => mvar_data%var(iLookMVAR%iLayerLiqFluxSoil)%dat       ! (soil only) ! liquid flux at layer interfaces at the end of the time step (m s-1)

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

 ! define upper boundary fluxes (m s-1)
 if(bc_upper==neumann) then
  if(nSnow==0) scalarRainPlusMelt = rainfall/iden_water + (scalarSfcMeltPond/dt)/iden_water    ! rainfall plus melt of the snow without a layer
  if(nSnow>0)  scalarRainPlusMelt = mvar_data%var(iLookMVAR%iLayerLiqFluxSnow)%dat(nSnow)  ! liquid water flux from the base of the snowpack (m s-1)
 endif
 !print*, 'scalarRainPlusMelt, rainfall/iden_water, (scalarSfcMeltPond/dt)/iden_water, scalarSfcMeltPond = ', &
 !         scalarRainPlusMelt, rainfall/iden_water, (scalarSfcMeltPond/dt)/iden_water, scalarSfcMeltPond

 ! compute initial fluxes
 if(iter==1)then
  ! compute the hydraulic conductivity for all layers
  call hydCond_all(mLayerMatricHead, & ! intent(in): matric head (m)
                   mLayerVolFracIce, & ! intent(in): volumetric fraction of ice (-)
                   mLayerVolFracLiq, & ! intent(in): volumetric fraction of liquid water (-)
                   err,cmessage)       ! intent(out): error control
  if(err/=0)then; message=trim(message)//trim(cmessage); return; endif
  ! compute liquid water fluxes at layer interfaces
  call iLayer_liq(mLayerMatricHead,      & ! intent(in):  matric head in each layer (m)
                  mLayerVolFracIce,      & ! intent(in):  volumetric ice content in each layer (-)
                  mLayerVolFracLiq,      & ! intent(in):  volumetric liquid water content (-)
                  iLayerInitLiqFluxSoil, & ! intent(out): liquid water flux at the soil layer interfaces (m s-1)
                  err,cmessage)            ! intent(out): error control
  if(err/=0)then; message=trim(message)//trim(cmessage); return; endif
 endif 

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
                  iLayerLiqFluxSoil,    & ! intent(out): liquid water flux at the soil layer interfaces (m s-1)
                  err,cmessage)           ! intent(out): error control
 if(err/=0)then; message=trim(message)//trim(cmessage); return; endif
 if(printflag) print*, 'residuals -- take1 ', rvec

 ! *****
 ! compute the tri-diagonal matrix
 call get_tridiag(mLayerMatricHeadIter,  & ! intent(in): matric head (m)
                  mLayerVolFracLiqIter,  & ! intent(in): volumetric fraction of liquid water (m)
                  d_m1,                  & ! intent(out): sub-diagonal vector
                  diag,                  & ! intent(out): diagonal vector
                  d_p1,                  & ! intent(out): super-diagonal vector
                  err,cmessage)            ! intent(out): error control
 if(err/=0)then; message=trim(message)//trim(cmessage); return; endif


 ! ***** compute the Jacobian using one-sided finite differences
 !print*, 'starting to compute Jacobian'
 !xtry=mLayerMatricHeadIter
 !xsav=xtry
 !h=EPS*abs(xsav)
 !where (h == 0.0) h=EPS
 !xph=xsav+h
 !h=xph-xsav
 !do ijac=1,5
 !print*, 'ijac = ', ijac
 ! xtry(ijac)=xph(ijac)
 ! call vlFrcLiqRes(xtry,                & ! intent(in): trial matric head (m)
 !                  mLayerVolFracIceIter,& ! intent(in): trial volumetric fraction of ice (-)
 !                  mLayerVolFracLiqIter,& ! intent(in): trial volumetric fraction of liquid water (-)
 !                  ftest,               & ! intent(out): residual vector (J m-3)
 !                  err,cmessage)          ! intent(out): error control
 ! if(err/=0)then; err=10; message=trim(message)//trim(cmessage); return; endif
 ! jmat(:,ijac)=(ftest(:)-rvec(:))/h(ijac)
 ! xtry(ijac)=xsav(ijac)
 !end do

 !print*, 'matric head Jacobian = '
 !do ijac=1,5
 ! write(*,'(i4,1x,5(e20.10,1x))') ijac, jmat(1:5,ijac)
 !end do

 !write(*,'(a,1x,5(e20.10,1x))') 'diag = ', diag(1:5)
 !write(*,'(a,1x,5(e20.10,1x))') 'd_m1 = ', d_m1(1:5)
 !write(*,'(a,1x,5(e20.10,1x))') 'd_p1 = ', d_p1(1:5)

 ! *****
 ! solve the tridiagonal system of equations -- returns mLayerMatricHeadNew
 call tridag(d_m1,diag,d_p1,-rvec,mLayerMatricHeadDiff,err,cmessage)
 if(err/=0)then
  write(*,'(50(e20.10,1x))') mLayerMatricHeadIter
  message=trim(message)//trim(cmessage); return
 endif

 ! test tridag
 if(printflag)then
  write(*,'(a)') 'test tridag: iLayer, d_m1(iLayer-1), diag(iLayer), d_p1(iLayer), mLayerMatricHeadDiff(iLayer), rvec(iLayer)'
  do iLayer=1,5
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

 ! update the fluxes
 

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

 ! check for an early return
 if(maxval(abs(mLayerMatricHeadDiff)) < epsilon(mLayerMatricHeadDiff))then
  ! save matric head
  mLayerMatricHeadNew = mLayerMatricHeadIter
  ! volumetric fraction of liquid water
  do iLayer=1,nSoil
   mLayerVolFracLiqNew(iLayer) = volFracLiq(mLayerMatricHeadNew(iLayer),vGn_alpha,theta_res,theta_sat,vGn_n,vGn_m)
  end do
  return
 endif

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
             iLayerLiqFluxSoil,       & ! intent(out): liquid water flux at the soil layer interfaces (m s-1)
             err,cmessage)              ! intent(out): error control
 if(err/=0)then; message=trim(message)//trim(cmessage); return; endif

 ! update the matric head increment
 mLayerMatricHeadDiff = mLayerMatricHeadNew - mLayerMatricHeadIter

 ! compute the boundary fluxes
 if(bc_upper==diriclet)&
  iLayerLiqFluxSoil(0)     = iLayerHydCond(0)*(upperBoundHead - mLayerMatricHeadNew(1))/(0.5_dp*mLayerDepth(1)) + iLayerHydCond(0)
 if(bc_lower==diriclet)&
  iLayerLiqFluxSoil(nSoil) = iLayerHydCond(nSoil)*(mLayerMatricHeadNew(nSoil) - lowerBoundHead)/(0.5_dp*mLayerDepth(nSoil)) + iLayerHydCond(nSoil)


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
  real(dp),dimension(size(mLayerMatricHeadTrial)) :: relIce    ! relative ice content in each layer soil layer (-)
  real(dp),dimension(size(mLayerMatricHeadTrial)) :: iceImpede ! ice impedence factor (-)
  real(dp)                      :: test                      ! scalar for testing
  ! initialize error control
  err=0; message="hydCond_all/"
 
  ! loop through layers
  do iLayer=1,nSoil
   ! check that the volumetric fraction of liquid water and ice is less than the porosity
   !if(mLayerVolFracIceTrial(iLayer) + mLayerVolFracLiqTrial(iLayer) > theta_sat)then
   ! write(message,'(a,i0,a,f20.15,a)') &
   !  trim(message)//"theta > theta_sat [iLayer=",iLayer,"; volFracIce=",mLayerVolFracIceTrial(iLayer),"]"
   ! err=20; return
   !endif
   ! compute the relative volumetric fraction of ice
   relIce(iLayer) = mLayerVolFracIceTrial(iLayer)/(theta_sat - theta_res)
   ! compute the ice impedence factor (-)
   iceImpede(iLayer) = 10._dp**(-f_impede*relIce(iLayer))
   ! compute the hydraulic conductivity for a given layer (m s-1)
   mLayerHydCond(iLayer) = hydCond(mLayerMatricHeadTrial(iLayer),k_soil,vGn_alpha,vGn_n,vGn_m) * iceImpede(iLayer)
   ! compute the hydraulic conductivity for layer interfaces (m s-1)
   if(iLayer>1) &  ! *** NOTE: use the geometric mean
    iLayerHydCond(iLayer-1) = (mLayerHydCond(iLayer-1) * mLayerHydCond(iLayer))**0.5_dp
  end do
  ! save the ice impedence factor in the top soil layer
  iceImpedeSurface = iceImpede(1)

  ! if diriclet, compute hydraulic conductivity at the boundaries (m s-1)
  if(bc_upper==diriclet) iLayerHydCond(0)     = hydCond(upperBoundHead,k_soil,vGn_alpha,vGn_n,vGn_m) * 10._dp**(f_impede*relIce(1))
  if(bc_lower==diriclet) iLayerHydCond(nSoil) = hydCond(lowerBoundHead,k_soil,vGn_alpha,vGn_n,vGn_m) * 10._dp**(f_impede*relIce(nSoil))

  ! if newmann, set hydraulic conductivity at boundaries to the layer conductivity (*** not used ***)
  if(bc_upper==neumann) iLayerHydCond(0)     = -huge(mLayerHydCond)
  if(bc_lower==neumann) iLayerHydCond(nSoil) = -huge(mLayerHydCond)

  end subroutine hydCond_all


  ! ************************************************************************************************
  ! internal subroutine: compute liquid fluxes at layer interfaces
  ! ************************************************************************************************
  subroutine iLayer_liq(mLayerMatricHeadTrial,  & ! intent(in):  matric head in each layer (m)
                        mLayerVolFracIceTrial,  & ! intent(in):  volumetric ice content in each layer (-)
                        mLayerVolFracLiqTrial,  & ! intent(in):  volumetric liquid water content (-)
                        iLayerLiqFlux,          & ! intent(out): liquid fluxes at layer interfaces (m s-1)
                        err,message)              ! intent(out): error control
  ! compute liquid flux at layer interfaces
  implicit none
  ! input
  real(dp),intent(in)           :: mLayerMatricHeadTrial(:)  ! matric head in each layer (m)
  real(dp),intent(in)           :: mLayerVolFracIceTrial(:)  ! volumetric ice content in each layer (-)
  real(dp),intent(in)           :: mLayerVolFracLiqTrial(:)  ! volumetric liquid water content in each layer (-)
  ! output
  real(dp),intent(out)          :: iLayerLiqFlux(0:)         ! liquid fluxes at layer interfaces (m s-1)
  integer(i4b),intent(out)      :: err                       ! error code
  character(*),intent(out)      :: message                   ! error message
  ! local variables
  character(LEN=256)            :: cmessage                  ! error message of downwind routine 
  real(dp)                      :: fracCap                   ! fraction of pore space filled with liquid water and ice
  real(dp)                      :: infUnsat                  ! infiltration over unsaturated areas (m s-1)
  real(dp)                      :: cflux                     ! capillary flux (m s-1)
  integer(i4b)                  :: iSoil                     ! index of soil layer
  real(dp)                      :: zWater                    ! depth to the water table (m)
  real(dp),parameter            :: zScale=0.1_dp             ! scaling factor (m)
  real(dp),parameter            :: kAnisotropic=10._dp       ! anisotropic factor for hydraulic conductivity (-)
  ! initialize error control
  err=0; message="iLayer_liq/"

  ! compute fluxes at the upper boundary -- positive downwards
  select case(bc_upper)
   case(neumann)
    ! compute the surface runoff (m s-1)
    fracCap = min((mLayerVolFracLiqTrial(1) + mLayerVolFracIceTrial(1))/theta_sat, 1._dp)  ! fraction of pore space filled with liquid water and ice
    !if(fracCap > 0.99_dp) fracCap = 1._dp
    infUnsat= min(k_soil*iceImpedeSurface,scalarRainPlusMelt)           ! infiltration over unsaturated areas (m s-1)
    !infUnsat= min(k_soil,scalarRainPlusMelt)        ! infiltration over unsaturated areas (m s-1)
    satArea = 1._dp - (1._dp - fracCap)**vic_bpar   ! saturated fraction (-)
    scalarSurfaceRunoff = satArea*scalarRainPlusMelt + (1._dp - satArea)*(scalarRainPlusMelt - infUnsat)
    !print*, 'scalarRainPlusMelt, fracCap, infUnsat, satArea, scalarSurfaceRunoff = ', &
    !         scalarRainPlusMelt, fracCap, infUnsat, satArea, scalarSurfaceRunoff
    ! compute the flux at the upper boundary
    iLayerLiqFlux(0) = scalarRainPlusMelt - scalarSurfaceRunoff
   case(diriclet)
    cflux = -iLayerHydCond(0)*(mLayerMatricHeadTrial(1) - upperBoundHead) / (mLayerDepth(1)*0.5_dp)
    iLayerLiqFlux(0) = cflux + iLayerHydCond(0)
  endselect

  ! compute fluxes at the lower boundary -- positive downwards
  select case(bc_lower)
   case(neumann)
    zWater = mLayerHeight(nSoil) - mLayerMatricHeadTrial(nSoil)
    iLayerLiqFlux(nSoil) = kAnisotropic*k_soil * exp(-zWater/zScale) 
    !iLayerLiqFlux(nSoil) = mLayerHydCond(nSoil)
   case(diriclet)
    cflux = -iLayerHydCond(nSoil)*(lowerBoundHead - mLayerMatricHeadTrial(nSoil)) / (mLayerDepth(nSoil)*0.5_dp)
    iLayerLiqFlux(nSoil) = cflux + iLayerHydCond(nSoil)
  endselect

  ! compute fluxes within the domain -- positive downwards
  do iSoil=1,nSoil-1 ! iSoil=0 is the upper boundary, so flux at iSoil=1 is bottom of the top layer
   ! compute the capillary flux (negative sign means positive downwards)
   cflux = -iLayerHydCond(iSoil)*(mLayerMatricHeadTrial(iSoil+1) - mLayerMatricHeadTrial(iSoil)) / &
                                 (mLayerHeight(iSoil+1) - mLayerHeight(iSoil))
   ! compute the total flux (add gravity flux, positive downwards)
   iLayerLiqFlux(iSoil) = cflux + iLayerHydCond(iSoil)
  end do ! looping through layers within the domain
  if(mLayerMatricHeadTrial(1) > 0._dp) print*, 'iLayerLiqFlux(0:1) = ', iLayerLiqFlux(0:1)

  end subroutine iLayer_liq


  ! ************************************************************************************************
  ! internal subroutine: compute the residual vector
  ! ************************************************************************************************
  subroutine vlFrcLiqRes(mLayerMatricHeadTrial, & ! intent(in): matric head (m)
                         mLayerVolFracIceTrial, & ! intent(in): volumetric fraction of ice (-)
                         mLayerVolFracLiqTrial, & ! intent(in): volumetric fraction of liquid water (-)
                         rvec,                  & ! intent(out): residual vector (-)
                         iLayerLiqFlux,         & ! intent(out): liquid fluxes at layer interfaces (m s-1)
                         err,message)             ! intent(out): error control
  ! compute residual vector for the volumetric fraction of liquid water
  implicit none
  ! input
  real(dp),intent(in)           :: mLayerMatricHeadTrial(:)  ! matric head in each layer (m)
  real(dp),intent(in)           :: mLayerVolFracIceTrial(:)  ! volumetric fraction of ice in each layer (-)
  real(dp),intent(in)           :: mLayerVolFracLiqTrial(:)  ! volumetric fraction of liquid water in each layer (-)
  ! output
  real(dp),intent(out)          :: rvec(:)                   ! residual vector (-)
  real(dp),intent(out)          :: iLayerLiqFlux(0:)         ! liquid fluxes at layer interfaces (m s-1)
  integer(i4b),intent(out)      :: err                       ! error code
  character(*),intent(out)      :: message                   ! error message
  ! local variables
  character(LEN=256)            :: cmessage                  ! error message of downwind routine 
  integer(i4b)                  :: iSoil                     ! index of soil layers
  ! initialize error control
  err=0; message="vlFrcLiqRes/"

  ! compute liquid water fluxes at layer interfaces
  call iLayer_liq(mLayerMatricHeadTrial,  & ! intent(in):  matric head in each layer (m)
                  mLayerVolFracIceTrial,  & ! intent(in):  volumetric ice content in each layer (-)
                  mLayerVolFracLiqTrial,  & ! intent(in):  volumetric liquid water content (-)
                  iLayerLiqFlux,          & ! intent(out): liquid water flux at the layer interfaces (m s-1)
                  err,cmessage)             ! intent(out): error control
  if(err/=0)then; message=trim(message)//trim(cmessage); return; endif

  ! compute the residual vector (-)
  do iSoil=1,nSoil
   rvec(iSoil) = mLayerVolFracLiqTrial(iSoil) - &
                 ( &
                   mLayerVolFracLiq(iSoil) + &
                   (-(iLayerInitLiqFluxSoil(iSoil) - iLayerInitLiqFluxSoil(iSoil-1)) + mLayerInitTranspire(iSoil))*(dt/mLayerDepth(iSoil))*wimplicit + &
                   (-(iLayerLiqFluxSoil(iSoil) - iLayerLiqFluxSoil(iSoil-1)) + mLayerTranspire(iSoil))*(dt/mLayerDepth(iSoil))*(1._dp - wimplicit) + &
                   (-(iden_ice/iden_water)*(mLayerVolFracIceTrial(iSoil)-mLayerVolFracIce(iSoil))) &
!                   (-(mLayervolFracLiqTrial(iSoil)/theta_sat)*specficStorage*(mLayerMatricHeadTrial(iSoil)-mLayerMatricHead(iSoil))) &
                 )
   !if(iSoil==1)then
    !print*, 'rvec = ', rvec(iSoil)
    !print*, 'iLayerLiqFlux(0:1) = ', iLayerLiqFlux(0:1)
    !print*, 'mLayerVolFracLiqTrial(iSoil) = ', mLayerVolFracLiqTrial(iSoil)
    !print*, 'mLayerVolFracLiq(iSoil) = ', mLayerVolFracLiq(iSoil)
    !print*, '-(iLayerInitLiqFluxSoil(iSoil) - iLayerInitLiqFluxSoil(iSoil-1)) = ', -(iLayerInitLiqFluxSoil(iSoil) - iLayerInitLiqFluxSoil(iSoil-1))
    !print*, '-(iLayerLiqFluxSoil(iSoil) - iLayerLiqFluxSoil(iSoil-1)) = ', -(iLayerLiqFluxSoil(iSoil) - iLayerLiqFluxSoil(iSoil-1))
    !print*, '(-(iLayerInitLiqFluxSoil(iSoil) - iLayerInitLiqFluxSoil(iSoil-1)) + mLayerInitTranspire(iSoil))*(dt/mLayerDepth(iSoil))*wimplicit = ',&
    !         (-(iLayerInitLiqFluxSoil(iSoil) - iLayerInitLiqFluxSoil(iSoil-1)) + mLayerInitTranspire(iSoil))*(dt/mLayerDepth(iSoil))*wimplicit
    !print*, '(-(iLayerLiqFluxSoil(iSoil) - iLayerLiqFluxSoil(iSoil-1)) + mLayerTranspire(iSoil))*(dt/mLayerDepth(iSoil))*(1._dp - wimplicit) = ', &
    !         (-(iLayerLiqFluxSoil(iSoil) - iLayerLiqFluxSoil(iSoil-1)) + mLayerTranspire(iSoil))*(dt/mLayerDepth(iSoil))*(1._dp - wimplicit) 
    !print*, '-(iden_ice/iden_water)*(mLayerVolFracIceTrial(iSoil)-mLayerVolFracIce(iSoil)) = ',&
    !         -(iden_ice/iden_water)*(mLayerVolFracIceTrial(iSoil)-mLayerVolFracIce(iSoil))
   !endif

   !if(mLayerVolFracLiqTrial(iSoil) > theta_sat - 0.01_dp)then
   ! print*, 'mLayerInitTranspire(iSoil), mLayerTranspire(iSoil) = ', mLayerInitTranspire(iSoil), mLayerTranspire(iSoil)
   !endif
   !print*, 'iSoil, mLayerMatricHeadTrial(iSoil), rvec(iSoil) = ', iSoil, mLayerMatricHeadTrial(iSoil), rvec(iSoil)
  end do  ! (looping through soil layers)

  end subroutine 


  ! ************************************************************************************************
  ! internal subroutine: assemble the tri-diagonal matrix
  ! ************************************************************************************************
  subroutine get_tridiag(mLayerMatricHeadTrial, & ! intent(in): matric head (m)
                         mLayerVolFracLiqTrial, & ! intent(in): volumetric fraction of liquid water (-)
                         d_m1,                  & ! intent(out): sub-diagonal vector
                         diag,                  & ! intent(out): diagonal vector
                         d_p1,                  & ! intent(out): super-diagonal vector
                         err,message)             ! intent(out): error control
  ! compute residual vector for the volumetric fraction of liquid water
  implicit none
  ! input
  real(dp),intent(in)           :: mLayerMatricHeadTrial(:)  ! matric head in each layer (m)
  real(dp),intent(in)           :: mLayerVolFracLiqTrial(:)  ! volumetric fraction of liquid water in each layer (-)
  ! output
  real(dp),intent(out)          :: d_m1(:)                   ! sub-diagonal vector (m-1) 
  real(dp),intent(out)          :: diag(:)                   ! diagonal vector (m-1) 
  real(dp),intent(out)          :: d_p1(:)                   ! super-diagonal vector (m-1) 
  integer(i4b),intent(out)      :: err                       ! error code
  character(*),intent(out)      :: message                   ! error message
  ! local variables
  character(LEN=256)            :: cmessage                  ! error message of downwind routine 
  real(dp),dimension(size(mLayerMatricHeadTrial)-1) :: dz_node  ! distance between the mid-point of a layer and a neighbouring layer (m)
  real(dp),dimension(size(mLayerMatricHeadTrial))   :: dt_dudz  ! the dt_dudz terms for the upper interface (s m-2)
  real(dp),dimension(size(mLayerMatricHeadTrial))   :: dt_dldz  ! the dt_dudz terms for the lower interface (s m-2)
  ! initialize error control
  err=0; message="get_tridiag/"

  ! compute the dt/dzdz terms for each layer(s m-2)
  dz_node(1:nSoil-1) = mLayerHeight(2:nSoil) - mLayerHeight(1:nSoil-1)
  dt_dudz(2:nSoil)   = (1._dp - wimplicit)*dt/(mLayerDepth(2:nSoil)*dz_node)    ! upper distance, valid 2..nSoil
  dt_dldz(1:nSoil-1) = (1._dp - wimplicit)*dt/(mLayerDepth(1:nSoil-1)*dz_node)  ! lower distance, valid 1..nSoil-1

  ! compute the dt/dzdz terms at the boundaries
  dt_dudz(1)     = (1._dp - wimplicit)*dt/(mLayerDepth(1)     * mLayerDepth(1)*0.5_dp)
  dt_dldz(nSoil) = (1._dp - wimplicit)*dt/(mLayerDepth(nSoil) * mLayerDepth(nSoil)*0.5_dp)

  ! loop through soil layers
  do iLayer=1,nSoil

   ! compute the derivative in the soil water characteristic (m-1)
   mLayerdTheta_dPsi(iLayer) = dTheta_dPsi(mLayerMatricHeadTrial(iLayer),vGn_alpha,theta_res,theta_sat,vGn_n,vGn_m)
   !if(mLayerMatricHeadTrial(iLayer) > 0._dp) print*, 'iLayer,mLayerdTheta_dPsi(iLayer) = ', iLayer, mLayerdTheta_dPsi(iLayer)

   ! compute the off-diagonal elements
   if(iLayer<nSoil) d_p1(iLayer)   = -dt_dldz(iLayer)*iLayerHydCond(iLayer)
   if(iLayer>1)     d_m1(iLayer-1) = -dt_dudz(iLayer)*iLayerHydCond(iLayer-1)

   ! compute the diagonal elements
   if(bc_lower==neumann .and. iLayer==nSoil)then  ! Neumann boundary conditions for the lower-most layer
    diag(iLayer) = mLayerdTheta_dPsi(iLayer) + dt_dudz(iLayer)*iLayerHydCond(iLayer-1)
   elseif(bc_upper==neumann .and. iLayer==1)then  ! Neumann boundary conditions for the upper-most layer
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
                    x,                       & ! intent(out): new matric head vector (m)
                    rvec,                    & ! intent(out): new residual vector (J m-3)
                    f,                       & ! intent(out): new function value (J m-3 J m-3)
                    mLayerVolFracLiqNew,     & ! intent(out): new volumetric fraction of liquid water (-)
                    iLayerLiqFlux,           & ! intent(out): liquid fluxes at layer interfaces (m s-1)
                    err,message)               ! intent(out): error control
  IMPLICIT NONE
  ! input variables
  REAL(DP), DIMENSION(:), INTENT(IN) :: xold,g
  REAL(DP), DIMENSION(:), INTENT(INOUT) :: p
  REAL(DP), INTENT(IN) :: stpmax,fold
  real(dp),intent(in)           :: mLayerVolFracLiqIter(:)  ! volumetric fraction of liquid water (-)
  real(dp),intent(in)           :: mLayerVolFracIceIter(:)  ! volumetric fraction of ice (-)
  ! output variables
  REAL(DP), DIMENSION(:), INTENT(OUT) :: x,rvec
  REAL(DP), INTENT(OUT) :: f
  real(dp),intent(out)          :: mLayerVolFracLiqNew(:)   ! volumetric fraction of liquid water (-)
  real(dp),intent(out)          :: iLayerLiqFlux(0:)        ! liquid fluxes at layer interfaces (m s-1)
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
   if(x(1) > 0._dp) x(1)=0._dp 
   ! compute matric head and volumetric fraction of liquid water
   do iLayer=1,nSoil
    mLayerVolFracLiqNew(iLayer) = volFracLiq(x(iLayer),vGn_alpha,theta_res,theta_sat,vGn_n,vGn_m)
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
                    iLayerLiqFlux,         & ! intent(out): liquid water flux at the soil layer interfaces (m s-1)
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

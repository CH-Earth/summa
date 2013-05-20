module snowHydrol_module
USE nrtype
implicit none
private
public::snowHydrol
contains

 ! ************************************************************************************************
 ! new subroutine: compute liquid water flux through the snowpack
 ! ************************************************************************************************
 subroutine snowHydrol(dt,                      & ! time step (seconds)
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
 real(dp)                      :: vDrainage1                 ! vertical drainage (m s-1)
 real(dp)                      :: dflw_dliq                  ! derivative in vertical drainage (m s-1)
 real(dp)                      :: dt_dz                      ! dt/dz (s m-1)
 real(dp)                      :: lin_error                  ! linearization error [-residual] (-)
 real(dp)                      :: increment                  ! iteration increment (-)
 real(dp),parameter            :: maxVolIceContent=0.7_dp    ! maximum volumetric ice content to store water (-)
 real(dp),parameter            :: dx=1.e-8_dp                ! finite difference increment (-)
 real(dp),parameter            :: atol=1.d-6                 ! absolute iteration tolerance (-)
 integer(i4b),parameter        :: maxiter=10                 ! maximum number of iterations
 integer(i4b)                  :: jiter                      ! internal iteration index
 logical(lgt)                  :: printflag                  ! flag to print crap to the screen
 ! initialize error control
 err=0; message="snowHydrol/"

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
   mLayerPoreSpace(iLayer)  = (iden_ice/iden_water) - mLayerVolFracIce(iLayer)
   ! compute the residual volumetric liquid water content (-)
   mLayerThetaResid(iLayer) = Fcapil*mLayerPoreSpace(iLayer) * multResid
   ! **(2)** compute the initial drainage flux
   if(mLayerVolFracIceIter(iLayer) > maxVolIceContent)then
    iLayerInitLiqFluxSnow(iLayer) = iLayerLiqFluxSnow(iLayer-1)
   else
    call mw_func( mLayerVolFracLiq(iLayer),mLayerThetaResid(iLayer),mLayerPoreSpace(iLayer), & ! input
                  iLayerInitLiqFluxSnow(iLayer) )                                              ! output
   endif
  end do  ! (looping through snow layers)
 endif  ! (if the first iteration)

 ! compute fluxes at the end of the step
 do iLayer=1,nSnow
  ! compute dt/dz
  dt_dz     = dt/mLayerDepth(iLayer)
  ! compute the liquid water equivalent associated with phase change
  phseChnge = (iden_ice/iden_water)*(mLayerVolFracIceIter(iLayer) - mLayerVolFracIce(iLayer))
  ! ** allow liquid water to pass through under very high density
  if(mLayerVolFracIce(iLayer) > maxVolIceContent)then ! NOTE: use start-of-step ice content, to avoid convergence problems
   iLayerLiqFluxSnow(iLayer)   = iLayerLiqFluxSnow(iLayer-1)
   mLayerVolFracLiqNew(iLayer) = mLayerVolFracLiq(iLayer) + &
                                 (          wimplicit *(iLayerInitLiqFluxSnow(iLayer-1) - iLayerInitLiqFluxSnow(iLayer)) +        &
                                   (1._dp - wimplicit)*(iLayerLiqFluxSnow(iLayer-1)     - iLayerLiqFluxSnow(iLayer)    ) )*dt_dz  &
                                      - phseChnge
  else
  ! ** iterate
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
                 vDrainage, dflw_dliq)
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
    ! test
    if(iLayer==1 .and. printflag)then
     write(*,'(a)')      'in snowHydrol, iter, jiter, iLayer, increment, volFracLiqTrial, volFracLiqNew, iLayerInitLiqFluxSnow(iLayer-1), vDrainage, phseChnge = '
     write(*,'(3(i4,1x),10(f20.10,1x))') iter, jiter, iLayer, increment, volFracLiqTrial, volFracLiqNew, iLayerInitLiqFluxSnow(iLayer-1), vDrainage, phseChnge
    endif
    ! check for convergence
    if(abs(increment) < atol) exit
    ! get ready for next iteration
    volFracLiqTrial = volFracLiqNew
   end do  ! (end iterations)
   ! save state
   mLayerVolFracLiqNew(iLayer) = volFracLiqNew
   ! get ready to process the next snow layer
   iLayerLiqFluxSnow(iLayer) = vDrainage + dflw_dliq*increment ! second term will be zero if converge completely
  endif  ! (if ice content is so high we need the direct pass through)
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

 end subroutine snowHydrol

end module snowHydrol_module

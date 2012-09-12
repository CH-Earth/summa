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
                       mLayerTempIter,          & ! layer temperature (K)
                       mLayerVolFracLiqIter,    & ! volumetric fraction of liquid water at the current iteration (-)
                       mLayerInfilFreezeNew,    & ! increase in volumetric ice content caused by freezing infiltrating flux (-)
                       mLayerVolFracLiqNew,     & ! volumetric fraction of liquid water at the next iteration (-)
                       err,message)
 USE multiconst,only:LH_fus                                                     ! latent heat of fusion (J kg-1)
 USE multiconst,only:iden_ice,iden_water                                        ! intrinsic density of ice and water (kg m-3)
 USE data_struc,only:model_decisions                                            ! model decision structure
 USE data_struc,only:mpar_data,forc_data,mvar_data,indx_data,ix_soil,ix_snow    ! data structures
 USE var_lookup,only:iLookPARAM,iLookFORCE,iLookMVAR,iLookINDEX,iLookDECISIONS  ! named variables for structure elements
 USE snow_utils_module,only:fracliquid      ! fractional liquid water content based on temperature (snow)
 ! compute energy fluxes within the domain
 implicit none
 ! dummy variables
 real(dp),intent(inout)        :: dt                         ! time step (seconds)
 real(dp),intent(in)           :: mLayerTempIter(:)          ! layer temperature (K)
 real(dp),intent(in)           :: mLayerVolFracLiqIter(:)    ! volumetric fraction of liquid water at the current iteration (-)
 real(dp),intent(out)          :: mLayerInfilFreezeNew(:)    ! increase in volumetric ice content caused by freezing infiltrating flux (-)
 real(dp),intent(out)          :: mLayerVolFracLiqNew(:)     ! volumetric fraction of liquid water at the next iteration (-)
 integer(i4b),intent(out)      :: err                        ! error code
 character(*),intent(out)      :: message                    ! error message
 ! local pointers to model parameters
 real(dp),pointer              :: Fcapil                     ! capillary retention as a fraction of the total pore volume (-)
 real(dp),pointer              :: k_snow                     ! hydraulic conductivity of snow (m s-1), 0.0055 = approx. 20 m/hr, from UEB
 real(dp),pointer              :: mw_exp                     ! exponent for meltwater flow (-)
 real(dp),pointer              :: snowfrz_scale              ! scaling parameter in the freezing curve for snow (K-1)
 ! local pointers to model forcing data
 real(dp),pointer              :: rainfall                   ! rainfall (kg m-2 s-1) 
 ! local pointers to model state variables 
 real(dp),pointer              :: mLayerDepth(:)             ! depth of the layer (m)
 real(dp),pointer              :: mLayerVolFracLiq(:)        ! volumetric fraction of liquid water in each snow layer (-)
 real(dp),pointer              :: mLayerRelSat(:)            ! relative saturation in each snow layer (-)
 ! local pointers to model diagnostic variables
 real(dp),pointer              :: mLayerVolFracAir(:)        ! volumetric fraction of air at the current iteration (-)
 real(dp),pointer              :: mLayerPoreSpace(:)         ! pore space in each snow layer at the current iteration (-)
 real(dp),pointer              :: iLayerLiquidWaterFlux(:)   ! vertical liquid water flux at layer interfaces (m s-1)
 ! local pointers to model index variables
 integer(i4b),pointer          :: layerType(:)               ! type of the layer (ix_soil or ix_snow)
 ! define local variables
 character(LEN=256)            :: cmessage                   ! error message of downwind routine
 integer(i4b)                  :: iLayer                     ! layer index
 integer(i4b)                  :: nSnow                      ! number of snow layers
 real(dp),parameter            :: cfrac=0.5_dp               ! critical fraction of liquid water to satisfy thermal requirements (-)
 real(dp)                      :: fracliq                    ! fraction of liquid water (-) 
 real(dp),pointer              :: ifcFlux                    ! vertical flux of liquid water at the layer interface (m s-1)
 real(dp)                      :: liqFlux                    ! unfrozen liquid water flux (m s-1)
 real(dp)                      :: vol_liq                    ! volumetric fraction of liquid water associated with infiltrating flux (-)
 real(dp)                      :: dzeqv                      ! equivalent depth of drainable pore space (m)
 real(dp)                      :: dt_dzeqv                   ! dt/dzeqv (s m-1)
 real(dp)                      :: relsatIter                 ! relative saturation of each layer for the current iteration (-)
 real(dp)                      :: relsatNew                  ! relative saturation of each layer at the next iteration (-)
 ! initialize error control
 err=0; message="liquidflow/"
 
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
 snowfrz_scale => mpar_data%var(iLookPARAM%snowfrz_scale)    ! scaling parameter in the freezing curve for snow (K-1)

 ! assign pointers to model forcing data
 rainfall => mvar_data%var(iLookMVAR%scalarRainfall)%dat(1)  ! computed rainfall rate (kg m-2 s-1)

 ! assign pointers to model state variables
 mLayerDepth           => mvar_data%var(iLookMVAR%mLayerDepth)%dat(1:nSnow)         ! depth of the layer (m)
 mLayerVolFracLiq      => mvar_data%var(iLookMVAR%mLayerVolFracLiq)%dat(1:nSnow)    ! volumetric fraction of liquid water in each snow layer (-)
 mLayerRelsat          => mvar_data%var(iLookMVAR%mLayerRelsat)%dat                 ! relative saturation in each snow layer 

 ! assign pointers to model diagnostic variables
 mLayerVolFracAir      => mvar_data%var(iLookMVAR%mLayerVolFracAir)%dat(1:nSnow)    ! volumetric fraction of air at the current iteration (-)
 mLayerPoreSpace       => mvar_data%var(iLookMVAR%mLayerPoreSpace)%dat              ! pore space in each snow layer at the current iteration (-)
 iLayerLiquidWaterFlux => mvar_data%var(iLookMVAR%iLayerLiquidWaterFlux)%dat        ! vertical liquid water flux at layer interfaces (m s-1)

 ! define the flux at the upper boundary (m s-1)
 iLayerLiquidWaterFlux(0) = rainfall/iden_water

 ! check the meltwater exponent is >=1
 if(mw_exp<1._dp)then; err=20; message=trim(message)//'meltwater exponent < 1'; return; endif

 do iLayer=1,nSnow
  ! *******************************************************************************************************
  ! part one: compute the change in volumetric ice content and the infiltrating flux
  ! *******************************************************************************************************
  if(iLayerLiquidWaterFlux(iLayer-1)>0._dp)then
   ! short-cut to the vertical flux of liquid water at the layer interface (m s-1) 
   ifcFlux =>iLayerLiquidWaterFlux(iLayer-1)
   ! estimate the volumetric fraction of liquid water associated with infiltrating flux (-)
   vol_liq = dt*ifcFlux/mLayerDepth(iLayer)
   ! compute the temperature when the fraction of liquid water is at cfrac
   crit_Tk = min(mLayerTempIter(iLayer),templiquid(cfrac,snowfrz_scale))
   ! compute the volumetric fraction of water required to raise temperatures to Tfreeze (-)
   critvol = ((mLayerVolFracIceIter(iLayer)*Cp_ice)/(LH_fus*iden_ice))*(crit_Tk - mLayerTempIter(iLayer))
   ! compute the fraction of the infiltrating flux that freezes (-)
   fracvol = min(vol_liq/critvol, 1._dp)
   ! compute the infiltrating flux of liquid water (m s-1)
   liqFlux = (1._dp - fracvol)*ifcFlux
 

 


   ! compute the infiltrating flux
   fracliq = fracliquid(mLayerTempIter(iLayer),snowfrz_scale)  ! fraction of liquid water (-)
   ifcFlux =>iLayerLiquidWaterFlux(iLayer-1)                   ! vertical flux of liquid water at the layer interface (m s-1)
   liqFlux = fracliq*ifcFlux                                   ! unfrozen liquid water flux (m s-1)
   vol_liq = dt*liqFlux/mLayerDepth(iLayer)                    ! volumetric fraction of liquid water associated with infiltrating flux (-)
   ! adjust time step if necessary
   if(vol_liq>0.1_dp*mLayerVolFracAir(iLayer))then
    err=-10; dt=0.1_dp*mLayerVolFracAir(iLayer)*mLayerDepth(iLayer)/liqFlux; return
    print*, 'adjusting time step...', dt, mLayerVolFracAir(iLayer)*mLayerDepth(iLayer), liqFlux 
   endif
   ! compute the change in volumetric ice content associated with freezing the infiltrating flux (-)
   mLayerInfilFreezeNew(iLayer) = (iden_water/iden_ice)*dt*(ifcFlux - liqFlux)/mLayerDepth(iLayer)  ! change in volumetric ice content (-)
   if(mLayerInfilFreezeNew(iLayer)>mLayerVolFracAir(iLayer)) mLayerInfilFreezeNew(iLayer)=mLayerVolFracAir(iLayer)
  else
   liqFlux = 0._dp
   vol_liq = 0._dp
   mLayerInfilFreezeNew(iLayer) = 0._dp
  endif
  ! *******************************************************************************************************
  ! part 2: compute fluid flow
  ! *******************************************************************************************************
  ! check if there is active drainage
  if(mLayerVolFracLiq(iLayer)+vol_liq > Fcapil*mLayerPoreSpace(iLayer))then 
   ! compute dt/dzeqv -- with dzeqv = equivalent depth of drainable pore space
   dzeqv    = (mLayerDepth(iLayer)*(mLayerPoreSpace(iLayer) - Fcapil)) ! equivalent depth of drainable pore space (m)
   dt_dzeqv = dt/dzeqv
   ! compute the initial estimate of the relative saturation -- note, use mLayerVolFracLiqIter
   if(mLayerVolFracLiqIter(iLayer)/mLayerPoreSpace(iLayer) > Fcapil) then
    relsatIter = (mLayerVolFracLiqIter(iLayer)/mLayerPoreSpace(iLayer) - Fcapil) / (1._dp - Fcapil)
   else
    relsatIter = mLayerRelsat(iLayer) + dt_dzeqv*liqFlux
   endif
   ! compute the relative saturation used at the next iteration
   call relsatSolve(.false.,dt_dzeqv,k_snow,mw_exp,mLayerRelsat(iLayer),liqFlux,relsatIter,relsatNew,err,cmessage)
   if(err/=0)then
    print*,trim(cmessage)
    call relsatSolve(.true.,dt_dzeqv,k_snow,mw_exp,mLayerRelsat(iLayer),liqFlux,relsatIter,relsatNew,err,cmessage)
    stop
   endif
   ! compute the flux draining from the layer
   iLayerLiquidWaterFlux(iLayer) = k_snow*relsatNew**mw_exp
   ! convert to volumetric liquid water content
   mLayerVolFracLiqNew(iLayer) = (Fcapil + relsatNew*(1._dp - Fcapil))*mLayerPoreSpace(iLayer)
  else
   ! no active drainage
   iLayerLiquidWaterFlux(iLayer) = 0._dp
   mLayerVolFracLiqNew(iLayer)   = mLayerVolFracLiq(iLayer) + vol_liq
  endif
 end do

 contains
  subroutine relsatSolve(idebug,dt_dz,k_snow,mw_exp,relsatStart,fluxIn,relsatIter,relsatNew,err,message)
  ! estimate relative saturation at the next iteration -- no need to go to convergence
  ! f = k_snow*relsatIter**mw_exp
  implicit none
  ! dummy variables
  logical(lgt),intent(in)       :: idebug                   ! debug flag
  real(dp),intent(in)           :: dt_dz                    ! time step / layer depth (s m-1)
  real(dp),intent(in)           :: k_snow                   ! hydraulic conductivity (m s-1)
  real(dp),intent(in)           :: mw_exp                   ! meltwater exponent (-)
  real(dp),intent(in)           :: relsatStart              ! relative saturation at the start of the step (-)
  real(dp),intent(in)           :: fluxIn                   ! meltwater flux input at the top of the layer (m s-1)
  real(dp),intent(in)           :: relsatIter               ! relative saturation at the current iteration (-)
  real(dp),intent(out)          :: relsatNew                ! new estimate of relative saturation
  integer(i4b),intent(out)      :: err                      ! error code
  character(*),intent(out)      :: message                  ! error message
  ! local variables
  integer(i4b)                  :: iter                     ! iteration index
  integer(i4b),parameter        :: maxiter=10               ! maximum number of iterations
  real(dp),parameter            :: converge=1.d-2           ! convergence criteria (-)
  real(dp)                      :: relsatTrial              ! trial value of relative saturation (-)
  real(dp)                      :: fluxOutTrial             ! meltwater flux out of the bottom of the layer (m s-1)
  real(dp)                      :: derivTrial               ! derivative in meltwater flux (m s-1)
  real(dp)                      :: residual                 ! residual (-)
  real(dp)                      :: delRelsat                ! iteration increment in relative saturation (-)
  real(dp)                      :: xmin,xmax                ! minimum and maximum values of relative saturation (-)
  ! initialize error control
  err=0; message='relsatSolve/'
  ! initialize bounds and relative saturation
  xmin=0._dp; xmax=1._dp
  relsatTrial  = relsatIter
  fluxOutTrial = k_snow*relsatTrial**mw_exp  ! compute the meltwater flux (m s-1)
  residual     = relsatTrial - (relsatStart + dt_dz*(FluxIn - fluxOutTrial))  ! compute the residual (-)
  ! ***** iterate
  do iter=1,maxiter
   ! compute the derivative of the meltwater flux (m s-1)
   derivTrial   = mw_exp*k_snow*relsatTrial**(mw_exp - 1._dp)
   ! compute the increment in relative saturation for the current iteration
   delRelsat   = -residual/(1._dp + dt_dz*derivTrial)
   ! compute the value of relative saturation for the next iteration
   relsatTrial  = relsatTrial + delRelsat
   ! update bounds
   if(derivTrial> 0._dp) xmin=relsatTrial
   if(derivTrial<=0._dp) xmax=relsatTrial
   ! use bi-section if outside bounds
   if(relsatTrial<xmin .or. relsatTrial>xmax) relsatTrial=0.5_dp*(xmin+xmax)
   if(idebug) write(*,'(i4,1x,3(e20.10,1x))') iter, residual, delRelsat, relsatTrial
   ! compute new flux (m s-1) and residual (-)
   fluxOutTrial = k_snow*relsatTrial**mw_exp
   residual     = relsatTrial - (relsatStart + dt_dz*(FluxIn - fluxOutTrial))
   ! exit if achieved desired tolerance
   if(residual<converge)then
    relsatNew=relsatTrial
    return
   endif
   ! check for non-convergence
   if(iter==maxiter)then; err=20; message=trim(message)//'failed to converge'; return; endif
  end do  ! iterating
  end subroutine relsatSolve

 end subroutine liquidflow

 ! ************************************************************************************************
 ! new subroutine: compute liquid water flux through the soil
 ! ************************************************************************************************
 subroutine masschange(dt,&                   ! time step (seconds)
                       mLayerMatricHeadIter,& ! matric head in each layer at the current iteration (m)
                       mLayerVolFracIceIter,& ! volumetric fraction of ice at the current iteration (-)
                       mLayerVolFracLiqIter,& ! volumetric fraction of liquid water at the current iteration (-)
                       mLayerMatricHeadNew, & ! matric head in each layer at the next iteration (m)
                       err,message)
 USE multiconst,only:iden_ice,iden_water                                        ! intrinsic density of ice and water (kg m-3)
 USE data_struc,only:model_decisions                                            ! model decision structure
 USE data_struc,only:mpar_data,forc_data,mvar_data,indx_data,ix_soil,ix_snow    ! data structures
 USE var_lookup,only:iLookPARAM,iLookFORCE,iLookMVAR,iLookINDEX,iLookDECISIONS  ! named variables for structure elements
 USE soil_utils_module,only:hydCond         ! compute hydraulic conductivity
 USE tridagSolv_module,only:tridag          ! solve tridiagonal system of equations
 ! compute energy fluxes within the domain
 implicit none
 ! dummy variables
 real(dp),intent(in)           :: dt                       ! time step (seconds)
 real(dp),intent(in)           :: mLayerMatricHeadIter(:)  ! matric head in each layer at the current iteration (m)
 real(dp),intent(in)           :: mLayerVolFracIceIter(:)  ! volumetric fraction of ice at the current iteration (-)
 real(dp),intent(in)           :: mLayerVolFracLiqIter(:)  ! volumetric fraction of liquid water at the current iteration (-)
 real(dp),intent(out)          :: mLayerMatricHeadNew(:)   ! matric head in each layer at the next iteration (m)
 integer(i4b),intent(out)      :: err                      ! error code
 character(*),intent(out)      :: message                  ! error message
 ! local pointers tp model parameters
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
 real(dp),pointer              :: snowmelt                 ! liquid water flux from the base of the snowpack (m s-1)
 ! local pointers to model variables that are constant over the simulation period 
 real(dp),pointer              :: vGn_m                    ! van Genutchen "m" parameter (-)
 ! local pointers to model state variables
 real(dp),pointer              :: mLayerVolFracIce(:)      ! volumetric fraction of ice at the start of the time step (-)
 real(dp),pointer              :: mLayerVolFracLiq(:)      ! volumetric fraction of liquid water at the start of the time step (-)
 real(dp),pointer              :: mLayerMatricHead(:)      ! matric head in each layer at the start of the time step (m)
 ! local pointers to variable short-cuts
 real(dp),pointer              :: volLatHt_fus             ! volumetric latent heat of fusion (J m-3 K-1)
 ! local pointers to model diagnostic variables
 real(dp),pointer              :: mLayerDepth(:)           ! depth of the layer (m)
 real(dp),pointer              :: mLayerHeight(:)          ! height of the layer mid-point (m)
 real(dp),pointer              :: mLayerdTheta_dPsi(:)     ! derivative in the soil water characteristic (m-1)
 ! local pointers to model index variables
 integer(i4b),pointer          :: layerType(:)             ! type of the layer (ix_soil or ix_snow)
 ! define local variables
 character(LEN=256)            :: cmessage                 ! error message of downwind routine
 integer(i4b)                  :: nSnow                    ! number of snow layers
 integer(i4b)                  :: nSoil                    ! number of soil layers
 integer(i4b)                  :: ibeg,iend                ! start and end indices of the soil layers in concatanated vector
 integer(i4b)                  :: iLayer                   ! layer index
 real(dp)                      :: q_top,q_bot              ! boundary fluxes (m s-1)
 real(dp),allocatable          :: mLayerHydCond(:)         ! hydraulic conductivity at layer mid-point (m s-1)
 real(dp),allocatable          :: iLayerHydCond(:)         ! hydraulic conductivity at layer interface (m s-1)
 real(dp),allocatable          :: relIce(:)                ! relative ice content in each layer (-)
 real(dp),allocatable          :: dz_node(:)               ! distance between the mid-point of a layer and a neighbouring layer (m)
 real(dp),allocatable          :: dt_dudz(:)               ! the dt_dudz terms for the upper interface (s m-2)
 real(dp),allocatable          :: dt_dldz(:)               ! the dt_dudz terms for the lower interface (s m-2)
 real(dp),allocatable          :: rv_temp(:)               ! common elements of the RHS vector (-)
 real(dp),allocatable          :: d_m1(:)                  ! sub-diagonal elements of the tridiagonal system (m-1)
 real(dp),allocatable          :: diag(:)                  ! diagonal elements of the tridiagonal system (m-1)
 real(dp),allocatable          :: d_p1(:)                  ! super-diagonal elements of the tridiagonal system (m-1)
 real(dp),allocatable          :: rvec(:)                  ! right-hand-side vector (-)
 integer(i4b),parameter        :: diriclet=1001            ! look-up value for diriclet boundary conditions
 integer(i4b),parameter        :: neumann=1002             ! look-up value for neumann boundary conditions
 integer(i4b)                  :: bc_upper,bc_lower        ! boundary condition index (diriclet or neumann)


 ! initialize error control
 err=0; message="masschange/"

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
 snowmelt          => mvar_data%var(iLookMVAR%iLayerLiquidWaterFlux)%dat(nSnow)  ! liquid water flux from the base of the snowpack (m s-1)

 ! assign pointers to model variables that are constant over the simulation period
 vGn_m             => mvar_data%var(iLookMVAR%scalarVGn_m)%dat(1)              ! van Genutchen "m" parameter (-)
 volLatHt_fus      => mvar_data%var(iLookMVAR%scalarvolLatHt_fus)%dat(1)       ! volumetric latent heat of fusion (J m-3)

 ! assign pointers to model state variables
 mLayerVolFracIce  => mvar_data%var(iLookMVAR%mLayerVolFracIce)%dat(ibeg:iend) ! volumetric fraction of ice in each layer (-)
 mLayerVolFracLiq  => mvar_data%var(iLookMVAR%mLayerVolFracLiq)%dat(ibeg:iend) ! volumetric fraction of liquid water in each layer (-)
 mLayerMatricHead  => mvar_data%var(iLookMVAR%mLayerMatricHead)%dat(ibeg:iend) ! matric head in each layer (m)

 ! assign pointers to model diagnostic variables
 mLayerDepth       => mvar_data%var(iLookMVAR%mLayerDepth)%dat(ibeg:iend)      ! depth of the layer (m)
 mLayerHeight      => mvar_data%var(iLookMVAR%mLayerHeight)%dat(ibeg:iend)     ! height of the layer mid-point (m)
 mLayerdTheta_dPsi => mvar_data%var(iLookMVAR%mLayerdTheta_dPsi)%dat ! (soil only) ! derivative in the soil water characteristic (m-1)


 ! allocate space for the dt/dzdz terms
 allocate(dz_node(nSoil-1),dt_dudz(nSoil),dt_dldz(nSoil),stat=err)
 if(err/=0)then; err=20; message=trim(message)//'problem allocating space for dt/dzdz vectors'; return; endif
 ! allocate space for hydraulic conductivity
 allocate(relIce(nSoil),mLayerHydCond(nSoil),iLayerHydCond(0:nSoil),stat=err)
 if(err/=0)then; err=20; message=trim(message)//'problem allocating space for hydraulic conductivity'; return; endif
 ! allocate space for temporary vectors
 allocate(rv_temp(nSoil),stat=err)
 if(err/=0)then; err=20; message=trim(message)//'problem allocating space for temporary vectors'; return; endif
 ! allocate space for the tri-diagonal vectors
 allocate(d_m1(nSoil-1),diag(nSoil),d_p1(nSoil-1),rvec(nSoil),stat=err)
 if(err/=0)then; err=20; message=trim(message)//'problem allocating space for tridag vectors'; return; endif

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
  if(nSnow==0) q_top = rainfall/iden_water
  if(nSnow>0) q_top = snowmelt
 endif
 if(bc_lower==neumann) q_bot = 0._dp
 !print*, 'rainfall = ', rainfall

 ! compute the hydraulic conductivity (m s-1)
 do iLayer=1,nSoil
  ! compute the relative volumetric fraction of ice
  relIce(iLayer) = mLayerVolFracIceIter(iLayer) / ( mLayerVolFracIceIter(iLayer) + (mLayerVolFracLiqIter(iLayer)-theta_res) )
  ! compute the hydraulic conductivity for a given layer (m s-1)
  mLayerHydCond(iLayer) = hydCond(mLayerMatricHeadIter(iLayer),k_soil,vGn_alpha,vGn_n,vGn_m) * 10._dp**(f_impede*relIce(iLayer)) 
  ! compute the hydraulic conductivity for layer interfaces (m s-1)
  if(iLayer>1) iLayerHydCond(iLayer-1) = (mLayerHydCond(iLayer-1) + mLayerHydCond(iLayer))/2._dp
 end do

 ! if diriclet, compute hydraulic conductivity at the boundaries (m s-1)
 if(bc_upper==diriclet) iLayerHydCond(0)     = hydCond(upperBoundHead,k_soil,vGn_alpha,vGn_n,vGn_m) * 10._dp**(f_impede*relIce(1))
 if(bc_lower==diriclet) iLayerHydCond(nSoil) = hydCond(lowerBoundHead,k_soil,vGn_alpha,vGn_n,vGn_m) * 10._dp**(f_impede*relIce(nSoil))

 ! if newmann, set hydraulic conductivity at boundaries to ridiculous values (not used)
 if(bc_upper==neumann) iLayerHydCond(0)       = huge(iLayerHydCond)
 if(bc_lower==neumann) iLayerHydCond(nSoil) = huge(iLayerHydCond)

 ! compute the dt/dzdz terms for each layer(s m-2)
 dz_node(1:nSoil-1) = mLayerHeight(2:nSoil) - mLayerHeight(1:nSoil-1)
 dt_dudz(2:nSoil)   = dt/(mLayerDepth(2:nSoil)*dz_node)    ! upper distance, valid 2..nSoil
 dt_dldz(1:nSoil-1) = dt/(mLayerDepth(1:nSoil-1)*dz_node)  ! lower distance, valid 1..nSoil-1

 ! compute the dt/dzdz terms at the boundaries
 dt_dudz(1)     = dt/(mLayerDepth(1)     * mLayerDepth(1))!*0.5_dp)
 dt_dldz(nSoil) = dt/(mLayerDepth(nSoil) * mLayerDepth(nSoil))!*0.5_dp)

 ! compute the common elements in the RHS vector (-)
 rv_temp = mLayerdTheta_dPsi*mLayerMatricHeadIter &
            - (iden_ice/iden_water)*(mLayerVolFracIceIter - mLayerVolFracIce) &
            - (mLayerVolFracLiqIter - mLayerVolFracLiq)

 ! assemble the tri-diagonal matrix
 do iLayer=1,nSoil
  ! compute the off-diagonal elements
  if(iLayer<nSoil) d_p1(iLayer)   = -dt_dldz(iLayer)*iLayerHydCond(iLayer)
  if(iLayer>1)       d_m1(iLayer-1) = -dt_dudz(iLayer)*iLayerHydCond(iLayer-1)
  ! ***** define the diagonal elements and the RHS vector
  if(iLayer==nSoil)then
   ! ** lower-most layer
   if(bc_lower==diriclet)then       ! (diriclet boundary condition)
    diag(iLayer) = mLayerdTheta_dPsi(iLayer) + dt_dldz(iLayer)*iLayerHydCond(iLayer) + dt_dudz(iLayer)*iLayerHydCond(iLayer-1)
    rvec(iLayer) = rv_temp(iLayer) + (dt/mLayerDepth(iLayer))*(iLayerHydCond(iLayer-1) - iLayerHydCond(iLayer)) + &
                    dt_dldz(iLayer)*iLayerHydCond(iLayer)*lowerBoundHead
   elseif(bc_lower==neumann)then    ! (Neumann boundary conditions)
    diag(iLayer) = mLayerdTheta_dPsi(iLayer) + dt_dudz(iLayer)*iLayerHydCond(iLayer-1)
    rvec(iLayer) = rv_temp(iLayer) + (dt/mLayerDepth(iLayer))*(iLayerHydCond(iLayer-1) - q_bot)
   else
    err=50; message=trim(message)//'cannot identify the type of boundary conditions' 
   endif
  elseif(iLayer==1)then 
   ! ** upper-most layer
   if(bc_upper==diriclet)then       ! (diriclet boundary condition)
    diag(iLayer) = mLayerdTheta_dPsi(iLayer) + dt_dldz(iLayer)*iLayerHydCond(iLayer) + dt_dudz(iLayer)*iLayerHydCond(iLayer-1)
    rvec(iLayer) = rv_temp(iLayer) + (dt/mLayerDepth(iLayer))*(iLayerHydCond(iLayer-1) - iLayerHydCond(iLayer)) + &
                    dt_dudz(iLayer)*iLayerHydCond(iLayer-1)*upperBoundHead
   elseif(bc_upper==neumann)then    ! (Neumann boundary conditions)
    diag(iLayer) = mLayerdTheta_dPsi(iLayer) + dt_dldz(iLayer)*iLayerHydCond(iLayer)
    rvec(iLayer) = rv_temp(iLayer) + (dt/mLayerDepth(iLayer))*(q_top - iLayerHydCond(iLayer))
   else
    err=50; message=trim(message)//'cannot identify the type of boundary conditions'
   endif
  else
   ! ** intermediate layers
   diag(iLayer) = mLayerdTheta_dPsi(iLayer) + dt_dudz(iLayer)*iLayerHydCond(iLayer-1) + dt_dldz(iLayer)*iLayerHydCond(iLayer)
   rvec(iLayer) = rv_temp(iLayer) + (dt/mLayerDepth(iLayer))*(iLayerHydCond(iLayer-1) - iLayerHydCond(iLayer))
  endif
 end do  ! (looping through layers)

 ! check the values
 !print*, 'tridiag'
 !do iLayer=1,nSoil
 ! if(iLayer==1)                   write(*,'(i4,1x,4(e20.8,1x))'), ilayer, 0._dp,          diag(iLayer), d_p1(iLayer), rvec(iLayer)
 ! if(iLayer>1 .and. iLayer<nSoil) write(*,'(i4,1x,4(e20.8,1x))'), ilayer, d_m1(iLayer),   diag(iLayer), d_p1(iLayer), rvec(iLayer)
 ! if(iLayer==nSoil)               write(*,'(i4,1x,4(e20.8,1x))'), ilayer, d_m1(iLayer-1), diag(iLayer), 0._dp,        rvec(iLayer)
 !end do

 ! solve the tridiagonal system of equations -- returns mLayerMatricHeadNew
 call tridag(d_m1,diag,d_p1,rvec,mLayerMatricHeadNew,err,cmessage)
 if(err/=0)then; message=trim(message)//trim(cmessage); return; endif

 ! compute the boundary fluxes
 if(bc_upper==diriclet)&
  q_top = iLayerHydCond(0)*(upperBoundHead - mLayerMatricHeadNew(1))/(0.5_dp*mLayerDepth(1)) + iLayerHydCond(0)
 if(bc_lower==diriclet)&
  q_bot = iLayerHydCond(nSoil)*(mLayerMatricHeadNew(nSoil) - lowerBoundHead)/(0.5_dp*mLayerDepth(nSoil)) + iLayerHydCond(nSoil)
 !print*, 'change in mass = ', q_top - q_bot

 ! deallocate vectors
 deallocate(dz_node,dt_dudz,dt_dldz,rv_temp,d_m1,diag,d_p1,rvec,stat=err)
 if(err/=0)then; err=30; message=trim(message)//'problem deallocating vectors for tridag solution'; return; endif

 end subroutine masschange

end module liquidflux_module

module energyflux_module
USE nrtype
! physical constants
USE multiconst,only:&
                    sigma,       & ! Stefan Boltzman constant      (W m-2 K-4)
                    Em_Sno,      & ! emissivity of snow            (-)
                    Cp_air,      & ! specific heat of air          (J kg-1 K-1)
                    LH_fus,      & ! latent heat of fusion         (J kg-1)
                    LH_vap,      & ! latent heat of vaporization   (J kg-1)
                    LH_sub,      & ! latent heat of sublimation    (J kg-1)
                    Tfreeze,     & ! freezing point of pure water  (K)
                    iden_air,    & ! intrinsic density of air      (kg m-3)
                    iden_ice,    & ! intrinsic density of ice      (kg m-3)
                    iden_water,  & ! intrinsic density of water    (kg m-3)
                    lambda_air,  & ! thermal conductivity of air   (J s-1 m-1)
                    lambda_ice,  & ! thermal conductivity of ice   (J s-1 m-1)
                    lambda_soil, & ! thermal conductivity of soil  (J s-1 m-1)
                    lambda_water   ! thermal conductivity of water (J s-1 m-1)
implicit none
private
public::diagn_evar
public::tempchange
public::phsechange
! local parameters
real(dp),parameter            :: RHsurf=1._dp             ! relative humidity of the surface (-)
! look-up values for the numerical method
integer(i4b)                  :: num_method               ! numerical method
integer(i4b),parameter        :: itertive=1001            ! named index for the iterative method
integer(i4b),parameter        :: non_iter=1002            ! named index for the non-iterative methof
integer(i4b),parameter        :: itersurf=1003            ! named index for the case where iterate only on the surface energy balance
! look-up values for method used to compute derivative
integer(i4b)                  :: fDerivMeth               ! index for the method used to calculate flux derivatives
integer(i4b),parameter        :: numerical=1001           ! look-up value for numerical solution
integer(i4b),parameter        :: analytical=1002          ! look-up value for analytical solution
! number of soil and snow layers
integer(i4b)                  :: nSoil                    ! number of soil layers
integer(i4b)                  :: nSnow                    ! number of snow layers
integer(i4b)                  :: nLayers                  ! total number of layers
contains


 ! ************************************************************************************************
 ! new subroutine: compute change in temperature over the time step
 ! ************************************************************************************************
 subroutine tempchange(dt,                   & ! intent(in): time step (seconds)
                       iter,                 & ! intent(in): current iteration count
                       mLayerTempIter,       & ! intent(in): trial temperature at the current iteration (K)
                       mLayerVolFracIceIter, & ! intent(in): volumetric fraction of ice at the current iteration (-)
                       mLayerVolFracLiqIter, & ! intent(in): volumetric fraction of liquid water at the current iteration (-)
                       mLayerMatricHeadIter, & ! intent(in): matric head at the current iteration (m)
                       mLayerTempDiffOld,    & ! intent(in): iteration increment for temperature at the last iteration (K)
                       mLayerTempDiff,       & ! intent(out): iteration increment for temperature (K)
                       mLayerTempNew,        & ! intent(out): new temperature (K)
                       mLayerVolFracIceNew,  & ! intent(out): new volumetric fraction of ice (-)
                       mLayerVolFracLiqNew,  & ! intent(out): new volumetric fraction of liquid water (-)
                       mLayerMatricHeadNew,  & ! intent(out): new matric head (m)
                       err,message)            ! intent(out): error control
 USE data_struc,only:mpar_data,forc_data,mvar_data,indx_data,ix_soil,ix_snow    ! data structures
 USE var_lookup,only:iLookPARAM,iLookFORCE,iLookMVAR,iLookINDEX                 ! named variables for structure elements
 ! compute change in temperature over the time step
 implicit none
 ! input
 real(dp),intent(in)           :: dt                       ! time step (seconds)
 integer(i4b),intent(in)       :: iter                     ! iteration count
 real(dp),intent(in)           :: mLayerTempIter(:)        ! trial temperature at the current iteration (K)
 real(dp),intent(in)           :: mLayerVolFracIceIter(:)  ! volumetric fraction of ice at the current iteration (-)
 real(dp),intent(in)           :: mLayerVolFracLiqIter(:)  ! volumetric fraction of liquid water at the current iteration (-)
 real(dp),intent(in)           :: mLayerMatricHeadIter(:)  ! matric head at the current iteration (m)
 real(dp),intent(in)           :: mLayerTempDiffOld(:)     ! iteration increment for temperature at the last iteration (K) 
 ! output
 real(dp),intent(out)          :: mLayerTempDiff(:)        ! iteration increment for temperature (K) 
 real(dp),intent(out)          :: mLayerTempNew(:)         ! new temperature (K)
 real(dp),intent(out)          :: mLayerVolFracIceNew(:)   ! new volumetric fraction of ice (-)
 real(dp),intent(out)          :: mLayerVolFracLiqNew(:)   ! new volumetric fraction of liquid water (-)
 real(dp),intent(out)          :: mLayerMatricHeadNew(:)   ! new matric head (m)
 integer(i4b),intent(out)      :: err                      ! error code
 character(*),intent(out)      :: message                  ! error message
 ! initialize error control
 err=0; message="tempchange/"

 ! define number of layers
 nSoil   = count(indx_data%var(iLookINDEX%layerType)%dat == ix_soil)  ! number of soil layers
 nSnow   = count(indx_data%var(iLookINDEX%layerType)%dat == ix_snow)  ! number of snow layers
 nLayers = indx_data%var(iLookINDEX%nLayers)%dat(1)                   ! total number of layers

 ! identify the numerical method
 select case(trim(model_decisions(iLookDECISIONS%num_method)%decision))
  case('itertive'); num_method=itertive  ! iterative
  case('non_iter'); num_method=non_iter  ! non-iterative
  case('itersurf'); num_method=itersurf  ! iterate only on the surface energy balance
  case default
   err=10; message=trim(message)//"unknown option for the numerical method"; return
 end select

 ! identify the method used to calculate flux derivatives
 select case(trim(model_decisions(iLookDECISIONS%fDerivMeth)%decision))
  case('numericl'); fDerivMeth=numerical
  case('analytic'); fDerivMeth=analytical
  case default
   err=10; message=trim(message)//"unknown method used to calculate flux derivatives [option="//trim(model_decisions(iLookDECISIONS%fDerivMeth)%decision)//"]"; return
 end select


 ! *****
 ! wrapper for the temperature change sub-routine...
 ! *************************************************
 ! used as a trick to both
 !  1) save typing; and
 !  2) "protect" variables through the use of the intent attribute
 call tempchange_muster(&
                        ! input variables from tempchange routine
                        dt,                                                & ! intent(in): time step (seconds)
                        iter,                                              & ! intent(in): current iteration count
                        mLayerTempIter,                                    & ! intent(in): trial temperature at the current iteration (K)
                        mLayerVolFracIceIter,                              & ! intent(in): volumetric fraction of ice at the current iteration (-)
                        mLayerVolFracLiqIter,                              & ! intent(in): volumetric fraction of liquid water at the current iteration (-)
                        mLayerMatricHeadIter,                              & ! intent(in): matric head at the current iteration (m)
                        mLayerTempDiffOld,                                 & ! intent(in): iteration increment for temperature at the last iteration (K)
                        ! index variables
                        indx_data%var(iLookINDEX%nLayers)%dat(1),          & ! intent(in): number of layers
                        indx_data%var(iLookINDEX%layerType)%dat,           & ! intent(in): layer type (ix_soil or ix_snow)
                        ! general model parameters
                        mpar_data%var(iLookPARAM%mheight),                 & ! intent(in): measurement height (m)
                        mpar_data%var(iLookPARAM%wimplicit),               & ! intent(in): weight assigned to start-of-step fluxes (-)
                        mpar_data%var(iLookPARAM%snowfrz_scale),           & ! intent(in): scaling parameter for the snow freezing curve (K-1)
                        mpar_data%var(iLookPARAM%lowerBoundTemp),          & ! intent(in): temperature of the lower boundary (K)
                        ! soil parameters
                        mpar_data%var(iLookPARAM%vGn_alpha),               & ! intent(in): van Genutchen "alpha" parameter (m-1)
                        mpar_data%var(iLookPARAM%vGn_n),                   & ! intent(in): van Genutchen "n" parameter (-)
                        mpar_data%var(iLookPARAM%theta_sat),               & ! intent(in): soil porosity (-)
                        mpar_data%var(iLookPARAM%theta_res),               & ! intent(in): soil residual volumetric water content (-)
                        ! vegetation parameters
                        mpar_data%var(iLookPARAM%LAI),                     & ! intent(in): leaf area index (m2 m-2)
                        mpar_data%var(iLookPARAM%minStomatalResist),       & ! intent(in): minimum stomatal resistance (s m-1)
                        mpar_data%var(iLookPARAM%plantWiltPsi),            & ! intent(in): critical matric head when stomatal resitance 2 x min (m)
                        mpar_data%var(iLookPARAM%plantWiltExp),            & ! intent(in): empirical exponent in plant wilting factor expression (-)
                        ! model variables that are constant over the simulation period
                        mvar_data%var(iLookMVAR%scalarVGn_m)%dat(1),       & ! intent(in): van Genutchen "m" parameter (-)
                        mvar_data%var(iLookMVAR%scalarKappa)%dat(1),       & ! intent(in): constant in the freezing curve function (m K-1)
                        mvar_data%var(iLookMVAR%mLayerRootDensity)%dat,    & ! intent(in): fraction of roots in each soil layer (-)
                        ! model forcing variables
                        forc_data%var(iLookFORCE%sw_down),                 & ! intent(in): downward shortwave radiation (W m-2)
                        forc_data%var(iLookFORCE%lw_down),                 & ! intent(in): downward longwave radiation (W m-2)
                        forc_data%var(iLookFORCE%airtemp),                 & ! intent(in): air temperature at 2 meter height (K)
                        forc_data%var(iLookFORCE%windspd),                 & ! intent(in): wind speed at 10 meter height (m s-1)
                        forc_data%var(iLookFORCE%airpres),                 & ! intent(in): air pressure at 2 meter height (Pa)
                        forc_data%var(iLookFORCE%spechum),                 & ! intent(in): specific humidity at 2 meter height (g g-1)
                        ! model state variables
                        mvar_data%var(iLookMVAR%scalarAlbedo)%dat(1),      & ! intent(in): surface albedo (-)
                        mvar_data%var(iLookMVAR%mLayerTemp)%dat,           & ! intent(in): temperature of each layer (K)
                        mvar_data%var(iLookMVAR%mLayerMatricHead)%dat,     & ! intent(in): matric head in each layer (m)
                        mvar_data%var(iLookMVAR%mLayerVolFracIce)%dat,     & ! intent(in): volumetric fraction of ice in each layer (-)
                        mvar_data%var(iLookMVAR%mLayerDepth)%dat,          & ! intent(in): depth of each layer (m)
                        ! model diagnostic variables (input)
                        mvar_data%var(iLookMVAR%mLayerHeight)%dat,         & ! intent(in): height at the mid-point of each layer (m)
                        mvar_data%var(iLookMVAR%iLayerHeight)%dat,         & ! intent(in): height at the interface of each layer (m)
                        mvar_data%var(iLookMVAR%mLayerThermalC)%dat,       & ! intent(in): thermal conductivity at the mid-point of each layer (W m-1 K-1)
                        mvar_data%var(iLookMVAR%iLayerThermalC)%dat,       & ! intent(in): thermal conductivity at the interface of each layer (W m-1 K-1)
                        mvar_data%var(iLookMVAR%mLayerVolHtCapBulk)%dat,   & ! intent(in): bulk volumetric heat capacity (J m-3 K-1)
                        ! model diagnostic variables (output)
                        mvar_data%var(iLookMVAR%mLayerTcrit)%dat,          & ! intent(out): critical soil temperature above which all water is unfrozen (K)
                        mvar_data%var(iLookMVAR%mLayerdTheta_dTk)%dat,     & ! intent(out): derivative in the freezing curve (K-1)
                        mvar_data%var(iLookMVAR%mLayerTranspireLim)%dat,   & ! intent(out): soil moist & veg limit on transpiration for each layer (-) 
                        mvar_data%var(iLookMVAR%mLayerInitTranspire)%dat,  & ! intent(out): transpiration loss from each soil layer at the start of the step (m s-1)
                        mvar_data%var(iLookMVAR%mLayerTranspire)%dat,      & ! intent(out): transpiration loss from each soil layer (m s-1)
                        mvar_data%var(iLookMVAR%iLayerInitNrgFlux)%dat,    & ! intent(out): energy flux at layer interfaces at the start of the time step (W m-2)
                        mvar_data%var(iLookMVAR%iLayerNrgFlux)%dat,        & ! intent(out): energy flux at layer interfaces at the end of the time step (W m-2)
                        ! diagnostic scalar variables (output)
                        mvar_data%var(iLookMVAR%scalarTranspireLim)%dat(1),& ! intent(out): aggregate soil moist & veg limit on transpiration, weighted by root density (-)
                        mvar_data%var(iLookMVAR%scalarPotentialET)%dat(1), & ! intent(out): potential ET (kg m-2 s-1)
                        mvar_data%var(iLookMVAR%scalarMassLiquid)%dat(1),  & ! intent(out): transpiration (kg m-2 s-1)
                        mvar_data%var(iLookMVAR%scalarMassSolid)%dat(1),   & ! intent(out): sublimation/frost (kg m-2 s-1)
                        mvar_data%var(iLookMVAR%scalarSenHeat)%dat(1),     & ! intent(out): sensible heat flux at the surface (W m-2)
                        mvar_data%var(iLookMVAR%scalarLatHeat)%dat(1),     & ! intent(out): latent heat flux at the surface (W m-2)
                        mvar_data%var(iLookMVAR%scalarExCoef)%dat(1),      & ! intent(out): turbulent exchange coefficient (-)
                        mvar_data%var(iLookMVAR%scalarExSen)%dat(1),       & ! intent(out): exchange factor for sensible heat (J m-2 s-1 K-1)
                        mvar_data%var(iLookMVAR%scalarExLat)%dat(1),       & ! intent(out): exchange factor for latent heat (J m-2 s-1)
                        ! output variables from tempchange subroutine
                        mLayerTempDiff,                                    & ! intent(out): iteration increment for temperature (K)
                        mLayerTempNew,                                     & ! intent(out): new temperature (K)
                        mLayerVolFracIceNew,                               & ! intent(out): new volumetric fraction of ice (-)
                        mLayerVolFracLiqNew,                               & ! intent(out): new volumetric fraction of liquid water (-)
                        mLayerMatricHeadNew,                               & ! intent(out): new matric head (m)
                        err,message)                                         ! intent(out): error control

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
 ! private subroutine: wrapper for the temperature change subroutine
 ! ************************************************************************************************
 subroutine tempchange_muster(&
                              ! input variables from tempchange routine
                              dt,                         & ! intent(in): time step (seconds)
                              iter,                       & ! intent(in): current iteration count
                              mLayerTempIter,             & ! intent(in): trial temperature at the current iteration (K)
                              mLayerVolFracIceIter,       & ! intent(in): volumetric fraction of ice at the current iteration (-)
                              mLayerVolFracLiqIter,       & ! intent(in): volumetric fraction of liquid water at the current iteration (-)
                              mLayerMatricHeadIter,       & ! intent(in): matric head at the current iteration (m)
                              mLayerTempDiffOld,          & ! intent(in): iteration increment for temperature at the last iteration (K)
                              ! index variables
                              nLayers,                    & ! intent(in): number of layers
                              layerType,                  & ! intent(in): layer type (ix_soil or ix_snow)
                              ! general model parameters
                              mheight,                    & ! intent(in): measurement height (m)
                              wimplicit,                  & ! intent(in): weight assigned to start-of-step fluxes (-)
                              snowfrz_scale,              & ! intent(in): scaling parameter for the snow freezing curve (K-1)
                              lowerBoundTemp,             & ! intent(in): temperature of the lower boundary (K)
                              ! soil parameters
                              vGn_alpha,                  & ! intent(in): van Genutchen "alpha" parameter (m-1)
                              vGn_n,                      & ! intent(in): van Genutchen "n" parameter (-)
                              theta_sat,                  & ! intent(in): soil porosity (-)
                              theta_res,                  & ! intent(in): soil residual volumetric water content (-)
                              ! vegetation parameters
                              LAI,                        & ! intent(in): leaf area index (m2 m-2)
                              minStomatalResist,          & ! intent(in): minimum stomatal resistance (s m-1)
                              plantWiltPsi,               & ! intent(in): critical matric head when stomatal resitance 2 x min (m)
                              plantWiltExp,               & ! intent(in): empirical exponent in plant wilting factor expression (-)
                              ! model variables that are constant over the simulation period
                              vGn_m,                      & ! intent(in): van Genutchen "m" parameter (-)
                              kappa,                      & ! intent(in): constant in the freezing curve function (m K-1)
                              mLayerRootDensity,          & ! intent(in): fraction of roots in each soil layer (-)
                              ! model forcing variables
                              sw_down,                    & ! intent(in): downward shortwave radiation (W m-2)
                              lw_down,                    & ! intent(in): downward longwave radiation (W m-2)
                              airtemp,                    & ! intent(in): air temperature at 2 meter height (K)
                              windspd,                    & ! intent(in): wind speed at 10 meter height (m s-1)
                              airpres,                    & ! intent(in): air pressure at 2 meter height (Pa)
                              spechum,                    & ! intent(in): specific humidity at 2 meter height (g g-1)
                              ! model state variables
                              scalarAlbedo,               & ! intent(in): surface albedo (-)
                              mLayerTemp,                 & ! intent(in): temperature of each layer (K)
                              mLayerMatricHead,           & ! intent(in): matric head in each layer (m)
                              mLayerVolFracIce,           & ! intent(in): volumetric fraction of ice in each layer (-)
                              mLayerDepth,                & ! intent(in): depth of each layer (m)
                              ! model diagnostic variables (input)
                              mLayerHeight,               & ! intent(in): height at the mid-point of each layer (m)
                              iLayerHeight,               & ! intent(in): height at the interface of each layer (m)
                              mLayerThermalC,             & ! intent(in): thermal conductivity at the mid-point of each layer (W m-1 K-1)
                              iLayerThermalC,             & ! intent(in): thermal conductivity at the interface of each layer (W m-1 K-1)
                              mLayerVolHtCapBulk,         & ! intent(in): bulk volumetric heat capacity (J m-3 K-1)
                              ! model diagnostic variables (output)
                              mLayerTcrit,                & ! intent(out): critical soil temperature above which all water is unfrozen (K)
                              mLayerdTheta_dTk,           & ! intent(out): derivative in the freezing curve (K-1)
                              mLayerTranspireLim,         & ! intent(out): soil moist & veg limit on transpiration for each layer (-) 
                              mLayerInitTranspire,        & ! intent(out): transpiration loss from each soil layer at the start of the step (m s-1)
                              mLayerTranspire,            & ! intent(out): transpiration loss from each soil layer (m s-1)
                              iLayerInitNrgFlux,          & ! intent(out): energy flux at layer interfaces at the start of the time step (W m-2)
                              iLayerNrgFlux,              & ! intent(out): energy flux at layer interfaces at the end of the time step (W m-2)
                              ! diagnostic scalar variables (output)
                              scalarTranspireLim,         & ! intent(out): aggregate soil moist & veg limit on transpiration, weighted by root density (-)
                              scalarPotentialET,          & ! intent(out): potential ET (kg m-2 s-1)
                              scalarMassLiquid,           & ! intent(out): transpiration (kg m-2 s-1)
                              scalarMassSolid,            & ! intent(out): sublimation/frost (kg m-2 s-1)
                              scalarSenHeat,              & ! intent(out): sensible heat flux at the surface (W m-2)
                              scalarLatHeat,              & ! intent(out): latent heat flux at the surface (W m-2)
                              scalarExCoef,               & ! intent(out): turbulent exchange coefficient (-)
                              scalarExSen,                & ! intent(out): exchange factor for sensible heat (J m-2 s-1 K-1)
                              scalarExLat,                & ! intent(out): exchange factor for latent heat (J m-2 s-1)
                              ! output variables from tempchange subroutine
                              mLayerTempDiff,             & ! intent(out): iteration increment for temperature (K)
                              mLayerTempNew,              & ! intent(out): new temperature (K)
                              mLayerVolFracIceNew,        & ! intent(out): new volumetric fraction of ice (-)
                              mLayerVolFracLiqNew,        & ! intent(out): new volumetric fraction of liquid water (-)
                              mLayerMatricHeadNew,        & ! intent(out): new matric head (m)
                              err,message)                  ! intent(out): error control
 ! compute change in temperature over the time step
 ! model decisions
 USE data_struc,only:model_decisions                        ! model decision structure
 USE data_struc,only:ix_soil,ix_snow                        ! names variables for snow and soil
 USE var_lookup,only:iLookDECISIONS                         ! named variables for elements of the decision structure
 ! utility modules
 USE snow_utils_module,only:dFracLiq_dTk                    ! differentiate the freezing curve w.r.t. temperature (snow)
 USE soil_utils_module,only:dTheta_dTk                      ! differentiate the freezing curve w.r.t. temperature (soil)
 USE conv_funcs_module,only:relhm2sphm                      ! compute specific humidity 
 USE tridagSolv_module,only:tridag                          ! solve tridiagonal system of equations
 implicit none
 ! input variables from the tempchange subroutine
 real(dp),intent(in)            :: dt                       ! time step (seconds)
 integer(i4b),intent(in)        :: iter                     ! iteration count
 real(dp),intent(in)            :: mLayerTempIter(:)        ! trial temperature at the current iteration (K)
 real(dp),intent(in)            :: mLayerVolFracIceIter(:)  ! volumetric fraction of ice at the current iteration (-)
 real(dp),intent(in)            :: mLayerVolFracLiqIter(:)  ! volumetric fraction of liquid water at the current iteration (-)
 real(dp),intent(in)            :: mLayerMatricHeadIter(:)  ! matric head at the current iteration (m)
 real(dp),intent(in)            :: mLayerTempDiffOld(:)     ! iteration increment for temperature at the last iteration (K) 
 ! model index variables
 integer(i4b),intent(in)        :: nLayers                  ! number of layers
 integer(i4b),intent(in)        :: layerType(:)             ! type of the layer (ix_soil or ix_snow)
 ! general model parameters
 real(dp),intent(in)            :: mheight                  ! measurement height (m)
 real(dp),intent(in)            :: wimplicit                ! weight assigned to start-of-step fluxes (-)
 real(dp),intent(in)            :: snowfrz_scale            ! scaling parameter for the snow freezing curve (K-1)
 real(dp),intent(in)            :: lowerBoundTemp           ! temperature of the lower boundary (K)
 ! soil parameters
 real(dp),intent(in)            :: vGn_alpha                ! van Genutchen "alpha" parameter
 real(dp),intent(in)            :: vGn_n                    ! van Genutchen "n" parameter
 real(dp),intent(in)            :: theta_sat                ! soil porosity (-)
 real(dp),intent(in)            :: theta_res                ! soil residual volumetric water content (-)
 ! vegetation parameters
 real(dp),intent(in)            :: LAI                      ! leaf area index (m2 m-2)
 real(dp),intent(in)            :: minStomatalResist        ! minimum stomatal resistance (s m-1)
 real(dp),intent(in)            :: plantWiltPsi             ! critical matric head when stomatal resitance 2 x min (m)
 real(dp),intent(in)            :: plantWiltExp             ! empirical exponent in plant wilting factor expression (-)
 ! derived model variables that are constant over the simulation period
 real(dp),intent(in)            :: vGn_m                    ! van Genutchen "m" parameter (-)
 real(dp),intent(in)            :: kappa                    ! constant in the freezing curve function (m K-1)
 real(dp),intent(in)            :: mLayerRootDensity(:)     ! fraction of roots in each soil layer (-)
 ! model forcing data
 real(dp),intent(in)            :: sw_down                  ! downward shortwave radiation (W m-2)
 real(dp),intent(in)            :: lw_down                  ! downward longwave radiation (W m-2)
 real(dp),intent(in)            :: airtemp                  ! air temperature at 2 meter height (K)
 real(dp),intent(in)            :: windspd                  ! wind speed at 10 meter height (m s-1)
 real(dp),intent(in)            :: airpres                  ! air pressure at 2 meter height (Pa)
 real(dp),intent(in)            :: spechum                  ! specific humidity at 2 meter height (g g-1)
 ! model state variables
 real(dp),intent(in)            :: scalarAlbedo             ! surface albedo (-)
 real(dp),intent(in)            :: mLayerTemp(:)            ! temperature of each layer (K)
 real(dp),intent(in)            :: mLayerMatricHead(:)      ! matric head in each layer (m)
 real(dp),intent(in)            :: mLayerVolFracIce(:)      ! volumetric fraction of ice in each layer (-)
 real(dp),intent(in)            :: mLayerDepth(:)           ! depth of each layer (m)
 ! model diagnostic variables (intent in)
 real(dp),intent(in)            :: mLayerHeight(:)          ! height at the mid-point of each layer (m)
 real(dp),intent(in)            :: iLayerHeight(0:)         ! height at the interface of each layer (m)
 real(dp),intent(in)            :: mLayerVolHtCapBulk(:)    ! bulk volumetric heat capacity (J m-3 K-1)
 real(dp),intent(in)            :: mLayerThermalC(:)        ! thermal conductivity at the mid-point of each layer (W m-1 K-1)
 real(dp),intent(in)            :: iLayerThermalC(0:)       ! thermal conductivity at the interface of each layer (W m-1 K-1)
 ! model diagnostic variables (intent out)
 real(dp),intent(out)           :: mLayerTcrit(:)           ! critical soil temperature above which all water is unfrozen (K)
 real(dp),intent(out)           :: mLayerdTheta_dTk(:)      ! derivative in the freezing curve (K-1)
 real(dp),intent(out)           :: mLayerTranspireLim(:)    ! moisture avail factor limiting transpiration in each layer (-)
 real(dp),intent(out)           :: mLayerInitTranspire(:)   ! transpiration loss from each soil layer at the start of the step (m s-1)
 real(dp),intent(out)           :: mLayerTranspire(:)       ! transpiration loss from each soil layer (m s-1)
 real(dp),intent(out)           :: iLayerInitNrgFlux(0:)    ! energy flux at layer interfaces at the start of the time step (W m-2)
 real(dp),intent(out)           :: iLayerNrgFlux(0:)        ! energy flux at layer interfaces at the end of the time step (W m-2)
 ! diagnostic scalar variables
 real(dp),intent(out)           :: scalarTranspireLim       ! aggregate soil moist & veg limit on transpiration, weighted by root density (-)
 real(dp),intent(out)           :: scalarPotentialET        ! potential ET (kg m-2 s-1)
 real(dp),intent(out)           :: scalarMassLiquid         ! evaporation/dew (kg m-2 s-1)
 real(dp),intent(out)           :: scalarMassSolid          ! sublimation/frost (kg m-2 s-1)
 real(dp),intent(out)           :: scalarSenHeat            ! sensible heat flux at the surface (W m-2)
 real(dp),intent(out)           :: scalarLatHeat            ! latent heat flux at the surface (W m-2)
 real(dp),intent(out)           :: scalarExCoef             ! turbulent exchange coefficient (-)
 real(dp),intent(out)           :: scalarExSen              ! exchange factor for sensible heat (J m-2 s-1 K-1)
 real(dp),intent(out)           :: scalarExLat              ! exchange factor for latent heat (J m-2 s-1)
 ! output variables from the tempchange subroutine
 real(dp),intent(out)           :: mLayerTempDiff(:)        ! iteration increment for temperature (K) 
 real(dp),intent(out)           :: mLayerTempNew(:)         ! new temperature (K)
 real(dp),intent(out)           :: mLayerVolFracIceNew(:)   ! new volumetric fraction of ice (-)
 real(dp),intent(out)           :: mLayerVolFracLiqNew(:)   ! new volumetric fraction of liquid water (-)
 real(dp),intent(out)           :: mLayerMatricHeadNew(:)   ! new matric head (m)
 integer(i4b),intent(out)       :: err                      ! error code
 character(*),intent(out)       :: message                  ! error message
 ! ---------------------------------------------------------------------------------------------------------------------------------
 ! define local variables
 ! ---------------------------------------------------------------------------------------------------------------------------------
 ! define general local variables
 character(LEN=256)             :: cmessage                 ! error message of downwind routine
 integer(i4b)                   :: iLayer                   ! index of model layers
 logical(lgt)                   :: printflag                ! .true. if print progress to the screen
 logical(lgt)                   :: fTranspire               ! .true. if computing transpiration
 real(dp)                       :: theta                    ! total volumetric water content (liquid plus ice)
 real(dp)                       :: critDiff                 ! temperature difference from critical temperature (K)
 real(dp),parameter             :: epsT=1.d-10              ! offset from Tcrit when re-setting iterations at the critical temperature (K)
 logical(lgt)                   :: computeDerivative        ! .true. if desire to compute the derivative
 ! define the local variables for the solution
 real(dp)                       :: totalSurfaceFlux         ! total surface flux (W m-2)
 real(dp)                       :: dTotalSurfaceFlux_dTemp  ! derivative in total surface flux w.r.t. temperature (W m-2 K-1)
 real(dp),dimension(0:size(mLayerTempIter)) :: dFlux_dTempAbove ! derivative in flux w.r.t. temperature in the layer above (J m-2 s-1 K-1)
 real(dp),dimension(0:size(mLayerTempIter)) :: dFlux_dTempBelow ! derivative in flux w.r.t. temperature in the layer below (J m-2 s-1 K-1)
 real(dp)                       :: nrg0,nrg1                ! energy content at the start of the time step / current iteration (J m-3)
 real(dp)                       :: flx0,flx1                ! fluxes at the start of the time step / current iteration (J m-3)
 real(dp)                       :: phse                     ! phase change term (J m-3)
 real(dp),dimension(size(mLayerTempIter))   :: rvec         ! residual vector (J m-3)
 real(dp)                                   :: wtim         ! weighted time (s-1)
 real(dp),dimension(size(mLayerTempIter)-1) :: d_m1         ! sub-diagonal elements of the tridiagonal system (J m-3 K-1)
 real(dp),dimension(size(mLayerTempIter))   :: diag         ! diagonal elements of the tridiagonal system (J m-3 K-1)
 real(dp),dimension(size(mLayerTempIter)-1) :: d_p1         ! super-diagonal elements of the tridiagonal system (J m-3 K-1)
 ! ---------------------------------------------------------------------------------------------------------------------------------
 ! initialize error control
 err=0; message="tempchange_muster/"

 ! initialize print flag
 printflag=.false.

 ! identify if there is a need to transpire
 fTranspire=.true.
 select case(trim(model_decisions(iLookDECISIONS%bound_cond)%decision))
  case('headflux'); fTranspire=.false.
  case('headhead'); fTranspire=.false.
 end select
 if(nSnow>0) fTranspire=.false.

 ! ***** compute fluxes at the surface
 call surfaceFlx(&
                 ! (model forcing variables)
                 airtemp,                 & ! intent(in): air temperature (K)
                 spechum,                 & ! intent(in): specific humidity (g g-1)
                 windspd,                 & ! intent(in): wind speed (m s-1)
                 airpres,                 & ! intent(in): air pressure (Pa)
                 sw_down,                 & ! intent(in): downwelling shortwave radiation (W m-2)   
                 lw_down,                 & ! intent(in): downwelling long wave radiation (W m-2)   
                 ! (model state variables)
                 scalarAlbedo,            & ! intent(in): surface albedo (-)
                 mLayerTempIter(1),       & ! intent(in): trial surface temperature at the current iteration (K)
                 mLayerMatricHeadIter,    & ! intent(in): trial matric head in all layers at the current iteration (m)
                 ! (control variables)
                 computeDerivative,       & ! intent(in): flag to compute the derivative
                 ! (diagnostic variables)
                 scalarExCoef,            & ! intent(out): surface-atmosphere exchange coeffcient (-)
                 scalarTranspireLim,      & ! intent(out): resistance to evaporation at the surface (-)
                 mLayerTranspireLim,      & ! intent(out): resistance to evaporation in each layer (-)
                 scalarExSen,             & ! intent(out): exchange factor for sensible heat (W m-2 K-1)
                 scalarExLat,             & ! intent(out): exchange factor for latent heat (W m-2)
                 scalarSenHeat,           & ! intent(out): sensible heat flux at the surface (W m-2)
                 scalarLatHeat,           & ! intent(out): latent heat flux at the surface (W m-2)
                 totalSurfaceFlux,        & ! intent(out): total surface flux (W m-2)
                 dTotalSurfaceFlux_dTemp, & ! intent(out): derivative in total surface flux w.r.t. temperature (W m-2 K-1)
                 err,cmessage)              ! intent(out): error control
 if(err/=0)then; message=trim(message)//trim(cmessage); return; endif

 ! ***** compute fluxes at layer interfaces and their derivatives (J m-2 s-1)
 call iLayer_nrg(&
                 ! (input)
                 mLayerDepth,            & ! intent(in): depth of each layer (m)
                 mLayerHeight,           & ! intent(in): height of layer mid-points (m)
                 iLayerThermalC,         & ! intent(in): thermal conductivity at layer interfaces (W m-1)
                 mLayerTempIter,         & ! intent(in): trial temperature at the current iteration (K)
                 totalSurfaceFlux,       & ! intent(in): total flux at the surface (W m-2)
                 dTotalSurfaceFlux_dTemp,& ! intent(in): derivative in total surface flux w.r.t. temperature (W m-2 K-1)
                 lowerBoundTemp,         & ! intent(in): temperature of the lower boundary (K)
                 computeDerivative,      & ! intent(in): flag to compute the derivative
                 ! (output)
                 iLayerNrgFlux,          & ! intent(out): energy flux at the layer interfaces (W m-2)
                 dFlux_dTempAbove,       & ! intent(out): derivatives in the flux w.r.t. temperature in the layer above (W m-2 K-1)
                 dFlux_dTempBelow,       & ! intent(out): derivatives in the flux w.r.t. temperature in the layer below (W m-2 K-1)
                 err,cmessage)             ! intent(out): error control
 if(err/=0)then; message=trim(message)//trim(cmessage); return; endif

 ! ***** compute the residual vector (J m-3)
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
 end do

 ! ***** compute the derivative in the freezing curve w.r.t. temperature (K-1)
 do iLayer=1,nLayers
  select case(layerType(iLayer))
   case(ix_snow) ! (snow layers)
    theta = mLayerVolFracIceIter(iLayer)*(iden_ice/iden_water) + mLayerVolFracLiqIter(iLayer)
    mLayerdTheta_dTk(iLayer) = dFracLiq_dTk(mLayerTempIter(iLayer),snowfrz_scale)*theta
   case(ix_soil) ! (soil layers)
    if(mLayerVolFracIceIter(iLayer)>0._dp)then
     mLayerdTheta_dTk(iLayer) = dTheta_dTk(mLayerTempIter(iLayer),theta_res,theta_sat,vGn_alpha,vGn_n,vGn_m)
    else
     mLayerdTheta_dTk(iLayer) = 0._dp
    endif
   case default
    err=40; message=trim(message)//"cannot identify the layer as snow or soil"; return
  endselect
 end do

 ! ***** assemble the tri-diagonal matrix
 wtim = (1._dp - wimplicit)*dt  ! weighted time
 diag = (wtim/mLayerDepth)*(-dFlux_dTempBelow(0:nLayers-1) + dFlux_dTempAbove(1:nLayers)) + mLayerdTheta_dTk*LH_fus*iden_water + mLayerVolHtCapBulk
 d_m1 = (wtim/mLayerDepth(1:nLayers-1))*(-dFlux_dTempAbove(1:nLayers-1) )
 d_p1 = (wtim/mLayerDepth(1:nLayers-1))*( dFlux_dTempBelow(1:nLayers-1) )

 ! ***** solve the tridiagonal system of equations -- returns mLayerTempDiff
 call tridag(d_m1,                    & ! intent(in): sub-diagonal elements of the tridiagonal system (J m-3 K-1)
             diag,                    & ! intent(in): diagonal elements of the tridiagonal system (J m-3 K-1)
             d_p1,                    & ! intent(in): super-diagonal elements of the tridiagonal system (J m-3 K-1)
             -rvec,                   & ! intent(in): residual vector (J m-3)
             mLayerTempDiff,          & ! intent(out): temperature increment (K)
             err,cmessage)              ! intent(out): error control
 if(err/=0)then; message=trim(message)//trim(cmessage); return; endif

 ! adjust del temperature in cases where snow temperature exceeds Tfreeze -- use bi-section
 if(nSnow>0)then
  do iLayer=1,nSnow
   if(mLayerTempIter(iLayer)+mLayerTempDiff(iLayer) > Tfreeze)then
    mLayerTempDiff(iLayer) = (Tfreeze-mLayerTempIter(iLayer))*0.5_dp
   endif
  end do
 endif

 ! adjust del temperature in cases where soil temperature crosses the critical temperature
 do iLayer=nSnow+1,nLayers
  ! get the difference from the critical temperature (K)
  critDiff = mLayerTcrit(iLayer-nSnow) - mLayerTempIter(iLayer)
  ! set temperature to Tcrit in cases where temperatures cross Tcrit
  if(critDiff > 0._dp)then  ! (mLayerTempIter < Tcrit)
   if(mLayerTempDiff(iLayer) > critDiff) mLayerTempDiff(iLayer) = critDiff + epsT 
  else                      ! (mLayerTempIter > Tcrit)
   if(mLayerTempDiff(iLayer) < critDiff) mLayerTempDiff(iLayer) = critDiff - epsT
  endif
 end do

 ! ***** update temperature
 mLayerTempNew = mLayerTempIter + mLayerTempDiff
 
 ! ***** compute phase change
 call phsechange(mLayerTempNew,       & ! intent(in): new temperature vector (K)
                 mLayerMatricHeadIter,& ! intent(in): matric head at the current iteration (m)
                 mLayerVolFracLiqIter,& ! intent(in): volumetric fraction of liquid water at the current iteration (-)
                 mLayerVolFracIceIter,& ! intent(in): volumetric fraction of ice at the current iteration (-)
                 mLayerMatricHeadNew, & ! intent(out): new matric head (m)
                 mLayerVolFracLiqNew, & ! intent(out): new volumetric fraction of liquid water (-)
                 mLayerVolFracIceNew, & ! intent(out): new volumetric fraction of ice (-)
                 err,cmessage)          ! intent(out): error control
 if(err/=0)then; message=trim(message)//trim(cmessage); return; endif

 ! ***** compute un-stressed ET (kg m-2 s-1)
 scalarPotentialET = iden_air * windspd * scalarExCoef * (spechum - relhm2sphm(RHsurf,airpres,mLayerTempNew(1)))

 ! ***** compute transpiration (m s-1)
 if(fTranspire)then
  ! ** compute actual transpiration from each soil layer (m s-1)
  mLayerTranspire(1:nSoil) =  (mLayerTranspireLim(1:nSoil)*mLayerRootDensity(1:nSoil)*scalarPotentialET)/iden_water
  ! ** compute total transpiration (kg m-2 s-1)
  scalarMassLiquid = ( wimplicit*(sum(mLayerInitTranspire)) + (1._dp - wimplicit)*sum(mLayerTranspire) )*iden_water
 else
  mLayerTranspire(1:nSoil) = 0._dp
  scalarMassLiquid = 0._dp
 endif

 ! ***** compute sublimation/frost (kg m-2 s-1)
 if(nSnow==0) scalarMassSolid  = 0._dp
 if(nsnow >0) scalarMassSolid  = scalarPotentialET

 ! update sensible and latent heat (W m-2) --> positive downwards
 scalarSenHeat = scalarExSen*(airtemp - mLayerTempNew(1))
 if(nSnow==0) scalarLatHeat = scalarMassLiquid*LH_vap
 if(nsnow >0) scalarLatHeat = scalarMassSolid*LH_sub

 ! ====================================================================================================================

 end subroutine tempchange_muster


 ! ************************************************************************************************
 ! private subroutine: compute surface energy flux and its derivative w.r.t. temperature
 ! ************************************************************************************************
 subroutine surfaceFlx(&
                       ! (model forcing variables)
                       airtemp,                 & ! intent(in): air temperature (K)
                       spechum,                 & ! intent(in): specific humidity (g g-1)
                       windspd,                 & ! intent(in): wind speed (m s-1)
                       airpres,                 & ! intent(in): air pressure (Pa)
                       sw_down,                 & ! intent(in): downwelling shortwave radiation (W m-2)   
                       lw_down,                 & ! intent(in): downwelling long wave radiation (W m-2)   
                       ! (model state variables)
                       scalarAlbedo,            & ! intent(in): surface albedo (-)
                       surfaceTempTrial,        & ! intent(in): trial surface temperature (K)
                       mLayerMatricHeadTrial,   & ! intent(in): trial matric head at the current iteration (m)
                       ! (control variables)
                       computeDerivative,       & ! intent(in): flag to compute the derivative
                       ! (diagnostic variables)
                       scalarExCoef,            & ! intent(out): surface-atmosphere exchange coeffcient (-)
                       scalarTranspireLim,      & ! intent(out): resistance to evaporation at the surface (-)
                       mLayerTranspireLim,      & ! intent(out): resistance to evaporation in each layer (-)
                       scalarExSen,             & ! intent(out): exchange factor for sensible heat (W m-2 K-1)
                       scalarExLat,             & ! intent(out): exchange factor for latent heat (W m-2)
                       scalarSenHeat,           & ! intent(out): sensible heat flux at the surface (W m-2)
                       scalarLatHeat,           & ! intent(out): latent heat flux at the surface (W m-2)
                       totalSurfaceFlux,        & ! intent(out): total surface flux (W m-2)
                       dTotalSurfaceFlux_dTemp, & ! intent(out): derivative in total surface flux w.r.t. temperature (W m-2 K-1)
                       err,message)               ! intent(out): error control
 ! -------------------------------------------------------------------------------------------------------
 USE conv_funcs_module,only:relhm2sphm            ! compute specific humidity
 implicit none
 ! input
 real(dp),intent(in)           :: airtemp                  ! air temperature (K)
 real(dp),intent(in)           :: spechum                  ! specific humidity (g g-1)
 real(dp),intent(in)           :: windspd                  ! wind speed (m s-1)
 real(dp),intent(in)           :: airpres                  ! air pressure (Pa)
 real(dp),intent(in)           :: sw_down                  ! downwelling shortwave radiation (W m-2)   
 real(dp),intent(in)           :: lw_down                  ! downwelling long wave radiation (W m-2)   
 real(dp),intent(in)           :: scalarAlbedo             ! surface albedo (-)
 real(dp),intent(in)           :: surfaceTempTrial         ! trial surface temperature (K)
 real(dp),intent(in)           :: mLayerMatricHeadTrial(:) ! trial matric head at the current iteration (m)
 logical(lgt),intent(in)       :: computeDerivative        ! flag to compute the derivative
 ! output
 real(dp),intent(out)          :: scalarExCoef             ! surface-atmosphere exchange coeffcient (-)
 real(dp),intent(out)          :: scalarTranspireLim       ! resistance to evaporation at the surface (-)
 real(dp),intent(out)          :: mLayerTranspireLim(:)    ! resistance to evaporation in each soil layer (-)
 real(dp),intent(out)          :: scalarExSen              ! exchange factor for sensible heat (W m-2 K-1)
 real(dp),intent(out)          :: scalarExLat              ! exchange factor for latent heat (W m-2)
 real(dp),intent(out)          :: scalarSenHeat            ! sensible heat flux at the surface (W m-2)
 real(dp),intent(out)          :: scalarLatHeat            ! latent heat flux at the surface (W m-2)
 real(dp),intent(out)          :: totalSurfaceFlux         ! total surface flux (W m-2)
 real(dp),intent(out)          :: dTotalSurfaceFlux_dTemp  ! derivative in total surface flux w.r.t. temperature (W m-2 K-1)
 integer(i4b),intent(out)      :: err                      ! error code
 character(*),intent(out)      :: message                  ! error message
 ! local variables
 character(len=256)            :: cmessage                 ! error message of downwind routine
 real(dp)                      :: dScalarExCoef_dTemp      ! derivative in surface-atmosphere exchange coeffcient w.r.t. temperature (K-1)
 real(dp)                      :: aerodynResist            ! aerodynamic resistance (s m-1)
 real(dp)                      :: Qh_temp                  ! "uncorrected" sensible heat flux (W m-2)
 real(dp)                      :: Qe_temp                  ! "uncorrected" latent heat flux (W m-2)
 real(dp)                      :: Qderiv                   ! derivative in specific humidity w.r.t. temperature (kg kg-1 K-1)
 real(dp)                      :: Qh_deriv                 ! derivative in sensible heat w.r.t. temperature (J m-2 s-1 K-1)
 real(dp)                      :: Qe_deriv                 ! derivative in latent heat w.r.t. temperature (J m-2 s-1 K-1)
 real(dp)                      :: LW_deriv                 ! derivative in longwave radiation w.r.t. temperature (J m-2 s-1 K-1)
 ! initialize error control
 err=0; message='surfaceFlx/'

 ! compute the surface exchange coefficients
 call exchCoefft(airtemp,                 & ! intent(in): air temperature (K)
                 windspd,                 & ! intent(in): wind speed (m s-1)
                 surfaceTempTrial,        & ! intent(in): trial surface temperature (K)
                 computeDerivative,       & ! intent(in): flag to compute the derivative
                 scalarExCoef,            & ! intent(out): surface-atmosphere exchange coeffcient (-)
                 dScalarExCoef_dTemp,     & ! intent(out): derivative in surface-atmosphere exchange coeffcient w.r.t. temperature (K-1)
                 err,cmessage)              ! intent(out): error control
 if(err/=0)then; message=trim(message)//trim(cmessage); return; endif

 ! compute aerodynamic resistance (s m-1)
 aerodynResist = 1._dp / (scalarExCoef*windspd)

 ! compute the ratio of actual:potential evapotranspiration
 call evapResist(aerodynResist,           & ! intent(in): aerodynamic resistance (s m-1)
                 mLayerMatricHeadTrial,   & ! intent(in): trial matric head in each layer (m)
                 scalarTranspireLim,      & ! intent(out): resistance to evaporation at the surface (-)
                 mLayerTranspireLim,      & ! intent(out): resistance to evaporation in each layer (-)
                 err,cmessage)              ! intent(out): error control
 if(err/=0)then; message=trim(message)//trim(cmessage); return; endif

 ! ***** compute the surface energy flux
 ! compute exchange factors for sensible and latent heat
 scalarExSen = Cp_air * iden_air * windspd * scalarExCoef                      ! J m-2 s-1 K-1
 if(nSnow >0) scalarExLat = LH_sub * iden_air * windspd                        ! J m-2 s-1
 if(nSnow==0) scalarExLat = LH_vap * iden_air * windspd  * scalarTranspireLim  ! J m-2 s-1
 ! compute sensible and latent heat at the current iteration (W m-2) --> positive downwards
 Qh_temp = scalarExSen*(airtemp - surfaceTempTrial)                            ! "uncorrected" sensible heat flux (W m-2)
 Qe_temp = scalarExSen*(spechum - relhm2sphm(RHsurf,airpres,surfaceTempTrial)) ! "uncorrected" latent heat flux (W m-2)
 scalarSenHeat = scalarExCoef * Qh_temp
 scalarLatHeat = scalarExCoef * Qe_temp
 ! compute the surface energy flux -- positive downwards
 totalSurfaceFlux = sw_down*(1._dp - scalarAlbedo) + lw_down - Em_Sno*sigma*surfaceTempTrial**4._dp + scalarSenHeat + scalarLatHeat

 ! ***** compute derivative in the surface energy flux
 if(computeDerivative)then
  ! compute the derivative in specific humidity w.r.t. temperature (kg kg-1 K-1)
  Qderiv = dsphum_dTk(RHsurf,airpres,surfaceTempTrial)
  ! compute the derivative in sensible and latent heat (J m-2 s-1 K-1) -- product rule
  Qh_deriv = dScalarExCoef_dTemp*Qh_temp - scalarExCoef*scalarExSen
  Qe_deriv = dScalarExCoef_dTemp*Qe_temp - scalarExCoef*scalarExLat*Qderiv
  ! compute the derivative in longwave radiation (J m-2 s-1 K-1)
  LW_deriv = -4._dp*Em_Sno*sigma*mLayerTempTrial**3._dp
  ! compute the total derivative
  dTotalSurfaceFlux_dTemp = Qh_deriv + Qe_deriv + LW_deriv
 else
  dTotalSurfaceFlux_dTemp = 1._dp
 endif

 end subroutine surfaceFlx



 ! ************************************************************************************************
 ! private subroutine: compute energy fluxes at layer interfaces, and their derivatives
 ! ************************************************************************************************
 subroutine iLayer_nrg(&
                       ! (input)
                       mLayerDepth,            & ! intent(in): depth of each layer (m)
                       mLayerHeight,           & ! intent(in): height of layer mid-points (m)
                       iLayerThermalC,         & ! intent(in): thermal conductivity at layer interfaces (W m-1)
                       mLayerTempTrial,        & ! intent(in): trial temperature at the current iteration (K)
                       totalSurfaceFlux,       & ! intent(in): total flux at the surface (W m-2)
                       dTotalSurfaceFlux_dTemp,& ! intent(in): derivative in total surface flux w.r.t. temperature (W m-2 K-1)
                       lowerBoundTemp,         & ! intent(in): temperature of the lower boundary (K)
                       computeDerivative,      & ! intent(in): flag to compute the derivative
                       ! (output)
                       iLayerNrgFlux,          & ! intent(out): energy flux at the layer interfaces (W m-2)
                       dFlux_dTempAbove,       & ! intent(out): derivatives in the flux w.r.t. temperature in the layer above (W m-2 K-1)
                       dFlux_dTempBelow,       & ! intent(out): derivatives in the flux w.r.t. temperature in the layer below (W m-2 K-1)
                       err,message)              ! intent(out): error control
 ! compute derivative in fluxes at layer interfaces w.r.t. temperature in the layer above and the layer below
 implicit none
 ! input
 real(dp),intent(in)           :: mLayerDepth(:)           ! intent(in): depth of each layer (m)
 real(dp),intent(in)           :: mLayerHeight(:)          ! intent(in): height of layer mid-points (m)
 real(dp),intent(in)           :: iLayerThermalC(0:)       ! thermal conductivity at layer interfaces (W m-1)
 real(dp),intent(in)           :: mLayerTempTrial(:)       ! trial temperature at the current iteration (K)
 real(dp),intent(in)           :: totalSurfaceFlux         ! total surface flux (W m-2)
 real(dp),intent(in)           :: dTotalSurfaceFlux_dTemp  ! derivative in total surface flux w.r.t. temperature (W m-2 K-1)
 real(dp),intent(in)           :: lowerBoundTemp           ! temperature of the lower boundary (K)
 logical(lgt),intent(in)       :: computeDerivative        ! flag to compute the derivative
 ! output
 real(dp),intent(out)          :: iLayerNrgFlux(0:)        ! energy flux at the layer interfaces (W m-2)
 real(dp),intent(out)          :: dFlux_dTempAbove(0:)     ! derivatives in the flux w.r.t. temperature in the layer above (J m-2 s-1 K-1)
 real(dp),intent(out)          :: dFlux_dTempBelow(0:)     ! derivatives in the flux w.r.t. temperature in the layer below (J m-2 s-1 K-1)
 integer(i4b),intent(out)      :: err                      ! error code
 character(*),intent(out)      :: message                  ! error message
 ! local variables
 integer(i4b)                  :: iLayer                   ! index of model layers
 ! initialize error control
 err=0; message='iLayer_nrg/'

 ! ***** compute fluxes at layer interfaces
 ! compute flux at the upper boundary -- positive downwards
 iLayerNrgFlux(0) = totalSurfaceFlux
 ! compute fluxes within the domain -- positive downwards
 do iLayer=1,nLayers-1
  iLayerNrgFlux(iLayer) = -iLayerThermalC(iLayer)*(mLayerTempIter(iLayer+1) - mLayerTempIter(iLayer)) / &
                                                  (mLayerHeight(iLayer+1) - mLayerHeight(iLayer))
 end do
 ! compute fluxes at the lower boundary
 iLayerNrgFlux(nLayers) = -iLayerThermalC(iLayer)*(lowerBoundTemp - mLayerTempIter(iLayer))/(mLayerDepth(iLayer)*0.5_dp)

 ! ***** compute the derivative in fluxes at layer interfaces w.r.t temperature in the layer above and the layer below
 if(computeDerivative)then

  do iLayer=0,nLayers  ! loop through interfaces
   ! ***** the upper boundary
   if(iLayer==0)then  ! (upper boundary)
    dFlux_dTempAbove(iLayer) = -huge(scalarExCoef)  ! don't expect this to be used, so deliberately set to a ridiculous value to cause problems
    dFlux_dTempBelow(iLayer) = dTotalSurfaceFlux_dTemp
   ! ***** the lower boundary
   elseif(iLayer==nLayers)then  ! (lower boundary)
    dFlux_dTempBelow(iLayer) = -huge(scalarExCoef)  ! don't expect this to be used, so deliberately set to a ridiculous value to cause problems
    dFlux_dTempAbove(iLayer) = iLayerThermalC(iLayer)/(mLayerDepth(iLayer)*0.5_dp)
   ! ***** internal layers
   else
    dFlux_dTempAbove(iLayer) =  iLayerThermalC(iLayer)/(mLayerHeight(iLayer+1) - mLayerHeight(iLayer))
    dFlux_dTempBelow(iLayer) = -iLayerThermalC(iLayer)/(mLayerHeight(iLayer+1) - mLayerHeight(iLayer))
   endif  ! type of layer (upper, internal, or lower)
  end do  ! (looping through layers)

 else  ! (not computing derivatives)
  dFlux_dTempAbove(:) = 1._dp
  dFlux_dTempBelow(:) = 1._dp
 endif 

 end subroutine iLayer_nrg



 ! ************************************************************************************************
 ! private subroutine: compute surface-atmosphere exchange coeffcient and its derivative w.r.t. temperature
 ! ************************************************************************************************
 subroutine exchCoefft(airtemp,               & ! intent(in): air temperature (K)
                       windspd,               & ! intent(in): wind speed (m s-1)
                       scalarTempTrial,       & ! intent(in): trial surface temperature (K)
                       computeDerivative,     & ! intent(in): flag to compute the derivative
                       exchangeCoefft,        & ! intent(out): surface-atmosphere exchange coeffcient (-)
                       dExchangeCoefft_dTemp, & ! intent(out): derivative in surface-atmosphere exchange coeffcient w.r.t. temperature (K-1)
                       err,message)             ! intent(out): error control
 USE snow_utils_module,only: bulkRichardson     ! compute the bulk Richardson number and its derivative w.r.t. temperature
 USE snow_utils_module,only: astability         ! compute surface exchange coefficient and its derivative w.r.t. temperature
 implicit none
 ! input
 real(dp),intent(in)           :: airtemp                  ! air temperature (K)
 real(dp),intent(in)           :: windspd                  ! wind speed (m s-1)
 real(dp),intent(in)           :: scalarTempTrial          ! trial surface temperature (K)
 logical(lgt),intent(in)       :: computeDerivative        ! flag to compute the derivative
 ! output
 real(dp),intent(out)          :: exchangeCoefft           ! surface-atmosphere exchange coeffcient (-) 
 real(dp),intent(out)          :: dExchangeCoefft_dTemp    ! derivative in surface-atmosphere exchange coeffcient w.r.t. temperature (K-1)
 integer(i4b),intent(out)      :: err                      ! error code
 character(*),intent(out)      :: message                  ! error message
 ! local variables
 real(dp)                      :: RiBulk                   ! bulk Richardson number (-)
 real(dp)                      :: dRiBulk_dTemp            ! derivative in the bulk Richardson number w.r.t. temperature (K-1) 

 ! initialize error control
 err=0; message='exchCoefft/'

 ! compute the bulk Richardson number and its derivative
 call bulkRichardson(airtemp,surfTemp,windspd,mheight,computeDerivative, & ! (input)
                     RiBulk,dRiBulk_dTemp,err,cmessage)                    ! (output)
 if(err/=0)then; err=10; message=trim(message)//trim(cmessage); return; endif

 ! compute surface-atmosphere exchange coeffcient and its derivative w.r.t. temperature
 call astability(RiBulk,dRiBulk_dTemp,computeDerivative, &           ! (input)
                 exchangeCoefft,dExchangeCoefft_dTemp,err,cmessage)  ! (output)
 if(err/=0)then; err=10; message=trim(message)//trim(cmessage); return; endif

 end subroutine exchCoefft


 ! ************************************************************************************************
 ! private subroutine: compute the resistance to evaporation (-)
 ! ************************************************************************************************
 subroutine evapResist(aerodynResist,        &  ! intent(in): aerodynamic resistance (s m-1)
                       mLayerMatricHeadIter, &  ! intent(in): matric head in each layer (m)
                       evapResistance,       &  ! intent(out): resistance to evaporation at the surface (-)
                       mLayerEvapResistance, &  ! intent(out): resistance to evaporation in each layer (-)
                       err,message)             ! intent(out): error control
 implicit none
 ! input
 real(dp),intent(in)           :: aerodynResist            ! aerodynamic resistance
 real(dp),intent(in)           :: mLayerMatricHeadIter(:)  ! matric head in each layer (m)
 ! output
 real(dp),intent(out)          :: evapResistance           ! resistance to evaporation at the surface (-) 
 real(dp),intent(out)          :: mLayerEvapResistance(:)  ! resistance to evaporation in each soil layer (-)
 integer(i4b),intent(out)      :: err                      ! error code
 character(*),intent(out)      :: message                  ! error message
 ! local variables
 real(dp)                      :: stomatalResist           ! stomatal resistance (s m-1)
 real(dp)                      :: plantWiltFactor          ! plant wilting factor (-) 
 integer(i4b)                  :: iLayer                   ! index of model layers
 ! initialize error control
 err=0; message='evapResist/'
 ! compute the ratio of actual:potential evapotranspiration
 if(nSnow>0)then
  evapResistance = 1._dp
 else
  ! ** compute the factor limiting evaporation for each soil layer (-)
  evapResistance = 0._dp  ! (initialize the weighted average)
  do iLayer=1,nSoil
   ! compute the stomatal resistance of the canopy (m s-1)
   plantWiltFactor = 1._dp + ( min(mLayerMatricHeadIter(iLayer),0._dp) / plantWiltPsi )**plantWiltExp
   stomatalResist  = (minStomatalResist/LAI)*plantWiltFactor
   ! compute the factor limiting evaporation for a given soil layer (-)
   mLayerEvapResistance(iLayer) = aerodynResist / (stomatalResist + aerodynResist)
   ! commpute the weighted average (weighted by root density)
   evapResistance = evapResistance + mLayerEvapResistance(iLayer)*mLayerRootDensity(iLayer)
  end do ! (looping through soil layers)
 endif  ! (if surface is snow-free)
 end function evapResist


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

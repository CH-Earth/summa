module ConvE2Temp_module
USE nrtype
implicit none
private
public::E2T_lookup
public::layerEthpy
public::layer_Temp
public::temp2ethpy
public::E2T_nosoil
! define the look-up table used to compute temperature based on enthalpy
integer(i4b),parameter            :: nlook=10001       ! number of elements in the lookup table
real(dp),dimension(nlook),public  :: E_lookup          ! enthalpy values (J kg-1)
real(dp),dimension(nlook),public  :: T_lookup          ! temperature values (K)
contains

 ! **********************************************************************************************************
 ! new subroutine: define a look-up table to compute specific enthalpy based on temperature, assuming no soil
 ! **********************************************************************************************************
 subroutine E2T_lookup(err,message)
 USE nr_utility_module,only:arth                       ! use to build vectors with regular increments
 USE spline_int_module,only:spline,splint              ! use for cubic spline interpolation
 USE multiconst,only:Tfreeze                           ! freezing point (K)
 USE data_struc,only:mpar_data                         ! model parameter structures (use snow_frz)
 USE var_lookup,only:iLookPARAM                        ! named variables to define structure element
 implicit none
 ! declare dummy variables
 integer(i4b),intent(out)      :: err                  ! error code
 character(*),intent(out)      :: message              ! error message
 ! define pointers to parameter structures
 real(dp),pointer              :: snow_frz             ! freezing curve parameter for snow (K-1)
 ! declare local variables
 character(len=128)            :: cmessage             ! error message in downwind routine
 real(dp),parameter            :: T_start=263.5_dp     ! start temperature value (K)
 real(dp),parameter            :: T_end  =Tfreeze      ! end temperature value (K)
 real(dp),parameter            :: E_start=-20000._dp   ! start enthalpy value (J kg-1)
 real(dp),parameter            :: E_end  =333500._dp   ! end enthalpy value (J kg-1)
 real(dp)                      :: T_incr,E_incr        ! temperature/enthalpy increments
 real(dp),dimension(nlook)     :: Tk                   ! initial temperature vector
 real(dp),dimension(nlook)     :: Ey                   ! initial enthalpy vector
 real(dp),parameter            :: soilWght=0._dp       ! weight applied to soil (kg m-3)
 real(dp),parameter            :: waterWght=1._dp      ! weight applied to total water (kg m-3) --- cancels out
 real(dp),dimension(nlook)     :: T2deriv              ! 2nd derivatives of the interpolating function at tabulated points
 integer(i4b)                  :: ilook                ! loop through lookup table
 ! initialize error control
 err=0; message="E2T_lookup/"
 ! assign pointers
 snow_frz => mpar_data%var(iLookPARAM%snow_frz)
 ! define initial temperature vector
 T_incr = (T_end - T_start) / real(nlook-1, kind(dp))  ! temperature increment
 Tk     = arth(T_start,T_incr,nlook)
 ! ***** compute specific enthalpy (NOTE: J m-3 --> J kg-1) *****
 do ilook=1,nlook
  Ey(ilook) = temp2ethpy(Tk(ilook),soilWght,waterWght,snow_frz)/waterWght  ! (J m-3 --> J kg-1)
 end do
 ! define the final enthalpy vector
 E_incr   = (E_end - E_start) / real(nlook-1, kind(dp))  ! enthalpy increment
 E_lookup = arth(E_start,E_incr,nlook)
 ! check that all values are covered
 if(E_start<Ey(1)) then; err=20; message="f-E2T_lookup/some desired enthalpy values outside computed range"; return; endif
 ! use cubic spline interpolation to obtain temperature values at the desired values of enthalpy
 call spline(Ey,Tk,1.e30_dp,1.e30_dp,T2deriv,err,cmessage)  ! get the second derivatives
 if(err/=0) then; message=trim(message)//trim(cmessage); return; endif
 do ilook=1,nlook
  call splint(Ey,Tk,T2deriv,E_lookup(ilook),T_lookup(ilook),err,cmessage)
  if(err/=0) then; message=trim(message)//trim(cmessage); return; endif
 end do
 end subroutine E2T_lookup

 ! **********************************************************************************************************
 ! new subroutine: loop through layers and compute total enthalpy based on temperature and mass
 ! **********************************************************************************************************
 subroutine layerEthpy(err,message)
 ! used to compute enthalpy based on temperature and total mass in layer (snow or soil)
 USE data_struc,only:ix_soil,ix_snow                    ! named variables to denote soil and snow
 USE data_struc,only:mpar_data,mvar_data,indx_data      ! model structure
 USE var_lookup,only:iLookPARAM,iLookMVAR,iLookINDEX    ! named variables that define structure elements
 implicit none
 ! declare dummy variables
 integer(i4b),intent(out)      :: err                   ! error code
 character(*),intent(out)      :: message               ! error message
 ! declare local parameters
 real(dp),pointer              :: soil_frz              ! freezing curve parameter for soil (K-1)
 real(dp),pointer              :: snow_frz              ! freezing curve parameter for snow (K-1)
 real(dp),pointer              :: soil_dens_bulk        ! bulk density of soil (kg m-3)
 ! declare local state variables
 real(dp),pointer,dimension(:) :: mlayerTemp(:)         ! temperature of each layer (K)
 real(dp),pointer,dimension(:) :: mlayerBulkDenWat(:)   ! bulk density of water in each layer (kg m-3)
 real(dp),pointer,dimension(:) :: mlayerEnthalpy(:)     ! enthalpy of each layer (J m-3)
 ! declare other local variables
 integer(i4b)                  :: ilayer                ! loop through model layers
 integer(i4b),pointer          :: nLayers               ! number of layers
 integer(i4b),pointer          :: layerType(:)          ! layer type (snow or soil)
 real(dp),pointer              :: fc_param              ! freezing curve parameter (K-1)
 real(dp)                      :: BulkDenSoil           ! bulk density of soil (kg m-3)
 ! initialize error control
 err=0; message="layerEthpy/"
 ! map local parameters onto the parameter data structure
 soil_frz         => mpar_data%var(iLookPARAM%soil_frz)             ! freezing curve parameter for soil (K-1)
 snow_frz         => mpar_data%var(iLookPARAM%snow_frz)             ! freezing curve parameter for snow (K-1)
 soil_dens_bulk   => mpar_data%var(iLookPARAM%soil_dens_bulk)       ! bulk density of soil (kg m-3)
 ! map local state variables onto the state data structure
 mLayerTemp       => mvar_data%var(iLookMVAR%mlayerTemp)%dat        ! temperature of each layer (K)
 mLayerBulkDenWat => mvar_data%var(iLookMVAR%mlayerBulkDenWat)%dat  ! bulk density of water in each layer (kg m-3)
 mLayerEnthalpy   => mvar_data%var(iLookMVAR%mLayerEnthalpy)%dat    ! enthalpy of each layer (J m-3)
 ! map state indices
 nLayers          => indx_data%var(iLookINDEX%nLayers)%dat(1)       ! number of layers
 layerType        => indx_data%var(iLookINDEX%layerType)%dat        ! layer type (snow or soil)
 ! loop through model layers
 do ilayer=1,nLayers
  ! identify the appropriate freezing curve parameter
  if(layerType(iLayer)==ix_soil) fc_param => soil_frz
  if(layerType(iLayer)==ix_snow) fc_param => snow_frz
  ! identify the soil density
  if(layerType(iLayer)==ix_soil) BulkDenSoil = soil_dens_bulk
  if(layerType(iLayer)==ix_snow) BulkDenSoil = 0._dp
  ! compute enthalpy
  mLayerEnthalpy(ilayer) = temp2ethpy(mlayerTemp(ilayer),BulkDenSoil,mlayerBulkDenWat(ilayer),fc_param)  ! (J m-3)
 end do  ! (looping through layers)
 end subroutine layerEthpy


 ! **********************************************************************************************************
 ! new subroutine: loop through layers and compute temperature based on enthalpy and mass
 ! **********************************************************************************************************
 subroutine layer_Temp(err,message)
 ! used to compute enthalpy based on temperature and total mass in layer (snow or soil)
 USE data_struc,only:ix_soil,ix_snow                    ! named variables to denote soil and snow
 USE data_struc,only:mpar_data,mvar_data,indx_data      ! model structure
 USE var_lookup,only:iLookPARAM,iLookMVAR,iLookINDEX    ! named variables that define structure elements
 implicit none
 ! declare dummy variables
 integer(i4b),intent(out)      :: err                   ! error code
 character(*),intent(out)      :: message               ! error message
 ! declare local parameters
 real(dp),pointer              :: soil_frz              ! freezing curve parameter for soil (K-1)
 real(dp),pointer              :: snow_frz              ! freezing curve parameter for snow (K-1)
 real(dp),pointer              :: soil_dens_bulk        ! bulk density of soil (kg m-3)
 ! declare local state variables
 real(dp),pointer,dimension(:) :: mlayerTemp(:)         ! temperature of each layer (K)
 real(dp),pointer,dimension(:) :: mlayerBulkDenWat(:)   ! bulk density of water in each layer (kg m-3)
 real(dp),pointer,dimension(:) :: mlayerEnthalpy(:)     ! enthalpy of each layer (J m-3)
 ! declare other local variables
 character(len=256)            :: cmessage              ! error message for downwind routine
 integer(i4b)                  :: ilayer                ! loop through model layers
 integer(i4b),pointer          :: nLayers               ! number of layers
 integer(i4b),pointer          :: layerType(:)          ! layer type (snow or soil)
 ! add temporary variables for testing
 !real(dp)                      :: test                  ! temporary variables for testing
 ! initialize error control
 err=0; message="layer_Temp/"
 ! map local parameters onto the parameter data structure
 soil_frz         => mpar_data%var(iLookPARAM%soil_frz)             ! freezing curve parameter for soil (K-1)
 snow_frz         => mpar_data%var(iLookPARAM%snow_frz)             ! freezing curve parameter for snow (K-1)
 soil_dens_bulk   => mpar_data%var(iLookPARAM%soil_dens_bulk)       ! bulk density of soil (kg m-3)
 ! map local state variables onto the state data structure
 mLayerTemp       => mvar_data%var(iLookMVAR%mlayerTemp)%dat        ! temperature of each layer (K)
 mLayerBulkDenWat => mvar_data%var(iLookMVAR%mlayerBulkDenWat)%dat  ! bulk density of water in each layer (kg m-3)
 mLayerEnthalpy   => mvar_data%var(iLookMVAR%mLayerEnthalpy)%dat    ! enthalpy of each layer (J m-3)
 ! map state indices
 nLayers          => indx_data%var(iLookINDEX%nLayers)%dat(1)       ! number of layers
 layerType        => indx_data%var(iLookINDEX%layerType)%dat        ! layer type (snow or soil)

 ! loop through model layers
 do ilayer=1,nLayers
  select case(layerType(iLayer))
   case(ix_soil); call Ethpy2Temp(mLayerEnthalpy(iLayer),soil_dens_bulk,mLayerBulkDenWat(iLayer),soil_frz,mLayerTemp(iLayer),err,cmessage)
   case(ix_snow); call E2T_nosoil(mLayerEnthalpy(iLayer),0._dp,mLayerBulkDenWat(iLayer),snow_frz,mLayerTemp(iLayer),err,cmessage)
   case default; err=20; message=trim(message)//"unknownLayerType"; return
  endselect
  if(err/=0)then; err=10; message=trim(message)//trim(cmessage); return; endif
 end do  ! (looping through layers)
 end subroutine layer_Temp


 ! **********************************************************************************************************
 ! new subroutine: compute temperature based on specific enthalpy -- appropriate when no dry mass, as in snow
 ! **********************************************************************************************************
 subroutine E2T_nosoil(Ey,BulkDenSoil,BulkDenWater,fc_param,Tk,err,message)
 ! compute temperature based on enthalpy -- appropriate when no dry mass, as in snow
 USE multiconst, only: Tfreeze, &                   ! freezing point of water (K)
                       Cp_soil,Cp_water,Cp_ice,&    ! specific heat of soil, water and ice (J kg-1 K-1)
                       LH_fus                       ! latent heat of fusion (J kg-1)
 implicit none
 ! declare dummy variables
 real(dp),intent(in)      :: Ey            ! total enthalpy (J m-3)
 real(dp),intent(in)      :: BulkDenSoil   ! bulk density of soil (kg m-3)
 real(dp),intent(in)      :: BulkDenWater  ! bulk density of water (kg m-3)
 real(dp),intent(in)      :: fc_param      ! freezing curve parameter (K-1)
 real(dp),intent(out)     :: Tk            ! initial temperature guess / final temperature value (K)
 integer(i4b),intent(out) :: err           ! error code
 character(*),intent(out) :: message       ! error message
 ! declare local variables
 real(dp),parameter       :: dx=-1.d-8     ! finite difference increment (J kg-1)
 real(dp),parameter       :: atol=1.d-12   ! convergence criteria (J kg-1)
 real(dp)                 :: E_spec        ! specific enthalpy (J kg-1)
 real(dp)                 :: E_incr        ! enthalpy increment
 integer(i4b)             :: niter=15      ! maximum number of iterations
 integer(i4b)             :: iter          ! iteration index
 integer(i4b)             :: i0            ! position in lookup table
 real(dp)                 :: Tg0,Tg1       ! trial temperatures (K)
 real(dp)                 :: Ht1           ! Total enthalpy, based on the trial temperatures (J m-3)
 real(dp)                 :: f0,f1         ! function evaluations (difference between enthalpy guesses)
 real(dp)                 :: dh            ! enthalpy derivative
 real(dp)                 :: dT            ! temperature increment
 ! initialize error control
 err=0; message="E2T_nosoil/"
 ! check that soil does not exist
 if(abs(BulkDenSoil)>epsilon(BulkDenSoil))then
  err=10; message=trim(message)//"soilExists"; return
 endif
 ! convert input of total enthalpy (J m-3) to total specific enthalpy (J kg-1)
 E_spec = Ey/BulkDenWater ! (NOTE: no soil)
 ! process cases below the limit (assume no fractional liquid water below limit)
 if(E_spec<E_lookup(1))then
  Tk = (E_spec - E_lookup(1))/Cp_ice + T_lookup(1)
  return
 endif
 ! get enthalpy increment
 E_incr = E_lookup(2) - E_lookup(1)
 ! get position in lookup table
 i0 = int( (E_spec - E_lookup(1)) / E_incr, kind(i4b) )
 ! get temperature guess
 Tg0 = T_lookup(i0)
 Tg1 = T_lookup(i0+1)
 ! compute function evaluations
 f0  = E_spec - E_lookup(i0)
 f1  = E_spec - E_lookup(i0+1)
 ! compute derivative
 dh  = (f0 - f1) / (Tg0 - Tg1)
 ! compute change in T
 dT  = f0/dh
 ! exit if already close enough
 if(dT<atol)then
  Tk = Tg0+dT
  return
 endif
 ! iterate a little
 do iter=1,niter
  ! save old function evaluation and temperature
  f0  = f1
  Tg0 = Tg1
  ! compute new value of tg
  Tg1 = Tg0+dT
  ! compute new value of specific enthalpy at Tg1 (J kg-1)
  Ht1 = temp2ethpy(Tg1,0._dp,1._dp,fc_param)
  ! compute new function evaluation
  f1  = E_spec - Ht1
  ! compute derivative
  dh  = (f0 - f1) / (Tg0 - Tg1)
  ! compute change in T
  dT  = f1/dh
  ! exit if converged
  if(dT<atol)then
   Tk = Tg1+dT
   return
  endif
  if(iter==niter)then; err=20; message=trim(message)//"failedToConverge"; return; endif
 enddo  ! (iteration loop
 end subroutine E2T_nosoil


 ! **********************************************************************************************************
 ! new subroutine: compute temperature based on enthalpy and mass
 ! **********************************************************************************************************
 subroutine Ethpy2Temp(Ey,BulkDenSoil,BulkDenWater,fc_param,Tk,err,message)
 ! compute temperature based on enthalpy and the bulk density of dry soil and total water
 USE multiconst, only: Tfreeze, &                   ! freezing point of water (K)
                       Cp_soil,Cp_water,Cp_ice,&    ! specific heat of soil, water and ice (J kg-1 K-1)
                       LH_fus                       ! latent heat of fusion (J kg-1)
 implicit none
 ! declare dummy variables
 real(dp),intent(in)      :: Ey            ! total enthalpy (J m-3)
 real(dp),intent(in)      :: BulkDenSoil   ! bulk density of soil (kg m-3)
 real(dp),intent(in)      :: BulkDenWater  ! bulk density of water (kg m-3)
 real(dp),intent(in)      :: fc_param      ! freezing curve parameter
 real(dp),intent(inout)   :: Tk            ! initial temperature guess / final temperature value (K)
 integer(i4b),intent(out) :: err           ! error code
 character(*),intent(out) :: message       ! error message
 ! declare local variables
 real(dp),parameter       :: dx=-1.d-8     ! finite difference increment (J kg-1)
 real(dp),parameter       :: atol=1.d-12   ! convergence criteria (J kg-1)
 real(dp)                 :: Tmin=-500._dp ! minimum temperature (K) -- avoid excessive extrapolation
 real(dp)                 :: Tmax= 900._dp ! maximum temperature (K) -- avoid excessive extrapolation
 integer(i4b)             :: niter=15      ! maximum number of iterations
 integer(i4b)             :: iter          ! iteration index
 real(dp)                 :: Tg0,Tg1       ! trial temperatures (K)
 real(dp)                 :: Ht0,Ht1       ! Total enthalpy, based on the trial temperatures (J m-3)
 real(dp)                 :: f0,f1         ! function evaluations (difference between enthalpy guesses)
 real(dp)                 :: dh            ! enthalpy derivative
 real(dp)                 :: dT            ! tenperature increment
 ! initialize error control
 err=0; message="Ethpy2Temp/"
 ! initialize temperature range
 Tmin=-500._dp; Tmax=900._dp
 ! begin iteration loop
 do iter=1,niter
  ! identify trial temperatures
  Tg0 = Tk
  Tg1 = Tg0 + dx
  ! compute total enthalpy, based on the trial temperature
  Ht0 = temp2ethpy(Tg0,BulkDenSoil,BulkDenWater,fc_param)
  Ht1 = temp2ethpy(Tg1,BulkDenSoil,BulkDenWater,fc_param)
  ! compute function evaluations
  f0  = Ey - ht0
  f1  = Ey - ht1
  ! compute derivative
  dh  = (f0 - f1) / dx
  ! compute change in T
  dT = f0/dh
  ! compute new value of T
  Tk = tg0 + f0/dh
  ! update bounds
  if(dT> 0._dp)tmin=tg0
  if(dT<=0._dp)tmax=tg0
  ! use bi-section if outside bounds
  if(Tk<tmin .or. Tk>tmax) Tk=0.5_dp*(tmax+tmin)
  ! exit if achieved tolerance
  if(abs(dT)<atol .or. abs(f0)<atol) exit
  ! check for failed convergence
  if(iter==niter)then; err=20; message=trim(message)//"failedToConverge"; return; endif
 enddo  ! (iteration loop
 end subroutine Ethpy2Temp


 ! **********************************************************************************************************
 ! new function: compute total enthalpy based on temperature and mass (J m-3)
 ! **********************************************************************************************************
 function temp2ethpy(Tk,BulkDenSoil,BulkDenWater,fc_param)
 ! used to compute enthalpy based on temperature and total mass in layer (snow or soil)
 USE multiconst, only: Tfreeze, &                   ! freezing point of water (K)
                       Cp_soil,Cp_water,Cp_ice,&    ! specific heat of soil, water and ice (J kg-1 K-1)
                       LH_fus                       ! latent heat of fusion (J kg-1)
 implicit none
 ! declare dummy variables
 real(dp),intent(in)  :: Tk            ! layer temperature (K)
 real(dp),intent(in)  :: BulkDenSoil   ! bulk density of soil (kg m-3)
 real(dp),intent(in)  :: BulkDenWater  ! bulk density of water (kg m-3)
 real(dp),intent(in)  :: fc_param      ! freezing curve parameter
 real(dp)             :: temp2ethpy    ! return value of the function, total specific enthalpy (J m-3)
 ! declare local variables
 real(dp)             :: frac_liq      ! fraction of liquid water
 real(dp)             :: enthTempSoil  ! temperature component of enthalpy for dry soil
 real(dp)             :: enthTempWater ! temperature component of enthalpy for total water (liquid and ice)
 real(dp)             :: enthMass      ! mass component of enthalpy
 ! compute the fraction of liquid water in the given layer
 frac_liq     = 1._dp / ( 1._dp + ( fc_param*( Tfreeze - min(Tk,Tfreeze) ) )**2._dp )
 ! compute the temperature component of enthalpy for the soil constituent (J m-3)
 enthTempSoil = Cp_soil*(Tk - Tfreeze) * BulkDenSoil
 ! compute the temperature component of enthalpy for total water (J m-3)
 if(Tk< Tfreeze) enthTempWater = BulkDenWater * Cp_ice*(Tk - Tfreeze) - (Cp_water - Cp_ice)*(atan(fc_param*(Tfreeze - Tk))/fc_param)
 if(Tk>=Tfreeze) enthTempWater = BulkDenWater * Cp_water*(Tk - Tfreeze)
 ! compute the mass component of enthalpy (J m-3) --> energy required to melt ice (note: enthalpy at Tfreeze is just mass component)
 enthMass     = LH_fus * frac_liq*BulkDenWater
 ! finally, compute the total enthalpy (J m-3)
 temp2ethpy   = enthTempSoil + enthTempWater + enthMass
 end function temp2ethpy

end module ConvE2Temp_module

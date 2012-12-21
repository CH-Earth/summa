module derivforce_module
USE nrtype
implicit none
private
public::derivforce
contains

 ! ************************************************************************************************
 ! new subroutine: compute derived forcing data
 ! ************************************************************************************************
 subroutine derivforce(err,message)
 USE nr_utility_module,only:erf                             ! provide access to the error function
 USE multiconst,only:Tfreeze                                ! freezing point of pure water (K)
 USE data_struc,only:mpar_data,forc_data,mvar_data          ! data structures
 USE var_lookup,only:iLookPARAM,iLookFORCE,iLookMVAR        ! named variables for structure elements
 USE conv_funcs_module,only:SPHM2RELHM,RELHM2SPHM,WETBULBTMP ! conversion functions
 USE snow_utils_module,only:fracliquid,templiquid           ! functions to compute temperature/liquid water
 ! compute derived forcing data variables
 implicit none
 ! dummy variables
 integer(i4b),intent(out)      :: err                       ! error code
 character(*),intent(out)      :: message                   ! error message
 ! local pointers to model parameters
 real(dp),pointer              :: Fabs_vis                  ! fraction radiation absorbed in visible part of spectrum (-)
 real(dp),pointer              :: minwind                   ! minimum windspeed (m s-1)
 real(dp),pointer              :: fc_param                  ! freeezing curve parameter for snow (K-1)
 real(dp),pointer              :: tempCritRain              ! critical temperature where precipitation is rain (K)
 real(dp),pointer              :: tempRangeTimestep         ! temperature range over the time step (K)
 ! local pointers to model forcing data
 real(dp),pointer              :: sw_down                   ! downward shortwave radiation (W m-2)
 real(dp),pointer              :: lw_down                   ! downward longwave radiation (W m-2)
 real(dp),pointer              :: airtemp                   ! air temperature at 2 meter height (K)
 real(dp),pointer              :: windspd                   ! wind speed at 10 meter height (m s-1)
 real(dp),pointer              :: airpres                   ! air pressure at 2 meter height (Pa)
 real(dp),pointer              :: spechum                   ! specific humidity at 2 meter height (g g-1)
 real(dp),pointer              :: pptrate                   ! precipitation rate (kg m-2 s-1)
 ! local pointers to model variables
 real(dp),pointer              :: swDownVis                 ! downwelling shortwave raadiation in the visible part of the spectrum (W m-2
 real(dp),pointer              :: swDownNir                 ! downwelling shortwave raadiation in the near-infrared part of the spectrum (W m-2)
 real(dp),pointer              :: twetbulb                  ! wet bulb temperature (K)
 real(dp),pointer              :: rainfall                  ! computed rainfall rate (kg m-2 s-1)
 real(dp),pointer              :: snowfall                  ! computed snowfall rate (kg m-2 s-1)
 real(dp),pointer              :: snowfallTemp              ! computed temperature of fresh snow (K)
 ! local variables
 real(dp)                      :: relhum                    ! relative humidity (-)
 real(dp)                      :: fracrain                  ! fraction of precipitation that falls as rain
 real(dp)                      :: maxFrozenSnowTemp         ! maximum temperature of snow when the snow is predominantely frozen (K)
 real(dp),parameter            :: unfrozenLiq=0.01_dp       ! unfrozen liquid water used to compute maxFrozenSnowTemp (-)
 real(dp),parameter            :: eps=epsilon(fracrain)     ! a number that is almost negligible
 real(dp)                      :: Tmin,Tmax                 ! minimum and maximum wet bulb temperature in the time step (K)
 ! initialize error control
 err=0; message="f-derivforce/"
 ! assign pointers to model parameters
 Fabs_vis          => mpar_data%var(iLookPARAM%Fabs_vis)            ! fraction radiation absorbed in visible part of spectrum (-)
 minwind           => mpar_data%var(iLookPARAM%minwind)             ! minimum windspeed (m s-1)
 fc_param          => mpar_data%var(iLookPARAM%snowfrz_scale)       ! freezing curve parameter for snow (K-1)
 tempCritRain      => mpar_data%var(iLookPARAM%tempCritRain)        ! critical temperature where precipitation is rain (K)
 tempRangeTimestep => mpar_data%var(iLookPARAM%tempRangeTimestep)   ! temperature range over the time step (K)
 ! assign pointers to model forcing data
 sw_down  => forc_data%var(iLookFORCE%sw_down)                      ! downward shortwave radiation (W m-2)
 lw_down  => forc_data%var(iLookFORCE%lw_down)                      ! downward longwave radiation (W m-2)
 airtemp  => forc_data%var(iLookFORCE%airtemp)                      ! air temperature at 2 meter height (K)
 windspd  => forc_data%var(iLookFORCE%windspd)                      ! wind speed at 10 meter height (m s-1)
 airpres  => forc_data%var(iLookFORCE%airpres)                      ! air pressure at 2 meter height (Pa)
 spechum  => forc_data%var(iLookFORCE%spechum)                      ! specific humidity at 2 meter height (g g-1)
 pptrate  => forc_data%var(iLookFORCE%pptrate)                      ! precipitation rate (kg m-2 s-1)
 ! assign pointers to model variables
 swDownVis    => mvar_data%var(iLookMVAR%scalarSwDownVis)%dat(1)    ! downwelling shortwave raadiation in the visible part of the spectrum (W m-2)
 swDownNir    => mvar_data%var(iLookMVAR%scalarSwDownNir)%dat(1)    ! downwelling shortwave raadiation in the near-infrared part of the spectrum (W m-2)
 twetbulb     => mvar_data%var(iLookMVAR%scalarTwetbulb)%dat(1)     ! wet bulb temperature (K)
 rainfall     => mvar_data%var(iLookMVAR%scalarRainfall)%dat(1)     ! computed rainfall rate (kg m-2 s-1)
 snowfall     => mvar_data%var(iLookMVAR%scalarSnowfall)%dat(1)     ! computed snowfall rate (kg m-2 s-1)
 snowfallTemp => mvar_data%var(iLookMVAR%scalarSnowfallTemp)%dat(1) ! computed temperature of fresh snow (K)
 ! compute shortwave radiation in the visible and near-infra-red part of the spectrum
 swDownVis = Fabs_vis*sw_down
 swDownNir = sw_down - swDownVis
 ! ensure wind speed is above a prescribed minimum value
 if(windspd < minwind) windspd=minwind
 ! compute relative humidity (-)
 relhum   = SPHM2RELHM(spechum, airpres, airtemp)
 ! if relative humidity exceeds saturation, then set relative and specific humidity to saturation
 if(relhum > 1._dp)then
  relhum  = 1._dp
  spechum = RELHM2SPHM(relhum, airpres, airtemp)
 endif
 ! compute wet bulb temperature
 twetbulb = WETBULBTMP(airtemp, relhum, airpres)
 ! ensure precipitation rate can be resolved by the data model
 if(pptrate<eps)then
  rainfall     = 0._dp
  snowfall     = 0._dp
  snowfallTemp = Tfreeze ! just so the value is populated
  return
 endif
 ! compute the maximum temperature of snow when the snow is predominantely frozen (K)
 maxFrozenSnowTemp = templiquid(unfrozenLiq,fc_param)
 ! compute fraction of rain and temperature of fresh snow
 Tmin = twetbulb - tempRangeTimestep/2._dp
 Tmax = twetbulb + tempRangeTimestep/2._dp
 if(Tmax < tempCritRain)then
  fracrain     = 0._dp
  snowfallTemp = twetbulb
 elseif(Tmin > tempCritRain)then
  fracrain     = 1._dp
  snowfallTemp = maxFrozenSnowTemp
 else
  fracrain     = (Tmax - tempCritRain)/(Tmax - Tmin)
  snowfallTemp = 0.5_dp*(Tmin + maxFrozenSnowTemp)
 endif
 ! ensure that snow falls at a temperature where all water 
 if(snowfallTemp < maxFrozenSnowTemp) snowfallTemp=maxFrozenSnowTemp
 ! ensure that snowfall temperature creates predominantely solid precipitation
 maxFrozenSnowTemp = templiquid(unfrozenLiq,fc_param)    ! snow temperature at fraction "unfrozenLiq" (K)
 snowfallTemp      = min(maxFrozenSnowTemp,snowfallTemp) ! snowfall temperature
 ! compute rainfall and snowfall
 rainfall = fracrain*pptrate
 snowfall = (1._dp - fracrain)*pptrate
 !print*, 'tempCritRain, tempRangeTimestep, pptrate, airtemp, rainfall, snowfall, twetbulb, relhum, snowfallTemp = '
 !print*, tempCritRain, tempRangeTimestep, pptrate, airtemp, rainfall, snowfall, twetbulb, relhum, snowfallTemp
 end subroutine derivforce

end module derivforce_module

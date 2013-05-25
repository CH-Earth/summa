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
 USE nr_utility_module,only:erf                              ! provide access to the error function
 USE multiconst,only:Tfreeze                                 ! freezing point of pure water (K)
 USE multiconst,only:secprhour                               ! number of seconds in an hour
 USE data_struc,only:data_step                               ! length of the data step (s)
 USE data_struc,only:time_data,forc_data                     ! forcing data structures
 USE data_struc,only:attr_data,mpar_data,mvar_data           ! model data structures
 USE var_lookup,only:iLookTIME,iLookATTR                     ! named variables for structure elements
 USE var_lookup,only:iLookPARAM,iLookFORCE,iLookMVAR         ! named variables for structure elements
 USE sunGeomtry_module,only:clrsky_rad                       ! compute cosine of the solar zenith angle
 USE conv_funcs_module,only:vapPress                         ! compute vapor pressure of air (Pa)
 USE conv_funcs_module,only:SPHM2RELHM,RELHM2SPHM,WETBULBTMP ! conversion functions
 USE snow_utils_module,only:fracliquid,templiquid            ! functions to compute temperature/liquid water
 ! compute derived forcing data variables
 implicit none
 ! dummy variables
 integer(i4b),intent(out)      :: err                        ! error code
 character(*),intent(out)      :: message                    ! error message
 ! variables for cosine of the solar zenith angle
 integer(i4b),pointer          :: im                         ! month
 integer(i4b),pointer          :: id                         ! day
 real(dp)                      :: ahour                      ! hour at start of time step
 real(dp)                      :: dataStep                   ! data step (hours)
 real(dp),parameter            :: slope=0._dp                ! terrain slope (assume flat)
 real(dp),parameter            :: azimuth=0._dp              ! terrain azimuth (assume zero)
 real(dp),pointer              :: latitude                   ! latitude (degrees north)
 real(dp)                      :: hri                        ! average radiation index over time step DT
 real(dp),pointer              :: cosZenith                  ! average cosine of the zenith angle over time step DT
 ! local pointers to model parameters
 real(dp),pointer              :: Frad_vis                   ! fraction radiation absorbed in visible part of spectrum (-)
 real(dp),pointer              :: Frad_direct                ! fraction direct radiation (-)
 real(dp),pointer              :: minwind                    ! minimum windspeed (m s-1)
 real(dp),pointer              :: fc_param                   ! freezing curve parameter for snow (K-1)
 real(dp),pointer              :: tempCritRain               ! critical temperature where precipitation is rain (K)
 real(dp),pointer              :: tempRangeTimestep          ! temperature range over the time step (K)
 real(dp),pointer              :: frozenPrecipMultip         ! frozen precipitation multiplier (-)
 real(dp),pointer              :: newSnowDenMin              ! minimum new snow density (kg m-3)
 real(dp),pointer              :: newSnowDenMult             ! multiplier for new snow density (kg m-3)
 real(dp),pointer              :: newSnowDenScal             ! scaling factor for new snow density (K)
 ! local pointers to model forcing data
 real(dp),pointer              :: SWRadAtm                   ! downward shortwave radiation (W m-2)
 real(dp),pointer              :: airtemp                    ! air temperature at 2 meter height (K)
 real(dp),pointer              :: windspd                    ! wind speed at 10 meter height (m s-1)
 real(dp),pointer              :: airpres                    ! air pressure at 2 meter height (Pa)
 real(dp),pointer              :: spechum                    ! specific humidity at 2 meter height (g g-1)
 real(dp),pointer              :: pptrate                    ! precipitation rate (kg m-2 s-1)
 ! local pointers to derived model forcing data
 real(dp),pointer              :: scalarO2air                ! atmospheric o2 concentration (Pa)
 real(dp),pointer              :: scalarCO2air               ! atmospheric co2 concentration (Pa)
 ! local pointers to model variables
 real(dp),pointer              :: spectralIncomingDirect(:)  ! downwelling direct shortwave radiation in each wave band (W m-2)
 real(dp),pointer              :: spectralIncomingDiffuse(:) ! downwelling diffuse shortwave radiation in each wave band (W m-2)
 real(dp),pointer              :: VPair                      ! vapor pressure of the air above the vegetation canopy (Pa)
 real(dp),pointer              :: twetbulb                   ! wet bulb temperature (K)
 real(dp),pointer              :: rainfall                   ! computed rainfall rate (kg m-2 s-1)
 real(dp),pointer              :: snowfall                   ! computed snowfall rate (kg m-2 s-1)
 real(dp),pointer              :: snowfallTemp               ! computed temperature of fresh snow (K)
 real(dp),pointer              :: newSnowDensity             ! computed density of fresh snow (kg m-3)
 ! local variables
 real(dp),parameter            :: valueMissing=-9999._dp     ! missing value
 real(dp),parameter            :: co2Factor=355.e-6_dp       ! empirical factor to obtain partial pressure of co2
 real(dp),parameter            :: o2Factor=0.209_dp          ! empirical factor to obtain partial pressure of o2
 real(dp)                      :: relhum                     ! relative humidity (-)
 real(dp)                      :: fracrain                   ! fraction of precipitation that falls as rain
 real(dp)                      :: maxFrozenSnowTemp          ! maximum temperature of snow when the snow is predominantely frozen (K)
 real(dp),parameter            :: unfrozenLiq=0.01_dp        ! unfrozen liquid water used to compute maxFrozenSnowTemp (-)
 real(dp),parameter            :: eps=epsilon(fracrain)      ! a number that is almost negligible
 real(dp)                      :: Tmin,Tmax                  ! minimum and maximum wet bulb temperature in the time step (K)
 ! initialize error control
 err=0; message="f-derivforce/"
 ! assign pointers to model parameters
 Frad_vis           => mpar_data%var(iLookPARAM%Frad_vis)           ! fraction radiation absorbed in visible part of spectrum (-)
 Frad_direct        => mpar_data%var(iLookPARAM%Frad_direct)        ! fraction direct radiation (-)
 minwind            => mpar_data%var(iLookPARAM%minwind)            ! minimum windspeed (m s-1)
 fc_param           => mpar_data%var(iLookPARAM%snowfrz_scale)      ! freezing curve parameter for snow (K-1)
 tempCritRain       => mpar_data%var(iLookPARAM%tempCritRain)       ! critical temperature where precipitation is rain (K)
 tempRangeTimestep  => mpar_data%var(iLookPARAM%tempRangeTimestep)  ! temperature range over the time step (K)
 frozenPrecipMultip => mpar_data%var(iLookPARAM%frozenPrecipMultip) ! frozen precipitation multiplier (-)
 newSnowDenMin      => mpar_data%var(iLookPARAM%newSnowDenMin)      ! minimum new snow density (kg m-3)
 newSnowDenMult     => mpar_data%var(iLookPARAM%newSnowDenMult)     ! multiplier for new snow density (kg m-3)
 newSnowDenScal     => mpar_data%var(iLookPARAM%newSnowDenScal)     ! scaling factor for new snow density (K)
 ! assign pointers to radiation geometry variables
 im        => time_data%var(iLookTIME%im)                           ! month
 id        => time_data%var(iLookTIME%id)                           ! day
 dataStep   = data_step/secprhour                                   ! time step (hours)
 ahour      = real(time_data%var(iLookTIME%ih),kind(dp)) - dataStep ! hour at start of time step
 latitude  => attr_data%var(iLookATTR%latitude)                     ! latitude (degrees north
 cosZenith => mvar_data%var(iLookMVAR%scalarCosZenith)%dat(1)       ! average cosine of the zenith angle over time step DT
 ! assign pointers to model forcing data
 SWRadAtm => forc_data%var(iLookFORCE%SWRadAtm)                     ! downward shortwave radiation (W m-2)
 airtemp  => forc_data%var(iLookFORCE%airtemp)                      ! air temperature at 2 meter height (K)
 windspd  => forc_data%var(iLookFORCE%windspd)                      ! wind speed at 10 meter height (m s-1)
 airpres  => forc_data%var(iLookFORCE%airpres)                      ! air pressure at 2 meter height (Pa)
 spechum  => forc_data%var(iLookFORCE%spechum)                      ! specific humidity at 2 meter height (g g-1)
 pptrate  => forc_data%var(iLookFORCE%pptrate)                      ! precipitation rate (kg m-2 s-1)
 ! assign pointers to derived model forcing data
 scalarO2air  => mvar_data%var(iLookMVAR%scalarO2air)%dat(1)        ! atmospheric o2 concentration (Pa)
 scalarCO2air => mvar_data%var(iLookMVAR%scalarCO2air)%dat(1)       ! atmospheric co2 concentration (Pa)
 ! assign pointers to radiation variables
 spectralIncomingDirect  => mvar_data%var(iLookMVAR%spectralIncomingDirect)%dat     ! downwelling direct shortwave radiation for each waveband (W m-2)
 spectralIncomingDiffuse => mvar_data%var(iLookMVAR%spectralIncomingDiffuse)%dat    ! downwelling diffuse shortwave radiation for each waveband (W m-2)
 if(size(spectralIncomingDirect) /= 2 .or. size(spectralIncomingDiffuse) /= 2)then
  err=20; message=trim(message)//'expect only two spectral classes for radiation'; return
 endif
 ! assign pointers to snow accumulation variables
 VPair          => mvar_data%var(iLookMVAR%scalarVPair)%dat(1)          ! vapor pressure of the air above the vegetation canopy (Pa)
 twetbulb       => mvar_data%var(iLookMVAR%scalarTwetbulb)%dat(1)       ! wet bulb temperature (K)
 rainfall       => mvar_data%var(iLookMVAR%scalarRainfall)%dat(1)       ! computed rainfall rate (kg m-2 s-1)
 snowfall       => mvar_data%var(iLookMVAR%scalarSnowfall)%dat(1)       ! computed snowfall rate (kg m-2 s-1)
 snowfallTemp   => mvar_data%var(iLookMVAR%scalarSnowfallTemp)%dat(1)   ! computed temperature of fresh snow (K)
 newSnowDensity => mvar_data%var(iLookMVAR%scalarNewSnowDensity)%dat(1) ! computed density of new snow (kg m-3) 

 ! compute the partial pressure of o2 and co2
 scalarCO2air = co2Factor * airpres  ! atmospheric co2 concentration (Pa)
 scalarO2air  = o2Factor * airpres   ! atmospheric o2 concentration (Pa)

 ! compute the cosine of the solar zenith angle
 call clrsky_rad(im,id,ahour,dataStep,   &  ! intent(in): time variables
                 slope,azimuth,latitude, &  ! intent(in): location variables
                 hri,cosZenith)             ! intent(out): cosine of the solar zenith angle
 ! ensure solar radiation is zero between sunset and sunrise
 if(cosZenith <= 0._dp) SWRadAtm = 0._dp
 ! compute direct shortwave radiation, in the visible and near-infra-red part of the spectrum
 spectralIncomingDirect(1) = SWRadAtm*Frad_direct*Frad_vis                         ! (direct vis)
 spectralIncomingDirect(2) = SWRadAtm*Frad_direct*(1._dp - Frad_vis)               ! (direct nir)
 ! compute diffuse shortwave radiation, in the visible and near-infra-red part of the spectrum
 spectralIncomingDiffuse(1) = SWRadAtm*(1._dp - Frad_direct)*Frad_vis              ! (diffuse vis)
 spectralIncomingDiffuse(2) = SWRadAtm*(1._dp - Frad_direct)*(1._dp - Frad_vis)    ! (diffuse nir)

 ! ensure wind speed is above a prescribed minimum value
 if(windspd < minwind) windspd=minwind

 ! compute relative humidity (-)
 relhum   = SPHM2RELHM(spechum, airpres, airtemp)
 ! if relative humidity exceeds saturation, then set relative and specific humidity to saturation
 if(relhum > 1._dp)then
  relhum  = 1._dp
  spechum = RELHM2SPHM(relhum, airpres, airtemp)
 endif

 ! compute vapor pressure of the air above the vegetation canopy (Pa)
 VPair = vapPress(spechum,airpres)
 !print*, 'VPair = ', VPair

 ! compute wet bulb temperature (K)
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
 snowfall = (1._dp - fracrain)*pptrate*frozenPrecipMultip
 !print*, 'tempCritRain, tempRangeTimestep, pptrate, airtemp, rainfall, snowfall, twetbulb, relhum, snowfallTemp = '
 !print*, tempCritRain, tempRangeTimestep, pptrate, airtemp, rainfall, snowfall, twetbulb, relhum, snowfallTemp

 ! compute density of new snow
 if(snowfall > tiny(fracrain))then
  newSnowDensity = newSnowDenMin + newSnowDenMult*exp((airtemp-Tfreeze)/newSnowDenScal)  ! new snow density (kg m-3)
 else
  newSnowDensity = valueMissing
  rainfall = rainfall + snowfall ! in most cases snowfall will be zero here
  snowfall = 0._dp
 endif

 end subroutine derivforce

end module derivforce_module

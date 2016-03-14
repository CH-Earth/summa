! SUMMA - Structure for Unifying Multiple Modeling Alternatives
! Copyright (C) 2014-2015 NCAR/RAL
!
! This file is part of SUMMA
!
! For more information see: http://www.ral.ucar.edu/projects/summa
!
! This program is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! This program is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License
! along with this program.  If not, see <http://www.gnu.org/licenses/>.

module derivforce_module
USE nrtype
implicit none
private
public::derivforce
contains


 ! ************************************************************************************************
 ! public subroutine derivforce: compute derived forcing data
 ! ************************************************************************************************
 subroutine derivforce(time_data,forc_data,attr_data,mpar_data,diag_data,flux_data,err,message)
 USE multiconst,only:Tfreeze                                 ! freezing point of pure water (K)
 USE multiconst,only:secprhour                               ! number of seconds in an hour
 USE multiconst,only:minprhour                               ! number of minutes in an hour
 USE data_struc,only:data_step                               ! length of the data step (s)
 USE data_struc,only:var_dlength                             ! data structure: x%var(:)%dat (dp)
 USE var_lookup,only:iLookTIME,iLookATTR                     ! named variables for structure elements
 USE var_lookup,only:iLookPARAM,iLookFORCE,iLookDIAG,iLookFLUX  ! named variables for structure elements
 USE sunGeomtry_module,only:clrsky_rad                       ! compute cosine of the solar zenith angle
 USE conv_funcs_module,only:vapPress                         ! compute vapor pressure of air (Pa)
 USE conv_funcs_module,only:SPHM2RELHM,RELHM2SPHM,WETBULBTMP ! conversion functions
 USE snow_utils_module,only:fracliquid,templiquid            ! functions to compute temperature/liquid water
 ! compute derived forcing data variables
 implicit none
 ! input variables
 integer(i4b),intent(in)         :: time_data(:)             ! vector of time data for a given time step
 real(dp),    intent(inout)      :: forc_data(:)             ! vector of forcing data for a given time step
 real(dp),    intent(in)         :: attr_data(:)             ! vector of model attributes
 real(dp),    intent(in)         :: mpar_data(:)             ! vector of model parameters
 ! output variables
 type(var_dlength),intent(inout) :: diag_data                ! data structure of model diagnostic variables for a local HRU
 type(var_dlength),intent(inout) :: flux_data                ! data structure of model fluxes for a local HRU
 integer(i4b),intent(out)        :: err                      ! error code
 character(*),intent(out)        :: message                  ! error message
 ! variables for cosine of the solar zenith angle
 real(dp)                        :: ahour                    ! hour at start of time step
 real(dp)                        :: dataStep                 ! data step (hours)
 real(dp),parameter              :: slope=0._dp              ! terrain slope (assume flat)
 real(dp),parameter              :: azimuth=0._dp            ! terrain azimuth (assume zero)
 real(dp)                        :: hri                      ! average radiation index over time step DT
 ! local variables
 real(dp),parameter              :: valueMissing=-9999._dp   ! missing value
 real(dp),parameter              :: co2Factor=355.e-6_dp     ! empirical factor to obtain partial pressure of co2
 real(dp),parameter              :: o2Factor=0.209_dp        ! empirical factor to obtain partial pressure of o2
 real(dp)                        :: relhum                   ! relative humidity (-)
 real(dp)                        :: fracrain                 ! fraction of precipitation that falls as rain
 real(dp)                        :: maxFrozenSnowTemp        ! maximum temperature of snow when the snow is predominantely frozen (K)
 real(dp),parameter              :: unfrozenLiq=0.01_dp      ! unfrozen liquid water used to compute maxFrozenSnowTemp (-)
 real(dp),parameter              :: eps=epsilon(fracrain)    ! a number that is almost negligible
 real(dp)                        :: Tmin,Tmax                ! minimum and maximum wet bulb temperature in the time step (K)
 ! ************************************************************************************************
 ! associate local variables with the information in the data structures
 associate(&
 ! model parameters
 Frad_vis                => mpar_data(iLookPARAM%Frad_vis)                        , & ! fraction radiation absorbed in visible part of spectrum (-)
 directScale             => mpar_data(iLookPARAM%directScale)                     , & ! scaling factor for fractional driect radiaion parameterization (-)
 Frad_direct             => mpar_data(iLookPARAM%Frad_direct)                     , & ! maximum fraction direct radiation (-)
 minwind                 => mpar_data(iLookPARAM%minwind)                         , & ! minimum windspeed (m s-1)
 fc_param                => mpar_data(iLookPARAM%snowfrz_scale)                   , & ! freezing curve parameter for snow (K-1)
 tempCritRain            => mpar_data(iLookPARAM%tempCritRain)                    , & ! critical temperature where precipitation is rain (K)
 tempRangeTimestep       => mpar_data(iLookPARAM%tempRangeTimestep)               , & ! temperature range over the time step (K)
 frozenPrecipMultip      => mpar_data(iLookPARAM%frozenPrecipMultip)              , & ! frozen precipitation multiplier (-)
 newSnowDenMin           => mpar_data(iLookPARAM%newSnowDenMin)                   , & ! minimum new snow density (kg m-3)
 newSnowDenMult          => mpar_data(iLookPARAM%newSnowDenMult)                  , & ! multiplier for new snow density (kg m-3)
 newSnowDenScal          => mpar_data(iLookPARAM%newSnowDenScal)                  , & ! scaling factor for new snow density (K)
 ! radiation geometry variables
 im                      => time_data(iLookTIME%im)                               , & ! month
 id                      => time_data(iLookTIME%id)                               , & ! day
 ih                      => time_data(iLookTIME%ih)                               , & ! hour
 imin                    => time_data(iLookTIME%imin)                             , & ! minute
 latitude                => attr_data(iLookATTR%latitude)                         , & ! latitude (degrees north)
 cosZenith               => diag_data%var(iLookDIAG%scalarCosZenith)%dat(1)       , & ! average cosine of the zenith angle over time step DT
 ! model forcing data
 SWRadAtm                => forc_data(iLookFORCE%SWRadAtm)                        , & ! downward shortwave radiation (W m-2)
 airtemp                 => forc_data(iLookFORCE%airtemp)                         , & ! air temperature at 2 meter height (K)
 windspd                 => forc_data(iLookFORCE%windspd)                         , & ! wind speed at 10 meter height (m s-1)
 airpres                 => forc_data(iLookFORCE%airpres)                         , & ! air pressure at 2 meter height (Pa)
 spechum                 => forc_data(iLookFORCE%spechum)                         , & ! specific humidity at 2 meter height (g g-1)
 pptrate                 => forc_data(iLookFORCE%pptrate)                         , & ! precipitation rate (kg m-2 s-1)
 ! derived model forcing data
 scalarO2air             => diag_data%var(iLookDIAG%scalarO2air)%dat(1)           , & ! atmospheric o2 concentration (Pa)
 scalarCO2air            => diag_data%var(iLookDIAG%scalarCO2air)%dat(1)          , & ! atmospheric co2 concentration (Pa)
 ! radiation variables
 scalarFractionDirect    => diag_data%var(iLookDIAG%scalarFractionDirect)%dat(1)  , & ! fraction of direct radiation (0-1)
 spectralIncomingDirect  => flux_data%var(iLookFLUX%spectralIncomingDirect)%dat   , & ! downwelling direct shortwave radiation for each waveband (W m-2)
 spectralIncomingDiffuse => flux_data%var(iLookFLUX%spectralIncomingDiffuse)%dat  , & ! downwelling diffuse shortwave radiation for each waveband (W m-2)
 ! snow accumulation variables
 rainfall                => flux_data%var(iLookFLUX%scalarRainfall)%dat(1)        , & ! computed rainfall rate (kg m-2 s-1)
 snowfall                => flux_data%var(iLookFLUX%scalarSnowfall)%dat(1)        , & ! computed snowfall rate (kg m-2 s-1)
 VPair                   => diag_data%var(iLookDIAG%scalarVPair)%dat(1)           , & ! vapor pressure of the air above the vegetation canopy (Pa)
 twetbulb                => diag_data%var(iLookDIAG%scalarTwetbulb)%dat(1)        , & ! wet bulb temperature (K)
 snowfallTemp            => diag_data%var(iLookDIAG%scalarSnowfallTemp)%dat(1)    , & ! computed temperature of fresh snow (K)
 newSnowDensity          => diag_data%var(iLookDIAG%scalarNewSnowDensity)%dat(1)    & ! computed density of new snow (kg m-3)
 ) ! (associating local variables with the information in the data structures)

 ! initialize error control
 err=0; message="derivforce/"

 ! check spectral dimension
 if(size(spectralIncomingDirect) /= 2 .or. size(spectralIncomingDiffuse) /= 2)then
  err=20; message=trim(message)//'expect two spectral classes for radiation'; return
 endif

 ! compute the partial pressure of o2 and co2
 scalarCO2air = co2Factor * airpres  ! atmospheric co2 concentration (Pa)
 scalarO2air  = o2Factor * airpres   ! atmospheric o2 concentration (Pa)

 ! compute the decimal hour at the start of the time step
 dataStep = data_step/secprhour  ! time step (hours)
 ahour    = real(ih,kind(dp)) + real(imin,kind(dp))/minprhour - data_step/secprhour  ! decimal hour (start of the step)

 ! compute the cosine of the solar zenith angle
 call clrsky_rad(im,id,ahour,dataStep,   &  ! intent(in): time variables
                 slope,azimuth,latitude, &  ! intent(in): location variables
                 hri,cosZenith)             ! intent(out): cosine of the solar zenith angle
 !write(*,'(a,1x,4(i2,1x),3(f9.3,1x))') 'im,id,ih,imin,ahour,dataStep,cosZenith = ', &
 !                                       im,id,ih,imin,ahour,dataStep,cosZenith
 ! check that we don't have considerable shortwave when the zenith angle is low
 ! NOTE: this is likely because the data are not in local time
 !if(cosZenith < epsilon(cosZenith) .and. SWRadAtm > 200._dp)then
 ! message=trim(message)//'SWRadAtm > 200 W m-2 when cos zenith angle is zero -- check that forcing data are in local time, '//&
 !                        'that the time stamp in forcing data is at the end of the data interval, and that the lat-lon '//&
 !                        'in the site characteristix file is correct'
 ! err=20; return
 !endif
 ! ensure solar radiation is zero between sunset and sunrise
 ! NOTE: also ensure that sw radiation is positive
 if(cosZenith <= 0._dp .or. SWRadAtm < 0._dp) SWRadAtm = 0._dp
 ! compute the fraction of direct radiation using the parameterization of Nijssen and Lettenmaier (1999)
 if(cosZenith > 0._dp)then
  scalarFractionDirect = Frad_direct*cosZenith/(cosZenith + directScale)
 else
  scalarFractionDirect = 0._dp
 endif
 ! compute direct shortwave radiation, in the visible and near-infra-red part of the spectrum
 spectralIncomingDirect(1) = SWRadAtm*scalarFractionDirect*Frad_vis                         ! (direct vis)
 spectralIncomingDirect(2) = SWRadAtm*scalarFractionDirect*(1._dp - Frad_vis)               ! (direct nir)
 ! compute diffuse shortwave radiation, in the visible and near-infra-red part of the spectrum
 spectralIncomingDiffuse(1) = SWRadAtm*(1._dp - scalarFractionDirect)*Frad_vis              ! (diffuse vis)
 spectralIncomingDiffuse(2) = SWRadAtm*(1._dp - scalarFractionDirect)*(1._dp - Frad_vis)    ! (diffuse nir)

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
 !write(*,'(a,1x,10(f20.10,1x))') 'Tmin, twetbulb, tempRangeTimestep, tempCritRain = ', &
 !                                 Tmin, twetbulb, tempRangeTimestep, tempCritRain

 ! ensure that snowfall temperature creates predominantely solid precipitation
 snowfallTemp      = min(maxFrozenSnowTemp,snowfallTemp) ! snowfall temperature

 ! ensure precipitation rate can be resolved by the data model
 if(pptrate<eps)then
  ! set rainfall and snowfall to zero
  rainfall     = 0._dp
  snowfall     = 0._dp
 else
  ! compute rainfall and snowfall
  rainfall = fracrain*pptrate
  snowfall = (1._dp - fracrain)*pptrate*frozenPrecipMultip
 endif

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

 ! end association of local variables with the information in the data structures
 end associate

 end subroutine derivforce


end module derivforce_module

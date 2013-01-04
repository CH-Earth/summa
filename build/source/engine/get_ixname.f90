module get_ixname_module
! used to get the index of a named variable
USE nrtype                                          ! variable types, etc.
implicit none
private
public::get_ixdecisions
public::get_ixTime
public::get_ixAttr
public::get_ixType
public::get_ixForce
public::get_ixParam
public::get_ixMvar
public::get_ixIndex
contains

 ! *******************************************************************************************************************
 ! new function: get the index of the named variables for the model decisions
 ! *******************************************************************************************************************
 function get_ixdecisions(varName)
 USE var_lookup,only:iLookDECISIONS                  ! indices of the named variables
 implicit none
 ! define dummy variables
 character(*), intent(in) :: varName                 ! variable name
 integer(i4b)             :: get_ixdecisions         ! index of the named variable
 ! define local variables
 integer(i4b), parameter  :: imiss = -999            ! missing value
 ! get the index of the named variables
 select case(trim(varName))
  case('soilCatTbl'      ); get_ixdecisions=iLookDECISIONS%soilCatTbl  ! ( 1) soil-category dateset
  case('vegeParTbl'      ); get_ixdecisions=iLookDECISIONS%vegeParTbl  ! ( 2) vegetation category dataset
  case('num_method'      ); get_ixdecisions=iLookDECISIONS%num_method  ! ( 3) choice of numerical method
  case('fDerivMeth'      ); get_ixdecisions=iLookDECISIONS%fDerivMeth  ! ( 4) choice of method to calculate flux derivatives
  case('f_Richards'      ); get_ixdecisions=iLookDECISIONS%f_Richards  ! ( 5) form of Richards' equation
  case('groundwatr'      ); get_ixdecisions=iLookDECISIONS%groundwatr  ! ( 6) choice of groundwater parameterization
  case('hc_profile'      ); get_ixdecisions=iLookDECISIONS%hc_profile  ! ( 7) choice of hydraulic conductivity profile
  case('bcUpprTdyn'      ); get_ixdecisions=iLookDECISIONS%bcUpprTdyn  ! ( 8) type of upper boundary condition for thermodynamics
  case('bcLowrTdyn'      ); get_ixdecisions=iLookDECISIONS%bcLowrTdyn  ! ( 9) type of lower boundary condition for thermodynamics
  case('bcUpprSoiH'      ); get_ixdecisions=iLookDECISIONS%bcUpprSoiH  ! (10) type of upper boundary condition for soil hydrology
  case('bcLowrSoiH'      ); get_ixdecisions=iLookDECISIONS%bcLowrSoiH  ! (11) type of lower boundary condition for soil hydrology
  case('astability'      ); get_ixdecisions=iLookDECISIONS%astability  ! (12) choice of stability function
  case('compaction'      ); get_ixdecisions=iLookDECISIONS%compaction  ! (13) choice of compaction routine
  case('thermlcond'      ); get_ixdecisions=iLookDECISIONS%thermlcond  ! (14) choice of thermal conductivity representation
  case('alb_method'      ); get_ixdecisions=iLookDECISIONS%alb_method  ! (15) choice of albedo representation
  case('subRouting'      ); get_ixdecisions=iLookDECISIONS%subRouting  ! (16) choice of method for sub-grid routing
  ! get to here if cannot find the variable
  case default
   get_ixdecisions = imiss
 endselect
 end function get_ixdecisions


 ! *******************************************************************************************************************
 ! new function: get the index of the named variables for the model time
 ! *******************************************************************************************************************
 function get_ixtime(varName)
 USE var_lookup,only:iLookTIME                       ! indices of the named variables
 implicit none
 ! define dummy variables
 character(*), intent(in) :: varName                 ! variable name
 integer(i4b)             :: get_ixtime              ! index of the named variable
 ! define local variables
 integer(i4b), parameter  :: imiss = -999            ! missing value
 ! get the index of the named variables
 select case(trim(varName))
  case('iyyy'            ); get_ixtime = iLookTIME%iyyy             ! year
  case('im'              ); get_ixtime = iLookTIME%im               ! month
  case('id'              ); get_ixtime = iLookTIME%id               ! day
  case('ih'              ); get_ixtime = iLookTIME%ih               ! hour
  case('imin'            ); get_ixtime = iLookTIME%imin             ! minute
  ! get to here if cannot find the variable
  case default
   get_ixtime = imiss
 endselect
 end function get_ixtime


 ! *******************************************************************************************************************
 ! new function: get the index of the named variables for the model forcing data
 ! *******************************************************************************************************************
 function get_ixforce(varName)
 USE var_lookup,only:iLookFORCE                      ! indices of the named variables
 implicit none
 ! define dummy variables
 character(*), intent(in) :: varName                 ! variable name
 integer(i4b)             :: get_ixforce             ! index of the named variable
 ! define local variables
 integer(i4b), parameter  :: imiss = -999            ! missing value
 ! get the index of the named variables
 select case(trim(varName))
  case('time'            ); get_ixforce = iLookFORCE%time             ! time since time reference       (s)
  case('pptrate'         ); get_ixforce = iLookFORCE%pptrate          ! precipitation rate              (kg m-2 s-1)
  case('airtemp'         ); get_ixforce = iLookFORCE%airtemp          ! air temperature                 (K)
  case('spechum'         ); get_ixforce = iLookFORCE%spechum          ! specific humidity               (g/g)
  case('windspd'         ); get_ixforce = iLookFORCE%windspd          ! windspeed                       (m/s)
  case('sw_down'         ); get_ixforce = iLookFORCE%sw_down          ! downwelling shortwave radiaiton (W m-2)
  case('lw_down'         ); get_ixforce = iLookFORCE%lw_down          ! downwelling longwave radiation  (W m-2)
  case('airpres'         ); get_ixforce = iLookFORCE%airpres          ! pressure                        (Pa)
  ! get to here if cannot find the variable
  case default
   get_ixforce = imiss
 endselect
 end function get_ixforce


 ! *******************************************************************************************************************
 ! new function: get the index of the named variables for the site characteristics
 ! *******************************************************************************************************************
 function get_ixAttr(varName)
 USE var_lookup,only:iLookATTR                       ! indices of the named variables
 implicit none
 ! define dummy variables
 character(*), intent(in) :: varName                 ! variable name
 integer(i4b)             :: get_ixAttr              ! index of the named variable
 ! define local variables
 integer(i4b), parameter  :: imiss = -999            ! missing value
 ! get the index of the named variables
 select case(trim(varName))
  case('latitude'   ); get_ixAttr = iLookATTR%latitude       ! latitude    (degrees north)
  case('longitude'  ); get_ixAttr = iLookATTR%longitude      ! longitude   (degrees east)
  case('elevation'  ); get_ixAttr = iLookATTR%elevation      ! elevation   (m)
  ! get to here if cannot find the variable
  case default
   get_ixAttr = imiss
 endselect
 end function get_ixAttr


 ! *******************************************************************************************************************
 ! new function: get the index of the named variables for the local classification of veg, soil, etc.
 ! *******************************************************************************************************************
 function get_ixType(varName)
 USE var_lookup,only:iLookTYPE                       ! indices of the named variables
 implicit none
 ! define dummy variables
 character(*), intent(in) :: varName                 ! variable name
 integer(i4b)             :: get_ixType              ! index of the named variable
 ! define local variables
 integer(i4b), parameter  :: imiss = -999            ! missing value
 ! get the index of the named variables
 select case(trim(varName))
  case('vegTypeIndex'   ); get_ixType = iLookTYPE%vegTypeIndex       ! index defining vegetation type
  case('soilTypeIndex'  ); get_ixType = iLookTYPE%soilTypeIndex      ! index defining soil type
  case('slopeTypeIndex' ); get_ixType = iLookTYPE%slopeTypeIndex     ! index defining slope
  ! get to here if cannot find the variable
  case default
   get_ixType = imiss
 endselect
 end function get_ixType


 ! *******************************************************************************************************************
 ! new function: get the index of the named variables for the model parameters
 ! *******************************************************************************************************************
 function get_ixparam(varName)
 USE var_lookup,only:iLookPARAM                      ! indices of the named variables
 implicit none
 ! define dummy variables
 character(*), intent(in) :: varName                 ! variable name
 integer(i4b)             :: get_ixparam             ! index of the named variable
 ! define local variables
 integer(i4b), parameter  :: imiss = -999            ! missing value
 ! get the index of the named variables
 select case(trim(varName))
  ! boundary conditions
  case('upperBoundHead'      ); get_ixparam = iLookPARAM%upperBoundHead       ! matric head of the upper boundary (m)
  case('lowerBoundHead'      ); get_ixparam = iLookPARAM%lowerBoundHead       ! matric head of the lower boundary (m)
  case('upperBoundTheta'     ); get_ixparam = iLookPARAM%upperBoundTheta      ! volumetric liquid water content at the upper boundary (-)
  case('lowerBoundTheta'     ); get_ixparam = iLookPARAM%lowerBoundTheta      ! volumetric liquid water content at the lower boundary (-)
  case('upperBoundTemp'      ); get_ixparam = iLookPARAM%upperBoundTemp       ! temperature of the upper boundary (K)
  case('lowerBoundTemp'      ); get_ixparam = iLookPARAM%lowerBoundTemp       ! temperature of the lower boundary (K)
  ! precipitation partitioning
  case('tempCritRain'        ); get_ixparam = iLookPARAM%tempCritRain         ! critical temperature where precipitation is rain (K)
  case('tempRangeTimestep'   ); get_ixparam = iLookPARAM%tempRangeTimestep    ! temperature range over the time step (K)
  ! freezing curve for snow
  case('snowfrz_scale'       ); get_ixparam = iLookPARAM%snowfrz_scale        ! scaling parameter for the freezing curve for snow (K-1)
  ! snow albedo
  case('snw_crit'            ); get_ixparam = iLookPARAM%snw_crit             ! critical mass necessary for albedo refreshment (kg m-2)
  case('alb_fresh'           ); get_ixparam = iLookPARAM%alb_fresh            ! fresh snow albedo (-)
  case('alb_dry'             ); get_ixparam = iLookPARAM%alb_dry              ! minimum snow albedo during winter (-)
  case('alb_wet'             ); get_ixparam = iLookPARAM%alb_wet              ! minimum snow albedo during spring (-)
  case('alb_decay'           ); get_ixparam = iLookPARAM%alb_decay            ! temporal decay factor for snow albedo (s-1)
  case('alb_scale'           ); get_ixparam = iLookPARAM%alb_scale            ! albedo scaling factor (s)
  case('soot_load'           ); get_ixparam = iLookPARAM%soot_load            ! temporal decay in snow albedo associated with the soot load (days-1)
  ! radiation transfer
  case('radExt_snow'         ); get_ixparam = iLookPARAM%radExt_snow          ! extinction coefficient for radiation penetration within the snowpack (m-1)
  case('Frad_vis'            ); get_ixparam = iLookPARAM%Frad_vis             ! fraction of radiation in the visible part of the spectrum (-)
  ! new snow density
  case('newSnowDenMin'       ); get_ixparam = iLookPARAM%newSnowDenMin        ! minimum new snow density (kg m-3)
  case('newSnowDenMult'      ); get_ixparam = iLookPARAM%newSnowDenMult       ! multiplier for new snow density (kg m-3)
  case('newSnowDenScal'      ); get_ixparam = iLookPARAM%newSnowDenScal       ! scaling factor for new snow density (K)
  ! snow compaction
  case('densScalGrowth'      ); get_ixparam = iLookPARAM%densScalGrowth       ! density scaling factor for grain growth (kg-1 m3)
  case('tempScalGrowth'      ); get_ixparam = iLookPARAM%tempScalGrowth       ! temperature scaling factor for grain growth (K-1)
  case('grainGrowthRate'     ); get_ixparam = iLookPARAM%grainGrowthRate      ! rate of grain growth (s-1)
  case('densScalOvrbdn'      ); get_ixparam = iLookPARAM%densScalOvrbdn       ! density scaling factor for overburden pressure (kg-1 m3)
  case('tempScalOvrbdn'      ); get_ixparam = iLookPARAM%tempScalOvrbdn       ! temperature scaling factor for overburden pressure (K-1)
  case('base_visc'           ); get_ixparam = iLookPARAM%base_visc            ! viscosity coefficient at T=T_frz and snow density=0  (kg s m-2)
  ! water flow through snow
  case('Fcapil'              ); get_ixparam = iLookPARAM%Fcapil               ! capillary retention as a fraction of the total pore volume (-)
  case('k_snow'              ); get_ixparam = iLookPARAM%k_snow               ! hydraulic conductivity of snow (m s-1), 0.0055 = approx. 20 m/hr, from UEB
  case('mw_exp'              ); get_ixparam = iLookPARAM%mw_exp               ! exponent for meltwater flow (-)
  ! turbulent heat fluxes
  case('mheight'             ); get_ixparam = iLookPARAM%mheight              ! measurement height (m)
  case('zon'                 ); get_ixparam = iLookPARAM%zon                  ! roughness length (m)
  case('c_star'              ); get_ixparam = iLookPARAM%c_star               ! parameter in Louis (1979) stability function
  case('bparam'              ); get_ixparam = iLookPARAM%bparam               ! parameter in Louis (1979) stability function
  case('Mahrt_m'             ); get_ixparam = iLookPARAM%Mahrt_m              ! the m parameter from the Mahrt (1987) stability function
  ! vegetation properties
  case('rootingDepth'        ); get_ixparam = iLookPARAM%rootingDepth         ! rooting depth (m)
  case('rootDistExp'         ); get_ixparam = iLookPARAM%rootDistExp          ! exponent for the vertical distriution of root density (-)
  case('minStomatalResist'   ); get_ixparam = iLookPARAM%minStomatalResist    ! minimum stomatal resistance (s m-1)
  case('maxStomatalResist'   ); get_ixparam = iLookPARAM%maxStomatalResist    ! maximum stomatal resistance (s m-1)
  case('plantWiltPsi'        ); get_ixparam = iLookPARAM%plantWiltPsi         ! critical matric head when stomatal resitance 2 x min (m)
  case('plantWiltExp'        ); get_ixparam = iLookPARAM%plantWiltExp         ! empirical exponent in plant wilting factor expression (-)
  case('critAquiferTranspire'); get_ixparam = iLookPARAM%critAquiferTranspire ! critical aquifer storage value when transpiration is limited (m)
  ! soil properties
  case('soilAlbedo'          ); get_ixparam = iLookPARAM%soilAlbedo           ! soil albedo (-)
  case('soil_dens_intr'      ); get_ixparam = iLookPARAM%soil_dens_intr       ! intrinsic soil density (kg m-3)
  case('frac_sand'           ); get_ixparam = iLookPARAM%frac_sand            ! fraction of sand (-)
  case('frac_silt'           ); get_ixparam = iLookPARAM%frac_silt            ! fraction of silt (-)
  case('frac_clay'           ); get_ixparam = iLookPARAM%frac_clay            ! fraction of clay (-)
  case('theta_sat'           ); get_ixparam = iLookPARAM%theta_sat            ! soil porosity (-)
  case('theta_res'           ); get_ixparam = iLookPARAM%theta_res            ! volumetric residual water content (-)
  case('vGn_alpha'           ); get_ixparam = iLookPARAM%vGn_alpha            ! van Genuchten "alpha" parameter (m-1)
  case('vGn_n'               ); get_ixparam = iLookPARAM%vGn_n                ! van Genuchten "n" parameter (-) 
  case('k_soil'              ); get_ixparam = iLookPARAM%k_soil               ! saturated hydraulic conductivity (m s-1)
  case('kAnisotropic'        ); get_ixparam = iLookPARAM%kAnisotropic         ! anisotropy factor for lateral hydraulic conductivity (-)
  case('zScale_TOPMODEL'     ); get_ixparam = iLookPARAM%zScale_TOPMODEL      ! scale factor for TOPMODEL-ish baseflow parameterization (m)
  case('compactedDepth'      ); get_ixparam = iLookPARAM%compactedDepth       ! depth where k_soil reaches the compacted value given by CH78 (m)
  case('bpar_VIC'            ); get_ixparam = iLookPARAM%bpar_VIC             ! b-parameter in the VIC surface runoff parameterization (-)
  case('specificYield'       ); get_ixparam = iLookPARAM%specificYield        ! specific yield (-)
  case('specificStorage'     ); get_ixparam = iLookPARAM%specificStorage      ! specific storage coefficient (m-1)
  case('aquiferScaleFactor'  ); get_ixparam = iLookPARAM%aquiferScaleFactor   ! scaling factor for aquifer storage in the big bucket (m)
  case('bucketBaseflowExp'   ); get_ixparam = iLookPARAM%bucketBaseflowExp    ! baseflow exponent for the big bucket (-)
  case('f_impede'            ); get_ixparam = iLookPARAM%f_impede             ! ice impedence factor (-)
  ! sub-grid routing
  case('routingGammaShape'   ); get_ixparam = iLookPARAM%routingGammaShape    ! shape parameter in Gamma distribution used for sub-grid routing (-)
  case('routingGammaScale'   ); get_ixparam = iLookPARAM%routingGammaScale    ! scale parameter in Gamma distribution used for sub-grid routing (s)
  ! algorithmic control parameters
  case('minwind'             ); get_ixparam = iLookPARAM%minwind              ! minimum wind speed (m s-1)
  case('minstep'             ); get_ixparam = iLookPARAM%minstep              ! minimum length of the time step
  case('maxstep'             ); get_ixparam = iLookPARAM%maxstep              ! maximum length of the time step
  case('wimplicit'           ); get_ixparam = iLookPARAM%wimplicit            ! weight assigned to start-of-step fluxes
  case('maxiter'             ); get_ixparam = iLookPARAM%maxiter              ! maximum number of iterations
  case('relConvTol_liquid'   ); get_ixparam = iLookPARAM%relConvTol_liquid    ! relative convergence tolerance for vol frac liq water (-)
  case('absConvTol_liquid'   ); get_ixparam = iLookPARAM%absConvTol_liquid    ! absolute convergence tolerance for vol frac liq water (-)
  case('relConvTol_matric'   ); get_ixparam = iLookPARAM%relConvTol_matric    ! relative convergence tolerance for matric head (-)
  case('absConvTol_matric'   ); get_ixparam = iLookPARAM%absConvTol_matric    ! absolute convergence tolerance for matric head (m)
  case('relConvTol_energy'   ); get_ixparam = iLookPARAM%relConvTol_energy    ! relative convergence tolerance for energy (-)
  case('absConvTol_energy'   ); get_ixparam = iLookPARAM%absConvTol_energy    ! absolute convergence tolerance for energy (J m-3)
  case('relConvTol_aquifr'   ); get_ixparam = iLookPARAM%relConvTol_aquifr    ! relative convergence tolerance for aquifer storage (-)
  case('absConvTol_aquifr'   ); get_ixparam = iLookPARAM%absConvTol_aquifr    ! absolute convergence tolerance for aquifer storage (m)
  case('zmin'                ); get_ixparam = iLookPARAM%zmin                 ! minimum layer depth (m)
  case('zmax'                ); get_ixparam = iLookPARAM%zmax                 ! maximum layer depth (m)
  ! get to here if cannot find the variable
  case default
   get_ixparam = imiss
 endselect
 end function get_ixparam


 ! *******************************************************************************************************************
 ! new function: get the index of the named variables for the model variables
 ! *******************************************************************************************************************
 function get_ixmvar(varName)
 USE var_lookup,only:iLookMVAR                       ! indices of the named variables
 implicit none
 ! define dummy variables
 character(*), intent(in) :: varName                 ! parameter name
 integer(i4b)             :: get_ixmvar              ! index of the named variable
 ! define local variables
 integer(i4b), parameter  :: imiss = -999            ! missing value
 ! get the index of the named variables
 select case(trim(varName))
  ! define timestep-average fluxes for a few key variables
  case('averageMassLiquid'           ); get_ixmvar = iLookMVAR%averageMassLiquid           ! evaporation or dew (kg m-2 s-1)
  case('averageMassSolid'            ); get_ixmvar = iLookMVAR%averageMassSolid            ! sublimation or frost (kg m-2 s-1)
  case('averageRainPlusMelt'         ); get_ixmvar = iLookMVAR%averageRainPlusMelt         ! rain plus melt, as input to soil before calculating surface runoff (m s-1)
  case('averageSurfaceRunoff'        ); get_ixmvar = iLookMVAR%averageSurfaceRunoff        ! surface runoff (m s-1)
  case('averageSoilInflux'           ); get_ixmvar = iLookMVAR%averageSoilInflux           ! influx of water at the top of the soil profile (m s-1)
  case('averageSoilBaseflow'         ); get_ixmvar = iLookMVAR%averageSoilBaseflow         ! total baseflow from throughout the soil profile (m s-1)
  case('averageSoilDrainage'         ); get_ixmvar = iLookMVAR%averageSoilDrainage         ! drainage from the bottom of the soil profile (m s-1)
  case('averageSoilEjection'         ); get_ixmvar = iLookMVAR%averageSoilEjection         ! ejected water from the soil matrix (m s-1)
  case('averageAquiferRecharge'      ); get_ixmvar = iLookMVAR%averageAquiferRecharge      ! recharge to the aquifer (m s-1)
  case('averageAquiferBaseflow'      ); get_ixmvar = iLookMVAR%averageAquiferBaseflow      ! baseflow from the aquifer (m s-1)
  ! NOAH-MP vegetation variables
  case('scalarCanopyHeight'          ); get_ixmvar = iLookMVAR%scalarCanopyHeight          ! height of the top of the canopy layer (m)
  case('scalarVegetationTemp'        ); get_ixmvar = iLookMVAR%scalarVegetationTemp        ! vegetation temperature (K)
  case('scalarRootZoneTemp'          ); get_ixmvar = iLookMVAR%scalarRootZoneTemp          ! average temperature of the root zone (K)
  case('scalarLAI'                   ); get_ixmvar = iLookMVAR%scalarLAI                   ! one-sided leaf area index (m2 m-2)
  case('scalarSAI'                   ); get_ixmvar = iLookMVAR%scalarSAI                   ! one-sided stem area index (m2 m-2)
  case('scalarEffectiveLAI'          ); get_ixmvar = iLookMVAR%scalarEffectiveLAI          ! effective leaf area index after burial by snow (m2 m-2)
  case('scalarEffectiveSAI'          ); get_ixmvar = iLookMVAR%scalarEffectiveSAI          ! effective stem area index after burial by snow(m2 m-2)
  case('scalarGrowingSeasonIndex'    ); get_ixmvar = iLookMVAR%scalarGrowingSeasonIndex    ! growing season index (0=off, 1=on)
  ! scalar variables
  case('scalarSwDownVis'             ); get_ixmvar = iLookMVAR%scalarSwDownVis             ! downwelling shortwave radiation in visible part of spectrum (W m-2)
  case('scalarSwDownNir'             ); get_ixmvar = iLookMVAR%scalarSwDownNir             ! downwelling shortwave radiation in near-infrared part of spectrum (W m-2)
  case('scalarTwetbulb'              ); get_ixmvar = iLookMVAR%scalarTwetbulb              ! wetbulb temperature (K)
  case('scalarRainfall'              ); get_ixmvar = iLookMVAR%scalarRainfall              ! computed rainfall rate (kg m-2 s-1)
  case('scalarSnowfall'              ); get_ixmvar = iLookMVAR%scalarSnowfall              ! computed snowfall rate (kg m-2 s-1)
  case('scalarSnowfallTemp'          ); get_ixmvar = iLookMVAR%scalarSnowfallTemp          ! temperature of fresh snow (K)
  case('scalarSurfaceTemp'           ); get_ixmvar = iLookMVAR%scalarSurfaceTemp           ! temperature of the top layer (K)
  case('scalarAlbedo'                ); get_ixmvar = iLookMVAR%scalarAlbedo                ! albedo of the surface, soil or snow (-)
  case('scalarSnowDepth'             ); get_ixmvar = iLookMVAR%scalarSnowDepth             ! total snow depth (m)
  case('scalarSWE'                   ); get_ixmvar = iLookMVAR%scalarSWE                   ! snow water equivalent (kg m-2)
  case('scalarSfcMeltPond'           ); get_ixmvar = iLookMVAR%scalarSfcMeltPond           ! ponded water caused by melt of the "snow without a layer" (kg m-2)
  case('scalarAquiferStorage'        ); get_ixmvar = iLookMVAR%scalarAquiferStorage        ! relative aquifer storage -- above bottom of the soil profile (m)
  case('scalarWaterTableDepth'       ); get_ixmvar = iLookMVAR%scalarWaterTableDepth       ! depth of the water table (m)
  case('scalarExCoef'                ); get_ixmvar = iLookMVAR%scalarExCoef                ! turbulent exchange coefficient (-)
  case('scalarExSen'                 ); get_ixmvar = iLookMVAR%scalarExSen                 ! exchange factor for sensible heat (W m-2 K-1)
  case('scalarExLat'                 ); get_ixmvar = iLookMVAR%scalarExLat                 ! exchange factor for latent heat (W m-2) 
  case('scalarSenHeat'               ); get_ixmvar = iLookMVAR%scalarSenHeat               ! sensible heat flux at the surface (W m-2)
  case('scalarLatHeat'               ); get_ixmvar = iLookMVAR%scalarLatHeat               ! latent heat flux at the surface (W m-2)
  case('scalarPotentialET'           ); get_ixmvar = iLookMVAR%scalarPotentialET           ! potential evapotranspiration (-)
  case('scalarTranspireLim'          ); get_ixmvar = iLookMVAR%scalarTranspireLim          ! aggregate soil moist & vegetation limit on transpiration (-)
  case('scalarTranspireLimAqfr'      ); get_ixmvar = iLookMVAR%scalarTranspireLimAqfr      ! soil moist & vegetation limit on transpiration in the aquifer (-)
  case('scalarMassLiquid'            ); get_ixmvar = iLookMVAR%scalarMassLiquid            ! evaporation or dew (kg m-2 s-1)
  case('scalarMassSolid'             ); get_ixmvar = iLookMVAR%scalarMassSolid             ! sublimation or frost (kg m-2 s-1)
  case('scalarHeatPrecip'            ); get_ixmvar = iLookMVAR%scalarHeatPrecip            ! sensible heat of precipitation (W m-2)
  case('scalarRainPlusMelt'          ); get_ixmvar = iLookMVAR%scalarRainPlusMelt          ! rain plus melt, as input to soil before calculating surface runoff (m s-1)
  case('scalarSurfaceRunoff'         ); get_ixmvar = iLookMVAR%scalarSurfaceRunoff         ! surface runoff (m s-1)
  case('scalarInitAquiferRecharge'   ); get_ixmvar = iLookMVAR%scalarInitAquiferRecharge   ! recharge to the aquifer -- at the start of the step (m s-1)
  case('scalarAquiferRecharge'       ); get_ixmvar = iLookMVAR%scalarAquiferRecharge       ! recharge to the aquifer (m s-1)
  case('scalarInitAquiferTranspire'  ); get_ixmvar = iLookMVAR%scalarInitAquiferTranspire  ! transpiration from the aquifer (m s-1)
  case('scalarAquiferTranspire'      ); get_ixmvar = iLookMVAR%scalarAquiferTranspire      ! transpiration from the aquifer (m s-1)
  case('scalarInitAquiferBaseflow'   ); get_ixmvar = iLookMVAR%scalarInitAquiferBaseflow   ! baseflow from the aquifer (m s-1)
  case('scalarAquiferBaseflow'       ); get_ixmvar = iLookMVAR%scalarAquiferBaseflow       ! baseflow from the aquifer (m s-1)
  case('scalarSoilInflux'            ); get_ixmvar = iLookMVAR%scalarSoilInflux            ! influx of water at the top of the soil profile (m s-1)
  case('scalarSoilBaseflow'          ); get_ixmvar = iLookMVAR%scalarSoilBaseflow          ! total baseflow from throughout the soil profile (m s-1)
  case('scalarSoilDrainage'          ); get_ixmvar = iLookMVAR%scalarSoilDrainage          ! drainage from the bottom of the soil profile (m s-1)
  case('scalarSoilEjection'          ); get_ixmvar = iLookMVAR%scalarSoilEjection          ! total water ejected from soil layers (m s-1)
  case('scalarSoilWatBalError'       ); get_ixmvar = iLookMVAR%scalarSoilWatBalError       ! error in the total soil water balance (kg m-2)
  case('scalarAquiferBalError'       ); get_ixmvar = iLookMVAR%scalarAquiferBalError       ! error in the aquifer water balance (kg m-2)
  case('scalarTotalSoilLiq'          ); get_ixmvar = iLookMVAR%scalarTotalSoilLiq          ! total mass of liquid water in the soil (kg m-2)
  case('scalarTotalSoilIce'          ); get_ixmvar = iLookMVAR%scalarTotalSoilIce          ! total mass of ice in the soil (kg m-2)
  ! variables at the mid-point of each layer -- domain geometry
  case('mLayerDepth'                 ); get_ixmvar = iLookMVAR%mLayerDepth                 ! depth of each layer (m)
  case('mLayerHeight'                ); get_ixmvar = iLookMVAR%mLayerHeight                ! height at the midpoint of each layer (m)
  case('mLayerRootDensity'           ); get_ixmvar = iLookMVAR%mLayerRootDensity           ! fraction of roots in each soil layer (-)
  ! variables at the mid-point of each layer -- coupled energy and mass
  case('mLayerTemp'                  ); get_ixmvar = iLookMVAR%mLayerTemp                  ! temperature of each layer (K)
  case('mLayerVolFracAir'            ); get_ixmvar = iLookMVAR%mLayerVolFracAir            ! volumetric fraction of air in each layer (-)
  case('mLayerVolFracIce'            ); get_ixmvar = iLookMVAR%mLayerVolFracIce            ! volumetric fraction of icein each layer (-)
  case('mLayerVolFracLiq'            ); get_ixmvar = iLookMVAR%mLayerVolFracLiq            ! volumetric fraction of liquid water in each layer (-)
  case('mLayerVolHtCapBulk'          ); get_ixmvar = iLookMVAR%mLayerVolHtCapBulk          ! volumetric heat capacity in each layer (J m-3 K-1)
  case('mLayerTcrit'                 ); get_ixmvar = iLookMVAR%mLayerTcrit                 ! critical soil temperature above which all water is unfrozen (K)
  case('mLayerdTheta_dTk'            ); get_ixmvar = iLookMVAR%mLayerdTheta_dTk            ! analytical derivative in the freezing curve (K-1)
  case('mLayerThermalC'              ); get_ixmvar = iLookMVAR%mLayerThermalC              ! thermal conductivity at the mid-point of each layer (W m-1 K-1)
  case('mLayerRadCondFlux'           ); get_ixmvar = iLookMVAR%mLayerRadCondFlux           ! temporal derivative in energy from radiative and conductive flux (J m-2 s-1)
  case('mLayerMeltFreeze'            ); get_ixmvar = iLookMVAR%mLayerMeltFreeze            ! melt/freeze in each layer (kg m-3 s-1)
  case('mLayerSatHydCond'            ); get_ixmvar = iLookMVAR%mLayerSatHydCond            ! saturated hydraulic conductivity in each layer (m s-1)
  case('mLayerMatricHead'            ); get_ixmvar = iLookMVAR%mLayerMatricHead            ! matric head of water in the soil (m)
  case('mLayerdTheta_dPsi'           ); get_ixmvar = iLookMVAR%mLayerdTheta_dPsi           ! analytical derivative in the soil water characteristic w.r.t. psi (m-1)
  case('mLayerdPsi_dTheta'           ); get_ixmvar = iLookMVAR%mLayerdPsi_dTheta           ! analytical derivative in the soil water characteristic w.r.t. theta (m)
  case('mLayerThetaResid'            ); get_ixmvar = iLookMVAR%mLayerThetaResid            ! residual volumetric water content in each snow layer (-)
  case('mLayerPoreSpace'             ); get_ixmvar = iLookMVAR%mLayerPoreSpace             ! total pore space in each snow layer (-)
  case('mLayerInfilFreeze'           ); get_ixmvar = iLookMVAR%mLayerInfilFreeze           ! volumetric ice content increase by freezing infiltrating flux (-)
  case('mLayerTranspireLim'          ); get_ixmvar = iLookMVAR%mLayerTranspireLim          ! moisture avail factor limiting transpiration in each layer (-)
  case('mLayerInitTranspire'         ); get_ixmvar = iLookMVAR%mLayerInitTranspire         ! transpiration loss from each soil layer at the start of the step (kg m-2 s-1)
  case('mLayerTranspire'             ); get_ixmvar = iLookMVAR%mLayerTranspire             ! transpiration loss from each soil layer (kg m-2 s-1)
  case('mLayerInitEjectWater'        ); get_ixmvar = iLookMVAR%mLayerInitEjectWater        ! water ejected from each soil layer at the start-of-step (m s-1)
  case('mLayerEjectWater'            ); get_ixmvar = iLookMVAR%mLayerEjectWater            ! water ejected from each soil layer (m s-1)
  case('mLayerInitBaseflow'          ); get_ixmvar = iLookMVAR%mLayerInitBaseflow          ! baseflow from each soil layer at the start of the time step (m s-1)
  case('mLayerBaseflow'              ); get_ixmvar = iLookMVAR%mLayerBaseflow              ! baseflow from each soil layer (m s-1)
  ! variables at the interface of each layer
  case('iLayerHeight'                ); get_ixmvar = iLookMVAR%iLayerHeight                ! height at the interface of each layer (m)
  case('iLayerThermalC'              ); get_ixmvar = iLookMVAR%iLayerThermalC              ! thermal conductivity at the interface of each layer (W m-1 K-1)
  case('iLayerInitNrgFlux'           ); get_ixmvar = iLookMVAR%iLayerInitNrgFlux           ! energy flux at layer interfaces at the start of the time step (W m-2)
  case('iLayerNrgFlux'               ); get_ixmvar = iLookMVAR%iLayerNrgFlux               ! energy flux at layer interfaces at the end of the time step (W m-2)
  case('iLayerSatHydCond'            ); get_ixmvar = iLookMVAR%iLayerSatHydCond            ! saturated hydraulic conductivity in each layer (m s-1)
  case('iLayerInitLiqFluxSnow'       ); get_ixmvar = iLookMVAR%iLayerInitLiqFluxSnow       ! liquid flux at snow layer interfaces at the start of the time step (m s-1)
  case('iLayerInitLiqFluxSoil'       ); get_ixmvar = iLookMVAR%iLayerInitLiqFluxSoil       ! liquid flux at soil layer interfaces at the start of the time step (m s-1)
  case('iLayerLiqFluxSnow'           ); get_ixmvar = iLookMVAR%iLayerLiqFluxSnow           ! liquid flux at snow layer interfaces at the end of the time step (m s-1)
  case('iLayerLiqFluxSoil'           ); get_ixmvar = iLookMVAR%iLayerLiqFluxSoil           ! liquid flux at soil layer interfaces at the end of the time step (m s-1)
  ! variables to compute runoff
  case('routingRunoffFuture'         ); get_ixmvar = iLookMVAR%routingRunoffFuture         ! runoff in future time steps (m s-1)
  case('routingFractionFuture'       ); get_ixmvar = iLookMVAR%routingFractionFuture       ! fraction of runoff in future time steps (-)
  case('averageInstantRunoff'        ); get_ixmvar = iLookMVAR%averageInstantRunoff        ! instantaneous runoff (m s-1)
  case('averageRoutedRunoff'         ); get_ixmvar = iLookMVAR%averageRoutedRunoff         ! routed runoff (m s-1)
  ! "short-cut" variables
  case('scalarExNeut'                ); get_ixmvar = iLookMVAR%scalarExNeut                ! exchange coefficient in neutral conditions (-)
  case('scalarBprime'                ); get_ixmvar = iLookMVAR%scalarBprime                ! stable b parameter in Louis (1979) stability function (-)
  case('scalarCparam'                ); get_ixmvar = iLookMVAR%scalarCparam                ! c parameter in Louis (1979) stability function 
  case('scalarVGn_m'                 ); get_ixmvar = iLookMVAR%scalarVGn_m                 ! van Genuchten "m" parameter (-) 
  case('scalarKappa'                 ); get_ixmvar = iLookMVAR%scalarKappa                 ! constant in the freezing curve function (m K-1)
  case('scalarVolHtCap_air'          ); get_ixmvar = iLookMVAR%scalarVolHtCap_air          ! volumetric heat capacity air         (J m-3 K-1)
  case('scalarVolHtCap_ice'          ); get_ixmvar = iLookMVAR%scalarVolHtCap_ice          ! volumetric heat capacity ice         (J m-3 K-1)
  case('scalarVolHtCap_soil'         ); get_ixmvar = iLookMVAR%scalarVolHtCap_soil         ! volumetric heat capacity dry soil    (J m-3 K-1)
  case('scalarVolHtCap_water'        ); get_ixmvar = iLookMVAR%scalarVolHtCap_water        ! volumetric heat capacity liquid wat  (J m-3 K-1)
  case('scalarLambda_drysoil'        ); get_ixmvar = iLookMVAR%scalarLambda_drysoil        ! thermal conductivity of dry soil     (W m-1)
  case('scalarLambda_wetsoil'        ); get_ixmvar = iLookMVAR%scalarLambda_wetsoil        ! thermal conductivity of wet soil     (W m-1)
  case('scalarVolLatHt_fus'          ); get_ixmvar = iLookMVAR%scalarVolLatHt_fus          ! volumetric latent heat of fusion     (J m-3)
  case('scalarAquiferRootFrac'       ); get_ixmvar = iLookMVAR%scalarAquiferRootFrac       ! fraction of roots below the soil profile (-)
  ! get to here if cannot find the variable
  case default
   get_ixmvar = imiss
 endselect
 end function get_ixmvar


 ! *******************************************************************************************************************
 ! new function: get the index of the named variables for the model indices
 ! *******************************************************************************************************************
 function get_ixindex(varName)
 USE var_lookup,only:iLookINDEX                      ! indices of the named variables
 implicit none
 ! define dummy variables
 character(*), intent(in) :: varName                 ! parameter name
 integer(i4b)             :: get_ixINDEX             ! index of the named variable
 ! define local variables
 integer(i4b), parameter  :: imiss = -999            ! missing value
 ! get the index of the named variables
 select case(trim(varName))
  case('nSnow'            ); get_ixindex = iLookINDEX%nSnow             ! number of snow layers
  case('nSoil'            ); get_ixindex = iLookINDEX%nSoil             ! number of soil layers
  case('nLayers'          ); get_ixindex = iLookINDEX%nLayers           ! total number of layers
  case('midSnowStartIndex'); get_ixindex = iLookINDEX%midSnowStartIndex ! start index of the midSnow vector for a given timestep
  case('midSoilStartIndex'); get_ixindex = iLookINDEX%midSoilStartIndex ! start index of the midSoil vector for a given timestep
  case('midTotoStartIndex'); get_ixindex = iLookINDEX%midTotoStartIndex ! start index of the midToto vector for a given timestep
  case('ifcSnowStartIndex'); get_ixindex = iLookINDEX%ifcSnowStartIndex ! start index of the ifcSnow vector for a given timestep
  case('ifcSoilStartIndex'); get_ixindex = iLookINDEX%ifcSoilStartIndex ! start index of the ifcSoil vector for a given timestep
  case('ifcTotoStartIndex'); get_ixindex = iLookINDEX%ifcTotoStartIndex ! start index of the ifcToto vector for a given timestep
  case('layerType'        ); get_ixindex = iLookINDEX%layerType         ! type of layer (soil or snow)
  ! get to here if cannot find the variable
  case default
   get_ixindex = imiss
 endselect
 end function get_ixindex

end module get_ixname_module

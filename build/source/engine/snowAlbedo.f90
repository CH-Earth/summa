module snowAlbedo_module
USE nrtype                        ! numerical recipes data types
USE mDecisions_module,only:  &    ! identify model options for snow albedo
 constantDecay,              &    ! constant decay in snow albedo (e.g., VIC, CLASS)
 variableDecay                    ! variable decay in snow albedo (e.g., BATS approach, with destructive metamorphism + soot content)
! -------------------------------------------------------------------------------------------------
implicit none
private
public::snowAlbedo
! dimensions
integer(i4b),parameter        :: nBands=2      ! number of spectral bands for shortwave radiation
contains

 ! ************************************************************************************************
 ! new subroutine: muster program to compute energy fluxes at vegetation and ground surfaces
 ! ************************************************************************************************
 subroutine snowAlbedo(&
                       ! input: model control
                       dt,                                    & ! intent(in): model time step (s)
                       snowPresence,                          & ! intent(in): logical flag to denote if snow is present
                       ixAlbedoMethod,                        & ! intent(in): index of the method used for snow albedo
                       ! input: model variables
                       snowfallRate,                          & ! intent(in): snowfall rate (kg m-2 s-1)
                       surfaceTemp,                           & ! intent(in): surface temperature (K)
                       cosZenith,                             & ! intent(in): cosine of the zenith angle
                       ! input: model parameters
                       Frad_vis,                              & ! intent(in): fraction of radiation in visible part of spectrum (-)
                       Frad_direct,                           & ! intent(in): fraction direct solar radiation (-)
                       albedoMax,                             & ! intent(in): maximum snow albedo for a single spectral band (-)
                       albedoMinWinter,                       & ! intent(in): minimum snow albedo during winter for a single spectral band (-)
                       albedoMinSpring,                       & ! intent(in): minimum snow albedo during spring for a single spectral band (-)
                       albedoMaxVisible,                      & ! intent(in): maximum snow albedo in the visible part of the spectrum (-)
                       albedoMinVisible,                      & ! intent(in): minimum snow albedo in the visible part of the spectrum (-)
                       albedoMaxNearIR,                       & ! intent(in): maximum snow albedo in the near infra-red part of the spectrum (-)
                       albedoMinNearIR,                       & ! intent(in): minimum snow albedo in the near infra-red part of the spectrum (-)
                       albedoDecayRate,                       & ! intent(in): albedo decay rate (s)
                       tempScalGrowth,                        & ! intent(in): temperature scaling factor for grain growth (K-1) 
                       albedoSootLoad,                        & ! intent(in): soot load factor (-)
                       albedoRefresh,                         & ! intent(in): critical mass necessary for albedo refreshment (kg m-2)
                       snowfrz_scale,                         & ! intent(in): scaling parameter for the freezing curve for snow (K-1)
                       ! input-output: snow albedo
                       spectralSnowAlbedoDiffuse,             & ! intent(inout): diffuse snow albedo in each spectral band (-)
                       spectralSnowAlbedoDirect,              & ! intent(inout): direct snow albedo in each spectral band (-)
                       scalarSnowAlbedo,                      & ! intent(inout): snow albedo for the entire spectral band (-)
                       ! output: error control
                       err,message)                             ! intent(out): error control
 ! --------------------------------------------------------------------------------------------------------------------------------------
 USE multiconst,only:Tfreeze                                    ! freezing point of pure water (K)
 USE snow_utils_module,only:fracliquid                          ! compute fraction of liquid water at a given temperature
 ! input: model control
 real(dp),intent(in)            :: dt                           ! model time step
 logical(lgt),intent(in)        :: snowPresence                 ! logical flag to denote if snow is present
 integer(i4b),intent(in)        :: ixAlbedoMethod               ! index of method used for snow albedo
 ! input: model variables
 real(dp),intent(in)            :: snowfallRate                 ! snowfall rate (kg m-2 s-1)
 real(dp),intent(in)            :: surfaceTemp                  ! surface temperature (K)
 real(dp),intent(in)            :: cosZenith                    ! cosine of the zenith angle 
 ! input: model parameters
 real(dp),intent(in)            :: Frad_vis                     ! fraction of radiation in visible part of spectrum (-)
 real(dp),intent(in)            :: Frad_direct                  ! fraction direct solar radiation (-)
 real(dp),intent(in)            :: albedoMax                    ! maximum snow albedo for a single spectral band (-)
 real(dp),intent(in)            :: albedoMinWinter              ! minimum snow albedo during winter for a single spectral band (-)
 real(dp),intent(in)            :: albedoMinSpring              ! minimum snow albedo during spring for a single spectral band (-)
 real(dp),intent(in)            :: albedoMaxVisible             ! maximum snow albedo in the visible part of the spectrum (-)
 real(dp),intent(in)            :: albedoMinVisible             ! minimum snow albedo in the visible part of the spectrum (-)
 real(dp),intent(in)            :: albedoMaxNearIR              ! maximum snow albedo in the near infra-red part of the spectrum (-)
 real(dp),intent(in)            :: albedoMinNearIR              ! minimum snow albedo in the near infra-red part of the spectrum (-)
 real(dp),intent(in)            :: albedoDecayRate              ! albedo decay rate (s)
 real(dp),intent(in)            :: tempScalGrowth               ! temperature scaling factor for grain growth (K-1) 
 real(dp),intent(in)            :: albedoSootLoad               ! soot load factor (-)
 real(dp),intent(in)            :: albedoRefresh                ! critical mass necessary for albedo refreshment (kg m-2)
 real(dp),intent(in)            :: snowfrz_scale                ! scaling parameter for the freezing curve for snow (K-1)
 ! input-output: snow albedo
 real(dp),intent(inout)         :: spectralSnowAlbedoDiffuse(:) ! diffuse snow albedo in each spectral band (-)
 real(dp),intent(inout)         :: spectralSnowAlbedoDirect(:)  ! direct snow albedo in each spectral band (-)
 real(dp),intent(inout)         :: scalarSnowAlbedo             ! snow albedo for the entire spectral band (-)
 ! output: error control
 integer(i4b),intent(out)       :: err                          ! error code
 character(*),intent(out)       :: message                      ! error message
 ! local variables
 integer(i4b),parameter         :: ixVisible=1                  ! named variable to define index in array of visible part of the spectrum
 integer(i4b),parameter         :: ixNearIR=2                   ! named variable to define index in array of near IR part of the spectrum
 real(dp),parameter             :: valueMissing=-9999._dp       ! missing value -- will cause problems if snow albedo is ever used for the non-snow case
 real(dp),parameter             :: slushExp=10._dp              ! "slush" exponent, to increase decay when snow is near Tfreeze
 real(dp),parameter             :: fractionLiqThresh=0.001_dp   ! threshold for the fraction of liquid water to switch to spring albedo minimum
 real(dp)                       :: fractionLiq                  ! fraction of liquid water (-)
 real(dp)                       :: age1,age2,age3               ! aging factors (-)
 real(dp)                       :: decayFactor                  ! albedo decay factor (-)
 real(dp)                       :: refreshFactor                ! albedo refreshment factor, representing albedo increase due to snowfall (-)
 real(dp)                       :: albedoMin                    ! minimum albedo -- depends if in winter or spring conditions (-)
 real(dp)                       :: fZen                         ! factor to modify albedo at low zenith angles (-)
 real(dp),parameter             :: bPar=2._dp                   ! empirical parameter in fZen
 ! initialize error control
 err=0; message='snowAlbedo/'

 ! return early if no snow
 if(.not. snowPresence)then
  scalarSnowAlbedo             = valueMissing
  spectralSnowAlbedoDirect(:)  = valueMissing
  spectralSnowAlbedoDiffuse(:) = valueMissing 
  return
 endif

 ! compute fractional increase in albedo associated with snowfall
 refreshFactor = dt*snowfallRate/albedoRefresh

 ! identify option for snow albedo
 select case(ixAlbedoMethod) 


  ! *** constant decay rate
  case(constantDecay)
   ! compute decay rate
   decayFactor = dt/albedoDecayRate
   ! compute minimum albedo
   fractionLiq = fracliquid(surfaceTemp,snowfrz_scale) ! fraction of liquid water
   if(scalarSnowAlbedo < albedoMinWinter .or. fractionLiq > fractionLiqThresh)then
    albedoMin = albedoMinSpring
   else
    albedoMin = albedoMinWinter
   endif
   ! compute average albedo
   call computeAlbedo(scalarSnowAlbedo,refreshFactor,decayFactor,albedoMax,albedoMin)
   ! assume albedo is the same in visible and near infra-red bands, and for direct and diffuse radiation
   spectralSnowAlbedoDiffuse(ixVisible) = scalarSnowAlbedo
   spectralSnowAlbedoDiffuse(ixNearIR)  = scalarSnowAlbedo
   spectralSnowAlbedoDirect(ixVisible)  = scalarSnowAlbedo
   spectralSnowAlbedoDirect(ixNearIR)   = scalarSnowAlbedo
   

  ! *** variable decay rate
  case(variableDecay)
   ! compute decay factor 
   age1 = exp(-tempScalGrowth*(Tfreeze - surfaceTemp ))  ! temperature dependence
   age2 = age1**slushExp                                 ! increase with liquid water
   age3 = albedoSootLoad                                 ! soot loading
   decayFactor = dt*(age1 + age2 + age3)/albedoDecayRate
   ! compute diffuse albedo for the different spectral bands
   call computeAlbedo(spectralSnowAlbedoDiffuse(ixVisible),refreshFactor,decayFactor,albedoMaxVisible,albedoMinVisible)
   call computeAlbedo(spectralSnowAlbedoDiffuse(ixNearIR), refreshFactor,decayFactor,albedoMaxNearIR, albedoMinNearIR)
   ! compute factor to modify direct albedo at low zenith angles
   if(cosZenith < 0.5_dp)then
    fZen = (1._dp/bPar)*( ((1._dp + bPar)/(1._dp + 2._dp*bPar*cosZenith)) - 1._dp)
   else
    fZen = 0._dp
   endif
   ! compute direct albedo
   spectralSnowAlbedoDirect(ixVisible) = spectralSnowAlbedoDiffuse(ixVisible) + 0.4_dp*fZen*(1._dp - spectralSnowAlbedoDiffuse(ixVisible))
   spectralSnowAlbedoDirect(ixNearIR)  = spectralSnowAlbedoDiffuse(ixNearIR)  + 0.4_dp*fZen*(1._dp - spectralSnowAlbedoDiffuse(ixNearIR))

   ! compute average albedo
   scalarSnowAlbedo = (        Frad_direct)*(Frad_vis*spectralSnowAlbedoDirect(ixVisible) + (1._dp - Frad_vis)*spectralSnowAlbedoDirect(ixNearIR) ) + &
                      (1._dp - Frad_direct)*(Frad_vis*spectralSnowAlbedoDirect(ixVisible) + (1._dp - Frad_vis)*spectralSnowAlbedoDirect(ixNearIR) )

  ! check that we identified the albedo option
  case default; err=20; message=trim(message)//'unable to identify option for snow albedo'; return

 end select  ! identify option for snow albedo

 ! check
 if(scalarSnowAlbedo < 0._dp)then; err=20; message=trim(message)//'unable to identify option for snow albedo'; return; endif

 end subroutine snowAlbedo

 ! ** private function
 ! compute change in albedo -- implicit solution
 subroutine computeAlbedo(snowAlbedo,refreshFactor,decayFactor,albedoMax,albedoMin)
 implicit none
 ! dummy variables
 real(dp),intent(inout)   :: snowAlbedo    ! snow albedo (-)
 real(dp),intent(in)      :: refreshFactor ! albedo refreshment factor (-)
 real(dp),intent(in)      :: decayFactor   ! albedo decay factor (-)
 real(dp),intent(in)      :: albedoMax     ! maximum albedo (-)
 real(dp),intent(in)      :: albedoMin     ! minimum albedo (-)
 ! local variables
 real(dp)                 :: albedoChange ! change in albedo over the time step (-)
 ! compute change in albedo
 albedoChange = refreshFactor*(albedoMax - snowAlbedo) - (decayFactor*(snowAlbedo - albedoMin)) / (1._dp + decayFactor) 
 snowAlbedo   = snowAlbedo + albedoChange
 if(snowAlbedo > albedoMax) snowAlbedo = albedoMax
 end subroutine computeAlbedo


end module snowAlbedo_module

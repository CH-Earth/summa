module vegSWavRad_module
USE nrtype
! -------------------------------------------------------------------------------------------------
implicit none
private
public::vegSWavRad
! dimensions
integer(i4b),parameter        :: nBands=2      ! number of spectral bands for shortwave radiation
! named variables
integer(i4b),parameter        :: ist     = 1   ! Surface type:  IST=1 => soil;  IST=2 => lake
! spatial indices
integer(i4b),parameter        :: iLoc    = 1   ! i-location
integer(i4b),parameter        :: jLoc    = 1   ! j-location
! algorithmic parameters
real(dp),parameter     :: missingValue=-9999._dp  ! missing value, used when diagnostic or state variables are undefined
real(dp),parameter     :: verySmall=1.e-6_dp   ! used as an additive constant to check if substantial difference among real numbers
real(dp),parameter     :: mpe=1.e-6_dp         ! prevents overflow error if division by zero 
real(dp),parameter     :: dx=1.e-6_dp          ! finite difference increment
! control
logical(lgt)           :: printflag            ! flag to turn on printing
contains

 ! ************************************************************************************************
 ! new subroutine: muster program to compute energy fluxes at vegetation and ground surfaces
 ! ************************************************************************************************
 subroutine vegSWavRad(&
                       ! input: model control
                       vegTypeIndex,                                       & ! intent(in): index of vegetation type
                       soilTypeIndex,                                      & ! intent(in): index of soil type
                       computeVegFlux,                                     & ! intent(in): logical flag to compute vegetation fluxes (.false. if veg buried by snow)
                       ! input: model variables
                       scalarCosZenith,                                    & ! intent(in): cosine of direct zenith angle (0-1)
                       spectralIncomingDirect,                             & ! intent(in): incoming direct solar radiation in each wave band (w m-2)
                       spectralIncomingDiffuse,                            & ! intent(in): incoming diffuse solar radiation in each wave band (w m-2)
                       spectralSnowAlbedoDirect,                           & ! intent(in): direct albedo of snow in each spectral band (-)
                       spectralSnowAlbedoDiffuse,                          & ! intent(in): diffuse albedo of snow in each spectral band (-)
                       scalarExposedLAI,                                   & ! intent(in): exposed leaf area index after burial by snow (m2 m-2)
                       scalarExposedSAI,                                   & ! intent(in): exposed stem area index after burial by snow (m2 m-2)
                       scalarVegFraction,                                  & ! intent(in): vegetation fraction (=1 forces no canopy gaps and open areas in radiation routine)
                       scalarCanopyWetFraction,                            & ! intent(in): fraction of lai, sai that is wetted (-)
                       scalarGroundSnowFraction,                           & ! intent(in): fraction of ground that is snow covered (-)
                       scalarVolFracLiqUpper,                              & ! intent(in): volumetric liquid water content in the upper-most soil layer (-)
                       scalarCanopyTempTrial,                              & ! intent(in): canopy temperature (K)
                       ! output
                       scalarBelowCanopySolar,                             & ! intent(out): radiation transmitted below the canopy (W m-2)
                       scalarCanopyAbsorbedSolar,                          & ! intent(out): radiation absorbed by the vegetation canopy (W m-2)
                       scalarGroundAbsorbedSolar,                          & ! intent(out): radiation absorbed by the ground (W m-2)
                       scalarCanopySunlitFraction,                         & ! intent(out): sunlit fraction of canopy (-)
                       scalarCanopySunlitLAI,                              & ! intent(out): sunlit leaf area (-)
                       scalarCanopyShadedLAI,                              & ! intent(out): shaded leaf area (-)
                       scalarCanopySunlitPAR,                              & ! intent(out): average absorbed par for sunlit leaves (w m-2)
                       scalarCanopyShadedPAR,                              & ! intent(out): average absorbed par for shaded leaves (w m-2)
                       err,message)                                          ! intent(out): error control
 ! Noah-MP modules
 USE NOAHMP_ROUTINES,only:twoStream                                          ! two-stream radiative transfer
 ! Noah vegetation tables
 USE NOAHMP_VEG_PARAMETERS, only: RHOS,RHOL                                  ! Noah-MP: stem and leaf reflectance for each wave band
 USE NOAHMP_VEG_PARAMETERS, only: TAUS,TAUL                                  ! Noah-MP: stem and leaf transmittance for each wave band
 ! input
 integer(i4b),intent(in)        :: vegTypeIndex                              ! vegetation type index
 integer(i4b),intent(in)        :: soilTypeIndex                             ! soil type index
 logical(lgt),intent(in)        :: computeVegFlux                            ! logical flag to compute vegetation fluxes (.false. if veg buried by snow)
 real(dp),intent(in)            :: scalarCosZenith                           ! cosine of the solar zenith angle (0-1)
 real(dp),intent(in)            :: spectralIncomingDirect(:)                 ! incoming direct solar radiation in each wave band (w m-2)
 real(dp),intent(in)            :: spectralIncomingDiffuse(:)                ! incoming diffuse solar radiation in each wave band (w m-2)
 real(dp),intent(in)            :: spectralSnowAlbedoDirect(:)               ! direct albedo of snow in each spectral band (-)
 real(dp),intent(in)            :: spectralSnowAlbedoDiffuse(:)              ! diffuse albedo of snow in each spectral band (-)
 real(dp),intent(in)            :: scalarExposedLAI                          ! exposed leaf area index after burial by snow (m2 m-2)
 real(dp),intent(in)            :: scalarExposedSAI                          ! exposed stem area index after burial by snow (m2 m-2)
 real(dp),intent(in)            :: scalarVegFraction                         ! vegetation fraction (=1 forces no canopy gaps and open areas in radiation routine)
 real(dp),intent(in)            :: scalarCanopyWetFraction                   ! fraction of canopy that is wet (-)
 real(dp),intent(in)            :: scalarGroundSnowFraction                  ! fraction of ground that is snow covered (-)
 real(dp),intent(in)            :: scalarVolFracLiqUpper                     ! volumetric liquid water content in the upper-most soil layer (-)
 real(dp),intent(in)            :: scalarCanopyTempTrial                     ! trial value of canopy temperature (K)
 ! output
 real(dp),intent(out)           :: scalarBelowCanopySolar                    ! radiation transmitted below the canopy (W m-2)
 real(dp),intent(out)           :: scalarCanopyAbsorbedSolar                 ! radiation absorbed by the vegetation canopy (W m-2)
 real(dp),intent(out)           :: scalarGroundAbsorbedSolar                 ! radiation absorbed by the ground (W m-2)
 real(dp),intent(out)           :: scalarCanopySunlitFraction                ! sunlit fraction of canopy (-)
 real(dp),intent(out)           :: scalarCanopySunlitLAI                     ! sunlit leaf area (-)
 real(dp),intent(out)           :: scalarCanopyShadedLAI                     ! shaded leaf area (-)
 real(dp),intent(out)           :: scalarCanopySunlitPAR                     ! average absorbed par for sunlit leaves (w m-2)
 real(dp),intent(out)           :: scalarCanopyShadedPAR                     ! average absorbed par for shaded leaves (w m-2)
 integer(i4b),intent(out)       :: err                                       ! error code
 character(*),intent(out)       :: message                                   ! error message
 ! -----------------------------------------------------------------------------------------------------------------------------------------------------------------
 ! local variables
 ! -----------------------------------------------------------------------------------------------------------------------------------------------------------------
 ! general
 integer(i4b),parameter                :: ixVisible=1                        ! index of the visible wave band
 integer(i4b)                          :: iBand                              ! index of wave band
 integer(i4b)                          :: ic                                 ! 0=unit incoming direct; 1=unit incoming diffuse
 ! vegetation properties
 real(dp)                              :: scalarExposedVAI                   ! one-sided leaf+stem area index (m2/m2)
 real(dp)                              :: weightLeaf                         ! fraction of exposed VAI that is leaf
 real(dp)                              :: weightStem                         ! fraction of exposed VAI that is stem
 real(dp),dimension(1:nBands)          :: spectralVegReflc                   ! leaf+stem reflectance (1:nbands)
 real(dp),dimension(1:nBands)          :: spectralVegTrans                   ! leaf+stem transmittance (1:nBands)
 ! albedo
 real(dp),dimension(1:nBands)          :: spectralAlbGndDirect               ! direct  albedo of underlying surface (1:nBands) (-)
 real(dp),dimension(1:nBands)          :: spectralAlbGndDiffuse              ! diffuse albedo of underlying surface (1:nBands) (-)
 ! output from two-stream -- direct-beam
 real(dp),dimension(1:nBands)          :: spectralCanopyAbsorbedDirect       ! flux abs by veg layer (per unit incoming flux), (1:nBands)
 real(dp),dimension(1:nBands)          :: spectralTotalReflectedDirect       ! flux refl above veg layer (per unit incoming flux), (1:nBands)
 real(dp),dimension(1:nBands)          :: spectralDirectBelowCanopyDirect    ! down dir flux below veg layer (per unit in flux), (1:nBands)
 real(dp),dimension(1:nBands)          :: spectralDiffuseBelowCanopyDirect   ! down dif flux below veg layer (per unit in flux), (1:nBands)
 real(dp),dimension(1:nBands)          :: spectralCanopyReflectedDirect      ! flux reflected by veg layer   (per unit incoming flux), (1:nBands)
 real(dp),dimension(1:nBands)          :: spectralGroundReflectedDirect      ! flux reflected by ground (per unit incoming flux), (1:nBands)
 ! output from two-stream -- diffuse
 real(dp),dimension(1:nBands)          :: spectralCanopyAbsorbedDiffuse      ! flux abs by veg layer (per unit incoming flux), (1:nBands)
 real(dp),dimension(1:nBands)          :: spectralTotalReflectedDiffuse      ! flux refl above veg layer (per unit incoming flux), (1:nBands)
 real(dp),dimension(1:nBands)          :: spectralDirectBelowCanopyDiffuse   ! down dir flux below veg layer (per unit in flux), (1:nBands)
 real(dp),dimension(1:nBands)          :: spectralDiffuseBelowCanopyDiffuse  ! down dif flux below veg layer (per unit in flux), (1:nBands)
 real(dp),dimension(1:nBands)          :: spectralCanopyReflectedDiffuse     ! flux reflected by veg layer   (per unit incoming flux), (1:nBands)
 real(dp),dimension(1:nBands)          :: spectralGroundReflectedDiffuse     ! flux reflected by ground (per unit incoming flux), (1:nBands)
 ! output from two-stream -- scalar variables
 real(dp)                              :: scalarGproj                        ! projected leaf+stem area in solar direction
 real(dp)                              :: scalarBetweenCanopyGapFraction     ! between canopy gap fraction for beam (-)
 real(dp)                              :: scalarWithinCanopyGapFraction      ! within canopy gap fraction for beam (-)
 ! radiation fluxes
 real(dp)                              :: ext                                ! optical depth of direct beam per unit leaf + stem area
 real(dp),dimension(1:nBands)          :: spectralBelowCanopyDirect          ! downward direct flux below veg layer (W m-2)
 real(dp),dimension(1:nBands)          :: spectralBelowCanopyDiffuse         ! downward diffuse flux below veg layer (W m-2)
 real(dp)                              :: scalarCanopyShadedFraction         ! shaded fraction of the canopy
 real(dp)                              :: fractionLAI                        ! fraction of vegetation that is leaves
 real(dp)                              :: visibleAbsDirect                   ! direct-beam radiation absorbed in the visible part of the spectrum (W m-2)
 real(dp)                              :: visibleAbsDiffuse                  ! diffuse radiation absorbed in the visible part of the spectrum (W m-2)
 ! -----------------------------------------------------------------------------------------------------------------------------------------------------------------
 ! initialize error control
 err=0; message='vegSWavRad/'

 ! compute the albedo of the ground surface
 call gndAlbedo(&
                ! input
                soilTypeIndex,                         & ! intent(in): index of soil type
                scalarGroundSnowFraction,              & ! intent(in): fraction of ground that is snow covered (-)
                scalarVolFracLiqUpper,                 & ! intent(in): volumetric liquid water content in upper-most soil layer (-)
                spectralSnowAlbedoDirect,              & ! intent(in): direct albedo of snow in each spectral band (-)
                spectralSnowAlbedoDiffuse,             & ! intent(in): diffuse albedo of snow in each spectral band (-)
                ! output
                spectralAlbGndDirect,                  & ! intent(out): direct  albedo of underlying surface (-)
                spectralAlbGndDiffuse,                 & ! intent(out): diffuse albedo of underlying surface (-)
                err,message)                             ! intent(out): error control

 ! initialize accumulated fluxes
 scalarBelowCanopySolar    = 0._dp  ! radiation transmitted below the canopy (W m-2)
 scalarCanopyAbsorbedSolar = 0._dp  ! radiation absorbed by the vegetation canopy (W m-2)
 scalarGroundAbsorbedSolar = 0._dp  ! radiation absorbed by the ground (W m-2)

 ! check for an early return (no radiation or no exposed canopy)
 if(.not.computeVegFlux .or. scalarCosZenith < tiny(scalarCosZenith))then
  ! set canopy radiation to zero
  scalarCanopySunlitFraction = 0._dp                ! sunlit fraction of canopy (-)
  scalarCanopySunlitLAI      = 0._dp                ! sunlit leaf area (-)
  scalarCanopyShadedLAI      = scalarExposedLAI     ! shaded leaf area (-)
  scalarCanopySunlitPAR      = 0._dp                ! average absorbed par for sunlit leaves (w m-2)
  scalarCanopyShadedPAR      = 0._dp                ! average absorbed par for shaded leaves (w m-2)
  ! compute below-canopy radiation
  do iBand=1,nBands
   ! (set below-canopy radiation to incoming radiation)
   if(scalarCosZenith > tiny(scalarCosZenith))then
    spectralBelowCanopyDirect(iBand)  = spectralIncomingDirect(iBand)
    spectralBelowCanopyDiffuse(iBand) = spectralIncomingDiffuse(iBand)
   else
    spectralBelowCanopyDirect(iBand)  = 0._dp 
    spectralBelowCanopyDiffuse(iBand) = 0._dp
   endif
   ! (accumulate radiation transmitted below the canopy)
   scalarBelowCanopySolar    = scalarBelowCanopySolar + &                                                  ! contribution from all previous wave bands
                               spectralBelowCanopyDirect(iBand) + spectralBelowCanopyDiffuse(iBand)        ! contribution from current wave band
   ! (accumulate radiation absorbed by the ground)
   scalarGroundAbsorbedSolar = scalarGroundAbsorbedSolar + &                                               ! contribution from all previous wave bands
                               spectralBelowCanopyDirect(iBand)*(1._dp - spectralAlbGndDirect(iBand)) + &  ! direct radiation from current wave band
                               spectralBelowCanopyDiffuse(iBand)*(1._dp - spectralAlbGndDiffuse(iBand))    ! diffuse radiation from current wave band
  end do  ! looping through wave bands
  return
 endif

 ! compute exposed leaf and stem area index
 scalarExposedVAI = scalarExposedLAI + scalarExposedSAI
 if(scalarExposedVAI < epsilon(scalarExposedVAI))then; err=20; message=trim(message)//'very small exposed vegetation area (covered with snow?)'; return; endif

 ! weight reflectance and transmittance by exposed leaf and stem area index
 weightLeaf       = scalarExposedLAI / scalarExposedVAI
 weightStem       = scalarExposedSAI / scalarExposedVAI
 do iBand = 1,nBands  ! loop through spectral bands
  spectralVegReflc(iBand) = RHOL(vegTypeIndex,iBand)*weightLeaf + RHOS(vegTypeIndex,iBand)*weightStem
  spectralVegTrans(iBand) = TAUL(vegTypeIndex,iBand)*weightLeaf + TAUS(vegTypeIndex,iBand)*weightStem
 end do  

 ! loop through wave bands
 do iBand=1,nBands

  ic = 0
  ! two-stream approximation for direct-beam radiation (from CLM/Noah-MP)
  call twoStream(&
                 ! input
                 iBand,                             & ! intent(in): waveband number
                 ic,                                & ! intent(in): 0=unit incoming direct; 1=unit incoming diffuse
                 vegTypeIndex,                      & ! intent(in): vegetation type
                 scalarCosZenith,                   & ! intent(in): cosine of direct zenith angle (0-1)
                 scalarExposedVAI,                  & ! intent(in): one-sided leaf+stem area index (m2/m2)
                 scalarCanopyWetFraction,           & ! intent(in): fraction of lai, sai that is wetted (-)
                 scalarCanopyTempTrial,             & ! intent(in): surface temperature (k)
                 spectralAlbGndDirect,              & ! intent(in): direct  albedo of underlying surface (1:nBands) (-)
                 spectralAlbGndDiffuse,             & ! intent(in): diffuse albedo of underlying surface (1:nBands) (-)
                 spectralVegReflc,                  & ! intent(in): leaf+stem reflectance (1:nbands)
                 spectralVegTrans,                  & ! intent(in): leaf+stem transmittance (1:nBands)
                 scalarVegFraction,                 & ! intent(in): vegetation fraction (=1 forces no canopy gaps and open areas in radiation routine)
                 ist,                               & ! intent(in): surface type
                 iLoc,jLoc,                         & ! intent(in): grid indices
                 ! output
                 spectralCanopyAbsorbedDirect,      & ! intent(out): flux abs by veg layer (per unit incoming flux), (1:nBands)
                 spectralTotalReflectedDirect,      & ! intent(out): flux refl above veg layer (per unit incoming flux), (1:nBands)
                 spectralDirectBelowCanopyDirect,   & ! intent(out): down dir flux below veg layer (per unit in flux), (1:nBands)
                 spectralDiffuseBelowCanopyDirect,  & ! intent(out): down dif flux below veg layer (per unit in flux), (1:nBands)
                 scalarGproj,                       & ! intent(out): projected leaf+stem area in solar direction
                 spectralCanopyReflectedDirect,     & ! intent(out): flux reflected by veg layer   (per unit incoming flux), (1:nBands)
                 spectralGroundReflectedDirect,     & ! intent(out): flux reflected by ground (per unit incoming flux), (1:nBands)
                 ! input-output
                 scalarBetweenCanopyGapFraction,    & ! intent(inout): between canopy gap fraction for beam (-)
                 scalarWithinCanopyGapFraction      ) ! intent(inout): within canopy gap fraction for beam (-)

  ic = 1
  ! two-stream approximation for diffuse radiation (from CLM/Noah-MP)
  call twoStream(&
                 ! input
                 iBand,                             & ! intent(in): waveband number
                 ic,                                & ! intent(in): 0=unit incoming direct; 1=unit incoming diffuse
                 vegTypeIndex,                      & ! intent(in): vegetation type
                 scalarCosZenith,                   & ! intent(in): cosine of direct zenith angle (0-1)
                 scalarExposedVAI,                  & ! intent(in): one-sided leaf+stem area index (m2/m2)
                 scalarCanopyWetFraction,           & ! intent(in): fraction of lai, sai that is wetted (-)
                 scalarCanopyTempTrial,             & ! intent(in): surface temperature (k)
                 spectralAlbGndDirect,              & ! intent(in): direct  albedo of underlying surface (1:nBands) (-)
                 spectralAlbGndDiffuse,             & ! intent(in): diffuse albedo of underlying surface (1:nBands) (-)
                 spectralVegReflc,                  & ! intent(in): leaf+stem reflectance (1:nbands)
                 spectralVegTrans,                  & ! intent(in): leaf+stem transmittance (1:nBands)
                 scalarVegFraction,                 & ! intent(in): vegetation fraction (=1 forces no canopy gaps and open areas in radiation routine)
                 ist,                               & ! intent(in): surface type
                 iLoc,jLoc,                         & ! intent(in): grid indices
                 ! output
                 spectralCanopyAbsorbedDiffuse,     & ! intent(out): flux abs by veg layer (per unit incoming flux), (1:nBands)
                 spectralTotalReflectedDiffuse,     & ! intent(out): flux refl above veg layer (per unit incoming flux), (1:nBands)
                 spectralDirectBelowCanopyDiffuse,  & ! intent(out): down dir flux below veg layer (per unit in flux), (1:nBands)
                 spectralDiffuseBelowCanopyDiffuse, & ! intent(out): down dif flux below veg layer (per unit in flux), (1:nBands)
                 scalarGproj,                       & ! intent(out): projected leaf+stem area in solar direction
                 spectralCanopyReflectedDiffuse,    & ! intent(out): flux reflected by veg layer   (per unit incoming flux), (1:nBands)
                 spectralGroundReflectedDiffuse,    & ! intent(out): flux reflected by ground (per unit incoming flux), (1:nBands)
                 ! input-output
                 scalarBetweenCanopyGapFraction,    & ! intent(inout): between canopy gap fraction for beam (-)
                 scalarWithinCanopyGapFraction      ) ! intent(inout): within canopy gap fraction for beam (-)

  ! compute below-canopy radiation
  spectralBelowCanopyDirect(iBand)  = spectralIncomingDirect(iBand)*spectralDirectBelowCanopyDirect(iBand)      ! direct radiation
  spectralBelowCanopyDiffuse(iBand) = spectralIncomingDirect(iBand)*spectralDiffuseBelowCanopyDirect(iBand) + & ! direct radiation transmitted as diffuse
                                      spectralIncomingDiffuse(iBand)*spectralDiffuseBelowCanopyDiffuse(iBand)   ! diffuse radiation transmitted as diffuse

  ! accumulate radiation transmitted below the canopy (W m-2)
  scalarBelowCanopySolar    = scalarBelowCanopySolar + &                                                  ! contribution from all previous wave bands
                              spectralBelowCanopyDirect(iBand) + spectralBelowCanopyDiffuse(iBand)        ! contribution from current wave band

  ! accumulate radiation absorbed by the vegetation canopy (W m-2)
  scalarCanopyAbsorbedSolar = scalarCanopyAbsorbedSolar + &                                               ! contribution from all previous wave bands
                              spectralIncomingDirect(iBand)*spectralCanopyAbsorbedDirect(iBand) + &       ! direct radiation from current wave band
                              spectralIncomingDiffuse(iBand)*spectralCanopyAbsorbedDiffuse(iBand)         ! diffuse radiation from current wave band

  ! accumulate radiation absorbed by the ground (W m-2)
  scalarGroundAbsorbedSolar = scalarGroundAbsorbedSolar + &                                               ! contribution from all previous wave bands
                              spectralBelowCanopyDirect(iBand)*(1._dp - spectralAlbGndDirect(iBand)) + &  ! direct radiation from current wave band
                              spectralBelowCanopyDiffuse(iBand)*(1._dp - spectralAlbGndDiffuse(iBand))    ! diffuse radiation from current wave band

 end do  ! (looping through wave bands)

 ! compute sunlit fraction of canopy (from CLM/Noah-MP)
 ext = scalarGproj/scalarCosZenith  ! optical depth of direct beam per unit leaf + stem area
 scalarCanopySunlitFraction = (1._dp - exp(-ext*scalarExposedVAI)) / max(ext*scalarExposedVAI,mpe)
 if(scalarCanopySunlitFraction < 0.01_dp) scalarCanopySunlitFraction = 0._dp

 ! compute sunlit and shaded LAI
 scalarCanopyShadedFraction = 1._dp - scalarCanopySunlitFraction
 scalarCanopySunlitLAI      = scalarExposedLAI*scalarCanopySunlitFraction
 scalarCanopyShadedLAI      = scalarExposedLAI*scalarCanopyShadedFraction

 ! compute PAR for sunlit and shaded leaves (from CLM/Noah-MP)
 fractionLAI       = scalarExposedLAI / max(scalarExposedVAI, mpe)
 visibleAbsDirect  = spectralIncomingDirect(ixVisible)*spectralCanopyAbsorbedDirect(ixVisible)
 visibleAbsDiffuse = spectralIncomingDiffuse(ixVisible)*spectralCanopyAbsorbedDiffuse(ixVisible)
 if(scalarCanopySunlitFraction > tiny(scalarCanopySunlitFraction))then
  scalarCanopySunlitPAR = (visibleAbsDirect + scalarCanopySunlitFraction*visibleAbsDiffuse) * fractionLAI / max(scalarCanopySunlitLAI, mpe)
  scalarCanopyShadedPAR = (                   scalarCanopyShadedFraction*visibleAbsDiffuse) * fractionLAI / max(scalarCanopyShadedLAI, mpe)
 else
  scalarCanopySunlitPAR = 0._dp
  scalarCanopyShadedPAR = (visibleAbsDirect + visibleAbsDiffuse) * fractionLAI / max(scalarCanopyShadedLAI, mpe)
 endif
 !print*, 'scalarCanopySunlitLAI, fractionLAI, visibleAbsDirect, visibleAbsDiffuse, scalarCanopySunlitPAR = ', &
 !         scalarCanopySunlitLAI, fractionLAI, visibleAbsDirect, visibleAbsDiffuse, scalarCanopySunlitPAR


 
 end subroutine vegSWavRad


 ! *************************************************************************************************************************************
 ! ***** PRIVATE SUBROUTINE: gndAlbedo: compute the albedo of the ground surface
 ! *************************************************************************************************************************************
 subroutine gndAlbedo(&
                      ! input
                      soilTypeIndex,                         & ! intent(in): index of soil type
                      scalarGroundSnowFraction,              & ! intent(in): fraction of ground that is snow covered (-)
                      scalarVolFracLiqUpper,                 & ! intent(in): volumetric liquid water content in upper-most soil layer (-)
                      spectralSnowAlbedoDirect,              & ! intent(in): direct albedo of snow in each spectral band (-)
                      spectralSnowAlbedoDiffuse,             & ! intent(in): diffuse albedo of snow in each spectral band (-)
                      ! output
                      spectralAlbGndDirect,                  & ! intent(out): direct  albedo of underlying surface (-)
                      spectralAlbGndDiffuse,                 & ! intent(out): diffuse albedo of underlying surface (-)
                      err,message)                             ! intent(out): error control
 ! --------------------------------------------------------------------------------------------------------------------------------------
 ! identify parameters for soil albedo
 USE NOAHMP_RAD_PARAMETERS, only: ALBSAT,ALBDRY  ! Noah-MP: saturated and dry soil albedos for each wave band
 ! --------------------------------------------------------------------------------------------------------------------------------------
 ! input: model control
 integer(i4b),intent(in)        :: soilTypeIndex                ! index of soil type
 real(dp),intent(in)            :: scalarGroundSnowFraction     ! fraction of ground that is snow covered (-)
 real(dp),intent(in)            :: scalarVolFracLiqUpper        ! volumetric liquid water content in upper-most soil layer (-)
 real(dp),intent(in)            :: spectralSnowAlbedoDirect(:)  ! direct albedo of snow in each spectral band (-)
 real(dp),intent(in)            :: spectralSnowAlbedoDiffuse(:) ! diffuse albedo of snow in each spectral band (-)
 ! output
 real(dp),intent(out)           :: spectralAlbGndDirect(:)      ! direct  albedo of underlying surface (-)
 real(dp),intent(out)           :: spectralAlbGndDiffuse(:)     ! diffuse albedo of underlying surface (-)
 integer(i4b),intent(out)       :: err                          ! error code
 character(*),intent(out)       :: message                      ! error message
 ! local variables
 integer(i4b)                   :: iBand                        ! index of spectral band
 real(dp)                       :: xInc                         ! soil water correction factor for soil albedo
 real(dp),dimension(1:nBands)   :: spectralSoilAlbedo           ! soil albedo in each spectral band
 ! initialize error control
 err=0; message='gndAlbedo/'

 ! compute soil albedo
 do iBand=1,nBands   ! loop through spectral bands
  xInc = max(0.11_dp - 0.40_dp*scalarVolFracLiqUpper, 0._dp)
  spectralSoilAlbedo(iBand)  = min(ALBSAT(soilTypeIndex,iBand)+xInc,ALBDRY(soilTypeIndex,iBand))
 end do  ! (looping through spectral bands)

 ! compute surface albedo (weighted combination of snow and soil)
 do iBand=1,nBands
  spectralAlbGndDirect(iBand)  = (1._dp - scalarGroundSnowFraction)*spectralSoilAlbedo(iBand)  + scalarGroundSnowFraction*spectralSnowAlbedoDirect(iBand)
  spectralAlbGndDiffuse(iBand) = (1._dp - scalarGroundSnowFraction)*spectralSoilAlbedo(iBand)  + scalarGroundSnowFraction*spectralSnowAlbedoDiffuse(iBand)
 end do  ! (looping through spectral bands)

 end subroutine gndAlbedo


end module vegSWavRad_module

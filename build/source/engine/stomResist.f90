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

module stomResist_module
USE nrtype
! physical constants
USE multiconst, only: R_gas    ! universal gas constant (J mol-1 K-1)
! look-up values for the stomatal resistance formulation
USE mDecisions_module,only:  &
 simpleResistance,           & ! simple resistance formulation
 BallBerryFlex,              & ! flexible Ball-Berry scheme           
 BallBerry,                  & ! Ball-Berry (from Noah-MP)
 Jarvis                        ! Jarvis (from Noah-MP)

implicit none
private
public::stomResist
! spatial indices
integer(i4b),parameter :: iLoc = 1   ! i-location
integer(i4b),parameter :: jLoc = 1   ! j-location
! algorithmic parameters
real(dp),parameter     :: missingValue=-9999._dp  ! missing value, used when diagnostic or state variables are undefined
real(dp),parameter     :: mpe=1.e-6_dp         ! prevents overflow error if division by zero

contains


 ! ************************************************************************************************
 ! public subroutine stomResist: compute stomatal resistance
 ! ************************************************************************************************
 subroutine stomResist(&
                       ! input: state and diagnostic variables
                       scalarVegetationTemp,                & ! intent(in): vegetation temperature (K)
                       scalarSatVP_VegTemp,                 & ! intent(in): saturation vapor pressure at vegetation temperature (Pa)
                       scalarVP_CanopyAir,                  & ! intent(in): canopy air vapor pressure (Pa)
                       ! input: data structures
                       type_data,                           & ! intent(in):    type of vegetation and soil
                       attr_data,                           & ! intent(in):    spatial attributes
                       forc_data,                           & ! intent(in):    model forcing data
                       mpar_data,                           & ! intent(in):    model parameters
                       model_decisions,                     & ! intent(in):    model decisions
                       ! input-output: data structures
                       mvar_data,                           & ! intent(inout): model variables for a local HRU
                       ! output: error control
                       err,message)                           ! intent(out): error control
 ! ------------------------------------------------------------------------------------------------------------------------------------------------------
 ! ------------------------------------------------------------------------------------------------------------------------------------------------------
 ! provide access to the derived types to define the data structures
 USE data_struc,only:&
                     var_i,            & ! data vector (i4b)
                     var_d,            & ! data vector (dp)
                     var_dlength,      & ! data vector with variable length dimension (dp)
                     model_options       ! defines the model decisions
 ! provide access to named variables defining elements in the data structures
 USE var_lookup,only:iLookTIME,iLookTYPE,iLookATTR,iLookFORCE,iLookPARAM,iLookMVAR,iLookBVAR,iLookINDEX  ! named variables for structure elements
 USE var_lookup,only:iLookDECISIONS                           ! named variables for elements of the decision structure
 ! ------------------------------------------------------------------------------------------------------------------------------------------------------
 ! input: state and diagnostic variables
 real(dp),intent(in)             :: scalarVegetationTemp      ! vegetation temperature (K)
 real(dp),intent(in)             :: scalarSatVP_VegTemp       ! saturation vapor pressure at vegetation temperature (Pa)
 real(dp),intent(in)             :: scalarVP_CanopyAir        ! canopy air vapor pressure (Pa)
 ! input: data structures
 type(var_i),intent(in)          :: type_data                 ! type of vegetation and soil
 type(var_d),intent(in)          :: attr_data                 ! spatial attributes
 type(var_d),intent(in)          :: forc_data                 ! model forcing data
 type(var_d),intent(in)          :: mpar_data                 ! model parameters
 type(model_options),intent(in)  :: model_decisions(:)        ! model decisions
 ! input-output: data structures
 type(var_dlength),intent(inout) :: mvar_data                 ! model variables for a local HRU
 ! output: error control
 integer(i4b),intent(out)        :: err                       ! error code
 character(*),intent(out)        :: message                   ! error message
 ! ------------------------------------------------------------------------------------------------------------------------------------------------------
 ! local variables
 character(LEN=256)              :: cmessage                  ! error message of downwind routine
 integer(i4b),parameter          :: ixSunlit=1                ! named variable for sunlit leaves
 integer(i4b),parameter          :: ixShaded=2                ! named variable for shaded leaves
 integer(i4b)                    :: iSunShade                 ! index defining sunlit or shaded leaves
 real(dp),pointer                :: absorbedPAR               ! absorbed PAR (W m-2)
 real(dp)                        :: scalarStomResist          ! stomatal resistance (s m-1)
 real(dp)                        :: scalarPhotosynthesis      ! photosynthesis (umol CO2 m-2 s-1)

 ! associate variables in the data structure
 associate(&

 ! input: model decisions
 ix_stomResist                   => model_decisions(iLookDECISIONS%stomResist)%iDecision,           & ! intent(in): [i4b] choice of function for stomatal resistance

 ! input: physical attributes
 vegTypeIndex                    => type_data%var(iLookTYPE%vegTypeIndex),                          & ! intent(in): [i4b] vegetation type index
 minStomatalResistance           => mpar_data%var(iLookPARAM%minStomatalResistance),                & ! intent(in): [dp] mimimum stomatal resistance (s m-1)

 ! input: forcing at the upper boundary
 airtemp                         => forc_data%var(iLookFORCE%airtemp),                              & ! intent(in): [dp] air temperature at some height above the surface (K)
 airpres                         => forc_data%var(iLookFORCE%airpres),                              & ! intent(in): [dp] air pressure at some height above the surface (Pa)
 scalarO2air                     => mvar_data%var(iLookMVAR%scalarO2air)%dat(1),                    & ! intent(in): [dp] atmospheric o2 concentration (Pa)
 scalarCO2air                    => mvar_data%var(iLookMVAR%scalarCO2air)%dat(1),                   & ! intent(in): [dp] atmospheric co2 concentration (Pa)
 scalarCanopySunlitPAR           => mvar_data%var(iLookMVAR%scalarCanopySunlitPAR)%dat(1),          & ! intent(in): [dp] average absorbed par for sunlit leaves (w m-2)
 scalarCanopyShadedPAR           => mvar_data%var(iLookMVAR%scalarCanopyShadedPAR)%dat(1),          & ! intent(in): [dp] average absorbed par for shaded leaves (w m-2)

 ! input: state and diagnostic variables
 scalarGrowingSeasonIndex        => mvar_data%var(iLookMVAR%scalarGrowingSeasonIndex)%dat(1),       & ! intent(in): [dp] growing season index (0=off, 1=on)
 scalarFoliageNitrogenFactor     => mvar_data%var(iLookMVAR%scalarFoliageNitrogenFactor)%dat(1),    & ! intent(in): [dp] foliage nitrogen concentration (1.0 = saturated)
 scalarTranspireLim              => mvar_data%var(iLookMVAR%scalarTranspireLim)%dat(1),             & ! intent(in): [dp] weighted average of the transpiration limiting factor (-)
 scalarLeafResistance            => mvar_data%var(iLookMVAR%scalarLeafResistance)%dat(1),           & ! intent(in): [dp] mean leaf boundary layer resistance per unit leaf area (s m-1)

 ! output: stomatal resistance and photosynthesis
 scalarStomResistSunlit          => mvar_data%var(iLookMVAR%scalarStomResistSunlit)%dat(1),         & ! intent(out): [dp] stomatal resistance for sunlit leaves (s m-1)
 scalarStomResistShaded          => mvar_data%var(iLookMVAR%scalarStomResistShaded)%dat(1),         & ! intent(out): [dp] stomatal resistance for shaded leaves (s m-1)
 scalarPhotosynthesisSunlit      => mvar_data%var(iLookMVAR%scalarPhotosynthesisSunlit)%dat(1),     & ! intent(out): [dp] sunlit photosynthesis (umolco2 m-2 s-1)
 scalarPhotosynthesisShaded      => mvar_data%var(iLookMVAR%scalarPhotosynthesisShaded)%dat(1)      & ! intent(out): [dp] shaded photosynthesis (umolco2 m-2 s-1)
 )
 ! ------------------------------------------------------------------------------------------------------------------------------------------------------
 ! initialize error control
 err=0; message="stomResist/"

 ! identify option for stomatal resistance
 select case(ix_stomResist)

  ! *******************************************************************************************************************************************

  ! simple resistance formulation
  case(simpleResistance)
   ! check that we don't divide by zero -- should be set to minimum of tiny in routine soilResist
   if(scalarTranspireLim < tiny(airpres))then; err=20; message=trim(message)//'soil moisture stress factor is < tiny -- this will cause problems'; return; endif
   ! compute stomatal resistance (assume equal for sunlit and shaded leaves)
   scalarStomResistSunlit = minStomatalResistance/scalarTranspireLim
   scalarStomResistShaded = scalarStomResistSunlit
   ! set photosynthesis to missing (not computed)
   scalarPhotosynthesisSunlit = missingValue
   scalarPhotosynthesisShaded = missingValue

  ! *******************************************************************************************************************************************

  ! flexible Ball-Berry
  case(BallBerryFlex)

   ! loop through sunlit and shaded leaves
   do iSunShade=1,2

    ! get appropriate value for PAR
    select case(iSunShade)
     case(ixSunlit); absorbedPAR => scalarCanopySunlitPAR       ! average absorbed par for sunlit leaves (w m-2)
     case(ixShaded); absorbedPAR => scalarCanopyShadedPAR       ! average absorbed par for shaded leaves (w m-2)
     case default; err=20; message=trim(message)//'unable to identify case for sunlit/shaded leaves'; return
    end select

    ! compute photosynthesis and stomatal resistance
    call stomResist_flex(&
                         ! input: state and diagnostic variables
                         scalarVegetationTemp,                & ! intent(in): vegetation temperature (K)
                         scalarSatVP_VegTemp,                 & ! intent(in): saturation vapor pressure at vegetation temperature (Pa)
                         scalarVP_CanopyAir,                  & ! intent(in): canopy air vapor pressure (Pa)
                         absorbedPAR,                         & ! intent(in): absorbed PAR (W m-2)
                         ! input: data structures
                         forc_data,                           & ! intent(in): model forcing data
                         mpar_data,                           & ! intent(in): model parameters
                         mvar_data,                           & ! intent(in): model variables for a local HRU
                         ! output:
                         scalarStomResist,                    & ! intent(out): stomatal resistance (s m-1)
                         scalarPhotosynthesis,                & ! intent(out): photosynthesis (umol CO2 m-2 s-1)
                         ! output: error control
                         err,cmessage)                          ! intent(out): error control
    if(err/=0)then; message=trim(message)//trim(cmessage); return; endif

    ! assign output variables
    select case(iSunShade)
     case(ixSunlit)
      scalarStomResistSunlit     = scalarStomResist 
      scalarPhotosynthesisSunlit = scalarPhotosynthesis
     case(ixShaded)
      scalarStomResistShaded     = scalarStomResist
      scalarPhotosynthesisShaded = scalarPhotosynthesis
     case default; err=20; message=trim(message)//'unable to identify case for sunlit/shaded leaves'; return
    end select

   end do  ! looping through sunlit and shaded leaves


  ! *******************************************************************************************************************************************

  ! compute stomatal resistance (wrapper around the Noah-MP routines)
  ! NOTE: canopy air vapor pressure is from the previous time step
  case(BallBerry,Jarvis)
   call stomResist_NoahMP(&
                          ! input (model decisions)
                          ix_stomResist,                     & ! intent(in): choice of function for stomatal resistance
                          ! input (local attributes)
                          vegTypeIndex,                      & ! intent(in): vegetation type index
                          iLoc, jLoc,                        & ! intent(in): spatial location indices
                          ! input (forcing)
                          airtemp,                           & ! intent(in): air temperature at some height above the surface (K)
                          airpres,                           & ! intent(in): air pressure at some height above the surface (Pa)
                          scalarO2air,                       & ! intent(in): atmospheric o2 concentration (Pa)
                          scalarCO2air,                      & ! intent(in): atmospheric co2 concentration (Pa)
                          scalarCanopySunlitPAR,             & ! intent(in): average absorbed par for sunlit leaves (w m-2)
                          scalarCanopyShadedPAR,             & ! intent(in): average absorbed par for shaded leaves (w m-2)
                          ! input (state and diagnostic variables)
                          scalarGrowingSeasonIndex,          & ! intent(in): growing season index (0=off, 1=on)
                          scalarFoliageNitrogenFactor,       & ! intent(in): foliage nitrogen concentration (1=saturated)
                          scalarTranspireLim,                & ! intent(in): weighted average of the soil moiture factor controlling stomatal resistance (-)
                          scalarLeafResistance,              & ! intent(in): leaf boundary layer resistance (s m-1)
                          scalarVegetationTemp,              & ! intent(in): temperature of the vegetation canopy (K)
                          scalarSatVP_VegTemp,               & ! intent(in): saturation vapor pressure at the temperature of the veg canopy (Pa)
                          scalarVP_CanopyAir,                & ! intent(in): canopy air vapor pressure (Pa)
                          ! output
                          scalarStomResistSunlit,            & ! intent(out): stomatal resistance for sunlit leaves (s m-1)
                          scalarStomResistShaded,            & ! intent(out): stomatal resistance for shaded leaves (s m-1)
                          scalarPhotosynthesisSunlit,        & ! intent(out): sunlit photosynthesis (umolco2 m-2 s-1)
                          scalarPhotosynthesisShaded,        & ! intent(out): shaded photosynthesis (umolco2 m-2 s-1)
                          err,cmessage                       ) ! intent(out): error control
   if(err/=0)then; message=trim(message)//trim(cmessage); return; endif

  ! *******************************************************************************************************************************************

  ! error check
  case default; err=20; message=trim(message)//'unable to identify option for stomatal resistance'; return

  ! *******************************************************************************************************************************************

 end select  ! (identifying option for stomatal resistance)

 ! print progress
 !write(*,'(a,1x,L1,1x,20(f16.8,1x))') 'ix_StomResist==BallBerryFlex, scalarPhotosynthesisSunlit, scalarPhotosynthesisShaded, scalarStomResistSunlit, scalarPhotosynthesisShaded = ', &
 !                                      ix_StomResist==BallBerryFlex, scalarPhotosynthesisSunlit, scalarPhotosynthesisShaded, scalarStomResistSunlit, scalarPhotosynthesisShaded
 !pause

 ! end association to variables in the data structures
 end associate

 end subroutine stomResist


 ! *******************************************************************************************************
 ! *******************************************************************************************************
 ! *** PRIVATE SUBROUTINES *******************************************************************************
 ! *******************************************************************************************************
 ! *******************************************************************************************************

 ! *******************************************************************************************************
 ! private subroutine stomResist_flex: flexible stomatal resistance routine to evaluate different options
 ! *******************************************************************************************************
 subroutine stomResist_flex(&
                            ! input: state and diagnostic variables
                            scalarVegetationTemp,                & ! intent(in): vegetation temperature (K)
                            scalarSatVP_VegTemp,                 & ! intent(in): saturation vapor pressure at vegetation temperature (Pa)
                            scalarVP_CanopyAir,                  & ! intent(in): canopy air vapor pressure (Pa)
                            absorbedPAR,                         & ! intent(in): absorbed PAR (W m-2)
                            ! input: data structures
                            forc_data,                           & ! intent(in): model forcing data
                            mpar_data,                           & ! intent(in): model parameters
                            mvar_data,                           & ! intent(in): model variables for a local HRU
                            ! output: stomatal resistance and photosynthesis
                            scalarStomResist,                    & ! intent(out): stomatal resistance (s m-1)
                            scalarPhotosynthesis,                & ! intent(out): photosynthesis (umol CO2 m-2 s-1)
                            ! output: error control
                            err,message)                           ! intent(out): error control

 ! ------------------------------------------------------------------------------------------------------------------------------------------------------
 ! ------------------------------------------------------------------------------------------------------------------------------------------------------
 ! provide access to the derived types to define the data structures
 USE data_struc,only:&
                     var_d,            & ! data vector (dp)
                     var_dlength         ! data vector with variable length dimension (dp)
 ! provide access to named variables defining elements in the data structures
 USE var_lookup,only:iLookTIME,iLookTYPE,iLookATTR,iLookFORCE,iLookPARAM,iLookMVAR,iLookBVAR,iLookINDEX  ! named variables for structure elements
 ! ------------------------------------------------------------------------------------------------------------------------------------------------------
 ! input: state and diagnostic variables
 real(dp),intent(in)             :: scalarVegetationTemp       ! vegetation temperature (K)
 real(dp),intent(in)             :: scalarSatVP_VegTemp        ! saturation vapor pressure at vegetation temperature (Pa)
 real(dp),intent(in)             :: scalarVP_CanopyAir         ! canopy air vapor pressure (Pa)
 real(dp),intent(in)             :: absorbedPAR                ! absorbed PAR (W m-2)
 ! input: data structures
 type(var_d),intent(in)          :: forc_data                  ! model forcing data
 type(var_d),intent(in)          :: mpar_data                  ! model parameters
 type(var_dlength),intent(in)    :: mvar_data                  ! model variables for a local HRU
 ! output: stomatal resistance and photosynthesis
 real(dp),intent(out)            :: scalarStomResist           ! stomatal resistance (s m-1)
 real(dp),intent(out)            :: scalarPhotosynthesis       ! photosynthesis (umol CO2 m-2 s-1)
 ! output: error control
 integer(i4b),intent(out)        :: err                        ! error code
 character(*),intent(out)        :: message                    ! error message
 ! ------------------------------------------------------------------------------------------------------------------------------------------------------
 ! local variables
 character(LEN=256)              :: cmessage                   ! error message of downwind routine
 real(dp),parameter              :: convToler=0.0001_dp        ! convergence tolerance (Pa)
 real(dp),parameter              :: umol_per_mol=1.e+6_dp      ! factor to relate umol to mol
 real(dp),parameter              :: o2scaleFactor=0.105_dp     ! scaling factor used to compute co2 compesation point (0.21/2)
 real(dp),parameter              :: h2o_co2__leafbl=1.37_dp    ! factor to represent the different diffusivities of h2o and co2 in the leaf boundary layer (-)
 real(dp),parameter              :: h2o_co2__stomPores=1.65_dp ! factor to represent the different diffusivities of h2o and co2 in the stomatal pores (-)
 real(dp),parameter              :: joule2umolConv=4.6_dp      ! conversion factor from joules to umol photons (umol J-1)
 real(dp),parameter              :: quantumYield=0.06_dp       ! quantam yield (mol co2 mol-1 photon)
 integer(i4b),parameter          :: nFactors=3                 ! number of limiting factors for assimilation (light, Rubisco, and export)
 integer(i4b),parameter          :: ixLight=1                  ! named variable for light-limited assimilation
 integer(i4b),parameter          :: ixRubi=2                   ! named variable for Rubisco-limited assimilation
 integer(i4b),parameter          :: ixExp=3                    ! named variable for export-limited assimilation
 integer(i4b)                    :: ixLimitVec(1),ixLimit      ! index of factor limiting assimilation
 real(dp)                        :: unitConv                   ! unit conversion factor (mol m-3, convert m s-1 --> mol H20 m-2 s-1)
 real(dp)                        :: xFac(nFactors)             ! temporary variable used to compute assimilation rate
 real(dp)                        :: xPSN(nFactors)             ! assimilation rate for different factors (light, Rubisco, and export)
 real(dp)                        :: x0,x1,x2,x3                ! temporary variables
 real(dp)                        :: leafConductance            ! leaf conductance (umol m-2 s-1)
 real(dp)                        :: Kc,Ko                      ! Michaelis-Menten constants for co2 and o2 (Pa)
 real(dp)                        :: vcmax                      ! maximum Rubisco carboxylation rate (umol m-2 s-1)
 real(dp)                        :: J                          ! electron transport rate (umol co2 m-2 s-1)
 real(dp)                        :: co2compPt                  ! co2 compensation point (Pa)
 real(dp)                        :: awb                        ! Michaelis-Menten control (Pa)
 real(dp)                        :: ciDiff                     ! difference between intercellular co2 concentraion and the compensation point (Pa)
 real(dp)                        :: cp2                        ! additional controls in light-limited assimilation (Pa)
 real(dp)                        :: ci,ci_old                  ! intercellular co2 partial pressure (Pa)
 real(dp)                        :: cs                         ! co2 partial pressure at leaf surface (Pa)
 real(dp)                        :: gs                         ! stomatal conductance (umol m-2 s-1)
 real(dp)                        :: psn                        ! leaf gross photosynthesis rate (umol co2 m-2 s-1)
 real(dp)                        :: fnf=0.6666667_dp           ! foliage nitrogen factor, fraction [0,1]
 real(dp)                        :: fHum                       ! humidity function, fraction [0,1] 
 real(dp)                        :: dA_dc                      ! derivative in photosynthesis w.r.t. intercellular co2 concentration 
 real(dp)                        :: dgs_dc                     ! derivative in stomatal conductance w.r.t. intercellular co2 concentration
 real(dp)                        :: dci_dc                     ! final derivative (-)
 real(dp)                        :: xInc                       ! iteration increment (Pa)
 integer(i4b)                    :: iter                       ! iteration index
 integer(i4b),parameter          :: maxiter=10                 ! maximum number of iterations
 ! ------------------------------------------------------------------------------------------------------------------------------------------------------
 ! associate variables in the data structure
 associate(&

 ! input: model parameters
 Kc25                            => mpar_data%var(iLookPARAM%Kc25),                                 & ! intent(in): [dp] Michaelis-Menten constant for CO2 at 25 degrees C (Pa)
 Ko25                            => mpar_data%var(iLookPARAM%Ko25),                                 & ! intent(in): [dp] Michaelis-Menten constant for O2 at 25 degrees C (Pa)
 vcmax25                         => mpar_data%var(iLookPARAM%vcmax25),                              & ! intent(in): [dp] potential carboxylation rate at 25 degrees C (umol co2 m-2 s-1)
 Kc_fac                          => mpar_data%var(iLookPARAM%Kc_fac),                               & ! intent(in): [dp] factor in the q10 function defining temperature controls on Kc (-)
 Ko_fac                          => mpar_data%var(iLookPARAM%Ko_fac),                               & ! intent(in): [dp] factor in the q10 function defining temperature controls on Ko (-)
 vcmax_fac                       => mpar_data%var(iLookPARAM%vcmax_fac),                            & ! intent(in): [dp] factor in the q10 function defining temperature controls on vcmax (-)
 hightemp_delS                   => mpar_data%var(iLookPARAM%hightemp_delS),                        & ! intent(in): [dp] entropy term in high temp inhibition function for vcmax (J K-1 mol-1)
 hightemp_delH                   => mpar_data%var(iLookPARAM%hightemp_delH),                        & ! intent(in): [dp] deactivation energy in high temp inhibition function for vcmax (J mol-1)
 cond2photo_slope                => mpar_data%var(iLookPARAM%cond2photo_slope),                     & ! intent(in): [dp] slope of conductance-photosynthesis relationship (-)
 minStomatalConductance          => mpar_data%var(iLookPARAM%minStomatalConductance),               & ! intent(in): [dp] mimimum stomatal conductance (umol H2O m-2 s-1)

 ! input: forcing at the upper boundary
 airtemp                         => forc_data%var(iLookFORCE%airtemp),                              & ! intent(in): [dp] air temperature at some height above the surface (K)
 airpres                         => forc_data%var(iLookFORCE%airpres),                              & ! intent(in): [dp] air pressure at some height above the surface (Pa)
 scalarO2air                     => mvar_data%var(iLookMVAR%scalarO2air)%dat(1),                    & ! intent(in): [dp] atmospheric o2 concentration (Pa)
 scalarCO2air                    => mvar_data%var(iLookMVAR%scalarCO2air)%dat(1),                   & ! intent(in): [dp] atmospheric co2 concentration (Pa)

 ! input: state and diagnostic variables
 scalarGrowingSeasonIndex        => mvar_data%var(iLookMVAR%scalarGrowingSeasonIndex)%dat(1),       & ! intent(in): [dp] growing season index (0=off, 1=on)
 scalarFoliageNitrogenFactor     => mvar_data%var(iLookMVAR%scalarFoliageNitrogenFactor)%dat(1),    & ! intent(in): [dp] foliage nitrogen concentration (1.0 = saturated)
 scalarTranspireLim              => mvar_data%var(iLookMVAR%scalarTranspireLim)%dat(1),             & ! intent(in): [dp] weighted average of the transpiration limiting factor (-)
 scalarLeafResistance            => mvar_data%var(iLookMVAR%scalarLeafResistance)%dat(1)            & ! intent(in): [dp] mean leaf boundary layer resistance per unit leaf area (s m-1)

 )
 ! ------------------------------------------------------------------------------------------------------------------------------------------------------
 ! initialize error control
 err=0; message="stomResist_flex/"

 ! define unit conversion (m s-1 --> mol m-2 s-1)
 ! NOTE: R_gas   = J K-1 Mol-1 (J = kg m2 s-2); airtemp = K; airpres = Pa (kg m-1 s-2)
 unitConv = airpres/(R_gas*airtemp)  ! mol m-3

 ! check there is light available for photosynthesis
 if(absorbedPAR < tiny(absorbedPAR) .or. scalarGrowingSeasonIndex < tiny(absorbedPAR))then
  scalarStomResist     = unitConv*umol_per_mol/(scalarTranspireLim*minStomatalConductance)
  scalarPhotosynthesis = 0._dp
  return
 endif

 ! define the foloage nitrogen factor (-)
 ! NOTE: need to implement something more sophisticated
 fnf=0.66666667_dp

 ! define the humidity function
 ! NOTE: need to implement an alternative
 fHum = min( max(0.25_dp, scalarVP_CanopyAir/scalarSatVP_VegTemp), 1._dp)

 ! compute the leaf conductance (umol m-2 s-1)
 leafConductance = umol_per_mol*unitConv/scalarLeafResistance  ! s m-1 --> umol m-2 s-1

 ! compute the Michaelis-Menten constants (Pa)
 Kc = Kc25*q10(Kc_fac,scalarVegetationTemp)
 Ko = Ko25*q10(Ko_fac,scalarVegetationTemp)

 ! compute maximum Rubisco carboxylation rate (umol co2 m-2 s-1)
 x0 = q10(vcmax_fac,scalarVegetationTemp)  ! temperature function
 x1 = fHigh(hightemp_delS,hightemp_delH,scalarVegetationTemp) ! high temperature inhibition function
 vcmax = vcmax25*fnf*scalarTranspireLim*x0/x1
 write(*,'(a,1x,20(f16.8,1x))') 'x0, x1, vcmax25, scalarVegetationTemp = ', x0, x1, vcmax25, scalarVegetationTemp-273.15_dp

 ! compute the electron transport rate (umol CO2 m-2 s-1)
 J = quantumYield*joule2umolConv*absorbedPAR  

 ! compute the co2 compensation point (Pa)
 co2compPt = (Kc/Ko)*scalarO2air*o2scaleFactor 

 ! compute the Michaelis-Menten controls (Pa)
 awb = Kc*(1._dp + scalarO2air/Ko)

 ! compute the additional controls in light-limited assimilation
 cp2 = co2compPt*2._dp

 ! define trial value of intercellular co2 (Pa)
 ci = 0.7_dp*scalarCO2air

 ! print progress
 write(*,'(a,1x,20(f16.8,1x))') 'fHum, vcmax, J = ', &
                                 fHum, vcmax, J

 ! ***
 ! iterate
 do iter=1,maxiter

  ! reset ci
  ci_old = ci

  ! compute the difference between intercellular co2 concentraion and the compensation point
  ciDiff = max(ci - co2compPt, 0.00001_dp)  ! small value to ensure convergence

  ! compute light-limited assimilation
  xFac(ixLight) = J/(ci + cp2)         ! umol co2 m-2 s-1 Pa-1
  xPSN(ixLight) = xFac(ixLight)*ciDiff ! umol co2 m-2 s-1

  ! compute Rubisco-limited assimilation
  xFac(ixRubi) = vcmax/(ci + awb)      ! umol co2 m-2 s-1 Pa-1
  xPSN(ixRubi) = xFac(ixRubi)*ciDiff   ! umol co2 m-2 s-1

  ! compute export limited assimilation
  xFac(ixExp) = 0.5_dp
  xPSN(ixExp) = xFac(ixExp)*vcmax      ! umol co2 m-2 s-1

  ! identify limiting factor
  ixLimitVec = minloc(xPSN)
  ixLimit = ixLimitVec(1)

  ! define photosynthesis





  x0  = xFac(ixLimit)
  psn = xPSN(ixLimit)





  ! define derivative in photosynthesis w.r.t. intercellular co2 concentration
  select case(ixLimit)
   case(ixLight); dA_dc = x0*x0*(co2compPt + cp2)/J
   case(ixRubi);  dA_dc = x0*x0*(co2compPt + awb)/vcmax
   case(ixExp);   dA_dc = 0._dp
   case default; message=trim(message)//'unable to identify limiting factor for photosynthesis'; err=20; return
  end select

  ! compute co2 concentration at leaf surface (Pa)
  x1 = h2o_co2__leafbl * airpres/leafConductance  ! Pa / (umol co2 m-2 s-1)
  cs = max(scalarCO2air - (x1 * psn), mpe)   ! Pa (avoid divide by zero)

  ! compute stomatal conductance (umol m-2 s-1)
  x2 = cond2photo_slope*airpres*fHum  ! Pa
  gs = x2*psn/cs + scalarTranspireLim*minStomatalConductance ! umol m-2 s-1

  ! compute intercellular co2 partial pressues (Pa)
  x3 = h2o_co2__stomPores * airpres  ! Pa
  ci = cs - x3*psn/gs  ! Pa


  ! compute derivatives
  if(ixLimit==ixLight .or. ixLimit==ixRubi)then
   dgs_dc = x2*dA_dc*(x1*psn/cs + 1._dp)/cs  ! derivative in stomatal conductance w.r.t. intercellular co2 concentration
   dci_dc = -x1*dA_dc + x3*(psn*dgs_dc/gs - dA_dc)/gs  ! final derivative
  else
   dci_dc = 0._dp  ! assimilation does not depend on intercellular co2 concentration
  endif

  ! compute iteration increment (Pa)
  xInc = (ci - ci_old)/(1._dp - dci_dc)

  ! update
  ci = ci_old + xInc

  ! print progress
  write(*,'(a,1x,2(i4,1x),20(f14.5,1x))') 'iter, ixLimit, xPSN, ci, co2compPt, cs, gs, xInc = ', &
                                           iter, ixLimit, xPSN, ci, co2compPt, cs, gs, xInc

  ! check for convergence
  if(abs(xInc) < convToler) exit
  if(iter==maxIter)then
   message=trim(message)//'did not converge in stomatal conductance iteration'
   err=20; return
  endif

 end do  ! iterating

 ! assign output variables
 scalarStomResist     = unitConv*umol_per_mol/gs  ! umol m-2 s-1 --> s/m
 scalarPhotosynthesis = psn

 end associate

 contains

  ! function for temperature dependence
  function q10(a,T)
  implicit none
  real(dp),intent(in) :: a              ! scale factor
  real(dp),intent(in) :: T              ! temperature (K)
  real(dp)            :: q10            ! temperature dependence (-)
  real(dp),parameter  :: Tmid=298.15_dp ! point where function is one (25 deg C)
  real(dp),parameter  :: Tscale=10._dp  ! scaling factor (K)
  q10 = a**((T - Tmid)/Tscale)
  end function q10

  ! function for high temperature inhibition
  function fHigh(delS,delH,T)
  implicit none
  real(dp),intent(in) :: delS     ! entropy term in high temp inhibition function (J K-1 mol-1) 
  real(dp),intent(in) :: delH     ! deactivation energy in high temp inhibition function (J mol-1)
  real(dp),intent(in) :: T        ! temperature (K)
  real(dp)            :: fHigh    ! high temperature inhibition (-)
  fHigh = 1._dp + exp( (delS*T - delH)/(R_gas*T) ) ! NOTE: R_gas = J K-1 mol-1
  end function fHigh

 end subroutine stomResist_flex


 ! *******************************************************************************************************
 ! private subroutine stomResist_NoahMP: use Noah-MP routines to compute stomatal resistance
 ! *******************************************************************************************************
 subroutine stomResist_NoahMP(&
                              ! input (model decisions)
                              ixStomResist,                        & ! intent(in): choice of function for stomatal resistance
                              ! input (local attributes)
                              vegTypeIndex,                        & ! intent(in): vegetation type index
                              iLoc, jLoc,                          & ! intent(in): spatial location indices
                              ! input (forcing)
                              airtemp,                             & ! intent(in): air temperature at some height above the surface (K)
                              airpres,                             & ! intent(in): air pressure at some height above the surface (Pa)
                              scalarO2air,                         & ! intent(in): atmospheric o2 concentration (Pa)
                              scalarCO2air,                        & ! intent(in): atmospheric co2 concentration (Pa)
                              scalarCanopySunlitPAR,               & ! intent(in): average absorbed par for sunlit leaves (w m-2)
                              scalarCanopyShadedPAR,               & ! intent(in): average absorbed par for shaded leaves (w m-2)
                              ! input (state and diagnostic variables)
                              scalarGrowingSeasonIndex,            & ! intent(in): growing season index (0=off, 1=on)
                              scalarFoliageNitrogenFactor,         & ! intent(in): foliage nitrogen concentration (1=saturated)
                              scalarTranspireLim,                  & ! intent(in): weighted average of the soil moiture factor controlling stomatal resistance (-)
                              scalarLeafResistance,                & ! intent(in): leaf boundary layer resistance (s m-1)
                              scalarVegetationTemp,                & ! intent(in): vegetation temperature (K)
                              scalarSatVP_VegTemp,                 & ! intent(in): saturation vapor pressure at vegetation temperature (Pa)
                              scalarVP_CanopyAir,                  & ! intent(in): canopy air vapor pressure (Pa)
                              ! output
                              scalarStomResistSunlit,              & ! intent(out): stomatal resistance for sunlit leaves (s m-1)
                              scalarStomResistShaded,              & ! intent(out): stomatal resistance for shaded leaves (s m-1)
                              scalarPhotosynthesisSunlit,          & ! intent(out): sunlit photosynthesis (umolco2 m-2 s-1)
                              scalarPhotosynthesisShaded,          & ! intent(out): shaded photosynthesis (umolco2 m-2 s-1)
                              err,message                          ) ! intent(out): error control
 ! -----------------------------------------------------------------------------------------------------------------------------------------
 ! Modified from Noah-MP
 ! Compute stomatal resistance and photosynthesis using either
 !  1) Ball-Berry
 !  2) Jarvis
 ! See Niu et al. JGR 2011 for more details
 USE mDecisions_module, only: BallBerry,Jarvis                ! options for the choice of function for stomatal resistance
 USE NOAHMP_ROUTINES,only:stomata                             ! compute canopy resistance based on Ball-Berry
 USE NOAHMP_ROUTINES,only:canres                              ! compute canopy resistance based Jarvis
 implicit none
 ! input (model decisions)
 integer(i4b),intent(in)       :: ixStomResist                ! choice of function for stomatal resistance
 ! input (local attributes)
 integer(i4b),intent(in)       :: vegTypeIndex                ! vegetation type index
 integer(i4b),intent(in)       :: iLoc, jLoc                  ! spatial location indices
 ! input (forcing)
 real(dp),intent(in)           :: airtemp                     ! measured air temperature at some height above the surface (K)
 real(dp),intent(in)           :: airpres                     ! measured air pressure at some height above the surface (Pa)
 real(dp),intent(in)           :: scalarO2air                 ! atmospheric o2 concentration (Pa)
 real(dp),intent(in)           :: scalarCO2air                ! atmospheric co2 concentration (Pa)
 real(dp),intent(in),target    :: scalarCanopySunlitPAR       ! average absorbed par for sunlit leaves (w m-2)
 real(dp),intent(in),target    :: scalarCanopyShadedPAR       ! average absorbed par for shaded leaves (w m-2)
 ! input (state and diagnostic variables)
 real(dp),intent(in)           :: scalarGrowingSeasonIndex    ! growing season index (0=off, 1=on)
 real(dp),intent(in)           :: scalarFoliageNitrogenFactor ! foliage nitrogen concentration (1=saturated)
 real(dp),intent(in)           :: scalarTranspireLim          ! weighted average of the soil moiture factor controlling stomatal resistance (-)
 real(dp),intent(in)           :: scalarLeafResistance        ! leaf boundary layer resistance (s m-1)
 real(dp),intent(in)           :: scalarVegetationTemp        ! vegetation temperature (K)
 real(dp),intent(in)           :: scalarSatVP_VegTemp         ! saturation vapor pressure at vegetation temperature (Pa)
 real(dp),intent(in)           :: scalarVP_CanopyAir          ! canopy air vapor pressure (Pa)
 ! output
 real(dp),intent(out)          :: scalarStomResistSunlit      ! stomatal resistance for sunlit leaves (s m-1)
 real(dp),intent(out)          :: scalarStomResistShaded      ! stomatal resistance for shaded leaves (s m-1)
 real(dp),intent(out)          :: scalarPhotosynthesisSunlit  ! sunlit photosynthesis (umolco2 m-2 s-1)
 real(dp),intent(out)          :: scalarPhotosynthesisShaded  ! sunlit photosynthesis (umolco2 m-2 s-1)
 integer(i4b),intent(out)      :: err                         ! error code
 character(*),intent(out)      :: message                     ! error message
 ! local variables
 integer(i4b),parameter        :: ixSunlit=1                  ! named variable for sunlit leaves
 integer(i4b),parameter        :: ixShaded=2                  ! named variable for shaded leaves
 integer(i4b)                  :: iSunShade                   ! index for sunlit/shaded leaves
 real(dp),pointer              :: PAR                         ! average absorbed PAR for sunlit/shaded leaves (w m-2)
 real(dp)                      :: scalarStomResist            ! stomatal resistance for sunlit/shaded leaves (s m-1)
 real(dp)                      :: scalarPhotosynthesis        ! photosynthesis for sunlit/shaded leaves (umolco2 m-2 s-1)
 ! initialize error control
 err=0; message='stomResist_NoahMP/'

 ! loop through sunlit and shaded leaves
 do iSunShade=1,2

  ! get appropriate value for PAR
  select case(iSunShade)
   case(ixSunlit); PAR => scalarCanopySunlitPAR               ! average absorbed par for sunlit leaves (w m-2)
   case(ixShaded); PAR => scalarCanopyShadedPAR               ! average absorbed par for shaded leaves (w m-2)
   case default; err=20; message=trim(message)//'unable to identify case for sunlit/shaded leaves'; return
  end select

  ! identify option for stomatal resistance
  select case(ixStomResist)

   ! Ball-Berry
   case(BallBerry)
   call stomata(&
                ! input
                vegTypeIndex,                       & ! intent(in): vegetation type index
                mpe,                                & ! intent(in): prevents overflow error if division by zero
                PAR,                                & ! intent(in): average absorbed par (w m-2)
                scalarFoliageNitrogenFactor,        & ! intent(in): foliage nitrogen concentration (1=saturated)
                iLoc, jLoc,                         & ! intent(in): spatial location indices
                scalarVegetationTemp,               & ! intent(in): vegetation temperature (K)
                scalarSatVP_VegTemp,                & ! intent(in): saturation vapor pressure at vegetation temperature (Pa)
                scalarVP_CanopyAir,                 & ! intent(in): canopy air vapor pressure (Pa)
                airtemp,                            & ! intent(in): air temperature at some height above the surface (K)
                airpres,                            & ! intent(in): air pressure at some height above the surface (Pa)
                scalarO2air,                        & ! intent(in): atmospheric o2 concentration (Pa)
                scalarCO2air,                       & ! intent(in): atmospheric co2 concentration (Pa)
                scalarGrowingSeasonIndex,           & ! intent(in): growing season index (0=off, 1=on)
                scalarTranspireLim,                 & ! intent(in): weighted average of the soil moiture factor controlling stomatal resistance (-)
                scalarLeafResistance,               & ! intent(in): leaf boundary layer resistance (s m-1)
                ! output
                scalarStomResist,                   & ! intent(out): stomatal resistance (s m-1)
                scalarPhotosynthesis                ) ! intent(out): photosynthesis (umolco2 m-2 s-1)

   ! Jarvis
   case(Jarvis)
   call canres(&
                ! input
                PAR,                                & ! intent(in): average absorbed par (w m-2)
                scalarVegetationTemp,               & ! intent(in): vegetation temperature (K)
                scalarTranspireLim,                 & ! intent(in): weighted average of the soil moiture factor controlling stomatal resistance (-)
                scalarVP_CanopyAir,                 & ! intent(in): canopy air vapor pressure (Pa)
                airpres,                            & ! intent(in): air pressure at some height above the surface (Pa)
                ! output
                scalarStomResist,                   & ! intent(out): stomatal resistance (s m-1)
                scalarPhotosynthesis,               & ! intent(out): photosynthesis (umolco2 m-2 s-1)
                ! location indices (input)
                iLoc, jLoc                          ) ! intent(in): spatial location indices

   ! check identified an option
   case default; err=20; message=trim(message)//'unable to identify case for stomatal resistance'; return

  end select  ! (selecting option for stomatal resistance)

  ! assign output variables
  select case(iSunShade)
   case(ixSunlit)
    scalarStomResistSunlit     = scalarStomResist
    scalarPhotosynthesisSunlit = scalarPhotosynthesis
   case(ixShaded)
    scalarStomResistShaded     = scalarStomResist
    scalarPhotosynthesisShaded = scalarPhotosynthesis
   case default; err=20; message=trim(message)//'unable to identify case for sunlit/shaded leaves'; return
  end select

 end do  ! (looping through sunlit and shaded leaves)

 end subroutine stomResist_NoahMP
 
 ! -- end private subroutines
 ! ------------------------------------------------------------------------------------------------------------
 ! ------------------------------------------------------------------------------------------------------------
 ! ------------------------------------------------------------------------------------------------------------

end module stomResist_module

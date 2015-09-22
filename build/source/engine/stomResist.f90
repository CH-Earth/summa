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
USE multiconst, only: Rgas     ! universal gas constant (J mol-1 K-1)
USE multiconst, only: Tfreeze  ! freezing point of pure water (K)
! conversion functions
USE conv_funcs_module,only:satVapPress   ! function to compute the saturated vapor pressure (Pa)
! look-up values for the stomatal resistance formulation
USE mDecisions_module,only:  &
 simpleResistance,           & ! simple resistance formulation
 BallBerryTest,              & ! flexible Ball-Berry scheme (testing)           
 BallBerryFlex,              & ! flexible Ball-Berry scheme           
 BallBerry,                  & ! Ball-Berry (from Noah-MP)
 Jarvis                        ! Jarvis (from Noah-MP)
! look-up values for the leaf temperature controls on photosynthesis + stomatal resistance
USE mDecisions_module,only:  &
 q10Func,                    & ! the q10 function used in CLM4 and Noah-MP
 Arrhenius                     ! the Arrhenius functions used in CLM5 and Cable 
! look-up values for the humidity controls on stomatal resistance
USE mDecisions_module,only:  &
 humidLeafSurface,           & ! humidity at the leaf surface [Bonan et al., 2011]
 scaledHyperbolic              ! scaled hyperbolic function [Leuning et al., 1995]
! look-up values for the electron transport function, dependence of photosynthesis on PAR
USE mDecisions_module,only:  &
 linear,                     & ! linear function used in CLM4 and Noah-MP
 linearJmax,                 & ! linear jmax function used in Cable [Wang et al., Ag Forest Met 1998, eq D5]
 quadraticJmax                 ! the quadratic Jmax function, used in SSiB and CLM5
! look-up values for the CO2 compensation point to calculate stomatal resistance
USE mDecisions_module,only:  &
 origBWB,                    & ! the original BWB function
 Leuning                       ! the Leuning function
implicit none
private
public::stomResist
! spatial indices
integer(i4b),parameter :: iLoc = 1   ! i-location
integer(i4b),parameter :: jLoc = 1   ! j-location
! conversion factors
real(dp),parameter     :: joule2umolConv=4.6_dp   ! conversion factor from joules to umol photons (umol J-1)
! algorithmic parameters
real(dp),parameter     :: missingValue=-9999._dp  ! missing value, used when diagnostic or state variables are undefined
real(dp),parameter     :: mpe=1.e-6_dp            ! prevents overflow error if division by zero
real(dp),parameter     :: dx=1.e-6_dp             ! finite difference increment

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
 ! provide access to the output directory
 USE summaFileManager,only:OUTPUT_PATH,OUTPUT_PREFIX         ! define output file
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
 real(dp)                        :: absorbedPAR               ! absorbed PAR (W m-2)
 real(dp)                        :: scalarStomResist          ! stomatal resistance (s m-1)
 real(dp)                        :: scalarPhotosynthesis      ! photosynthesis (umol CO2 m-2 s-1)
 ! ------------------------------------------------------------------------------------------------------------------------------------------------------
 ! local variables for testing
 integer(i4b),parameter          :: ixUnit=66                 ! file unit
 integer(i4b)                    :: iTrial                    ! index of model trial
 integer(i4b),parameter          :: nTrial=100                ! number of model trials
 integer(i4b)                    :: iExp                      ! index of model experiment
 integer(i4b),parameter          :: nExp=3                    ! number of model experiments
 integer(i4b),parameter          :: ixTemp=1                  ! named variable defining perturbations in leaf temperature
 integer(i4b),parameter          :: ixVPD=2                   ! named variable defining perturbations in vapor pressure deficit
 integer(i4b),parameter          :: ixPAR=3                   ! named variable defining perturbations in PAR
 real(dp)                        :: xTrial                    ! normalized trial value
 real(dp)                        :: leafTemp                  ! leaf temperature (K)
 real(dp)                        :: vpd                       ! vapor pressure deficit (Pa)
 real(dp)                        :: par                       ! absorbed PAR (W m-2)
 real(dp)                        :: SatVP_leafTemp            ! saturated vapor pressure at the leaf temperature (Pa)
 real(dp)                        :: derivNotUsed              ! derivative in saturated vapor pressure w.r.t. temperature (not used)
 real(dp)                        :: vpCanair                  ! vapor pressure of the canopy air space (Pa)
 real(dp)                        :: unitConv                  ! unit conversion factor (mol m-3, convert m s-1 --> mol H20 m-2 s-1)
 character(LEN=256)              :: filename                  ! filename
 character(LEN=4),parameter      :: cTemp='Temp'              ! named variable for temperature
 character(LEN=4),parameter      :: cVPD='VPD '               ! named variable for vapor pressure deficit
 character(LEN=4),parameter      :: cPAR='PAR '               ! named variable for PAR
 ! ------------------------------------------------------------------------------------------------------------------------------------------------------

 ! associate variables in the data structure
 associate(&

 ! input: model decisions
 ix_stomResist                   => model_decisions(iLookDECISIONS%stomResist)%iDecision,           & ! intent(in): [i4b] choice of function for stomatal resistance

 ! input: physical attributes
 vegTypeIndex                    => type_data%var(iLookTYPE%vegTypeIndex),                          & ! intent(in): [i4b] vegetation type index
 minStomatalResistance           => mpar_data%var(iLookPARAM%minStomatalResistance),                & ! intent(in): [dp] mimimum stomatal resistance (s m-1)
 vcmax25                         => mpar_data%var(iLookPARAM%vcmax25),                              & ! intent(in): [dp] potential carboxylation rate at 25 degrees C (umol co2 m-2 s-1)

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
     case(ixSunlit); absorbedPAR = scalarCanopySunlitPAR       ! average absorbed par for sunlit leaves (w m-2)
     case(ixShaded); absorbedPAR = scalarCanopyShadedPAR       ! average absorbed par for shaded leaves (w m-2)
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
                         model_decisions,                     & ! intent(in): model decisions
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
  case(BallBerryTest)
   ! synthetic test to explore basic parameter sensitivity

   ! set vcmax25 to 40
   vcmax25 = 40._dp

   ! set air pressure to the standard atmosphere
   airpres = 101325._dp

   ! define partial pressure of o2 and co2
   scalarO2air  = 0.209_dp    * airpres ! from Bonan et al. JGR 2011
   scalarCO2air = 0.000379_dp * airpres ! from Bonan et al. JGR 2011

   ! set the growing season index to 1
   scalarGrowingSeasonIndex = 1._dp

   ! set the leaf boundary layer conductance to 0.05 m s-1
   scalarLeafResistance = 1._dp / 0.05_dp  

   ! set the soil moisture control
   scalarTranspireLim = 1._dp

   ! define unit conversion (m s-1 --> mol m-2 s-1)
   ! NOTE: Rgas   = J K-1 Mol-1 (J = kg m2 s-2); airtemp = K; airpres = Pa (kg m-1 s-2)
   unitConv = airpres/(Rgas*airtemp)  ! mol m-3

   ! loop through experiments
   do iExp=1,nExp

    ! define default values
    leafTemp =   25._dp + Tfreeze
    vpd      =    0._dp
    par      = 2000._dp / joule2umolConv 

    ! define output file
    select case(iExp)
     case(ixTemp); filename=trim(OUTPUT_PATH)//trim(OUTPUT_PREFIX)//'_stomatalResistance.'//trim(cTemp)//'.txt'
     case(ixVPD);  filename=trim(OUTPUT_PATH)//trim(OUTPUT_PREFIX)//'_stomatalResistance.'//trim(cVPD)//'.txt'
     case(ixPAR);  filename=trim(OUTPUT_PATH)//trim(OUTPUT_PREFIX)//'_stomatalResistance.'//trim(cPAR)//'.txt'
     case default; err=20; message=trim(message)//'unable to identify case for forcing variable'; return
    end select

    ! open up output file
    open(unit=ixUnit,file=trim(filename))

    ! loop through trial values
    do iTrial=1,nTrial

     ! define a number between zero and one
     xTrial = real(iTrial-1, kind(1._dp)) / real(nTrial-1, kind(1._dp))

     ! get forcing variables for a given experiment
     select case(iExp)
      case(ixTemp); leafTemp = xTrial*39._dp + 0.5_dp + Tfreeze
      case(ixVPD);  vpd      = xTrial*3000._dp
      case(ixPAR);  par      = xTrial*500._dp
      case default; err=20; message=trim(message)//'unable to identify case for forcing variable'; return
     end select

     ! compute saturated vapor pressure at the leaf temperature
     call satVapPress(leafTemp-Tfreeze, SatVP_leafTemp, derivNotUsed)

     ! compute the vapor pressure of the canopy air space
     vpCanair = SatVP_leafTemp - vpd

     ! compute photosynthesis and stomatal resistance
     call stomResist_flex(&
                          ! input: state and diagnostic variables
                          leafTemp,                            & ! intent(in): vegetation temperature (K)
                          SatVP_leafTemp,                      & ! intent(in): saturation vapor pressure at vegetation temperature (Pa)
                          vpCanair,                            & ! intent(in): canopy air vapor pressure (Pa)
                          par,                                 & ! intent(in): absorbed PAR (W m-2)
                          ! input: data structures
                          forc_data,                           & ! intent(in): model forcing data
                          mpar_data,                           & ! intent(in): model parameters
                          mvar_data,                           & ! intent(in): model variables for a local HRU
                          model_decisions,                     & ! intent(in): model decisions
                          ! output:
                          scalarStomResist,                    & ! intent(out): stomatal resistance (s m-1)
                          scalarPhotosynthesis,                & ! intent(out): photosynthesis (umol CO2 m-2 s-1)
                          ! output: error control
                          err,cmessage)                          ! intent(out): error control
     if(err/=0)then; message=trim(message)//trim(cmessage); return; endif

     ! print progress
     !write(*,'(10(f12.5,1x))') leafTemp-Tfreeze, vpd, par, unitConv/scalarStomResist, scalarPhotosynthesis
     write(ixUnit,'(10(f12.5,1x))') leafTemp-Tfreeze, vpd, par, unitConv/scalarStomResist, scalarPhotosynthesis

    end do  ! looping through trial values

    ! close file unit
    close(ixUnit)

   end do  ! looping through experiments

   stop 'completed test simulation'

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
                            model_decisions,                     & ! intent(in): model decisions
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
                     var_dlength,      & ! data vector with variable length dimension (dp)
                     model_options       ! defines the model decisions
 ! provide access to named variables defining elements in the data structures
 USE var_lookup,only:iLookTIME,iLookTYPE,iLookATTR,iLookFORCE,iLookPARAM,iLookMVAR,iLookBVAR,iLookINDEX  ! named variables for structure elements
 USE var_lookup,only:iLookDECISIONS                           ! named variables for elements of the decision structure
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
 type(model_options),intent(in)  :: model_decisions(:)         ! model decisions
 ! output: stomatal resistance and photosynthesis
 real(dp),intent(out)            :: scalarStomResist           ! stomatal resistance (s m-1)
 real(dp),intent(out)            :: scalarPhotosynthesis       ! photosynthesis (umol CO2 m-2 s-1)
 ! output: error control
 integer(i4b),intent(out)        :: err                        ! error code
 character(*),intent(out)        :: message                    ! error message
 ! ------------------------------------------------------------------------------------------------------------------------------------------------------
 ! general local variables
 character(LEN=256)              :: cmessage                   ! error message of downwind routine
 logical(lgt),parameter          :: testDerivs=.false.         ! flag to test the derivatives
 real(dp)                        :: unitConv                   ! unit conversion factor (mol m-3, convert m s-1 --> mol H20 m-2 s-1)
 real(dp)                        :: rlb                        ! leaf boundary layer rersistance (umol-1 m2 s)
 real(dp)                        :: x0,x1,x2,x3,x4             ! temporary variables
 real(dp)                        :: co2compPt                  ! co2 compensation point (Pa)
 real(dp)                        :: fHum                       ! humidity function, fraction [0,1] 
 ! ------------------------------------------------------------------------------------------------------------------------------------------------------
 ! fixed parameters
 integer(i4b),parameter          :: maxiter=10                 ! maximum number of iterations
 integer(i4b),parameter          :: maxiter_noahMP=3           ! maximum number of iterations for Noah-MP
 real(dp),parameter              :: convToler=0.0001_dp        ! convergence tolerance (Pa)
 real(dp),parameter              :: umol_per_mol=1.e+6_dp      ! factor to relate umol to mol
 real(dp),parameter              :: o2scaleFactor=0.105_dp     ! scaling factor used to compute co2 compesation point (0.21/2)
 real(dp),parameter              :: h2o_co2__leafbl=1.37_dp    ! factor to represent the different diffusivities of h2o and co2 in the leaf boundary layer (-)
 real(dp),parameter              :: h2o_co2__stomPores=1.65_dp ! factor to represent the different diffusivities of h2o and co2 in the stomatal pores (-)
 real(dp),parameter              :: Tref=298.16_dp             ! reference temperature (25 deg C)
 real(dp),parameter              :: Tscale=10._dp              ! scaling factor in q10 function (K)
 real(dp),parameter              :: c_ps2=0.7_dp               ! curvature factor for electron transport (-)
 ! ------------------------------------------------------------------------------------------------------------------------------------------------------
 ! photosynthesis
 real(dp)                        :: Kc,Ko                      ! Michaelis-Menten constants for co2 and o2 (Pa)
 real(dp)                        :: vcmax                      ! maximum Rubisco carboxylation rate (umol m-2 s-1)
 real(dp)                        :: jmax25                     ! maximum electron transport rate at 25 deg C (umol m-2 s-1)
 real(dp)                        :: jmax                       ! maximum electron transport rate (umol m-2 s-1)
 real(dp)                        :: aQuad                      ! the quadratic coefficient in the quadratic equation 
 real(dp)                        :: bQuad                      ! the linear coefficient in the quadratic equation
 real(dp)                        :: cQuad                      ! the constant in the quadratic equation
 real(dp)                        :: bSign                      ! sign of the linear coeffcient
 real(dp)                        :: xTemp                      ! temporary variable in the quadratic equation
 real(dp)                        :: qQuad                      ! the "q" term in the quadratic equation
 real(dp)                        :: root1,root2                ! roots of the quadratic function
 real(dp)                        :: Js                         ! scaled electron transport rate (umol co2 m-2 s-1)
 real(dp)                        :: I_ps2                      ! PAR absorbed by PS2 (umol photon m-2 s-1)
 real(dp)                        :: awb                        ! Michaelis-Menten control (Pa)
 real(dp)                        :: cp2                        ! additional controls in light-limited assimilation (Pa)
 real(dp)                        :: psn                        ! leaf gross photosynthesis rate (umol co2 m-2 s-1)
 real(dp)                        :: dA_dc                      ! derivative in photosynthesis w.r.t. intercellular co2 concentration 
 ! ------------------------------------------------------------------------------------------------------------------------------------------------------
 ! stomatal resistance
 real(dp)                        :: gMin                       ! scaled minimum conductance (umol m-2 s-1)
 real(dp)                        :: cs                         ! co2 partial pressure at leaf surface (Pa)
 real(dp)                        :: csx                        ! control of co2 partial pressure at leaf surface on stomatal conductance (Pa)
 real(dp)                        :: g0                         ! stomatal conductance in the absence of humidity controls (umol m-2 s-1)
 real(dp)                        :: ci,ci_old                  ! intercellular co2 partial pressure (Pa)
 real(dp)                        :: rs                         ! stomatal resistance (umol-1 m2 s)
 real(dp)                        :: dg0_dc                     ! derivative in g0 w.r.t intercellular co2 concentration (umol m-2 s-1 Pa-1)
 real(dp)                        :: drs_dc                     ! derivative in stomatal resistance w.r.t. intercellular co2 concentration
 real(dp)                        :: dci_dc                     ! final derivative (-)
 ! ------------------------------------------------------------------------------------------------------------------------------------------------------
 ! iterative solution
 integer(i4b),parameter          :: ixSumma=1001               ! named variable for the Noah-MP solution (fixed point iteration, max 3 iterations)
 integer(i4b),parameter          :: ixNoahMP=1002              ! named variable for the Noah-MP solution (fixed point iteration, max 3 iterations)
 integer(i4b)                    :: ixSolution=ixSumma         ! define solution method
 real(dp)                        :: func1,func2                ! functions for numerical derivative calculation
 real(dp)                        :: cMin,cMax                  ! solution brackets
 real(dp)                        :: xInc                       ! iteration increment (Pa)
 integer(i4b)                    :: iter                       ! iteration index
 ! ------------------------------------------------------------------------------------------------------------------------------------------------------
 ! associate variables in the data structure
 associate(&

 ! input: model decisions
 ix_bbTempFunc                   => model_decisions(iLookDECISIONS%bbTempFunc)%iDecision,           & ! intent(in): [i4b] leaf temperature controls on photosynthesis + stomatal resistance
 ix_bbHumdFunc                   => model_decisions(iLookDECISIONS%bbHumdFunc)%iDecision,           & ! intent(in): [i4b] humidity controls on stomatal resistance 
 ix_bbElecFunc                   => model_decisions(iLookDECISIONS%bbElecFunc)%iDecision,           & ! intent(in): [i4b] dependence of photosynthesis on PAR
 ix_bbCO2point                   => model_decisions(iLookDECISIONS%bbCO2point)%iDecision,           & ! intent(in): [i4b] use of CO2 compensation point to calculate stomatal resistance

 ! input: model parameters
 Kc25                            => mpar_data%var(iLookPARAM%Kc25),                                 & ! intent(in): [dp] Michaelis-Menten constant for CO2 at 25 degrees C (umol mol-1)
 Ko25                            => mpar_data%var(iLookPARAM%Ko25),                                 & ! intent(in): [dp] Michaelis-Menten constant for O2 at 25 degrees C (mol mol-1)
 Kc_qFac                         => mpar_data%var(iLookPARAM%Kc_qFac),                              & ! intent(in): [dp] factor in the q10 function defining temperature controls on Kc (-)
 Ko_qFac                         => mpar_data%var(iLookPARAM%Ko_qFac),                              & ! intent(in): [dp] factor in the q10 function defining temperature controls on Ko (-)
 kc_Ha                           => mpar_data%var(iLookPARAM%kc_Ha),                                & ! intent(in): [dp] activation energy for the Michaelis-Menten constant for CO2 (J mol-1)
 ko_Ha                           => mpar_data%var(iLookPARAM%ko_Ha),                                & ! intent(in): [dp] activation energy for the Michaelis-Menten constant for O2 (J mol-1)
 vcmax25                         => mpar_data%var(iLookPARAM%vcmax25),                              & ! intent(in): [dp] potential carboxylation rate at 25 degrees C (umol co2 m-2 s-1)
 vcmax_qFac                      => mpar_data%var(iLookPARAM%vcmax_qFac),                           & ! intent(in): [dp] factor in the q10 function defining temperature controls on vcmax (-)
 vcmax_Ha                        => mpar_data%var(iLookPARAM%vcmax_Ha),                             & ! intent(in): [dp] activation energy in the vcmax function (J mol-1)
 vcmax_Hd                        => mpar_data%var(iLookPARAM%vcmax_Hd),                             & ! intent(in): [dp] deactivation energy in the vcmax function (J mol-1)
 vcmax_Sv                        => mpar_data%var(iLookPARAM%vcmax_Sv),                             & ! intent(in): [dp] entropy term in the vcmax function (J mol-1 K-1)
 jmax25_scale                    => mpar_data%var(iLookPARAM%jmax25_scale),                         & ! intent(in): [dp] scaling factor to relate jmax25 to vcmax25 (-)
 jmax_Ha                         => mpar_data%var(iLookPARAM%jmax_Ha),                              & ! intent(in): [dp] activation energy in the jmax function (J mol-1)
 jmax_Hd                         => mpar_data%var(iLookPARAM%jmax_Hd),                              & ! intent(in): [dp] deactivation energy in the jmax function (J mol-1)
 jmax_Sv                         => mpar_data%var(iLookPARAM%jmax_Sv),                              & ! intent(in): [dp] entropy term in the jmax function (J mol-1 K-1)
 fractionJ                       => mpar_data%var(iLookPARAM%fractionJ),                            & ! intent(in): [dp] fraction of light lost by other than the chloroplast lamellae (-)
 quantamYield                    => mpar_data%var(iLookPARAM%quantamYield),                         & ! intent(in): [dp] quantam yield (mol e mol-1 q)
 vpScaleFactor                   => mpar_data%var(iLookPARAM%vpScaleFactor),                        & ! intent(in): [dp] vapor pressure scaling factor in stomatal conductance function (Pa)
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

 !print*, '**'

 ! *****
 ! * preliminaries...
 ! ******************

 ! define unit conversion (m s-1 --> mol m-2 s-1)
 ! NOTE: Rgas   = J K-1 Mol-1 (J = kg m2 s-2); airtemp = K; airpres = Pa (kg m-1 s-2)
 unitConv = airpres/(Rgas*airtemp)  ! mol m-3

 ! check there is light available for photosynthesis
 if(absorbedPAR < tiny(absorbedPAR) .or. scalarGrowingSeasonIndex < tiny(absorbedPAR))then
  scalarStomResist     = unitConv*umol_per_mol/(scalarTranspireLim*minStomatalConductance)
  scalarPhotosynthesis = 0._dp
  return
 endif

 ! define the scaled minimum conductance
 gMin = scalarTranspireLim*minStomatalConductance

 ! compute the leaf conductance (umol m-2 s-1)
 rlb = scalarLeafResistance/(umol_per_mol*unitConv)  ! s m-1 --> umol-1 m2 s

 ! compute the maximum electron transport rate at 25 deg C (umol m-2 s-1)
 jmax25 = jmax25_scale * vcmax25

 ! *****
 ! * compute temperature controls on stomatal conductance...
 ! *********************************************************

 ! identify the temperature function
 select case(ix_bbTempFunc)

  ! q10 function used in CLM4 and Noah-MP
  case(q10Func)
   ! compute the Michaelis-Menten constants (Pa)
   Kc = airpres*(Kc25/umol_per_mol)*q10(Kc_qFac,scalarVegetationTemp,Tref,Tscale)  ! umol mol-1 --> mol mol-1 --> Pa
   Ko = airpres*Ko25*q10(Ko_qFac,scalarVegetationTemp,Tref,Tscale)  ! mol mol-1 --> Pa
   ! compute maximum Rubisco carboxylation rate (umol co2 m-2 s-1)
   x0 = q10(vcmax_qFac,scalarVegetationTemp,Tref,Tscale)  ! temperature function
   x1 = fHigh(vcmax_Hd,vcmax_Sv,scalarVegetationTemp) ! high temperature inhibition function
   vcmax = vcmax25*x0/x1

  ! Arrhenius function used in CLM5 and Cable
  case(Arrhenius)
   ! compute the Michaelis-Menten constants (Pa)
   Kc = airpres*(Kc25/umol_per_mol)*fT(kc_Ha,scalarVegetationTemp,Tref)  ! umol mol-1 --> mol mol-1 --> Pa
   Ko = airpres*Ko25*fT(ko_Ha,scalarVegetationTemp,Tref)  ! mol mol-1 --> Pa
   ! compute maximum Rubisco carboxylation rate (umol co2 m-2 s-1)
   x0 = fT(vcmax_Ha,scalarVegetationTemp,Tref)
   x1 = fhigh(vcmax_Hd,vcmax_Sv,Tref) / fhigh(vcmax_Hd,vcmax_Sv,scalarVegetationTemp)
   vcmax = vcmax25*x0*x1
   ! compute the maximum electron transport rate (umol electron m-2 s-1)
   x0 = fT(jmax_Ha,scalarVegetationTemp,Tref)
   x1 = fhigh(jmax_Hd,jmax_Sv,Tref) / fhigh(jmax_Hd,jmax_Sv,scalarVegetationTemp)
   jmax = jmax25*x0*x1

  ! check found an appropriate option
  case default; err=20; message=trim(message)//'unable to find option for leaf temperature controls on stomatal conductance'; return

 end select  ! temperature controls

 ! *****
 ! * compute electron transport controls on stomatal conductance...
 ! ****************************************************************

 ! identify the electron transport function
 select case(ix_bbElecFunc)

  ! linear model, as used in CLM4 and Noah-MP
  case(linear); Js = quantamYield*joule2umolConv*absorbedPAR

  ! linear function of qmax, as used in Cable [Wang et al., Ag Forest Met 1998, eq D5]
  case(linearJmax);
   x0 = quantamYield*joule2umolConv*absorbedPAR
   x1 = x0*jmax / (x0 + 2.1_dp*jmax)
   Js = x1/4._dp ! scaled electron transport

  ! quadraric function of jmax, as used in CLM5 (Bonan et al., JGR 2011, Table B2)
  case(quadraticJmax)
   !  PAR absorbed by PS2 (umol photon m-2 s-1)
   I_ps2 = 0.5_dp*(1._dp - fractionJ) * joule2umolConv*absorbedPAR   ! Farquar (1980), eq 8: PAR absorbed by PS2 (umol photon m-2 s-1)
   ! define coefficients in the quadratic equation
   aQuad = c_ps2            ! quadratic coefficient = cuurvature factor for electron transport
   bQuad = -(I_ps2 + jmax)  ! linear coefficient
   cQuad =  I_ps2 * jmax    ! free term
   ! compute the q term (NOTE: bQuad is always positive)
   bSign = abs(bQuad)/bQuad
   xTemp = bQuad*bQuad - 4._dp *aQuad*cQuad
   qQuad = -0.5_dp * (bQuad + bSign*sqrt(xTemp))
   ! compute roots
   root1 = qQuad / aQuad
   root2 = cQuad / qQuad
   ! select minimum root, required to ensure J=0 when par=0
   ! NOTE: Wittig et al. select the first root, which is the max in all cases I tried
   Js = min(root1,root2) / 4._dp  ! scaled J

  ! check found an appropriate option
  case default; err=20; message=trim(message)//'unable to find option for electron transport controls on stomatal conductance'; return

 end select  ! electron transport controls

 ! *****
 ! * define additional controls on stomatal conductance...
 ! ****************************************************************

 ! define the humidity function
 select case(ix_bbHumdFunc)
  case(humidLeafSurface); fHum = min( max(0.25_dp, scalarVP_CanopyAir/scalarSatVP_VegTemp), 1._dp)
  case(scaledHyperbolic); fHum = (scalarSatVP_VegTemp - scalarVP_CanopyAir)/vpScaleFactor
  case default; err=20; message=trim(message)//'unable to identify humidity control on stomatal conductance'; return
 end select

 ! compute the co2 compensation point (Pa)
 co2compPt = (Kc/Ko)*scalarO2air*o2scaleFactor 

 ! compute the Michaelis-Menten controls (Pa)
 awb = Kc*(1._dp + scalarO2air/Ko)

 ! compute the additional controls in light-limited assimilation
 cp2 = co2compPt*2._dp

 ! define trial value of intercellular co2 (Pa)
 ci = 0.7_dp*scalarCO2air

 ! initialize brackets for the solution
 cMin = 0._dp
 cMax = scalarCO2air

 ! *********************************************************************************************************************************
 ! *********************************************************************************************************************************
 ! *********************************************************************************************************************************
 ! *********************************************************************************************************************************
 ! *********************************************************************************************************************************
 ! *********************************************************************************************************************************

 !print *, '**'
 !print *, '**'

 ! ***
 ! iterate
 do iter=1,maxiter

  ! reset ci
  ci_old = ci

  ! *****
  ! * compute photosynthesis and stomatal resistance...
  ! ***************************************************

  ! compute gross photosynthesis [follow Farquar (Planta, 1980), as implemented in CLM4 and Noah-MP]
  call photosynthesis(.true., ci, co2compPt, awb, cp2, vcmax, Js, psn, dA_dc)

  ! compute co2 concentration at leaf surface (Pa)
  x1 = h2o_co2__leafbl * airpres * rlb  ! Pa / (umol co2 m-2 s-1)
  cs = max(scalarCO2air - (x1 * psn), mpe)   ! Pa (avoid divide by zero)

  ! compute control of the compensation point on stomatal conductance
  if(ix_bbCO2point == origBWB)then
   csx = cs
  else
   csx = cs - co2compPt
  endif

  ! compute conductance in the absence of humidity
  g0     = cond2photo_slope*airpres*psn/csx 
  dg0_dc = cond2photo_slope*airpres*dA_dc*(x1*psn/cs + 1._dp)/csx

  ! use quadratic function to compute stomatal resistance
  call quadResist(.true.,ix_bbHumdFunc,rlb,fHum,gMin,g0,dg0_dc,rs,drs_dc)

  ! compute intercellular co2 partial pressues (Pa)
  x2 = h2o_co2__stomPores * airpres  ! Pa
  ci = max(cs - x2*psn*rs, 0._dp)    ! Pa

  ! final derivative
  dci_dc = -x1*dA_dc - x2*(psn*drs_dc + rs*dA_dc)

  ! test derivatives
  if(testDerivs)then
   func1 = testFunc(ci_old,    cond2photo_slope, airpres, scalarCO2air, ix_bbHumdFunc, ix_bbCO2point)
   func2 = testFunc(ci_old+dx, cond2photo_slope, airpres, scalarCO2air, ix_bbHumdFunc, ix_bbCO2point)
   write(*,'(a,1x,20(e20.10,1x))') '(func2 - func1)/dx, dci_dc = ', &
                                    (func2 - func1)/dx, dci_dc
  endif  ! if testing the derivatives

  ! *****
  ! * iterative solution...
  ! ***********************

  ! CLM4 and Noah-MP use fixed point iteration, continuing the next iteration from this point
  ! here we try and improve matters by calculating the derivatives

  ! update the brackets for the solution
  if(ci_old > ci)then
   cMax = ci_old
  else
   cMin = ci_old
  endif

  ! case for Noah-MP (use fixed point iteration)
  if(ixSolution==ixNoahMP)then
   if(iter==maxiter_NoahMP) exit ! exit after a specified number of iterations (normally 3)
   cycle  ! fixed-point iteration
  endif

  ! compute iteration increment (Pa)
  xInc = (ci - ci_old)/(1._dp - dci_dc)

  ! update
  ci = max(ci_old + xInc, 0._dp)

  ! ensure that we stay within brackets
  if(ci > cMax .or. ci < cMin)then
   ci = 0.5_dp * (cMin + cMax)
  endif

  ! print progress
  !write(*,'(a,1x,i4,1x,20(f12.7,1x))') 'iter, cMin, cMax, psn, ci_old, ci, scalarCO2air, rs, xInc = ', &
  !                                      iter, cMin, cMax, psn, ci_old, ci, scalarCO2air, rs, xInc

  ! check for convergence
  if(abs(xInc) < convToler) exit
  if(iter==maxIter)then
   message=trim(message)//'did not converge in stomatal conductance iteration'
   err=20; return
  endif

 end do  ! iterating
 !pause 'iterating'

 ! assign output variables
 scalarStomResist     = unitConv*umol_per_mol*rs  ! umol m-2 s-1 --> s/m
 scalarPhotosynthesis = psn

 end associate

 contains

  ! ******************************************************
  ! ******************************************************

  ! internal function used to test derivatives
  function testFunc(ci, cond2photo_slope, airpres, scalarCO2air, ix_bbHumdFunc, ix_bbCO2point)
  real(dp),intent(in)     :: ci, cond2photo_slope, airpres, scalarCO2air
  integer(i4b),intent(in) :: ix_bbHumdFunc, ix_bbCO2point
  real(dp)                :: testFunc
  real(dp)                :: unUsed

  ! compute gross photosynthesis [follow Farquar (Planta, 1980), as implemented in CLM4 and Noah-MP]
  call photosynthesis(.false., ci, co2compPt, awb, cp2, vcmax, Js, psn, unUsed)

  ! compute co2 concentration at leaf surface (Pa)
  x1 = h2o_co2__leafbl * airpres * rlb  ! Pa / (umol co2 m-2 s-1)
  if(ix_bbCO2point == origBWB)then
   cs = max(scalarCO2air - (x1 * psn), mpe)   ! Pa (avoid divide by zero)
  else
   cs = max(scalarCO2air - (x1 * psn), mpe) - co2compPt
  endif

  ! compute conductance in the absence of humidity
  g0     = cond2photo_slope*airpres*psn/cs
  dg0_dc = 0._dp 

  ! use quadratic function to compute stomatal resistance
  call quadResist(.false.,ix_bbHumdFunc,rlb,fHum,gMin,g0,dg0_dc,rs,unUsed)

  ! compute intercellular co2 partial pressues (Pa)
  x2 = h2o_co2__stomPores * airpres  ! Pa
  testFunc = max(cs - x2*psn*rs, 0._dp)    ! Pa

  end function testFunc

 end subroutine stomResist_flex

 ! *******************************************************************************************************
 ! private subroutine photosynthesis: compute gross photosynthesis
 ! *******************************************************************************************************
 subroutine photosynthesis(desireDeriv, ci, co2compPt, awb, cp2, vcmax, Js, psn, dA_dc)
 implicit none
 ! dummy variables
 logical(lgt),intent(in) :: desireDeriv   ! .true. if the derivative is desired
 real(dp),intent(in)     :: ci            ! intercellular co2 concentration (Pa)
 real(dp),intent(in)     :: co2compPt     ! co2 compensation point (Pa)
 real(dp),intent(in)     :: awb           ! Michaelis-Menten control (Pa)
 real(dp),intent(in)     :: cp2           ! additional controls in light-limited assimilation (Pa)
 real(dp),intent(in)     :: vcmax         ! maximum Rubisco carboxylation rate (umol co2 m-2 s-1)
 real(dp),intent(in)     :: Js            ! scaled electron transport rate (umol co2 m-2 s-1)
 real(dp),intent(out)    :: psn           ! leaf gross photosynthesis rate (umol co2 m-2 s-1)
 real(dp),intent(out)    :: dA_dc         ! derivative in photosynthesis w.r.t. intercellular co2 concentration (umol co2 m-2 s-1 Pa-1)
 ! local variables
 logical(lgt),parameter          :: colimitation=.true.         ! colimitation, as described by Collatz et al. 1991, Sellers et al., 1996, and Bonan et al., 2011
 integer(i4b),parameter          :: nFactors=3                  ! number of limiting factors for assimilation (light, Rubisco, and export)
 integer(i4b),parameter          :: ixRubi=1                    ! named variable for Rubisco-limited assimilation
 integer(i4b),parameter          :: ixLight=2                   ! named variable for light-limited assimilation
 integer(i4b),parameter          :: ixExport=3                  ! named variable for export-limited assimilation
 integer(i4b)                    :: ixLimitVec(1),ixLimit       ! index of factor limiting assimilation
 real(dp)                        :: xFac(nFactors)              ! temporary variable used to compute assimilation rate
 real(dp)                        :: xPSN(nFactors)              ! assimilation rate for different factors (light, Rubisco, and export)
 real(dp)                        :: ciDiff                      ! difference between intercellular co2 and the co2 compensation point
 real(dp)                        :: ciDer                       ! factor to account for constainted intercellular co2 in calculating derivatives
 real(dp)                        :: x0                          ! temporary variable
 real(dp)                        :: xsPSN                       ! intermediate smoothed photosynthesis
 real(dp)                        :: dAc_dc,dAj_dc,dAe_dc,dAi_dc ! derivatives in assimilation w.r.t. intercellular co2 concentration
 real(dp),parameter              :: theta_cj=0.98_dp            ! coupling coefficient (see Sellers et al., 1996 [eq C6]; Bonan et al., 2011 [Table B1])
 real(dp),parameter              :: theta_ie=0.95_dp            ! coupling coefficient (see Sellers et al., 1996 [eq C6]; Bonan et al., 2011 [Table B1])
 ! ------------------------------------------------------------
 ! this method follows Farquar (Planta, 1980), as implemented in CLM4 and Noah-MP

 ! compute the difference between intercellular co2 concentraion and the compensation point
 ciDiff = max(0._dp, ci - co2compPt)

 ! impose constraints (NOTE: derivative is zero if constraints are imposed)
 if(ci < co2compPt)then; ciDer = 0._dp; else; ciDer = 1._dp; endif

 ! compute Rubisco-limited assimilation
 xFac(ixRubi) = vcmax/(ci + awb)      ! umol co2 m-2 s-1 Pa-1
 xPSN(ixRubi) = xFac(ixRubi)*ciDiff   ! umol co2 m-2 s-1

 ! compute light-limited assimilation
 xFac(ixLight) = Js/(ci + cp2)        ! umol co2 m-2 s-1 Pa-1
 xPSN(ixLight) = xFac(ixLight)*ciDiff ! umol co2 m-2 s-1

 ! compute export limited assimilation
 xFac(ixExport) = 0.5_dp
 xPSN(ixExport) = xFac(ixExport)*vcmax   ! umol co2 m-2 s-1

 !write(*,'(a,1x,10(f12.5,1x))') 'xPSN, vcmax, ciDiff, ci, awb = ', xPSN, vcmax, ciDiff, ci, awb 

 ! if not colimitation (standard CLM4 and NoahMP) 
 if(.not.colimitation)then

  ! identify limiting factor
  ixLimitVec = minloc(xPSN)
  ixLimit = ixLimitVec(1)

  ! define photosynthesis
  x0  = xFac(ixLimit)
  psn = xPSN(ixLimit)

  ! if derivatives are desired
  if(desireDeriv)then

   ! compute derivatives in assimilation (no colimitation)
   select case(ixLimit)
    case(ixRubi);   dA_dc = x0*ciDer - ciDiff*x0*x0/vcmax  ! Rubisco-limited assimilation
    case(ixLight);  dA_dc = x0*ciDer - ciDiff*x0*x0/Js     ! light-limited assimilation
    case(ixExport); dA_dc = 0._dp                          ! export-limited assimilation
   endselect

  ! derivatives are not desired
  else
   dA_dc = 0._dp
  endif

 ! * colimitation (Collatz et al., 1991; Sellers et al., 1996; Bonan et al., 2011)
 else

  ! compute derivatives for individual terms
  if(desireDeriv)then
   dAc_dc = xFac(ixRubi)*ciDer - ciDiff*xFac(ixRubi)*xFac(ixRubi)/vcmax
   dAj_dc = xFac(ixLight)*ciDer - ciDiff*xFac(ixLight)*xFac(ixLight)/Js
   dAe_dc = 0._dp
  else
   dAc_dc = 0._dp
   dAj_dc = 0._dp
   dAe_dc = 0._dp
  endif

  ! smooth Rubisco-limitation and light limitation
  call quadSmooth(desireDeriv, xPSN(ixRubi), xPSN(ixLight), theta_cj, dAc_dc, dAj_dc, xsPSN, dAi_dc)

  ! smooth intermediate-limitation and export limitation
  call quadSmooth(desireDeriv, xsPSN, xPSN(ixExport), theta_ie, dAi_dc, dAe_dc, psn, dA_dc)

 endif  ! colimitation

 end subroutine photosynthesis

 ! *******************************************************************************************************
 ! private subroutine quadResist: compute stomatal resistance
 ! *******************************************************************************************************

 ! use quadratic function to compute stomatal resistance

 ! this method follows CLM4, described most fully in Oleson et al. (NCAR Tech. Note, 2010)
 ! details are also provided in Sellers et al., part 1 (J. Climate, 1996) and Bonan et al. (JGR 2011)

 ! stomatal conductance can be given as
 !     1/rs = m * (A/cs) * (es/ei) * Patm + b * beta   ! see Bonan et al. (2011) for inclusion of beta in the 2nd term
 ! here es is the (unknown) vapor pressure at the leaf surface

 ! the photosynthesis (computed above) assumes equality in co2 gradients between the atmosphere and the leaf surface,
 !  and between the leaf surface and the leaf interior, as
 !     A = (ca - cs)/(1.37*rb*Patm) = (cs - ci)/(1.65*rs*Patm)
 !  which requires that
 !     (ea - ei)/(rb + rs) = (ea - es)/rb = (es - ei)/rb
 ! where ea is the vapor pressure in the vegetation canopy, ei is the saturated vapor pressure at the leaf temperature,
 !  and then
 !     es = (ea*rs + ei*rb) / (rb + rs)
 ! more details are in Appendix C of Sellers et al. (J. Climate 1996) and Oleson et al. (NCAR Tech. Note, 2010)

 ! stomatal resistance is the larger of two roots in the quadratic equation defined in Oleson et al. (2010)

 ! -----------------------------------------------------------------------------------------------------------------
 subroutine quadResist(desireDeriv,ix_bbHumdFunc,rlb,fHum,gMin,g0,dg0_dc,rs,drs_dc)
 implicit none
 ! dummy variables
 logical(lgt),intent(in) :: desireDeriv   ! flag to denote if the derivative is desired
 integer(i4b),intent(in) :: ix_bbHumdFunc ! option for humidity control on stomatal resistance 
 real(dp),intent(in)     :: rlb           ! leaf boundary layer resistance (umol-1 m2 s)
 real(dp),intent(in)     :: fHum          ! scaled humidity function (-)
 real(dp),intent(in)     :: gMin          ! scaled minimum stomatal consuctance (umol m-2 s-1)
 real(dp),intent(in)     :: g0            ! stomatal conductance in the absence of humidity controls (umol m-2 s-1)
 real(dp),intent(in)     :: dg0_dc        ! derivative in g0 w.r.t intercellular co2 concentration (umol m-2 s-1 Pa-1)
 real(dp),intent(out)    :: rs            ! stomatal resistance ((umol-1 m2 s)
 real(dp),intent(out)    :: drs_dc        ! derivaive in rs w.r.t intercellular co2 concentration (umol-1 m2 s Pa-1)
 ! local variables
 real(dp)                :: aQuad,bQuad,cQuad ! coefficients in the quadratic function
 real(dp)                :: bSign,xTemp,qQuad ! q term in the quadratic
 real(dp)                :: root1,root2       ! roots of the quadratic
 real(dp)                :: dxT_dc,dqq_dc     ! derivatives in the q term

 ! define terms for the quadratic function
 select case(ix_bbHumdFunc)

  ! original Ball-Berry
  case(humidLeafSurface)
   aQuad = g0*fHum + gMin
   bQuad = (g0 + gMin)*rlb - 1._dp
   cQuad = -rlb

  ! Leuning 1995
  case(scaledHyperbolic)
   aQuad =  g0 + gMin*(1._dp + fHum)
   bQuad = (g0 + gMin)*rlb - fHum - 1._dp
   cQuad = -rlb

 endselect

 ! compute the q term in the quadratic
 bSign = abs(bQuad)/bQuad
 xTemp = bQuad*bQuad - 4._dp *aQuad*cQuad
 qquad = -0.5_dp * (bQuad + bSign*sqrt(xTemp))

 ! compute roots
 root1 = qQuad / aQuad
 root2 = cQuad / qQuad
 rs = max(root1,root2)

 ! compute derivatives
 if(desireDeriv)then

  ! compute derivatives in qquad w.r.t. ci
  select case(ix_bbHumdFunc)
   case(humidLeafSurface); dXt_dc = dg0_dc*(rlb*bQuad*2._dp - fHum*cQuad*4._dp)
   case(scaledHyperbolic); dXt_dc = dg0_dc*(rlb*bQuad*2._dp - cQuad*4._dp)
  end select
  dqq_dc = -0.5_dp * (rlb*dg0_dc + bSign*dXt_dc*0.5_dp / sqrt(xTemp) )

  ! compute derivatives in rs
  if(root1 > root2)then
   select case(ix_bbHumdFunc)
    case(humidLeafSurface); drs_dc = (dqq_dc - root1*fHum*dg0_dc)/aQuad 
    case(scaledHyperbolic); drs_dc = (dqq_dc - root1*dg0_dc)/aQuad
   end select
  else
   drs_dc = -root2*dqq_dc/qQuad
  endif

 ! derivatives not desired
 else
  drs_dc = 0._dp
 endif

 end subroutine quadResist

 ! *****
 ! * quadratic smoother...
 ! ***********************

 subroutine quadSmooth(desireDeriv, x1, x2, xsFac, dx1_dc, dx2_dc, xs, dxs_dc)
 implicit none
 ! dummy variables
 logical(lgt),intent(in) :: desireDeriv       ! flag to denote if a derivative is desired
 real(dp),intent(in)     :: x1,x2             ! variables to be smoothed
 real(dp),intent(in)     :: xsFac             ! smoothing factor
 real(dp),intent(in)     :: dx1_dc,dx2_dc     ! derivatives in variables w.r.t. something important
 real(dp),intent(out)    :: xs                ! smoothed variable
 real(dp),intent(out)    :: dxs_dc            ! derivative w.r.t. something important
 ! local variables
 real(dp)                :: aQuad,bQuad,cQuad ! coefficients in the quadratic function
 real(dp)                :: bSign,xTemp,qQuad ! q term in the quadratic
 real(dp)                :: root1,root2       ! roots of the quadratic
 real(dp)                :: dbq_dc,dcq_dc     ! derivatives in quadratic coefficients
 real(dp)                :: dxT_dc,dqq_dc     ! derivatives in the q term

 ! uses the quadratic of the form
 !  xsFac*xs^2 - (x1 + x2)*xs + x1*x2 = 0
 ! to smooth variables x1 and x2

 ! define the terms in the quadratic
 aQuad = xsFac
 bQuad = -(x1 + x2)
 cQuad = x1*x2

 ! compute the q term in the quadratic
 bSign = abs(bQuad)/bQuad
 xTemp = bQuad*bQuad - 4._dp *aQuad*cQuad
 qquad = -0.5_dp * (bQuad + bSign*sqrt(xTemp))

 ! compute roots
 root1 = qQuad / aQuad
 root2 = cQuad / qQuad
 xs    = min(root1,root2)

 ! compute derivatives
 if(desireDeriv)then

  ! compute derivatives for the terms in the quadratic
  dbq_dc = -(dx1_dc + dx2_dc)
  dcq_dc = x1*dx2_dc + x2*dx1_dc

  ! compute derivatives for xTemp
  dxT_dc = 2._dp*(bQuad*dbq_dc) - 4._dp*aQuad*dcq_dc
  dqq_dc = -0.5_dp * (dbq_dc + bsign*dxT_dc/(2._dp*sqrt(xTemp)))

  ! compute derivatives in the desired root
  if(root1 < root2)then
   dxs_dc = dqq_dc/aQuad
  else
   dxs_dc = (dcq_dc - root2*dqq_dc)/qQuad
  endif

 ! derivatives not required
 else
  dxs_dc = 0._dp
 endif

 end subroutine quadSmooth


 ! *****
 ! * temperature functions...
 ! **************************

 ! q10 function for temperature dependence
 function q10(a,T,Tmid,Tscale)
 implicit none
 real(dp),intent(in) :: a              ! scale factor
 real(dp),intent(in) :: T              ! temperature (K)
 real(dp),intent(in) :: Tmid           ! point where function is one (25 deg C)
 real(dp),intent(in) :: Tscale         ! scaling factor (K)
 real(dp)            :: q10            ! temperature dependence (-)
 q10 = a**((T - Tmid)/Tscale)
 end function q10

 ! Arrhenius function for temperature dependence
 function fT(delH,T,Tref)
 implicit none
 real(dp),intent(in) :: delH     ! activation energy in temperature function (J mol-1)
 real(dp),intent(in) :: T        ! temperature (K)
 real(dp),intent(in) :: Tref     ! reference temperature (K)
 real(dp)            :: fT       ! temperature dependence (-)
 fT = exp((delH/(Tref*Rgas))*(1._dp - Tref/T))  ! NOTE: Rgas = J K-1 mol-1
 end function fT

 ! function for high temperature inhibition
 function fHigh(delH,delS,T)
 implicit none
 real(dp),intent(in) :: delH     ! deactivation energy in high temp inhibition function (J mol-1)
 real(dp),intent(in) :: delS     ! entropy term in high temp inhibition function (J K-1 mol-1) 
 real(dp),intent(in) :: T        ! temperature (K)
 real(dp)            :: fHigh    ! high temperature inhibition (-)
 fHigh = 1._dp + exp( (delS*T - delH)/(Rgas*T) ) ! NOTE: Rgas = J K-1 mol-1
 end function fHigh


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


 ! *******************************************************************************************************
 ! private subroutine quadratic: solve the equation ax2 + bx + c = 0
 ! *******************************************************************************************************
 subroutine quadratic(a, b, c, r1, r2)
 implicit none
 real(dp)  :: a, b, c  ! terms in the quadratic equation
 real(dp)  :: r1, r2   ! roots
 real(dp)  :: q

 if (b >= 0._dp)then
  q = -0.5_dp * (b + sqrt(b*b - 4._dp*a*c))
 else
  q = -0.5_dp * (b - sqrt(b*b - 4._dp*a*c))
 endif

 r1 = q / a
 if (q /= 0._dp)then
  r2 = c / q
 else
  r2 = 1.e36_dp
 endif

 end subroutine quadratic


 
 ! -- end private subroutines
 ! ------------------------------------------------------------------------------------------------------------
 ! ------------------------------------------------------------------------------------------------------------
 ! ------------------------------------------------------------------------------------------------------------

end module stomResist_module

! SUMMA - Structure for Unifying Multiple Modeling Alternatives
! Copyright (C) 2014-2020 NCAR/RAL; University of Saskatchewan; University of Washington
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

module t2enthalpyAddPrime_module

! constants
USE multiconst, only: gravity, &                          ! gravitational acceleration (m s-1)
                      Tfreeze, &                          ! freezing point of water (K)
                      Cp_soil,Cp_water,Cp_ice,Cp_air,&    ! specific heat of soil, water and ice (J kg-1 K-1)
                      iden_water,iden_ice,iden_air,&      ! intrinsic density of water and ice (kg m-3)
                      LH_fus                              ! latent heat of fusion (J kg-1)

! data types
USE nrtype
USE data_types,only:var_iLength                    ! var(:)%dat(:)
USE data_types,only:var_dLength                    ! var(:)%dat(:)
USE data_types,only:zLookup                        ! z(:)%var(:)%lookup(:)

! indices within parameter structure
USE var_lookup,only:iLookPARAM                     ! named variables to define structure element
USE var_lookup,only:iLookINDEX                     ! named variables to define structure element
USE var_lookup,only:iLookLOOKUP                    ! named variables to define structure element
USE var_lookup,only:iLookDIAG       ! named variables for structure elements

! data dimensions
USE var_lookup,only:maxvarLookup                   ! maximum number of variables in the lookup tables

! domain types
USE globalData,only:iname_cas                      ! named variables for canopy air space
USE globalData,only:iname_veg                      ! named variables for vegetation canopy
USE globalData,only:iname_snow                     ! named variables for snow
USE globalData,only:iname_soil                     ! named variables for soil
USE globalData,only:iname_aquifer                  ! named variables for the aquifer

! named variables to describe the state variable type
USE globalData,only:iname_nrgCanair                ! named variable defining the energy of the canopy air space
USE globalData,only:iname_nrgCanopy                ! named variable defining the energy of the vegetation canopy
USE globalData,only:iname_nrgLayer                 ! named variable defining the energy state variable for snow+soil layers

! missing values
USE globalData,only:integerMissing                 ! missing integer
USE globalData,only:realMissing                    ! missing real number

! privacy
implicit none
private
public::t2enthalpyPrime

! define the look-up table used to compute temperature based on enthalpy
contains


! ************************************************************************************************************************
! public subroutine t2enthalpyPrime: compute enthalpy prime from temperature and total water content
! ************************************************************************************************************************
subroutine t2enthalpyPrime(&
                      ! input: data structures
                      diag_data,                         & ! intent(in):   model diagnostic variables for a local HRU
                      mpar_data,                         & ! intent(in):   parameter data structure
                      indx_data,                         & ! intent(in):   model indices
                      lookup_data,                       & ! intent(in):   lookup table data structure
                      ! input: state variables for the vegetation canopy
                      scalarCanairTempPrime,             & ! intent(in):   prime value of canopy air temperature (K)
                      scalarCanopyTempTrial,             & ! intent(in):   trial value of canopy temperature (K)
                      scalarCanopyWatTrial,              & ! intent(in):   trial value of canopy total water (kg m-2)
                      scalarCanopyTempPrime,             & ! intent(in):   prime value of canopy temperature (K)
                      scalarCanopyWatPrime,              & ! intent(in):   prime value of canopy total water (kg m-2)
                      ! input: variables for the snow-soil domain
                      mLayerTempTrial,                   & ! intent(in):   trial vector of layer temperature (K)
                      mLayerVolFracWatTrial,             & ! intent(in):   trial vector of volumetric total water content (-)
                      mLayerMatricHeadTrial,             & ! intent(in):   trial vector of total water matric potential (m)
                      mLayerTempPrime,                   & ! intent(in):   prime vector of layer temperature (K)
                      mLayerVolFracWatPrime,             & ! intent(in):   prime vector of volumetric total water content (-)
                      mLayerMatricHeadPrime,             & ! intent(in):   prime vector of total water matric potential (m)
                      ! input: pre-computed derivatives
                      dVolTot_dPsi0,                     & ! intent(in):   derivative in total water content w.r.t. total water matric potential (m-1)
                      d2VolTot_dPsi02,                   & ! intent(in):   second derivative in total water content w.r.t. total water matric potential (m-2)
                      ! output: enthalpy prime and derivatives
                      scalarCanairEnthalpyPrime,         & ! intent(out):  prime enthalpy of the canopy air space (J m-3)
                      scalarCanopyEnthalpyPrime,         & ! intent(out):  prime enthalpy of the vegetation canopy (J m-3)
                      mLayerEnthalpyPrime,               & ! intent(out):  prime enthalpy of each snow+soil layer (J m-3)
                      dCanEnthalpyPrime_dTk,             & ! intent(out):  derivatives in prime canopy enthalpy w.r.t. temperature
                      dCanEnthalpyPrime_dWat,            & ! intent(out):  derivatives in prime canopy enthalpy w.r.t. water state
                      dEnthalpyPrime_dTk,                & ! intent(out):  derivatives in prime layer enthalpy w.r.t. temperature
                      dEnthalpyPrime_dWat,               & ! intent(out):  derivatives in prime layer enthalpy w.r.t. water state
                      dCanEnthalpyPrime_dTkPrime,        & ! intent(out):  derivatives in prime canopy enthalpy w.r.t. prime temperature
                      dCanEnthalpyPrime_dWatPrime,       & ! intent(out):  derivatives in prime canopy enthalpy w.r.t. prime water state
                      dEnthalpyPrime_dTkPrime,           & ! intent(out):  derivatives in prime layer enthalpy w.r.t. prime temperature
                      dEnthalpyPrime_dWatPrime,          & ! intent(out):  derivatives in prime layer enthalpy w.r.t. prime water state
                      ! output: error control
                      err,message)                         ! intent(out): error control
  ! -------------------------------------------------------------------------------------------------------------------------
  ! downwind routines
  USE soil_utils_module,only:crit_soilT     ! compute critical temperature below which ice exists
  USE soil_utils_module,only:volFracLiq     ! compute volumetric fraction of liquid water
  USE spline_int_module,only:splint         ! use for cubic spline interpolation
  implicit none
  ! delare dummy variables
  ! -------------------------------------------------------------------------------------------------------------------------
  ! input: data structures
  type(var_dlength),intent(in)     :: diag_data                    ! diagnostic variables for a local HRU
  type(var_dlength),intent(in)     :: mpar_data                    ! model parameters
  type(var_ilength),intent(in)     :: indx_data                    ! model indices
  type(zLookup),intent(in)         :: lookup_data                  ! lookup tables
  ! input: state variables for the vegetation canopy
  real(rkind),intent(in)           :: scalarCanairTempPrime        ! prime value of canopy air temperature (K)
  real(rkind),intent(in)           :: scalarCanopyTempTrial        ! trial value of canopy temperature (K)
  real(rkind),intent(in)           :: scalarCanopyWatTrial         ! trial value of canopy total water (kg m-2)
  real(rkind),intent(in)           :: scalarCanopyTempPrime        ! prime value of canopy temperature (K)
  real(rkind),intent(in)           :: scalarCanopyWatPrime         ! prime value of canopy total water (kg m-2)
  ! input: variables for the snow-soil domain
  real(rkind),intent(in)           :: mLayerTempTrial(:)           ! trial vector of layer temperature (K)
  real(rkind),intent(in)           :: mLayerVolFracWatTrial(:)     ! trial vector of volumetric total water content (-)
  real(rkind),intent(in)           :: mLayerMatricHeadTrial(:)     ! trial vector of total water matric potential (m)
  real(rkind),intent(in)           :: mLayerTempPrime(:)           ! prime vector of layer temperature (K)
  real(rkind),intent(in)           :: mLayerVolFracWatPrime(:)     ! prime vector of volumetric total water content (-)
  real(rkind),intent(in)           :: mLayerMatricHeadPrime(:)     ! prime vector of total water matric potential (m)
  ! input: pre-computed derivatives
  real(rkind),intent(in)           :: dVolTot_dPsi0(:)             ! derivative in total water content w.r.t. total water matric potential (m-1)
  real(rkind),intent(in)           :: d2VolTot_dPsi02(:)           ! second derivative in total water content w.r.t. total water matric potential (m-2)
  ! output: enthalpy prime
  real(rkind),intent(out)          :: scalarCanairEnthalpyPrime    ! prime enthalpy of the canopy air space (J m-3)
  real(rkind),intent(out)          :: scalarCanopyEnthalpyPrime    ! prime enthalpy of the vegetation canopy (J m-3)
  real(rkind),intent(out)          :: mLayerEnthalpyPrime(:)       ! prime enthalpy of each snow+soil layer (J m-3)
  ! output: derivatives
  real(rkind),intent(out)          :: dCanEnthalpyPrime_dTk        ! derivatives in prime canopy enthalpy w.r.t. temperature
  real(rkind),intent(out)          :: dCanEnthalpyPrime_dWat       ! derivatives in prime canopy enthalpy w.r.t. water state
  real(rkind),intent(out)          :: dEnthalpyPrime_dTk(:)        ! derivatives in prime layer enthalpy w.r.t. temperature
  real(rkind),intent(out)          :: dEnthalpyPrime_dWat(:)       ! derivatives in prime layer enthalpy w.r.t. water state  
  real(rkind),intent(out)          :: dCanEnthalpyPrime_dTkPrime   ! derivatives in prime canopy enthalpy w.r.t. priem temperature
  real(rkind),intent(out)          :: dCanEnthalpyPrime_dWatPrime  ! derivatives in prime canopy enthalpy w.r.t. prime water state
  real(rkind),intent(out)          :: dEnthalpyPrime_dTkPrime(:)   ! derivatives in prime layer enthalpy w.r.t. prime temperature
  real(rkind),intent(out)          :: dEnthalpyPrime_dWatPrime(:)  ! derivatives in prime layer enthalpy w.r.t. prime water state

  ! output: error control
  integer(i4b),intent(out)         :: err                          ! error code
  character(*),intent(out)         :: message                      ! error message
  ! -------------------------------------------------------------------------------------------------------------------------
  ! declare local variables
  character(len=128)               :: cmessage                  ! error message in downwind routine
  integer(i4b)                     :: iState                    ! index of model state variable
  integer(i4b)                     :: iLayer                    ! index of model layer
  integer(i4b)                     :: ixFullVector              ! index within full state vector
  integer(i4b)                     :: ixDomainType              ! name of a given model domain
  integer(i4b)                     :: ixControlIndex            ! index within a given model domain
  real(rkind)                      :: vGn_m                     ! van Genuchten "m" parameter (-)
  real(rkind)                      :: Tcrit                     ! temperature where all water is unfrozen (K)
  real(rkind)                      :: volFracWat                ! volumetric fraction of total water, liquid+ice (-)
  real(rkind)                      :: diffT                     ! temperature difference from Tfreeze
  real(rkind)                      :: integral                  ! integral of snow freezing curve
  real(rkind)                      :: dTcrit_dPsi0              ! derivative of temperature where all water is unfrozen (K) with matric head
  real(rkind)                      :: d_integral_dTk            ! derivative of integral with temperature
  real(rkind)                      :: d2_integral_dTk2          ! second derivative of integral with temperature
  real(rkind)                      :: dE                        ! derivative of enthalpy with temperature at layer temperature
  real(rkind)                      :: dEcrit                    ! derivative of enthalpy with temperature at critical temperature
  ! enthalpy
  real(rkind)                      :: enthVegP                   ! prime enthalpy of the vegetation (J m-3)
  real(rkind)                      :: enthSoilP                  ! prime enthalpy of soil particles (J m-3)
  real(rkind)                      :: enthMixP                   ! prime enthalpy of the mixed region, liquid+ice (J m-3)
  real(rkind)                      :: enthLiqP                   ! prime enthalpy of the liquid region (J m-3)
  real(rkind)                      :: enthIceP                   ! prime enthalpy of the ice region (J m-3)
  real(rkind)                      :: enthAirP                   ! prime enthalpy of air (J m-3)
  real(rkind)                      :: enthTemp                   ! enthalpy at the temperature of the control volume (J m-3)
  real(rkind)                      :: enthTcrit                  ! enthalpy at the critical temperature where all water is unfrozen (J m-3)
  real(rkind)                      :: enthPhaseP                 ! prime enthalpy associated with phase change (J m-3)
  real(rkind)                      :: enthWaterP                 ! prime enthalpy of total water (J m-3)
  ! enthalpy derivatives
  real(rkind)                      :: dEnthVegP_dTk              ! derivative of prime enthalpy of the vegetation with temperature
  real(rkind)                      :: dEnthSoilP_dTk             ! derivative of prime enthalpy of the soil with temperature
  real(rkind)                      :: dEnthLiqP_dTk              ! derivative of prime enthalpy of the liquid with temperature
  real(rkind)                      :: dEnthIceP_dTk              ! derivative of prime enthalpy of the ice with temperature
  real(rkind)                      :: dEnthAirP_dTk              ! derivative of prime enthalpy of the air with temperature
  real(rkind)                      :: dEnthWaterP_dTk            ! derivative of prime enthalpy of the total water with temperature
  real(rkind)                      :: dEnthVegP_dWat             ! derivative of prime enthalpy of the vegetation with water state
  real(rkind)                      :: dEnthSoilP_dWat            ! derivative of prime enthalpy of the soil with water state
  real(rkind)                      :: dEnthLiqP_dWat             ! derivative of prime enthalpy of the liquid with water state
  real(rkind)                      :: dEnthIceP_dWat             ! derivative of prime enthalpy of the ice with water state
  real(rkind)                      :: dEnthAirP_dWat             ! derivative of prime enthalpy of the air with water state
  real(rkind)                      :: dEnthWaterP_dWat           ! derivative of prime enthalpy of the total water with water state
  real(rkind)                      :: dEnthVegP_dTkP             ! derivative of prime enthalpy of the vegetation with prime temperature
  real(rkind)                      :: dEnthSoilP_dTkP            ! derivative of prime enthalpy of the soil with prime temperature
  real(rkind)                      :: dEnthLiqP_dTkP             ! derivative of prime enthalpy of the liquid with prime temperature
  real(rkind)                      :: dEnthIceP_dTkP             ! derivative of prime enthalpy of the ice with prime temperature
  real(rkind)                      :: dEnthAirP_dTkP             ! derivative of prime enthalpy of the air with prime temperature
  real(rkind)                      :: dEnthWaterP_dTkP           ! derivative of prime enthalpy of the total water with prime temperature
  real(rkind)                      :: dEnthVegP_dWatP            ! derivative of prime enthalpy of the vegetation with prime water state
  real(rkind)                      :: dEnthSoilP_dWatP           ! derivative of prime enthalpy of the soil with prime water state
  real(rkind)                      :: dEnthLiqP_dWatP            ! derivative of prime enthalpy of the liquid with prime water state
  real(rkind)                      :: dEnthIceP_dWatP            ! derivative of prime enthalpy of the ice with prime water state
  real(rkind)                      :: dEnthAirP_dWatP            ! derivative of prime enthalpy of the air with prime water state
  real(rkind)                      :: dEnthWaterP_dWatP          ! derivative of prime enthalpy of the total water with prime water state

  ! --------------------------------------------------------------------------------------------------------------------------------
  ! make association with variables in the data structures
  generalVars: associate(&
    ! number of model layers, and layer type
    nSnow                   => indx_data%var(iLookINDEX%nSnow)%dat(1)                 ,& ! intent(in):  [i4b]    total number of snow layers
    nSoil                   => indx_data%var(iLookINDEX%nSoil)%dat(1)                 ,& ! intent(in):  [i4b]    total number of soil layers
    nLayers                 => indx_data%var(iLookINDEX%nLayers)%dat(1)               ,& ! intent(in):  [i4b]    total number of snow and soil layers
    ! mapping between the full state vector and the state subset
    ixMapFull2Subset        => indx_data%var(iLookINDEX%ixMapFull2Subset)%dat         ,& ! intent(in):  [i4b(:)] list of indices in the state subset for each state in the full state vector
    ixMapSubset2Full        => indx_data%var(iLookINDEX%ixMapSubset2Full)%dat         ,& ! intent(in):  [i4b(:)] [state subset] list of indices of the full state vector in the state subset
    ! type of domain, type of state variable, and index of control volume within domain
    ixDomainType_subset     => indx_data%var(iLookINDEX%ixDomainType_subset)%dat      ,& ! intent(in):  [i4b(:)] [state subset] id of domain for desired model state variables
    ixControlVolume         => indx_data%var(iLookINDEX%ixControlVolume)%dat          ,& ! intent(in):  [i4b(:)] index of the control volume for different domains (veg, snow, soil)
    ixStateType             => indx_data%var(iLookINDEX%ixStateType)%dat               & ! intent(in):  [i4b(:)] indices defining the type of the state (iname_nrgLayer...)
    ) ! end associate statement
    ! --------------------------------------------------------------------------------------------------------------------------------

    ! initialize error control
    err=0; message="t2enthalpyPrime/"

    ! loop through model state variables
    do iState=1,size(ixMapSubset2Full)

      ! -----
      ! - compute indices...
      ! --------------------

      ! get domain type, and index of the control volume within the domain
      ixFullVector   = ixMapSubset2Full(iState)       ! index within full state vector
      ixDomainType   = ixDomainType_subset(iState)    ! named variables defining the domain (iname_cas, iname_veg, etc.)
      ixControlIndex = ixControlVolume(ixFullVector)  ! index within a given domain

      ! check an energy state
      if(ixStateType(ixFullVector)==iname_nrgCanair .or. ixStateType(ixFullVector)==iname_nrgCanopy .or. ixStateType(ixFullVector)==iname_nrgLayer)then

        ! get the layer index
        select case(ixDomainType)
          case(iname_cas);     iLayer = integerMissing
          case(iname_veg);     iLayer = integerMissing
          case(iname_snow);    iLayer = ixControlIndex
          case(iname_soil);    iLayer = ixControlIndex + nSnow
          case(iname_aquifer); cycle ! aquifer: do nothing (no thermodynamics in the aquifer)
          case default; err=20; message=trim(message)//'expect case to be iname_cas, iname_veg, iname_snow, iname_soil, iname_aquifer'; return
        end select

        ! identify domain
        select case(ixDomainType)
          case(iname_cas)
            scalarCanairEnthalpyPrime = Cp_air * iden_air * scalarCanairTempPrime

          case(iname_veg)
            ! association to necessary variables for vegetation
            vegVars: associate(&
              canopyDepth             => diag_data%var(iLookDIAG%scalarCanopyDepth)%dat(1),         & ! intent(in): [dp]      canopy depth (m)
              specificHeatVeg         => mpar_data%var(iLookPARAM%specificHeatVeg)%dat(1),          & ! intent(in): specific heat of vegetation (J kg-1 K-1)
              maxMassVegetation       => mpar_data%var(iLookPARAM%maxMassVegetation)%dat(1),        & ! intent(in): maximum mass of vegetation (kg m-2)
              snowfrz_scale           => mpar_data%var(iLookPARAM%snowfrz_scale)%dat(1)   & ! intent(in):  [dp] scaling parameter for the snow freezing curve (K-1)
              )

              diffT = scalarCanopyTempTrial - Tfreeze
              enthVegP = specificHeatVeg * maxMassVegetation * scalarCanopyTempPrime / canopyDepth
              ! enthalpy prime derivatives
              dEnthVegP_dTk   = specificHeatVeg * maxMassVegetation / canopyDepth
              dEnthVegP_dWat  = 0._rkind
              dEnthVegP_dTkP  = 0._rkind
              dEnthVegP_dWatP = 0._rkind

              if(diffT>=0._rkind)then
                enthLiqP = Cp_water * ( scalarCanopyWatPrime * diffT + scalarCanopyWatTrial * scalarCanopyTempPrime )/ canopyDepth
                enthIceP = 0._rkind
                ! enthalpy prime derivatives
                dEnthLiqP_dTk   = Cp_water * scalarCanopyWatPrime / canopyDepth
                dEnthLiqP_dWat  = Cp_water * scalarCanopyTempPrime / canopyDepth
                dEnthLiqP_dTkP  = Cp_water * scalarCanopyWatTrial / canopyDepth
                dEnthLiqP_dWatP = Cp_water * diffT / canopyDepth
                dEnthIceP_dTk   = 0._rkind
                dEnthIceP_dWat  = 0._rkind
                dEnthIceP_dTkP  = 0._rkind
                dEnthIceP_dWatP = 0._rkind
              else
                integral = (1._rkind/snowfrz_scale) * atan(snowfrz_scale * diffT)
                d_integral_dTk = 1._rkind / (1._rkind + (snowfrz_scale * diffT)**2_i4b) ! Note: VolFracLiq = d_integral_dTk*VolFracWat
                enthLiqP = Cp_water * ( scalarCanopyWatPrime*integral + scalarCanopyWatTrial*scalarCanopyTempPrime*d_integral_dTk )/ canopyDepth
                enthIceP = Cp_ice * ( scalarCanopyWatPrime*( diffT - integral ) + scalarCanopyWatTrial*scalarCanopyTempPrime*( 1._rkind - d_integral_dTk ) ) / canopyDepth
                ! derivatives
                d2_integral_dTk2 = -2._rkind*snowfrz_scale*diffT / (1._rkind + (snowfrz_scale * diffT)**2_i4b)**2_i4b
                ! enthalpy prime derivatives
                dEnthLiqP_dTk   = Cp_water * ( scalarCanopyWatPrime*d_integral_dTk + scalarCanopyWatTrial*scalarCanopyTempPrime*d2_integral_dTk2 ) / canopyDepth
                dEnthLiqP_dWat  = Cp_water * scalarCanopyTempPrime * d_integral_dTk / canopyDepth
                dEnthLiqP_dTkP  = Cp_water *  scalarCanopyWatTrial * d_integral_dTk / canopyDepth
                dEnthLiqP_dWatP = Cp_water * integral / canopyDepth
                dEnthIceP_dTk   = Cp_ice * ( scalarCanopyWatPrime*( 1._rkind - d_integral_dTk ) - scalarCanopyWatTrial*scalarCanopyTempPrime*d2_integral_dTk2 ) / canopyDepth
                dEnthIceP_dWat  = Cp_ice * ( scalarCanopyTempPrime*( 1._rkind - d_integral_dTk ) ) / canopyDepth
                dEnthIceP_dTkP  = Cp_ice * scalarCanopyWatTrial*( 1._rkind - d_integral_dTk ) / canopyDepth
                dEnthIceP_dWatP = Cp_ice * ( diffT - integral ) / canopyDepth
              endif

              scalarCanopyEnthalpyPrime = enthVegP + enthLiqP + enthIceP
              dCanEnthalpyPrime_dTk       = dEnthVegP_dTk + dEnthLiqP_dTk + dEnthIceP_dTk
              dCanEnthalpyPrime_dWat      = dEnthVegP_dWat + dEnthLiqP_dWat + dEnthIceP_dWat
              dCanEnthalpyPrime_dTkPrime  = dEnthVegP_dTkP + dEnthLiqP_dTkP + dEnthIceP_dTkP
              dCanEnthalpyPrime_dWatPrime = dEnthVegP_dWatP + dEnthLiqP_dWatP + dEnthIceP_dWatP

            end associate vegVars

          case(iname_snow)

            ! association to necessary variables for snow
            snowVars: associate(&
              snowfrz_scale           => mpar_data%var(iLookPARAM%snowfrz_scale)%dat(1)   & ! intent(in):  [dp] scaling parameter for the snow freezing curve (K-1)
              )

              diffT = mLayerTempTrial(iLayer) - Tfreeze  ! diffT<0._rkind because snow is frozen
              integral = (1._rkind/snowfrz_scale) * atan(snowfrz_scale * diffT)
              d_integral_dTk = 1._rkind / (1._rkind + (snowfrz_scale * diffT)**2_i4b) ! Note: VolFracLiq = d_integral_dTk*VolFracWat
              enthLiqP = iden_water * Cp_water * ( mLayerVolFracWatPrime(iLayer)*integral + mLayerVolFracWatTrial(iLayer)*mLayerTempPrime(iLayer)*d_integral_dTk )
              enthIceP = iden_water * Cp_ice * ( mLayerVolFracWatPrime(iLayer)*( diffT - integral ) + mLayerVolFracWatTrial(iLayer)*mLayerTempPrime(iLayer)*( 1._rkind - d_integral_dTk ) )
              enthAirP = iden_air * Cp_air * ( mLayerTempPrime(iLayer)  - mLayerVolFracWatPrime(iLayer) * ( (iden_water/iden_ice)*( diffT - integral ) + integral )  &
                                              - mLayerVolFracWatTrial(iLayer)*mLayerTempPrime(iLayer)* ( (iden_water/iden_ice)*( 1._rkind - d_integral_dTk ) + d_integral_dTk ) )
              ! derivatives
              d2_integral_dTk2 = -2._rkind*snowfrz_scale*diffT / (1._rkind + (snowfrz_scale * diffT)**2_i4b)**2_i4b
              ! enthalpy prime derivatives
              dEnthLiqP_dTk   = iden_water * Cp_water * ( mLayerVolFracWatPrime(iLayer)*d_integral_dTk + mLayerVolFracWatTrial(iLayer)*mLayerTempPrime(iLayer)*d2_integral_dTk2 )
              dEnthLiqP_dWat  = iden_water * Cp_water * mLayerTempPrime(iLayer) * d_integral_dTk
              dEnthLiqP_dTkP  = iden_water * Cp_water * mLayerVolFracWatTrial(iLayer) * d_integral_dTk
              dEnthLiqP_dWatP = iden_water * Cp_water * integral
              dEnthIceP_dTk   = iden_water * Cp_ice * ( mLayerVolFracWatPrime(iLayer)*( 1._rkind - d_integral_dTk ) - mLayerVolFracWatTrial(iLayer)*mLayerTempPrime(iLayer)*d2_integral_dTk2 )
              dEnthIceP_dWat  = iden_water * Cp_ice * ( mLayerTempPrime(iLayer)*( 1._rkind - d_integral_dTk ) )
              dEnthIceP_dTkP  = iden_water * Cp_ice * mLayerVolFracWatTrial(iLayer)*( 1._rkind - d_integral_dTk )
              dEnthIceP_dWatP = iden_water * Cp_ice * ( diffT - integral )
              dEnthAirP_dTk   =  iden_air  * Cp_air * ( -mLayerVolFracWatPrime(iLayer) * ( (iden_water/iden_ice)*( 1._rkind - d_integral_dTk ) + d_integral_dTk )  &
                                                       - mLayerVolFracWatTrial(iLayer)*mLayerTempPrime(iLayer) * (-(iden_water/iden_ice)*d2_integral_dTk2 + d2_integral_dTk2 ) )
              dEnthAirP_dWat  = -iden_air  * Cp_air * ( mLayerTempPrime(iLayer) * (iden_water/iden_ice)*( 1._rkind - d_integral_dTk ) + d_integral_dTk )
              dEnthAirP_dTkP  =  iden_air  * Cp_air * ( 1._rkind - mLayerVolFracWatTrial(iLayer)* ( (iden_water/iden_ice)*( 1._rkind - d_integral_dTk ) + d_integral_dTk ) )
              dEnthAirP_dWatP =- iden_air  * Cp_air * ( (iden_water/iden_ice)*( diffT - integral ) + integral )

              mLayerEnthalpyPrime(iLayer)      = enthLiqP + enthIceP + enthAirP
              dEnthalpyPrime_dTk(iLayer)       = dEnthLiqP_dTk + dEnthIceP_dTk + dEnthAirP_dTk
              dEnthalpyPrime_dWat(iLayer)      = dEnthLiqP_dWat + dEnthIceP_dWat + dEnthAirP_dWat
              dEnthalpyPrime_dTkPrime(iLayer)  = dEnthLiqP_dTkP + dEnthIceP_dTkP + dEnthAirP_dTkP
              dEnthalpyPrime_dWatPrime(iLayer) = dEnthLiqP_dWatP + dEnthIceP_dWatP + dEnthAirP_dWatP

            end associate snowVars

          case(iname_soil)

            ! make association to variables in the data structures...
            soilVars: associate(&

              ! associate model parameters
              soil_dens_intr => mpar_data%var(iLookPARAM%soil_dens_intr)%dat(ixControlIndex)      , & ! intrinsic soil density             (kg m-3)
              theta_sat      => mpar_data%var(iLookPARAM%theta_sat)%dat(ixControlIndex)           , & ! soil porosity                      (-)
              theta_res      => mpar_data%var(iLookPARAM%theta_res)%dat(ixControlIndex)           , & ! volumetric residual water content  (-)
              vGn_alpha      => mpar_data%var(iLookPARAM%vGn_alpha)%dat(ixControlIndex)           , & ! van Genuchten "alpha" parameter    (m-1)
              vGn_n          => mpar_data%var(iLookPARAM%vGn_n)%dat(ixControlIndex)               , & ! van Genuchten "n" parameter        (-)

              ! associate values in the lookup table
              Tk            => lookup_data%z(ixControlIndex)%var(iLookLOOKUP%temperature)%lookup  , & ! temperature (K)
              Ey            => lookup_data%z(ixControlIndex)%var(iLookLOOKUP%enthalpy)%lookup     , & ! enthalpy (J m-3)
              E2            => lookup_data%z(ixControlIndex)%var(iLookLOOKUP%deriv2)%lookup         & ! second derivative of the interpolating function

              ) ! end associate statement

              ! diagnostic variables
              vGn_m    = 1._rkind - 1._rkind/vGn_n
              Tcrit    = crit_soilT( mLayerMatricHeadTrial(ixControlIndex) )
              volFracWat = volFracLiq(mLayerMatricHeadTrial(ixControlIndex),vGn_alpha,theta_res,theta_sat,vGn_n,vGn_m)
              diffT = mLayerTempTrial(iLayer) - Tfreeze
              dTcrit_dPsi0 = 0._rkind
              if (mLayerMatricHeadTrial(ixControlIndex)<0._rkind) dTcrit_dPsi0 = gravity * Tfreeze/LH_fus

              ! *** compute enthalpy prime of water for unfrozen conditions
              if(mlayerTempTrial(iLayer)>=Tcrit)then
                enthWaterP = iden_water * Cp_water * ( mLayerMatricHeadPrime(ixControlIndex)*dVolTot_dPsi0(ixControlIndex)*diffT + volFracWat*mLayerTempPrime(iLayer) )  ! valid for temperatures below freezing also
                ! enthalpy derivatives
                dEnthWaterP_dTk   = iden_water * Cp_water *  mLayerMatricHeadPrime(ixControlIndex) * dVolTot_dPsi0(ixControlIndex) ! dVolFracWat_dWat = dVolTot_dPsi0(ixControlIndex)
                dEnthWaterP_dWat  = iden_water * Cp_water * ( mLayerMatricHeadPrime(ixControlIndex)*d2VolTot_dPsi02(ixControlIndex)*diffT &
                                                             + dVolTot_dPsi0(ixControlIndex)*mLayerTempPrime(iLayer) )
                dEnthWaterP_dTkP  = iden_water * Cp_water * volFracWat
                dEnthWaterP_dWatP = iden_water * Cp_water * dVolTot_dPsi0(ixControlIndex)*diffT

              ! *** compute enthalpy prime of water for frozen conditions
              else
                ! calculate enthalpy prime at the temperature (cubic spline interpolation)
                call splint(Tk,Ey,E2,mlayerTempTrial(iLayer),enthTemp,dE,err,cmessage)
                if(err/=0) then; message=trim(message)//trim(cmessage); return; end if
                
                ! calculate enthalpy prime at the critical temperature (cubic spline interpolation)
                call splint(Tk,Ey,E2,Tcrit,enthTcrit,dEcrit,err,cmessage)
                if(err/=0) then; message=trim(message)//trim(cmessage); return; end if

                ! calculate the enthalpy prime of water
                enthMixP   = mlayerTempPrime(iLayer)*dE - mLayerMatricHeadPrime(ixControlIndex)*dEcrit ! enthalpy prime of the liquid+ice mix
                enthLiqP   = iden_water * Cp_water * mLayerMatricHeadPrime(ixControlIndex) * ( dVolTot_dPsi0(ixControlIndex) * (Tcrit - Tfreeze) + volFracWat * dTcrit_dPsi0 )
                enthWaterP = enthMixP + enthLiqP
                ! enthalpy derivatives
                dEnthWaterP_dTk = 0._rkind ! mlayerTempPrime(iLayer)*dE_dTk - mLayerMatricHeadPrime(ixControlIndex)*dEcrit_dTk, but derivatives of dE and dEcrit are zero
                dEnthWaterP_dWat =  0._rkind ! -mLayerMatricHeadPrime(ixControlIndex)*dEcrit_dTk*dTcrit_dPsi0, but derivatives of dEcrit are zero
                dEnthWaterP_dTkP = dE
                dEnthWaterP_dWatP =  -dEcrit 
              endif ! (if frozen conditions)

              ! *** compute the enthalpy of soil
              enthSoilP = soil_dens_intr * Cp_soil * (1._rkind - theta_sat)*mlayerTempPrime(iLayer)
              ! enthalpy derivatives
              dEnthSoilP_dTk   = 0._rkind
              dEnthSoilP_dWat  = 0._rkind
              dEnthSoilP_dTkP  = soil_dens_intr * Cp_soil * (1._rkind - theta_sat)
              dEnthSoilP_dWatP = 0._rkind

              ! *** compute the enthalpy of air
              enthAirP = iden_air * Cp_air * ( mLayerMatricHeadPrime(ixControlIndex)*dVolTot_dPsi0(ixControlIndex)*diffT + (1._rkind - theta_sat - volFracWat)*mlayerTempPrime(iLayer) )
              ! enthalpy derivatives
              dEnthAirP_dTk   = iden_air * Cp_air * mLayerMatricHeadPrime(ixControlIndex) * dVolTot_dPsi0(ixControlIndex)
              dEnthAirP_dWat  = iden_air * Cp_air * ( mLayerMatricHeadPrime(ixControlIndex)*d2VolTot_dPsi02(ixControlIndex)*diffT - dVolTot_dPsi0(ixControlIndex)*mlayerTempPrime(iLayer) )
              dEnthAirP_dTkP  = iden_air * Cp_air * (1._rkind - theta_sat - volFracWat)
              dEnthAirP_dWatP = iden_air * Cp_air * dVolTot_dPsi0(ixControlIndex) * diffT 

              ! *** compute the total enthalpy (J m-3)
              mLayerEnthalpyPrime(iLayer)      = enthWaterP + enthSoilP + enthAirP
              dEnthalpyPrime_dTk(iLayer)       = dEnthWaterP_dTk + dEnthSoilP_dTk + dEnthAirP_dTk
              dEnthalpyPrime_dWat(iLayer)      = dEnthWaterP_dWat + dEnthSoilP_dWat + dEnthAirP_dWat
              dEnthalpyPrime_dTkPrime(iLayer)  = dEnthWaterP_dTkP + dEnthSoilP_dTkP + dEnthAirP_dTkP
              dEnthalpyPrime_dWatPrime(iLayer) = dEnthWaterP_dWatP + dEnthSoilP_dWatP + dEnthAirP_dWatP

            end associate soilVars

          ! -----
          ! - checks...
          ! -----------
          case(iname_aquifer); cycle ! aquifer: do nothing
          case default; err=20; message=trim(message)//'expect case to be iname_cas, iname_veg, iname_snow, iname_soil, iname_aquifer'; return
        end select

      end if  ! if an energy layer
    end do  ! looping through state variables

  end associate generalVars

end subroutine t2enthalpyPrime

end module t2enthalpyAddPrime_module

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

module enthalpyTempAddPrime_module

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

! indices within parameter structure
USE var_lookup,only:iLookPARAM                     ! named variables to define structure element
USE var_lookup,only:iLookINDEX                     ! named variables to define structure element
USE var_lookup,only:iLookLOOKUP                    ! named variables to define structure element
USE var_lookup,only:iLookDIAG                      ! named variables for structure elements

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
public::t2enthalpyPrime

contains


! ************************************************************************************************************************
! public subroutine t2enthalpyPrime: compute enthalpy prime from temperature and total water content
! ************************************************************************************************************************
subroutine t2enthalpyPrime(&
                      ! input: data structures
                      diag_data,                         & ! intent(in):   model diagnostic variables for a local HRU
                      mpar_data,                         & ! intent(in):   parameter data structure
                      indx_data,                         & ! intent(in):   model indices
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
                      ! output: enthalpy prime
                      scalarCanairEnthalpyPrime,         & ! intent(out):  prime enthalpy of the canopy air space (J m-3)
                      scalarCanopyEnthalpyPrime,         & ! intent(out):  prime enthalpy of the vegetation canopy (J m-3)
                      mLayerEnthalpyPrime,               & ! intent(out):  prime enthalpy of each snow+soil layer (J m-3)
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
  type(var_dlength),intent(in)     :: diag_data                 ! diagnostic variables for a local HRU
  type(var_dlength),intent(in)     :: mpar_data                 ! model parameters
  type(var_ilength),intent(in)     :: indx_data                 ! model indices
  ! input: state variables for the vegetation canopy
  real(rkind),intent(in)           :: scalarCanairTempPrime     ! prime value of canopy air temperature (K)
  real(rkind),intent(in)           :: scalarCanopyTempTrial     ! trial value of canopy temperature (K)
  real(rkind),intent(in)           :: scalarCanopyWatTrial      ! trial value of canopy total water (kg m-2)
  real(rkind),intent(in)           :: scalarCanopyTempPrime     ! prime value of canopy temperature (K)
  real(rkind),intent(in)           :: scalarCanopyWatPrime      ! prime value of canopy total water (kg m-2)
  ! input: variables for the snow-soil domain
  real(rkind),intent(in)           :: mLayerTempTrial(:)        ! trial vector of layer temperature (K)
  real(rkind),intent(in)           :: mLayerVolFracWatTrial(:)  ! trial vector of volumetric total water content (-)
  real(rkind),intent(in)           :: mLayerMatricHeadTrial(:)  ! trial vector of total water matric potential (m)
  real(rkind),intent(in)           :: mLayerTempPrime(:)        ! prime vector of layer temperature (K)
  real(rkind),intent(in)           :: mLayerVolFracWatPrime(:)  ! prime vector of volumetric total water content (-)
  real(rkind),intent(in)           :: mLayerMatricHeadPrime(:)  ! prime vector of total water matric potential (m)
  ! input: pre-computed derivatives
  real(rkind),intent(in)           :: dVolTot_dPsi0(:)          ! derivative in total water content w.r.t. total water matric potential (m-1)
   ! output: enthalpy prime
  real(rkind),intent(out)          :: scalarCanairEnthalpyPrime ! prime enthalpy of the canopy air space (J m-3)
  real(rkind),intent(out)          :: scalarCanopyEnthalpyPrime ! prime enthalpy of the vegetation canopy (J m-3)
  real(rkind),intent(out)          :: mLayerEnthalpyPrime(:)    ! prime enthalpy of each snow+soil layer (J m-3)
  ! output: error control
  integer(i4b),intent(out)         :: err                       ! error code
  character(*),intent(out)         :: message                   ! error message
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
  real(rkind)                      :: d_integral_psiLiq_dTk     ! derivative with temperature of integral of soil mLayerPsiLiq from Tfreeze to layer temperature
  real(rkind)                      :: xConst                    ! constant in the freezing curve function (m K-1)
  real(rkind)                      :: mLayerPsiLiq              ! liquid water matric potential (m)
  ! enthalpy
  real(rkind)                      :: enthVegP                  ! prime enthalpy of the vegetation (J m-3)
  real(rkind)                      :: enthSoilP                 ! prime enthalpy of soil particles (J m-3)
  real(rkind)                      :: enthLiqP                  ! prime enthalpy of the liquid region (J m-3)
  real(rkind)                      :: enthIceP                  ! prime enthalpy of the ice region (J m-3)
  real(rkind)                      :: enthAirP                  ! prime enthalpy of air (J m-3)
  real(rkind)                      :: enthPhaseP                ! prime enthalpy associated with phase change (J m-3)
  real(rkind)                      :: enthWaterP                ! prime enthalpy of total water (J m-3)
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
              canopyDepth       => diag_data%var(iLookDIAG%scalarCanopyDepth)%dat(1),   & ! canopy depth                                   (m)
              specificHeatVeg   => mpar_data%var(iLookPARAM%specificHeatVeg)%dat(1),    & ! specific heat of vegetation                    (J kg-1 K-1)
              maxMassVegetation => mpar_data%var(iLookPARAM%maxMassVegetation)%dat(1),  & ! maximum mass of vegetation                     (kg m-2)
              snowfrz_scale     => mpar_data%var(iLookPARAM%snowfrz_scale)%dat(1)       & ! scaling parameter for the snow freezing curve  (K-1)
              )

              diffT = scalarCanopyTempTrial - Tfreeze
              enthVegP = specificHeatVeg * maxMassVegetation * scalarCanopyTempPrime / canopyDepth

              if(diffT>=0._rkind)then
                enthLiqP = Cp_water * ( scalarCanopyWatPrime * diffT + scalarCanopyWatTrial * scalarCanopyTempPrime )/ canopyDepth
                enthIceP = 0._rkind
              else
                integral = (1._rkind/snowfrz_scale) * atan(snowfrz_scale * diffT)
                d_integral_dTk = 1._rkind / (1._rkind + (snowfrz_scale * diffT)**2_i4b) ! Note: volFracLiq = d_integral_dTk*volFracWat
                enthLiqP = Cp_water * ( scalarCanopyWatPrime*integral + scalarCanopyWatTrial*scalarCanopyTempPrime*d_integral_dTk )/ canopyDepth
                enthIceP = Cp_ice * ( scalarCanopyWatPrime*( diffT - integral ) + scalarCanopyWatTrial*scalarCanopyTempPrime*( 1._rkind - d_integral_dTk ) ) / canopyDepth
              endif

              scalarCanopyEnthalpyPrime = enthVegP + enthLiqP + enthIceP

            end associate vegVars

          case(iname_snow)

            ! association to necessary variables for snow
            snowVars: associate(&
              snowfrz_scale => mpar_data%var(iLookPARAM%snowfrz_scale)%dat(1)   & ! scaling parameter for the snow freezing curve (K-1)
              )

              diffT = mLayerTempTrial(iLayer) - Tfreeze  ! diffT<0._rkind because snow is frozen
              integral = (1._rkind/snowfrz_scale) * atan(snowfrz_scale * diffT)
              d_integral_dTk = 1._rkind / (1._rkind + (snowfrz_scale * diffT)**2_i4b) ! Note: volFracLiq = d_integral_dTk*volFracWat
              enthLiqP = iden_water * Cp_water * ( mLayerVolFracWatPrime(iLayer)*integral &
                                                  + mLayerVolFracWatTrial(iLayer)*mLayerTempPrime(iLayer)*d_integral_dTk )
              enthIceP = iden_water * Cp_ice * ( mLayerVolFracWatPrime(iLayer)*( diffT - integral ) &
                                                  + mLayerVolFracWatTrial(iLayer)*mLayerTempPrime(iLayer)*( 1._rkind - d_integral_dTk ) )
              enthAirP = iden_air * Cp_air * ( mLayerTempPrime(iLayer)  - mLayerVolFracWatPrime(iLayer) * ( (iden_water/iden_ice)*( diffT - integral ) + integral )  &
                                              - mLayerVolFracWatTrial(iLayer)*mLayerTempPrime(iLayer)*( (iden_water/iden_ice)*( 1._rkind - d_integral_dTk ) + d_integral_dTk ) )

              mLayerEnthalpyPrime(iLayer)      = enthLiqP + enthIceP + enthAirP

            end associate snowVars

          case(iname_soil)
            ! make association to variables for soil
            soilVars: associate(&
              soil_dens_intr => mpar_data%var(iLookPARAM%soil_dens_intr)%dat(ixControlIndex)      , & ! intrinsic soil density             (kg m-3)
              theta_sat      => mpar_data%var(iLookPARAM%theta_sat)%dat(ixControlIndex)           , & ! soil porosity                      (-)
              theta_res      => mpar_data%var(iLookPARAM%theta_res)%dat(ixControlIndex)           , & ! volumetric residual water content  (-)
              vGn_alpha      => mpar_data%var(iLookPARAM%vGn_alpha)%dat(ixControlIndex)           , & ! van Genuchten "alpha" parameter    (m-1)
              vGn_n          => mpar_data%var(iLookPARAM%vGn_n)%dat(ixControlIndex)                 & ! van Genuchten "n" parameter        (-)
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
                enthWaterP = iden_water * Cp_water * ( mLayerMatricHeadPrime(ixControlIndex)*dVolTot_dPsi0(ixControlIndex)*diffT &
                                                      + volFracWat*mLayerTempPrime(iLayer) )  ! valid for temperatures below freezing also

              ! *** compute enthalpy prime of water for frozen conditions
              else
                ! NOTE: mLayerPsiLiq is the liquid water matric potential from the Clapeyron equation, used to separate the total water into liquid water and ice
                !       mLayerPsiLiq is DIFFERENT from the liquid water matric potential used in the flux calculations
                xConst        = LH_fus/(gravity*Tfreeze)        ! m K-1 (NOTE: J = kg m2 s-2)
                mLayerPsiLiq  = xConst*diffT   ! liquid water matric potential from the Clapeyron eqution
                ! NOTE: integral_psiLiq = diffT * ( (theta_sat - theta_res)*gauss_hg_T + theta_res )
                d_integral_psiLiq_dTk = volFracLiq(mLayerPsiLiq,vGn_alpha,theta_res,theta_sat,vGn_n,vGn_m)

                enthLiqP = iden_water * Cp_water * mLayerTempPrime(iLayer)*d_integral_psiLiq_dTk
                enthIceP = iden_ice * Cp_ice * ( mLayerMatricHeadPrime(ixControlIndex)*dVolTot_dPsi0(ixControlIndex)*diffT &
                                               + volFracWat*mLayerTempPrime(iLayer) ) - iden_ice * Cp_ice *  mLayerTempPrime(iLayer)* d_integral_psiLiq_dTk
                enthWaterP = enthIceP + enthLiqP

              endif ! (if frozen conditions)

              enthSoilP = soil_dens_intr * Cp_soil * (1._rkind - theta_sat)*mlayerTempPrime(iLayer)
              enthAirP = iden_air * Cp_air * ( mLayerMatricHeadPrime(ixControlIndex)*dVolTot_dPsi0(ixControlIndex)*diffT &
                                               + (1._rkind - theta_sat - volFracWat)*mlayerTempPrime(iLayer) )

              mLayerEnthalpyPrime(iLayer) = enthWaterP + enthSoilP + enthAirP

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

end module enthalpyTempAddPrime_module

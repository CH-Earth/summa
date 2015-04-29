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

module ctranspire_module
USE nrtype
implicit none
private
public::ctranspire
contains


 ! ************************************************************************************************
 ! public subroutine ctranspire: compute energy flux
 ! ************************************************************************************************
 subroutine ctranspire(err,message)
 USE multiconst,only:ave_slp,&  ! mean sea level pressure     (Pa)
                     gascnst,&  ! gas constant for dry air    (Pa K-1 m3 kg-1)
                     Cp_air,&   ! specific heat of air        (J kg-1 K-1)
                     gravity,&  ! acceleration of gravity     (m s-2)
                     LH_sub,&   ! latent heat of sublimation  (J kg-1)
                     Em_Sno,&   ! emissivity of snow          (-)
                     sigma      ! Stefan Boltzman constant    (W m-2 K-4)
 USE conv_funcs_module,only:relhm2sphm                                          ! compute specific humidity
 USE snow_utils_module,only:astability                                          ! compute atmospheric stability
 USE data_struc,only:mpar_data,forc_data,mvar_data,indx_data,ix_soil,ix_snow    ! data structures
 USE var_lookup,only:iLookPARAM,iLookFORCE,iLookMVAR,iLookINDEX                 ! named variables for structure elements
 ! compute energy fluxes within the domain
 implicit none
 ! dummy variables
 integer(i4b),intent(out)      :: err                  ! error code
 character(*),intent(out)      :: message              ! error message
 ! local pointers
 real(dp),pointer              :: mheight              ! measurement height (m)
 real(dp),pointer              :: Fabs_vis             ! fraction of radiation in visible part of spectrum (-)
 real(dp),pointer              :: rad_ext              ! extinction coefficient for penetrating sw radiation (m-1)
 real(dp),pointer              :: pptrate              ! precipitation rate (kg m-2 s-1)
 real(dp),pointer              :: sw_down              ! downward shortwave radiation (W m-2)
 real(dp),pointer              :: lw_down              ! downward longwave radiation (W m-2)
 real(dp),pointer              :: airtemp              ! air temperature at 2 meter height (K)
 real(dp),pointer              :: windspd              ! wind speed at 10 meter height (m s-1)
 real(dp),pointer              :: airpres              ! air pressure at 2 meter height (Pa)
 real(dp),pointer              :: spechum              ! specific humidity at 2 meter height (g g-1)
 real(dp),pointer              :: scalarAlbedo         ! surface albedo (-)
 real(dp),pointer              :: scalarSenHeat        ! sensible heat flux at the surface (W m-2)
 real(dp),pointer              :: scalarLatHeat        ! latent heat flux at the surface (W m-2)
 real(dp),pointer              :: mLayerTemp(:)        ! temperature of each layer (K)
 real(dp),pointer              :: iLayerThermalC(:)    ! thermal conductivity at the interface of each layer (W m-1 K-1)
 real(dp),pointer              :: mLayerHeight(:)      ! height at the mid-point of each layer (m)
 real(dp),pointer              :: iLayerHeight(:)      ! height at the interface of each layer (m)
 real(dp),pointer              :: iLayerRadCondFlux(:) ! energy flux at the interface of each layer (W m-2)
 real(dp),pointer              :: mLayerRadCondFlux(:) ! change in energy in each layer (J m-2 s-1)
 integer(i4b),pointer          :: nLayers              ! number of layers
 ! define local variables
 character(LEN=256)            :: cmessage             ! error message of downwind routine
 real(dp)                      :: ptcorr               ! potential temperature correction
 real(dp)                      :: T_surf               ! surface temperature (corrected)
 real(dp)                      :: T__air               ! air temperature (corrected)
 real(dp)                      :: airden               ! air density (kg m-3)
 real(dp)                      :: vhtcap               ! vol heat capacity (J m-3 K-1)
 real(dp)                      :: Q_surf               ! specific humidity at the surface (g/g)
 real(dp)                      :: T_grad               ! potential temperature gradient between air and surface
 real(dp)                      :: T_mean               ! mean temperature between air and surface
 real(dp)                      :: RiBulk               ! bulk Richardson number
 real(dp)                      :: ExCoef               ! turbulent exchange coefficient
 real(dp)                      :: LW__up               ! upward longwave radiation, positive downward (W m-2)
 real(dp)                      :: sw__net              ! net shortwave radiation (W m-2)
 real(dp)                      :: sw__vis              ! shortwave radiation in the visible part of the spectrum
 real(dp)                      :: sw__nir              ! shortwave radiation in the near-infra-red part of the spectrum
 real(dp)                      :: depth_layr           ! depth at the layer interface (m)
 real(dp)                      :: condv_flux           ! conductive heat flux (W m-2)
 real(dp)                      :: swrad_flux           ! radiative heat flux (W m-2)
 integer(i4b)                  :: iLayer               ! index of model layer
 ! initialize error control
 err=0; message="f-energyflux/"
 ! assign pointers to model parameters
 mheight           => mpar_data%var(iLookPARAM%mheight)               ! measurement height (m)
 Fabs_vis          => mpar_data%var(iLookPARAM%Fabs_vis)              ! fraction of radiation in visible part of spectrum (-)
 rad_ext           => mpar_data%var(iLookPARAM%rad_ext)               ! extinction coefficient for penetrating sw radiation (m-1)
 ! assign pointers to model forcing variables
 pptrate           => forc_data%var(iLookFORCE%pptrate)               ! precipitation rate (kg m-2 s-1)
 sw_down           => forc_data%var(iLookFORCE%sw_down)               ! downward shortwave radiation (W m-2)
 lw_down           => forc_data%var(iLookFORCE%lw_down)               ! downward longwave radiation (W m-2)
 airtemp           => forc_data%var(iLookFORCE%airtemp)               ! air temperature at 2 meter height (K)
 windspd           => forc_data%var(iLookFORCE%windspd)               ! wind speed at 10 meter height (m s-1)
 airpres           => forc_data%var(iLookFORCE%airpres)               ! air pressure at 2 meter height (Pa)
 spechum           => forc_data%var(iLookFORCE%spechum)               ! specific humidity at 2 meter height (g g-1)
 ! assign pointers to model variables
 scalarAlbedo      => mvar_data%var(iLookMVAR%scalarAlbedo)%dat(1)    ! surface albedo (-)
 scalarSenHeat     => mvar_data%var(iLookMVAR%scalarSenHeat)%dat(1)   ! sensible heat flux at the surface (W m-2)
 scalarLatHeat     => mvar_data%var(iLookMVAR%scalarLatHeat)%dat(1)   ! latent heat flux at the surface (W m-2)
 mLayerTemp        => mvar_data%var(iLookMVAR%mLayerTemp)%dat         ! temperature of each layer (K)
 iLayerThermalC    => mvar_data%var(iLookMVAR%iLayerThermalC)%dat     ! thermal conductivity at the interface of each layer (W m-1 K-1)
 mLayerHeight      => mvar_data%var(iLookMVAR%mLayerHeight)%dat       ! height at the mid-point of each layer (m)
 iLayerHeight      => mvar_data%var(iLookMVAR%iLayerHeight)%dat       ! height at the interface of each layer (m)
 iLayerRadCondFlux => mvar_data%var(iLookMVAR%iLayerRadCondFlux)%dat  ! energy flux at the interface of each layer (W m-2)
 mLayerRadCondFlux => mvar_data%var(iLookMVAR%mLayerRadCondFlux)%dat  ! change in energy in each layer (J m-2 s-1)
 ! assign pointers to index variables
 nLayers         => indx_data%var(iLookINDEX%nLayers)%dat(1)      ! number of layers


 CanInterceptFlux  ! canopy interception flux
 fracAreaLeafWet   ! fractional area of a leaf that collects water (0.25)
 lai               ! leaf area index
 sai               ! stem area index

 ! mean maximum canopy storage
 canopyMeanMaxStorGlobal = 0.1_dp  ! mean maximum storage of canopy water per unit (LAI+SAI)
 canopyMeanMaxStorLocal  = canopyStorGlobal*(lai + sai) ! mean maximum canopy storage scaled by local (lai + sai) values

 ! coefficient of variation in maximum canopy storage per unit (LAI+SAI)
 canopyCvarMaxStorage    = 1.0_dp

 ! define the fractional area of a leaf that collects water
 canopyFracCollectH20    = 0.25_dp

 ! define the time constant for unloading/drip
 dripTimeRain = 1.5d+5  ! s
 dripTimeSnow = 1.5d+6  ! s

 ! define parameters of the probability distribution describing maximum canopy storage
 canopyMaxPdisZeta   = sqrt( log(1._dp) + canopyCvarMaxStorage**2._dp)
 canopyMaxPdisLambda = log(canopyMeanMaxStorLocal) - 0.5_dp*(canopyMaxPdisZeta**2._dp)

 ! define canopy storage
 canopyStorage = 0.1_dp  ! state variable

 ! compute dry (transpiring) area of the canopy and the wet (evaporating) area of the canopy
 erf_arg       = (log(canopyStorage) - canopyMaxPdisLambda) / (sqrt(2._dp)*canopyMaxPdisZeta)
 canopyFracDry = 0.5_dp*(1._dp - erf(erf_arg)) * canopyFracCollectH20
 canopyFracWet = 1._dp - canopyFracDry

 ! compute throughfall (kg m-2 s-1)
 canopyThruRain = canopyFracWet * rainfall
 canopyThruSnow = canopyFracWet * snowfall

 ! compute canopy drip (kg m-2 s-1)
 canopyDripRain = canopyStorage/dripTimeRain
 canopyDripSnow = canopyStorage/dripTimeSnow

 ! compute canopy evaporation/sublimation (kg m-2 s-1)
 canopyEvap     =

 ! define the canopy water balance (kg m-2 s-1)
 delCanopyWat = rainfall - canopyThruRain - canopyDripRain

(alog(stor)-alam)/zeta
fdry = 0.5d*erfc(z_dm/sqrt(2.d))

 ! canopy interception flux (kg m-2 s-1)

 canopyInterceptFlux = rainfall * fracAreaLeafWet * ( 1._dp - exp(-0.5_dp*(lai + sai)) )
 ! canopy throughfall
 canopyThrufallFlux  =

 ! compute potential temperature
 ptcorr = (ave_slp/airpres)**(gascnst/Cp_air) ! potential temperature correction
 T_surf = mLayerTemp(nLayers) * ptcorr
 T__air = airtemp * ptcorr
 ! compute air density and volumetric heat capacity -- independent of snow temp
 airden = airpres / (airtemp*gascnst)         ! air density (kg m-3)
 vhtcap = Cp_air * airden                     ! vol heat capacity (J m-3 K-1)
 ! compute specific humidity at the surface (assume saturated)
 Q_surf = relhm2sphm(1._dp,airpres,mLayerTemp(nLayers))
 ! compute stability corrections
 T_grad = T__air - T_surf
 T_mean = 0.5_dp * (T__air + T_surf)
 RiBulk = (T_grad/T_mean) * ((gravity*mheight)/(windspd*windspd))
 call astability(RiBulk,ExCoef,err,cmessage)
 if(err/=0)then; err=10; message=trim(message)//trim(cmessage); return; endif
 ! compute sensible and latent heat
 scalarSenHeat = vhtcap * windspd * ExCoef * (T__air - T_surf)                 ! W m-2
 scalarLatHeat = LH_sub * airden * windspd * ExCoef * (spechum - Q_surf)       ! W m-2
 ! compute upward longwave radiation -- note sign convention is positive downward
 LW__up = -Em_Sno * sigma * mLayerTemp(nLayers)**4._dp                         ! W m-2
 ! separate absorbed shortwave radiation into visible and near infra-red components
 sw__net = sw_down * (1._dp - scalarAlbedo)
 sw__vis = Fabs_vis * sw__net  ! sw radiation in the visible part of the spectrum
 sw__nir = sw__net - sw__vis   ! sw radiation in the near infra-red part of the spectrum
 ! compute flux at the top of the domain
 iLayerRadCondFlux(nLayers) = scalarSenHeat + scalarLatHeat + lw_down + LW__up  ! W m-2
 ! prescribe zero flux at the bottom of the domain
 iLayerRadCondFlux(iLayer) = 0._dp
 ! compute the energy fluxes at layer interfaces (W m-2) -- the **top** of each layer
 !   NOTE: fluxes are positive downward
 do iLayer=1,nLayers
  ! flux at top of the domain is computed earlier
  if(ilayer<nLayers)then
   ! compute conductive flux (W m-2)
   condv_flux = iLayerThermalC(iLayer) * ( (mLayerTemp(iLayer+1)   - mLayerTemp(iLayer)) /  &
                                           (mLayerHeight(iLayer+1) - mLayerHeight(iLayer))   )
   ! compute shortwave radiative flux (W m-2)
   if(mLayerHeight(iLayer+1)>0._dp)then ! snow is above iLayerHeight(iLayer)
    depth_layr = iLayerHeight(nLayers) - iLayerHeight(iLayer) ! depth of the top of the layer
    swrad_flux = sw__vis * exp(-rad_ext*depth_layr)           ! flux at the top of the layer
   else
    swrad_flux = 0._dp  ! soil is above iLayerHeight(iLayer)
   endif
   ! compute total energy flux (W m-2)
   iLayerRadCondFlux(iLayer) = condv_flux + swrad_flux
  endif
  ! compute temporal change in layer energy (J m-2 s-1)
  mLayerRadCondFlux(iLayer)  = iLayerRadCondFlux(iLayer) - iLayerRadCondFlux(iLayer-1)
 end do  ! looping through layers
 end subroutine ctranspire


end module ctranspire_module

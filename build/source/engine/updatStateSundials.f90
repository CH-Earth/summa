
module updatStateSundials_module
USE nrtype
! physical constants
USE multiconst,only:&
                    Tfreeze,     & ! freezing point of pure water  (K)
                    iden_air,    & ! intrinsic density of air      (kg m-3)
                    iden_ice,    & ! intrinsic density of ice      (kg m-3)
                    iden_water,  & ! intrinsic density of water    (kg m-3)
                    gravity,     & ! gravitational acceleteration  (m s-2)
                    LH_fus         ! latent heat of fusion         (J kg-1)
implicit none
private
public::updateSnowSundials
public::updateSoilSundials
public::updateSoilSundials2
public::updateVegSundials

real(rkind),parameter     :: verySmall=1e-14_rkind ! a very small number (used to avoid divide by zero)

contains
 
 
  ! *************************************************************************************************************
 ! public subroutine updateVegSundials: compute phase change impacts on volumetric liquid water and ice
 ! Input: Theta * canopyDepth * iden_water
 ! Outputs: VolFracLiq * canopyDepth * iden_water    and  VolFracIce * canopyDepth * iden_ice
 ! *************************************************************************************************************
 subroutine updateVegSundials(&   
                       ! input
                       Temp            ,& ! intent(in): temperature (K)
                       Theta           ,& ! intent(in): volume fraction of total water (-)
                       snowfrz_scale   ,& ! intent(in): scaling parameter for the snow freezing curve (K-1)
                       TempPrime       ,& ! intent(in): temperature (K)
                       ThetaPrime      ,& ! intent(in): volume fraction of total water (-)
                       ! output
                       VolFracLiq      ,& ! intent(out): volumetric fraction of liquid water (-)
                       VolFracIce      ,& ! intent(out): volumetric fraction of ice (-)
                       VolFracLiqPrime ,& ! intent(out): volumetric fraction of liquid water (-)
                       VolFracIcePrime ,& ! intent(out): volumetric fraction of ice (-)
                       fLiq                  ,& ! intent(out): fraction of liquid water (-)
                       err,message)        ! intent(out): error control
 ! utility routines
 USE snow_utils_module,only:fracliquid     ! compute volumetric fraction of liquid water
 USE snow_utils_module,only:dFracLiq_dTk   ! differentiate the freezing curve w.r.t. temperature (snow)
 implicit none
 ! input variables
 real(rkind),intent(in)           :: Temp           ! temperature (K)
 real(rkind),intent(in)           :: Theta          ! volume fraction of total water (-)
 real(rkind),intent(in)           :: snowfrz_scale        ! scaling parameter for the snow freezing curve (K-1)
 real(rkind),intent(in)           :: TempPrime           ! temperature (K)
 real(rkind),intent(in)           :: ThetaPrime          ! volume fraction of total water (-)
 ! output variables
 real(rkind),intent(out)          :: VolFracLiq     ! volumetric fraction of liquid water (-)
 real(rkind),intent(out)          :: VolFracIce     ! volumetric fraction of ice (-)
 real(rkind),intent(out)          :: VolFracLiqPrime     ! volumetric fraction of liquid water (-)
 real(rkind),intent(out)          :: VolFracIcePrime     ! volumetric fraction of ice (-)
 real(rkind),intent(out)          :: fLiq                 ! fraction of liquid water (-)
 ! error control
 integer(i4b),intent(out)      :: err                  ! error code
 character(*),intent(out)      :: message              ! error message
 ! initialize error control
 err=0; message="updateVegSundials/"

 ! compute the volumetric fraction of liquid water and ice (-)
 fLiq = fracliquid(Temp,snowfrz_scale)
 VolFracLiq = fLiq*Theta
 VolFracIce = (1._rkind - fLiq)*Theta
 VolFracLiqPrime = fLiq * ThetaPrime + dFracLiq_dTk(Temp,snowfrz_scale) * Theta * TempPrime
 VolFracIcePrime = ( ThetaPrime - VolFracLiqPrime )
 

 end subroutine updateVegSundials
 
 
   ! *************************************************************************************************************
 ! public subroutine updateSnowSundials: compute phase change impacts on volumetric liquid water and ice
 ! *************************************************************************************************************
 subroutine updateSnowSundials(&   
                       ! input
                       mLayerTemp            ,& ! intent(in): temperature (K)
                       mLayerTheta           ,& ! intent(in): volume fraction of total water (-)
                       snowfrz_scale         ,& ! intent(in): scaling parameter for the snow freezing curve (K-1)
                       mLayerTempPrime       ,& ! intent(in): temperature (K)
                       mLayerThetaPrime      ,& ! intent(in): volume fraction of total water (-)
                       ! output
                       mLayerVolFracLiq      ,& ! intent(out): volumetric fraction of liquid water (-)
                       mLayerVolFracIce      ,& ! intent(out): volumetric fraction of ice (-)
                       mLayerVolFracLiqPrime ,& ! intent(out): volumetric fraction of liquid water (-)
                       mLayerVolFracIcePrime ,& ! intent(out): volumetric fraction of ice (-)
                       fLiq                  ,& ! intent(out): fraction of liquid water (-)
                       err,message)        ! intent(out): error control
 ! utility routines
 USE snow_utils_module,only:fracliquid     ! compute volumetric fraction of liquid water
 USE snow_utils_module,only:dFracLiq_dTk   ! differentiate the freezing curve w.r.t. temperature (snow)
 implicit none
 ! input variables
 real(rkind),intent(in)           :: mLayerTemp           ! temperature (K)
 real(rkind),intent(in)           :: mLayerTheta          ! volume fraction of total water (-)
 real(rkind),intent(in)           :: snowfrz_scale        ! scaling parameter for the snow freezing curve (K-1)
 real(rkind),intent(in)           :: mLayerTempPrime           ! temperature (K)
 real(rkind),intent(in)           :: mLayerThetaPrime          ! volume fraction of total water (-)
 ! output variables
 real(rkind),intent(out)          :: mLayerVolFracLiq     ! volumetric fraction of liquid water (-)
 real(rkind),intent(out)          :: mLayerVolFracIce     ! volumetric fraction of ice (-)
 real(rkind),intent(out)          :: mLayerVolFracLiqPrime     ! volumetric fraction of liquid water (-)
 real(rkind),intent(out)          :: mLayerVolFracIcePrime     ! volumetric fraction of ice (-)
 real(rkind),intent(out)          :: fLiq                 ! fraction of liquid water (-)
 ! error control
 integer(i4b),intent(out)      :: err                  ! error code
 character(*),intent(out)      :: message              ! error message
 ! initialize error control
 err=0; message="updateSnowSundials/"

 ! compute the volumetric fraction of liquid water and ice (-)
 fLiq = fracliquid(mLayerTemp,snowfrz_scale)
 mLayerVolFracLiq = fLiq*mLayerTheta
 mLayerVolFracIce = (1._rkind - fLiq)*mLayerTheta*(iden_water/iden_ice)
 mLayerVolFracLiqPrime = fLiq * mLayerThetaPrime + dFracLiq_dTk(mLayerTemp,snowfrz_scale) * mLayerTheta * mLayerTempPrime
 mLayerVolFracIcePrime = ( mLayerThetaPrime - mLayerVolFracLiqPrime ) * (iden_water/iden_ice)
 

 end subroutine updateSnowSundials

 ! *************************************************************************************************************
 ! public subroutine updateSoilSundials: compute phase change impacts on matric head and volumetric liquid water and ice
 ! *************************************************************************************************************
 subroutine updateSoilSundials(&
                       ! input
                       dt_cur,                &
                       mLayerTemp            ,& ! intent(in): temperature vector (K)
                       mLayerMatricHead      ,& ! intent(in): matric head (m)
                       mLayerMatricHeadPrev  ,&
                       mLayerVolFracWatPrev  ,&
                       mLayerTempPrime       ,& ! intent(in)
                       mLayerMatricHeadPrime, & ! intent(in)
                       vGn_alpha             ,& ! intent(in): van Genutchen "alpha" parameter
                       vGn_n                 ,& ! intent(in): van Genutchen "n" parameter
                       theta_sat             ,& ! intent(in): soil porosity (-)
                       theta_res             ,& ! intent(in): soil residual volumetric water content (-)
                       vGn_m                 ,& ! intent(in): van Genutchen "m" parameter (-)
                       ! output
                       mLayerVolFracWat ,& ! intent(out): volumetric fraction of total water (-)
                       mLayerVolFracLiq ,& ! intent(out): volumetric fraction of liquid water (-)
                       mLayerVolFracIce ,& ! intent(out): volumetric fraction of ice (-)
                       mLayerVolFracWatPrime ,& ! intent(out): volumetric fraction of total water (-)
                       mLayerVolFracLiqPrime ,& ! intent(out): volumetric fraction of liquid water (-)
                       mLayerVolFracIcePrime ,& ! intent(out): volumetric fraction of ice (-)
                       err,message)        ! intent(out): error control
 ! utility routines
 USE soil_utils_module,only:volFracLiq     ! compute volumetric fraction of liquid water based on matric head
 USE soil_utils_module,only:matricHead     ! compute the matric head based on volumetric liquid water content
 USE soil_utils_module,only:dTheta_dPsi 
 implicit none
 ! input variables
 real(rkind),intent(in)           :: dt_cur
 real(rkind),intent(in)           :: mLayerTemp           ! estimate of temperature (K)
 real(rkind),intent(in)           :: mLayerMatricHead     ! matric head (m)
 real(rkind),intent(in)           :: mLayerMatricHeadPrev     ! matric head (m)
 real(rkind),intent(in)           :: mLayerVolFracWatPrev
 real(rkind),intent(in)           :: mLayerTempPrime
 real(rkind),intent(in)           :: mLayerMatricHeadPrime ! matric head (m)
 real(rkind),intent(in)           :: vGn_alpha            ! van Genutchen "alpha" parameter
 real(rkind),intent(in)           :: vGn_n                ! van Genutchen "n" parameter
 real(rkind),intent(in)           :: theta_sat            ! soil porosity (-)
 real(rkind),intent(in)           :: theta_res            ! soil residual volumetric water content (-)
 real(rkind),intent(in)           :: vGn_m                ! van Genutchen "m" parameter (-)
 ! output variables
 real(rkind),intent(out)          :: mLayerVolFracWat     ! fractional volume of total water (-)
 real(rkind),intent(out)          :: mLayerVolFracLiq     ! volumetric fraction of liquid water (-)
 real(rkind),intent(out)          :: mLayerVolFracIce     ! volumetric fraction of ice (-)
 real(rkind),intent(out)          :: mLayerVolFracWatPrime     ! fractional volume of total water (-)
 real(rkind),intent(out)          :: mLayerVolFracLiqPrime     ! volumetric fraction of liquid water (-)
 real(rkind),intent(out)          :: mLayerVolFracIcePrime     ! volumetric fraction of ice (-)
 integer(i4b),intent(out)      :: err                  ! error code
 character(*),intent(out)      :: message              ! error message
 ! define local variables
 real(rkind)                      :: TcSoil               ! critical soil temperature when all water is unfrozen (K)
 real(rkind)                      :: xConst               ! constant in the freezing curve function (m K-1)
 real(rkind)                      :: mLayerPsiLiq         ! liquid water matric potential (m)
 ! initialize error control
 err=0; message="updateSoilSundials/"
 
 ! compute fractional **volume** of total water (liquid plus ice)
  mLayerVolFracWat = volFracLiq(mLayerMatricHead,vGn_alpha,theta_res,theta_sat,vGn_n,vGn_m)
 ! mLayerVolFracWatPrime = dTheta_dPsi(mLayerMatricHead,vGn_alpha,theta_res,theta_sat,vGn_n,vGn_m) * mLayerMatricHeadPrime
  if( abs(mLayerMatricHead - mLayerMatricHeadPrev) < verySmall )then
      mLayerVolFracWatPrime =  (mLayerVolFracWat - mLayerVolFracWatPrev) / dt_cur
  else
      mLayerVolFracWatPrime = mLayerMatricHeadPrime * (mLayerVolFracWat - mLayerVolFracWatPrev) / (mLayerMatricHead - mLayerMatricHeadPrev)
  endif

  


 
 if(mLayerVolFracWat > theta_sat)then; err=20; message=trim(message)//'volume of liquid and ice exceeds porosity'; return; end if
 
 ! compute the critical soil temperature where all water is unfrozen (K)
 ! (eq 17 in Dall'Amico 2011)
 TcSoil = Tfreeze + min(mLayerMatricHead,0._rkind)*gravity*Tfreeze/LH_fus  ! (NOTE: J = kg m2 s-2, so LH_fus is in units of m2 s-2)

 ! *** compute volumetric fraction of liquid water and ice for partially frozen soil
 if(mLayerTemp < TcSoil)then ! (check if soil temperature is less than the critical temperature)
  ! - volumetric liquid water content (-)
  ! NOTE: mLayerPsiLiq is the liquid water matric potential from the Clapeyron equation, used to separate the total water into liquid water and ice
  !       mLayerPsiLiq is DIFFERENT from the liquid water matric potential used in the flux calculations
  xConst           = LH_fus/(gravity*Tfreeze)        ! m K-1 (NOTE: J = kg m2 s-2)
  mLayerPsiLiq     = xConst*(mLayerTemp - Tfreeze)   ! liquid water matric potential from the Clapeyron eqution
  mLayerVolFracLiq = volFracLiq(mLayerPsiLiq,vGn_alpha,theta_res,theta_sat,vGn_n,vGn_m)
  if(mLayerPsiLiq<0._rkind)then
   mLayerVolFracLiqPrime = dTheta_dPsi(mLayerPsiLiq,vGn_alpha,theta_res,theta_sat,vGn_n,vGn_m) * xConst * mLayerTempPrime
  else
   mLayerVolFracLiqPrime = 0._rkind
  endif

  ! - volumetric ice content (-)
  mLayerVolFracIce = mLayerVolFracWat - mLayerVolFracLiq
  mLayerVolFracIcePrime = mLayerVolFracWatPrime - mLayerVolFracLiqPrime
  

 ! *** compute volumetric fraction of liquid water and ice for unfrozen soil
 else
 
  ! all water is unfrozen
  mLayerVolFracLiq = mLayerVolFracWat
  mLayerVolFracIce = 0._rkind
  
  mLayerVolFracLiqPrime = mLayerVolFracWatPrime
  mLayerVolFracIcePrime = 0._rkind

 end if  ! (check if soil is partially frozen)


 end subroutine updateSoilSundials

 ! *************************************************************************************************************
 ! public subroutine updateSoilSundials: compute phase change impacts on matric head and volumetric liquid water and ice
 ! *************************************************************************************************************
 subroutine updateSoilSundials2(&
                       ! input
                       mLayerTemp            ,& ! intent(in): temperature vector (K)
                       mLayerMatricHead      ,& ! intent(in): matric head (m)
                       mLayerTempPrime       ,& ! intent(in)
                       mLayerMatricHeadPrime, & ! intent(in)
                       vGn_alpha             ,& ! intent(in): van Genutchen "alpha" parameter
                       vGn_n                 ,& ! intent(in): van Genutchen "n" parameter
                       theta_sat             ,& ! intent(in): soil porosity (-)
                       theta_res             ,& ! intent(in): soil residual volumetric water content (-)
                       vGn_m                 ,& ! intent(in): van Genutchen "m" parameter (-)
                       ! output
                       mLayerVolFracWat ,& ! intent(out): volumetric fraction of total water (-)
                       mLayerVolFracLiq ,& ! intent(out): volumetric fraction of liquid water (-)
                       mLayerVolFracIce ,& ! intent(out): volumetric fraction of ice (-)
                       mLayerVolFracWatPrime ,& ! intent(out): volumetric fraction of total water (-)
                       mLayerVolFracLiqPrime ,& ! intent(out): volumetric fraction of liquid water (-)
                       mLayerVolFracIcePrime ,& ! intent(out): volumetric fraction of ice (-)
                       err,message)        ! intent(out): error control
 ! utility routines
 USE soil_utils_module,only:volFracLiq     ! compute volumetric fraction of liquid water based on matric head
 USE soil_utils_module,only:matricHead     ! compute the matric head based on volumetric liquid water content
 USE soil_utils_module,only:dTheta_dPsi 
 implicit none
 ! input variables
 real(rkind),intent(in)           :: mLayerTemp           ! estimate of temperature (K)
 real(rkind),intent(in)           :: mLayerMatricHead     ! matric head (m)
 real(rkind),intent(in)           :: mLayerTempPrime
 real(rkind),intent(in)           :: mLayerMatricHeadPrime ! matric head (m)
 real(rkind),intent(in)           :: vGn_alpha            ! van Genutchen "alpha" parameter
 real(rkind),intent(in)           :: vGn_n                ! van Genutchen "n" parameter
 real(rkind),intent(in)           :: theta_sat            ! soil porosity (-)
 real(rkind),intent(in)           :: theta_res            ! soil residual volumetric water content (-)
 real(rkind),intent(in)           :: vGn_m                ! van Genutchen "m" parameter (-)
 ! output variables
 real(rkind),intent(out)          :: mLayerVolFracWat     ! fractional volume of total water (-)
 real(rkind),intent(out)          :: mLayerVolFracLiq     ! volumetric fraction of liquid water (-)
 real(rkind),intent(out)          :: mLayerVolFracIce     ! volumetric fraction of ice (-)
 real(rkind),intent(out)          :: mLayerVolFracWatPrime     ! fractional volume of total water (-)
 real(rkind),intent(out)          :: mLayerVolFracLiqPrime     ! volumetric fraction of liquid water (-)
 real(rkind),intent(out)          :: mLayerVolFracIcePrime     ! volumetric fraction of ice (-)
 integer(i4b),intent(out)      :: err                  ! error code
 character(*),intent(out)      :: message              ! error message
 ! define local variables
 real(rkind)                      :: TcSoil               ! critical soil temperature when all water is unfrozen (K)
 real(rkind)                      :: xConst               ! constant in the freezing curve function (m K-1)
 real(rkind)                      :: mLayerPsiLiq         ! liquid water matric potential (m)
 ! initialize error control
 err=0; message="updateSoilSundials2/"
 
 ! compute fractional **volume** of total water (liquid plus ice)
 mLayerVolFracWat = volFracLiq(mLayerMatricHead,vGn_alpha,theta_res,theta_sat,vGn_n,vGn_m)
 mLayerVolFracWatPrime = dTheta_dPsi(mLayerMatricHead,vGn_alpha,theta_res,theta_sat,vGn_n,vGn_m) * mLayerMatricHeadPrime
 
 if(mLayerVolFracWat > theta_sat)then; err=20; message=trim(message)//'volume of liquid and ice exceeds porosity'; return; end if
 

 ! compute the critical soil temperature where all water is unfrozen (K)
 ! (eq 17 in Dall'Amico 2011)
 TcSoil = Tfreeze + min(mLayerMatricHead,0._rkind)*gravity*Tfreeze/LH_fus  ! (NOTE: J = kg m2 s-2, so LH_fus is in units of m2 s-2)

 ! *** compute volumetric fraction of liquid water and ice for partially frozen soil
 if(mLayerTemp < TcSoil)then ! (check if soil temperature is less than the critical temperature)
  ! - volumetric liquid water content (-)
  ! NOTE: mLayerPsiLiq is the liquid water matric potential from the Clapeyron equation, used to separate the total water into liquid water and ice
  !       mLayerPsiLiq is DIFFERENT from the liquid water matric potential used in the flux calculations
  xConst           = LH_fus/(gravity*Tfreeze)        ! m K-1 (NOTE: J = kg m2 s-2)
  mLayerPsiLiq     = xConst*(mLayerTemp - Tfreeze)   ! liquid water matric potential from the Clapeyron eqution
  mLayerVolFracLiq = volFracLiq(mLayerPsiLiq,vGn_alpha,theta_res,theta_sat,vGn_n,vGn_m)
  if(mLayerPsiLiq<0._rkind)then
   mLayerVolFracLiqPrime = dTheta_dPsi(mLayerPsiLiq,vGn_alpha,theta_res,theta_sat,vGn_n,vGn_m) * xConst * mLayerTempPrime
  else
   mLayerVolFracLiqPrime = 0._rkind
  endif

  ! - volumetric ice content (-)
  mLayerVolFracIce = mLayerVolFracWat - mLayerVolFracLiq
  mLayerVolFracIcePrime = mLayerVolFracWatPrime - mLayerVolFracLiqPrime
  

 ! *** compute volumetric fraction of liquid water and ice for unfrozen soil
 else
 
  ! all water is unfrozen
  mLayerVolFracLiq = mLayerVolFracWat
  mLayerVolFracIce = 0._rkind
  
  mLayerVolFracLiqPrime = mLayerVolFracWatPrime
  mLayerVolFracIcePrime = 0._rkind

 end if  ! (check if soil is partially frozen)


 end subroutine updateSoilSundials2
end module updatStateSundials_module

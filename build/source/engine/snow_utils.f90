module snow_utils_module
USE nrtype
implicit none
private
public::fracliquid
public::templiquid
public::dFracLiq_dTk
public::tcond_snow
public::bulkRichardson
public::astability
contains

 ! ***********************************************************************************************************
 ! new function: compute fraction of liquid water
 ! ***********************************************************************************************************
 function fracliquid(Tk,fc_param)
 USE multiconst,only:Tfreeze
 implicit none
 real(dp),intent(in) :: Tk         ! temperature (K)
 real(dp),intent(in) :: fc_param   ! freezing curve parameter (K-1)
 real(dp)            :: fracliquid ! fraction of liquid water (-)
 ! compute fraction of liquid water (-)
 fracliquid = 1._dp / ( 1._dp + (fc_param*( Tfreeze - min(Tk,Tfreeze) ))**2._dp )
 end function fracliquid


 ! ***********************************************************************************************************
 ! new function: invert the fraction of liquid water function
 ! ***********************************************************************************************************
 function templiquid(fracliquid,fc_param)
 USE multiconst,only:Tfreeze
 implicit none
 real(dp),intent(in) :: fracliquid ! fraction of liquid water (-)
 real(dp),intent(in) :: fc_param   ! freezing curve parameter (K-1)
 real(dp)            :: templiquid ! temperature (K)
 ! compute temperature based on the fraction of liquid water (K)
 templiquid = Tfreeze - ((1._dp/fracliquid - 1._dp)/fc_param**2._dp)**(0.5_dp)
 end function templiquid


 ! ***********************************************************************************************************
 ! new function: differentiate the freezing curve
 ! ***********************************************************************************************************
 function dFracLiq_dTk(Tk,fc_param)
 USE multiconst,only:Tfreeze
 implicit none
 ! dummies
 real(dp),intent(in) :: Tk           ! temperature (K)
 real(dp),intent(in) :: fc_param     ! freezing curve parameter (K-1)
 real(dp)            :: dFracLiq_dTk ! differentiate the freezing curve (K-1)
 ! locals
 real(dp)            :: Tdep         ! temperature depression (K)
 real(dp)            :: Tdim         ! dimensionless temperature (-)
 ! compute local variables (just to make things more efficient)
 Tdep = Tfreeze - min(Tk,Tfreeze)
 Tdim = fc_param*Tdep
 ! differentiate the freezing curve w.r.t temperature
 dFracLiq_dTk = (fc_param*2._dp*Tdim) / ( ( 1._dp + Tdim**2._dp)**2._dp )
 end function dFracLiq_dTk


 ! ***********************************************************************************************************
 ! new subroutine: compute thermal conductivity of snow
 ! ***********************************************************************************************************
 subroutine tcond_snow(BulkDenIce,thermlcond,err,message)
 USE multiconst,only:lambda_air,lambda_ice  ! thermal conductivity of air and ice
 USE data_struc,only:model_decisions        ! model decision structure
 USE var_lookup,only:iLookDECISIONS         ! named variables for elements of the decision structure
 implicit none
 real(dp),intent(in)      :: BulkDenIce     ! bulk density of ice (kg m-3)
 real(dp),intent(out)     :: thermlcond     ! thermal conductivity of snow (W m-1 K-1)
 integer(i4b),intent(out) :: err            ! error code
 character(*),intent(out) :: message        ! error message
 ! initialize error control
 err=0; message="f-tcond_snow/"
 ! compute thermal conductivity of snow
 select case(trim(model_decisions(iLookDECISIONS%thermlcond)%decision))
  case('tyen1965'); thermlcond = 3.217d-6 * BulkDenIce**2._dp               ! Yen (1965)
  case('melr1977'); thermlcond = 2.576d-6 * BulkDenIce**2._dp + 7.4d-2      ! Mellor (1977)
  case('jrdn1991'); thermlcond = lambda_air + (7.75d-5*BulkDenIce + 1.105d-6*(BulkDenIce**2._dp)) &
                     * (lambda_ice-lambda_air)                              ! Jordan (1991)
  case('smnv2000'); thermlcond = 0.35_dp                                    ! Smirnova et al. (2000)
  case default
   err=10; message=trim(message)//"unknownOption"; return
 end select
 end subroutine tcond_snow

 ! ***********************************************************************************************************
 ! new function: compute bulk Richardson number
 ! *********************************************************************************************************** 
 subroutine bulkRichardson(airtemp,surfTemp,windspd,mheight,computeDerivative, & ! (input)
                           RiBulk,dRiBulk_dTemp,err,message)                     ! (output)
 USE multiconst,only:gravity         ! acceleration of gravity     (m s-2)
 implicit none
 ! input
 real(dp),intent(in)           :: airtemp       ! air temperature (K)
 real(dp),intent(in)           :: surfTemp      ! surface temperature (K)
 real(dp),intent(in)           :: windspd       ! wind speed (m s-1)
 real(dp),intent(in)           :: mheight       ! measurement height (m)
 logical(lgt),intent(in)       :: computeDerivative ! flag to compute the derivative
 ! output
 real(dp),intent(out)          :: RiBulk        ! bulk Richardson number (-)
 real(dp),intent(out)          :: dRiBulk_dTemp ! derivative in the bulk Richardson number w.r.t. temperature (K-1)
 integer(i4b),intent(out)      :: err           ! error code
 character(*),intent(out)      :: message       ! error message
 ! local variables
 real(dp)                      :: T_grad        ! gradient in temperature between the atmosphere and surface (K)
 real(dp)                      :: T_mean        ! mean of the atmosphere and surface temperature (K)
 real(dp)                      :: RiMult        ! dimensionless scaling factor (-)
 ! initialize error control
 err=0; message='bulkRichardson/'
 ! compute local variables
 T_grad = airtemp - surfTemp
 T_mean = 0.5_dp*(airtemp + surfTemp)
 RiMult = (gravity*mheight)/(windspd*windspd)
 ! compute the Richardson number
 RiBulk = (T_grad/T_mean) * RiMult
 ! compute the derivative in the Richardson number
 if(computeDerivative)then
  dRiBulk_dTemp = (-1._dp/T_mean + T_grad*(-0.5_dp)*T_mean**(-2._dp) ) * RiMult
 else
  dRiBulk_dTemp = 1._dp
 endif
 end subroutine bulkRichardson


 ! ***********************************************************************************************************
 ! new subroutine: compute exchange coefficient for turbulent heat transfer
 ! *********************************************************************************************************** 
 subroutine astability(RiBulk,dRiBulk_dTemp,computeDerivative, & ! input
                       ExCoef,dExCoef_dTemp,err,message)         ! output
 USE data_struc,only:mvar_data,mpar_data,model_decisions ! model structures
 USE var_lookup,only:iLookDECISIONS,iLookPARAM,iLookMVAR ! named variables for structure elements
 implicit none
 ! declare input variables
 real(dp),intent(in)      :: RiBulk            ! bulk Richardson number (-)
 real(dp),intent(in)      :: dRiBulk_dTemp     ! derivative in the bulk Richardson number w.r.t. temperature (K-1)
 logical(lgt),intent(in)  :: computeDerivative ! flag to compute the derivative
 ! declare output variables
 real(dp),intent(out)     :: ExCoef            ! surface-atmosphere exchange coefficient (-)
 real(dp),intent(out)     :: dExCoef_dTemp     ! derivative in surface-atmosphere exchange coefficient w.r.t. temperature (K-1)
 integer(i4b),intent(out) :: err               ! error code
 character(*),intent(out) :: message           ! error message
 ! declare local pointers
 real(dp),pointer         :: ExNeut            ! exchange coefficient in neutral conditions
 real(dp),pointer         :: bprime            ! parameter in the Louis (1979) stability function
 real(dp),pointer         :: Mahrt_m           ! parameter in the Mahrt (1987) stability function
 ! initialize error control
 err=0; message="astability/"
 ! assign pointers
 ExNeut  => mvar_data%var(iLookMVAR%scalarExNeut)%dat(1)
 bprime  => mvar_data%var(iLookMVAR%scalarBprime)%dat(1)
 Mahrt_m => mpar_data%var(iLookPARAM%Mahrt_m)

 ! set derivative to one if not computing it
 if(.not.computeDerivative) dExCoef_dTemp = 1._dp

 ! ***** process unstable cases
 if(RiBulk<0._dp)then
  ! compute surface-atmosphere exchange coefficient (-)
  ExCoef = ExNeut * (1._dp - 16._dp*RiBulk)**0.5_dp
  ! compute derivative in surface-atmosphere exchange coefficient w.r.t. temperature (K-1)
  if(computeDerivative)&
   dExCoef_dTemp = dRiBulk_dTemp * (-16._dp) * 0.5_dp*(1._dp - 16._dp*RiBulk)**(-0.5_dp) * ExNeut
  return
 endif

 ! ***** process stable cases
 select case(trim(model_decisions(iLookDECISIONS%astability)%decision))

  ! ("standard" stability correction, a la Anderson 1976)
  case('standard')
   ! compute surface-atmosphere exchange coefficient (-)   
   if(RiBulk<=0.2_dp) ExCoef = ExNeut * (1._dp - 5._dp*RiBulk)**2._dp
   if(RiBulk> 0.2_dp) ExCoef = 0._dp
   ! compute derivative in surface-atmosphere exchange coefficient w.r.t. temperature (K-1)
   if(computeDerivative)then
    if(RiBulk<=0.2_dp) dExCoef_dTemp = dRiBulk_dTemp * (-5._dp) * 2._dp*(1._dp - 5._dp*RiBulk) * ExNeut 
    if(RiBulk> 0.2_dp) dExCoef_dTemp = 0._dp
   endif

  ! (Louis 1979)
  case('louisinv')
   ! compute surface-atmosphere exchange coefficient (-)
   ExCoef = ExNeut / ( (1._dp + bprime*RiBulk)**2._dp )
   ! compute derivative in surface-atmosphere exchange coefficient w.r.t. temperature (K-1)
   if(computeDerivative)&
    dExCoef_dTemp = dRiBulk_dTemp * bprime * (-2._dp)*(1._dp + bprime*RiBulk)**(-3._dp) * ExNeut

  ! (Mahrt 1987)
  case('mahrtexp')
   ! compute surface-atmosphere exchange coefficient (-)
   ExCoef = ExNeut * exp(-Mahrt_m * RiBulk)
   ! compute derivative in surface-atmosphere exchange coefficient w.r.t. temperature (K-1)
   if(computeDerivative)&
    dExCoef_dTemp = dRiBulk_dTemp * (-Mahrt_m) * exp(-Mahrt_m * RiBulk) * ExNeut

  ! (return error if the stability correction method is not found)
  case default
   err=10; message=trim(message)//"optionNotFound"; return

 endselect
 end subroutine astability 

end module snow_utils_module

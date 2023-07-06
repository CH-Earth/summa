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

module tempAdjust_module

! data types
USE nrtype

! derived types to define the data structures
USE data_types,only:&
                    var_d,        & ! data vector (rkind)
                    var_dlength     ! data vector with variable length dimension (rkind)

! named variables defining elements in the data structures
USE var_lookup,only:iLookPARAM,iLookPROG,iLookDIAG  ! named variables for structure elements

! physical constants
USE multiconst,only:Tfreeze         ! freezing point of pure water (K)
USE multiconst,only:LH_fus          ! latent heat of fusion (J kg-1)
USE multiconst,only:Cp_ice          ! specific heat of ice (J kg-1 K-1)
USE multiconst,only:Cp_water        ! specific heat of liquid water (J kg-1 K-1)
USE multiconst,only:iden_water      ! intrinsic density of water (kg m-3)

! privacy
implicit none
private
public::tempAdjust

contains


 ! ************************************************************************************************
 ! public subroutine tempAdjust: compute change in snow stored on the vegetation canopy
 ! ************************************************************************************************
 subroutine tempAdjust(&
                       ! input: derived parameters
                       canopyDepth,                 & ! intent(in):    canopy depth (m)
                       ! input: derivatives
                       dCanopyIce_dWat,             & ! intent(in):    derivative of canopy ice with canopy water    
                       dCanopyIce_dTk,              & ! intent(in):    derivative of canopy ice with canopy temperature 
                       ! input/output: data structures
                       mpar_data,                   & ! intent(in):    model parameters
                       prog_data,                   & ! intent(inout): model prognostic variables for a local HRU
                       diag_data,                   & ! intent(inout): model diagnostic variables for a local HRU
                       ! output: error control
                       err,message)                   ! intent(out): error control
 ! ------------------------------------------------------------------------------------------------
 ! utility routines
 USE snow_utils_module,only:fracliquid     ! compute fraction of liquid water
 USE snow_utils_module,only:dFracLiq_dTk   ! differentiate the freezing curve w.r.t. temperature (snow)
 USE snow_utils_module,only:d2FracLiq_dTk2 ! differentiate the freezing curve twice w.r.t. temperature (snow)
 implicit none
 ! ------------------------------------------------------------------------------------------------
 ! input: derived parameters
 real(rkind),intent(in)          :: canopyDepth              ! depth of the vegetation canopy (m)
 ! input: derivatives
 real(rkind),intent(in)          :: dCanopyIce_dWat          ! derivative of canopy ice with canopy water
 real(rkind),intent(in)          :: dCanopyIce_dTk           ! derivative of canopy ice with canopy temperature 
 ! input/output: data structures  
 type(var_dlength),intent(in)    :: mpar_data                ! model parameters
 type(var_dlength),intent(inout) :: prog_data                ! model prognostic variables for a local HRU
 type(var_dlength),intent(inout) :: diag_data                ! model diagnostic variables for a local HRU
 ! output: error control  
 integer(i4b),intent(out)        :: err                      ! error code
 character(*),intent(out)        :: message                  ! error message
 ! ------------------------------------------------------------------------------------------------
 ! local variables for canopy thermodynamics
 integer(i4b)                    :: iTry                     ! trial index
 integer(i4b)                    :: iter                     ! iteration index
 integer(i4b),parameter          :: maxiter=100              ! maximum number of iterations
 real(rkind)                     :: fLiq                     ! fraction of liquid water (-)
 real(rkind)                     :: tempMin,tempMax          ! solution constraints for temperature (K)
 real(rkind)                     :: nrgMeltFreeze            ! energy required to melt-freeze the water to the current canopy temperature (J m-3)
 real(rkind)                     :: scalarCanopyWat          ! total canopy water (kg m-2)
 real(rkind)                     :: scalarCanopyIceOld       ! canopy ice content after melt-freeze to the initial temperature (kg m-2)
 real(rkind),parameter           :: resNrgToler=0.1_rkind    ! tolerance for the energy residual (J m-3)
 real(rkind)                     :: f1,f2,x1,x2,fTry,xTry,fDer,xInc ! iteration variables
 logical(lgt)                    :: fBis                     ! .true. if bisection
 ! local variables for computing derivatives
 real(rkind)                     :: dCp_dWat                 ! derivative of heat capacity with canopy water    
 real(rkind)                     :: dCp_dTk                  ! derivative of heat capacity with canopy temperature 
 real(rkind)                     :: dxTry_dWat,dx2_dWat,dtempMin_dWat,dtempMax_dWat ! derivative of iteration temperature and brackets w.r.t. canopy water   
 real(rkind)                     :: dxTry_dTk, dx2_dTk, dtempMin_dTk, dtempMax_dTk  ! derivative of iteration temperature and brackets w.r.t. canopy temperature
 real(rkind)                     :: df1_dWat,dfTry_dWat,dfDer_dWat ! derivative of iteration variables w.r.t. canopy water
 real(rkind)                     :: df1_dTk, dfTry_dTk, dfDer_dTk  ! derivative of iteration variables w.r.t. canopy temperature
 ! -------------------------------------------------------------------------------------------------------------------------------
 ! initialize error control
 err=0; message='tempAdjust/'
 ! ------------------------------------------------------------------------------------------------
 ! associate variables in the data structure
 associate(&
 ! model parameters for canopy thermodynamics (input)
 snowfrz_scale             => mpar_data%var(iLookPARAM%snowfrz_scale)%dat(1),              & ! intent(in):    [dp] scaling factor for snow freezing curve (K)
 specificHeatVeg           => mpar_data%var(iLookPARAM%specificHeatVeg)%dat(1),            & ! intent(in):    [dp] specific heat of vegetation mass (J kg-1 K-1)
 maxMassVegetation         => mpar_data%var(iLookPARAM%maxMassVegetation)%dat(1),          & ! intent(in):    [dp] maximum mass of vegetation (full foliage) (kg m-2)
 ! state variables (input/output)
 scalarCanopyLiq           => prog_data%var(iLookPROG%scalarCanopyLiq)%dat(1),             & ! intent(inout): [dp] mass of liquid water on the vegetation canopy (kg m-2)
 scalarCanopyIce           => prog_data%var(iLookPROG%scalarCanopyIce)%dat(1),             & ! intent(inout): [dp] mass of ice on the vegetation canopy (kg m-2)
 scalarCanopyTemp          => prog_data%var(iLookPROG%scalarCanopyTemp)%dat(1),            & ! intent(inout): [dp] temperature of the vegetation canopy (K)
 ! diagnostic variables (output)
 scalarBulkVolHeatCapVeg   => diag_data%var(iLookDIAG%scalarBulkVolHeatCapVeg)%dat(1),     & ! intent(out):   [dp] volumetric heat capacity of the vegetation (J m-3 K-1)
 ! canopy derivatives from adjusting the canopy temperature
 dTkCanopyAdj_dTkCanopy    => diag_data%var(iLookDIAG%dTkCanopyAdj_dTkCanopy)%dat(1),      & ! intent(in):    [dp] derivative in the adjusted temperature w.r.t. original temperature
 dTkCanopyAdj_dCanWat      => diag_data%var(iLookDIAG%dTkCanopyAdj_dCanWat)%dat(1)         & ! intent(in):    [dp] derivative in the adjusted temperature w.r.t. canopy water
 ! output: derivatives
 )  ! associate variables in the data structures
 ! -----------------------------------------------------------------------------------------------------------------------------------------------------
 ! ** preliminaries

 ! compute the total canopy water (state variable: will not change)
 scalarCanopyWat = scalarCanopyLiq + scalarCanopyIce
 
 ! compute the fraction of liquid water associated with the canopy temperature
 fLiq = fracliquid(scalarCanopyTemp,snowfrz_scale)

 ! compute the new volumetric ice content
 ! NOTE: new value; iterations will adjust this value for consistency with temperature
 scalarCanopyIceOld = (1._rkind - fLiq)*scalarCanopyWat

 ! compute volumetric heat capacity of vegetation (J m-3 K-1)
 scalarBulkVolHeatCapVeg = specificHeatVeg*maxMassVegetation/canopyDepth + & ! vegetation component
                           Cp_water*scalarCanopyLiq/canopyDepth          + & ! liquid water component
                           Cp_ice*scalarCanopyIce/canopyDepth                ! ice component

 ! heat capacity derivatives
 dCp_dWat = ( -Cp_ice*( fLiq-1._rkind ) + Cp_water*fLiq )/canopyDepth 
 if(scalarCanopyTemp < Tfreeze)then
     dCp_dTk = (-Cp_ice + Cp_water) * fLiq * scalarCanopyLiq/canopyDepth ! no derivative in air
 else
     dCp_dTk = 0._rkind
 endif

 ! compute the energy required to melt-freeze the water to the current canopy temperature (J m-3)
 nrgMeltFreeze = LH_fus*(scalarCanopyIceOld - scalarCanopyIce)/canopyDepth

 ! -----------------------------------------------------------------------------------------------------------------------------------------------------

 ! ** get ready for iterating

 ! compute initial function and derivatives
 x1   = scalarCanopyTemp
 f1   = nrgMeltFreeze
 df1_dWat = LH_fus*((1._rkind - fLiq) - dCanopyIce_dWat )/canopyDepth
 df1_dTk  = LH_fus*( -dFracLiq_dTk(scalarCanopyTemp,snowfrz_scale) - dCanopyIce_dTk )/canopyDepth
 fDer = resNrgDer(x1,scalarBulkVolHeatCapVeg,snowfrz_scale)
 dfDer_dWat = dCp_dWat + dFracLiq_dTk(x1,snowfrz_scale)*LH_fus/canopyDepth
 dfDer_dTk  = dCp_dTk + scalarCanopyWat*d2FracLiq_dTk2(x1,snowfrz_scale)*LH_fus/canopyDepth

 ! compute new function based on newton step from the first function
 x2 = x1 + f1 / fDer
 f2 = resNrgFunc(x2,scalarCanopyTemp,scalarBulkVolHeatCapVeg,snowfrz_scale)
 ! compute derivatives
 dx2_dWat = df1_dWat/fDer - f1/(dfDer_dWat**2_i4b)
 dx2_dTk  = 1._rkind + df1_dTk/fDer - f1/(dfDer_dTk**2_i4b)

 ! ensure that we bracket the root and recompute x2 and derivatives if not
 if(f1*f2 > 0._rkind)then
   xInc = f1 / fDer
   x2   = 1._rkind
   dx2_dWat = 0._rkind
   dx2_dTk  = 0._rkind
   do iter=1,maxiter
     ! successively expand limit in order to bracket the root
     x2 = x1 + sign(x2,xInc)*2._rkind
     f2 = resNrgFunc(x2,scalarCanopyTemp,scalarBulkVolHeatCapVeg,snowfrz_scale)
     ! compute derivatives
     dx2_dWat = 0._rkind
     dx2_dTk  = 1._rkind
     if(f1*f2 < 0._rkind)exit
     ! check that we bracketed the root (should get here in just a couple of expansions)
     if(iter==maxiter)then
      message=trim(message)//'unable to bracket the root'
      err=20; return
     end if
   end do ! trying to bracket the root
 end if  ! first check that we bracketed the root

 ! define initial constraints and derivatives
 if(x1 < x2)then
   tempMin = x1
   tempMax = x2
   dtempMin_dWat = 0._rkind
   dtempMin_dTk  = 1._rkind
   dtempMax_dWat = dx2_dWat
   dtempMax_dTk  = dx2_dTk
 else
   tempMin = x2
   tempMax = x1
   dtempMin_dWat = dx2_dWat
   dtempMin_dTk  = dx2_dTk
   dtempMax_dWat = 0._rkind
   dtempMax_dTk  = 1._rkind
 end if

 ! get starting trial
 xInc = huge(1._rkind)
 xTry = 0.5_rkind*(x1 + x2)
 fTry = resNrgFunc(xTry,scalarCanopyTemp,scalarBulkVolHeatCapVeg,snowfrz_scale)
 fDer = resNrgDer(xTry,scalarBulkVolHeatCapVeg,snowfrz_scale)

 ! compute derivatives
 dxTry_dWat = 0.5_rkind*(dtempMin_dWat + dtempMax_dWat)
 dxTry_dTk  = 0.5_rkind*(dtempMin_dTk  + dtempMax_dTk)
 ! xIce = (1._rkind - fracliquid(xTry,snowfrz_scale))*scalarCanopyWat
 ! fTry = -scalarBulkVolHeatCapVeg*(xTry - scalarCanopyTemp) + LH_fus*(xIce - scalarCanopyIce)/canopyDepth
 dfTry_dWat = -dCp_dWat*(xTry - scalarCanopyTemp) - scalarBulkVolHeatCapVeg*dxTry_dWat &
             + LH_fus* ( 1._rkind - fracliquid(xTry,snowfrz_scale) &
             - dFracLiq_dTk(xTry,snowfrz_scale)*scalarCanopyWat*dxTry_dWat - dCanopyIce_dWat )/canopyDepth
 dfTry_dTk  = -dCp_dTk*(xTry - scalarCanopyTemp) + scalarBulkVolHeatCapVeg*dxTry_dTk &
             + LH_fus* ( -dFracLiq_dTk(xTry,snowfrz_scale)*scalarCanopyWat*dxTry_dTk - dCanopyIce_dTk )/canopyDepth
 ! fDer = scalarBulkVolHeatCapVeg + scalarCanopyWat*dFracLiq_dTk(xTry,snowfrz_scale)*LH_fus/canopyDepth
 dfDer_dWat = dCp_dWat + ( dFracLiq_dTk(xTry,snowfrz_scale) &
             + scalarCanopyWat*d2FracLiq_dTk2(xTry,snowfrz_scale)*dxTry_dWat )*LH_fus/canopyDepth
 dfDer_dTk  = dCp_dTk + scalarCanopyWat*d2FracLiq_dTk2(xTry,snowfrz_scale)*dxTry_dTk*LH_fus/canopyDepth

 ! -----------------------------------------------------------------------------------------------------------------------------------------------------
 ! iterate
 do iter=1,maxiter
   if(xTry <= tempMin .or. xTry >= tempMax)then ! bisect if out of range
     xTry = 0.5_rkind*(tempMin + tempMax)  ! new value
     fBis = .true.
     ! compute derivatives
     dxTry_dWat = 0.5_rkind*(dtempMin_dWat + dtempMax_dWat)
     dxTry_dTk  = 0.5_rkind*(dtempMin_dTk  + dtempMax_dTk)
   else ! value in range; use the newton step
     xInc = fTry/fDer
     xTry = xTry + xInc
     fBis = .false.
     ! compute derivatives
     dxTry_dWat = dxTry_dWat+ dfTry_dWat/fDer - fTry/(dfDer_dWat**2_i4b)
     dxTry_dTk  = dxTry_dTk + dfTry_dTk/fDer - fTry/(dfDer_dTk**2_i4b)
   end if  ! (switch between bi-section and newton)

   ! compute new function and derivatives
   fTry = resNrgFunc(xTry,scalarCanopyTemp,scalarBulkVolHeatCapVeg,snowfrz_scale)
   fDer = resNrgDer(xTry,scalarBulkVolHeatCapVeg,snowfrz_scale)  
   dfTry_dWat = -dCp_dWat*(xTry - scalarCanopyTemp) - scalarBulkVolHeatCapVeg*dxTry_dWat &
               + LH_fus* ( 1._rkind - fracliquid(xTry,snowfrz_scale) &
               - dFracLiq_dTk(xTry,snowfrz_scale)*scalarCanopyWat*dxTry_dWat - dCanopyIce_dWat )/canopyDepth
   dfTry_dTk  = -dCp_dTk*(xTry - scalarCanopyTemp) + scalarBulkVolHeatCapVeg*dxTry_dTk &
               + LH_fus* ( -dFracLiq_dTk(xTry,snowfrz_scale)*scalarCanopyWat*dxTry_dTk - dCanopyIce_dTk )/canopyDepth
   dfDer_dWat = dCp_dWat + ( dFracLiq_dTk(xTry,snowfrz_scale) &
               + scalarCanopyWat*d2FracLiq_dTk2(xTry,snowfrz_scale)*dxTry_dWat )*LH_fus/canopyDepth
   dfDer_dTk  = dCp_dTk + scalarCanopyWat*d2FracLiq_dTk2(xTry,snowfrz_scale)*dxTry_dTk*LH_fus/canopyDepth

   ! update limits and derivatives
   if(fTry < 0._rkind)then
     tempMax = min(xTry,tempMax)
     if(xTry<tempMax)then
       dtempMax_dWat = dxTry_dWat
       dtempMax_dTk  = dxTry_dTk
     endif
   else
     tempMin = max(tempMin,xTry)
     if(xTry>tempMin)then
       dtempMin_dWat = dxTry_dWat
       dtempMin_dTk  = dxTry_dTk
     endif
   end if

   ! check the functions at the limits (should be of opposing sign)
   !f1 = resNrgFunc(tempMax,scalarCanopyTemp,scalarBulkVolHeatCapVeg,snowfrz_scale)
   !f2 = resNrgFunc(tempMin,scalarCanopyTemp,scalarBulkVolHeatCapVeg,snowfrz_scale)

   ! check convergence
   if(abs(fTry) < resNrgToler) exit

   ! check non-convergence
   if(iter==maxiter)then
     ! (print out a 1-d x-section)
     do iTry=1,maxiter
       xTry = 1.0_rkind*real(iTry,kind(1._rkind))/real(maxiter,kind(1._rkind)) + 272.5_rkind
       fTry = resNrgFunc(xTry,scalarCanopyTemp,scalarBulkVolHeatCapVeg,snowfrz_scale)
       write(*,'(a,1x,i4,1x,e20.10,1x,4(f20.10,1x))') 'iTry, fTry, xTry = ', iTry, fTry, xTry
     end do
     ! (return with error)
     message=trim(message)//'unable to converge'
     err=20; return
   end if
 end do  ! iterating
 ! -----------------------------------------------------------------------------------------------------------------------------------------------------

 ! update state variables
 scalarCanopyTemp = xTry
 scalarCanopyIce  = (1._rkind - fracliquid(xTry,snowfrz_scale))*scalarCanopyWat
 scalarCanopyLiq  = scalarCanopyWat - scalarCanopyIce

 ! update derivatives
 dTkCanopyAdj_dCanWat   = dxTry_dWat
 dTkCanopyAdj_dTkCanopy = dxTry_dTk

 ! update bulk heat capacity
 scalarBulkVolHeatCapVeg = specificHeatVeg*maxMassVegetation/canopyDepth + & ! vegetation component
                           Cp_water*scalarCanopyLiq/canopyDepth          + & ! liquid water component
                           Cp_ice*scalarCanopyIce/canopyDepth                ! ice component

 ! end association to variables in the data structure
 end associate

 contains


  ! ************************************************************************************************
  ! internal function resNrgFunc: calculate the residual in energy (J m-3)
  ! ************************************************************************************************
  function resNrgFunc(xTemp,xTemp0,scalarBulkVolHeatCapVeg,snowfrz_scale)
  !
  implicit none
  real(rkind),intent(in) :: xTemp              ! temperature (K)
  real(rkind),intent(in) :: xTemp0             ! initial temperature (K)
  real(rkind),intent(in) :: scalarBulkVolHeatCapVeg  ! volumetric heat capacity of veg (J m-3 K-1)
  real(rkind),intent(in) :: snowfrz_scale      ! scaling factor in freezing curve (K-1)
  real(rkind)            :: xIce               ! canopy ice content (kg m-2)
  real(rkind)            :: resNrgFunc         ! residual in energy (J m-3)
  xIce       = (1._rkind - fracliquid(xTemp,snowfrz_scale))*scalarCanopyWat
  resNrgFunc = -scalarBulkVolHeatCapVeg*(xTemp - xTemp0) + LH_fus*(xIce - scalarCanopyIceOld)/canopyDepth + nrgMeltFreeze
  return
  end function resNrgFunc

  ! ************************************************************************************************
  ! internal function resNrgDer: calculate the derivative (J m-3 K-1)
  ! ************************************************************************************************
  function resNrgDer(xTemp,scalarBulkVolHeatCapVeg,snowfrz_scale)
  implicit none
  real(rkind),intent(in) :: xTemp              ! temperature (K)
  real(rkind),intent(in) :: scalarBulkVolHeatCapVeg  ! volumetric heat capacity of veg (J m-3 K-1)
  real(rkind),intent(in) :: snowfrz_scale      ! scaling factor in freezing curve (K-1)
  real(rkind)            :: dW_dT              ! derivative in canopy ice content w.r.t. temperature (kg m-2 K-1)
  real(rkind)            :: resNrgDer          ! derivative (J m-3 K-1)
  dW_dT     = -scalarCanopyWat*dFracLiq_dTk(xTemp,snowfrz_scale)
  resNrgDer = scalarBulkVolHeatCapVeg - dW_dT*LH_fus/canopyDepth
  return
  end function resNrgDer


 end subroutine tempAdjust


end module tempAdjust_module

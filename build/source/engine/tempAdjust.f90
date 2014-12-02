module tempAdjust_module
! data types
USE nrtype
! physical constants
USE multiconst,only:Tfreeze         ! freezing point of pure water (K)
USE multiconst,only:LH_fus          ! latent heat of fusion (J kg-1)
USE multiconst,only:Cp_ice          ! specific heat of ice (J kg-1 K-1)
USE multiconst,only:Cp_water        ! specific heat of liquid water (J kg-1 K-1)
USE multiconst,only:iden_water      ! intrinsic density of water (kg m-3)
implicit none
private
public::tempAdjust

contains

 ! ************************************************************************************************
 ! new subroutine: compute change in snow stored on the vegetation canopy
 ! ************************************************************************************************
 subroutine tempAdjust(&
                       ! input: derived parameters
                       canopyDepth,                 & ! intent(in): canopy depth (m)
                       ! input/output: data structures
                       mpar_data,                   & ! intent(in):    model parameters
                       mvar_data,                   & ! intent(inout): model variables for a local HRU
                       ! output: error control
                       err,message)                   ! intent(out): error control
 ! ------------------------------------------------------------------------------------------------
 ! provide access to the derived types to define the data structures
 USE data_struc,only:&
                     var_d,              & ! data vector (dp)
                     var_dlength           ! data vector with variable length dimension (dp)
 ! provide access to named variables defining elements in the data structures
 USE var_lookup,only:iLookPARAM,iLookMVAR  ! named variables for structure elements
 ! utility routines
 USE snow_utils_module,only:fracliquid     ! compute fraction of liquid water
 USE snow_utils_module,only:dFracLiq_dTk   ! differentiate the freezing curve w.r.t. temperature (snow)
 implicit none
 ! ------------------------------------------------------------------------------------------------
 ! input: derived parameters
 real(dp),intent(in)             :: canopyDepth         ! depth of the vegetation canopy (m)
 ! input/output: data structures
 type(var_d),intent(in)          :: mpar_data           ! model parameters
 type(var_dlength),intent(inout) :: mvar_data           ! model variables for a local HRU
 ! output: error control
 integer(i4b),intent(out)        :: err                 ! error code
 character(*),intent(out)        :: message             ! error message
 ! -------------------------------------------------------------------------------------------------------------------------------
 ! variables in the data structures
 ! input: model parameters for canopy thermodynamics
 real(dp)                      :: snowfrz_scale              ! scaling factor for snow freezing curve (K)
 real(dp)                      :: specificHeatVeg            ! specific heat of vegetation mass (J kg-1 K-1)
 real(dp)                      :: maxMassVegetation          ! maximum mass of vegetation (full foliage) (kg m-2)
 ! input-output: state variables
 real(dp)                      :: scalarCanopyLiq            ! mass of liquid water on the vegetation canopy (kg m-2)
 real(dp)                      :: scalarCanopyIce            ! mass of ice on the vegetation canopy (kg m-2)
 real(dp)                      :: scalarCanopyTemp           ! temperature of the vegetation canopy (K)
 ! output: diagnostic variables
 real(dp)                      :: scalarBulkVolHeatCapVeg    ! bulk volumetric heat capacity of vegetation (J m-3 K-1)
 ! ------------------------------------------------------------------------------------------------
 ! local variables for canopy thermodynamics
 integer(i4b)                  :: iTry                       ! trial index
 integer(i4b)                  :: iter                       ! iteration index
 integer(i4b),parameter        :: maxiter=100                ! maximum number of iterations
 real(dp),parameter            :: dx=1.e-6_dp                ! finite difference increment (used to test derivatives)
 real(dp)                      :: fLiq                       ! fraction of liquid water (-)
 real(dp)                      :: dW_dT                      ! derivative in canopy ice content w.r.t canopy temperature (kg m-2 K-1)
 real(dp)                      :: tempMin,tempMax            ! solution constraints for temperature (K)
 real(dp)                      :: nrgMeltFreeze              ! energy required to melt-freeze the water to the current canopy temperature (J m-3)
 real(dp)                      :: scalarCanopyWat            ! total canopy water (kg m-2)
 real(dp)                      :: scalarCanopyIceOld         ! canopy ice content after melt-freeze to the initial temperature (kg m-2)
 real(dp)                      :: scalarCanopyIceIter        ! trial value for canopy ice content (kg m-2)
 real(dp)                      :: scalarCanopyTempIter       ! trial value for canopy temperature (K)
 real(dp)                      :: scalarCanopyTempTrial      ! trial value for canopy temperature, before update (K)
 real(dp)                      :: resNrg,resNrgOld           ! energy residual (J m-3)
 real(dp),parameter            :: resNrgToler=0.1_dp         ! tolerance for the energy residual (J m-3)
 real(dp)                      :: delTemp                    ! iteration increment for temperature
 real(dp)                      :: xLambda                    ! scaling factor for the iteration increment (-)
 real(dp)                      :: adjTemp                    ! adjusted iteration increment (K)
 real(dp)                      :: f1,f2,x1,x2,fTry,xTry,fDer,xInc
 logical(lgt) :: fBis  ! .true. if bisection
 real(dp) :: term1, term2
 ! -------------------------------------------------------------------------------------------------------------------------------
 ! initialize error control
 err=0; message='tempAdjust/'
 ! ------------------------------------------------------------------------------------------------
 ! associate variables in the data structure
 associate(&

 ! model parameters for canopy thermodynamics (input)
 snowfrz_scale             => mpar_data%var(iLookPARAM%snowfrz_scale),                     & ! intent(in): [dp] scaling factor for snow freezing curve (K)
 specificHeatVeg           => mpar_data%var(iLookPARAM%specificHeatVeg),                   & ! intent(in): [dp] specific heat of vegetation mass (J kg-1 K-1)
 maxMassVegetation         => mpar_data%var(iLookPARAM%maxMassVegetation),                 & ! intent(in): [dp] maximum mass of vegetation (full foliage) (kg m-2)

 ! state variables (input/output)
 scalarCanopyLiq           => mvar_data%var(iLookMVAR%scalarCanopyLiq)%dat(1),             & ! intent(inout): [dp] mass of liquid water on the vegetation canopy (kg m-2)
 scalarCanopyIce           => mvar_data%var(iLookMVAR%scalarCanopyIce)%dat(1),             & ! intent(inout): [dp] mass of ice on the vegetation canopy (kg m-2)
 scalarCanopyTemp          => mvar_data%var(iLookMVAR%scalarCanopyTemp)%dat(1),            & ! intent(inout): [dp] temperature of the vegetation canopy (K)
 
 ! diagnostic variables (output)
 scalarBulkVolHeatCapVeg   => mvar_data%var(iLookMVAR%scalarBulkVolHeatCapVeg)%dat(1)      & ! intent(out): [dp] volumetric heat capacity of the vegetation (J m-3 K-1)

 )  ! associate variables in the data structures
 ! -----------------------------------------------------------------------------------------------------------------------------------------------------

 ! ** preliminaries

 ! compute the total canopy water (state variable: will not change)
 scalarCanopyWat = scalarCanopyLiq + scalarCanopyIce
 !write(*,'(a,1x,3(f20.10,1x))') 'scalarCanopyWat, scalarCanopyLiq, scalarCanopyIce = ', scalarCanopyWat, scalarCanopyLiq, scalarCanopyIce

 ! compute the fraction of liquid water associated with the canopy temperature
 fLiq = fracliquid(scalarCanopyTemp,snowfrz_scale)

 ! compute the new volumetric ice content
 ! NOTE: new value; iterations will adjust this value for consistency with temperature
 scalarCanopyIceOld = (1._dp - fLiq)*scalarCanopyWat

 ! compute volumetric heat capacity of vegetation (J m-3 K-1)
 scalarBulkVolHeatCapVeg = specificHeatVeg*maxMassVegetation/canopyDepth + & ! vegetation component
                           Cp_water*scalarCanopyLiq/canopyDepth          + & ! liquid water component
                           Cp_ice*scalarCanopyIce/canopyDepth                ! ice component

 ! compute the energy required to melt-freeze the water to the current canopy temperature (J m-3)
 nrgMeltFreeze = LH_fus*(scalarCanopyIceOld - scalarCanopyIce)/canopyDepth

 ! -----------------------------------------------------------------------------------------------------------------------------------------------------

 ! ** get ready for iterating

 ! compute initial function and derivative
 x1   = scalarCanopyTemp
 f1   = nrgMeltFreeze
 fDer = resNrgDer(x1,scalarBulkVolHeatCapVeg,snowfrz_scale)

 ! compute new function based on newton step from the first function
 x2 = x1 + f1 / fDer
 f2 = resNrgFunc(x2,scalarCanopyTemp,scalarBulkVolHeatCapVeg,snowfrz_scale)
 !print*, 'x1, x2 = ', x1, x2
 !print*, 'f1, f2 = ', f1, f2

 ! ensure that we bracket the root
 if(f1*f2 > 0._dp)then
  xInc = f1 / fDer
  x2   = 1._dp
  do iter=1,maxiter
   ! successively expand limit in order to bracket the root
   x2 = x1 + sign(x2,xInc)*2._dp
   f2 = resNrgFunc(x2,scalarCanopyTemp,scalarBulkVolHeatCapVeg,snowfrz_scale)
   if(f1*f2 < 0._dp)exit
   ! check that we bracketed the root
   ! (should get here in just a couple of expansions)
   if(iter==maxiter)then
    message=trim(message)//'unable to bracket the root'
    err=20; return
   endif
  end do ! trying to bracket the root
 endif  ! first check that we bracketed the root
 !print*, 'x1, x2 = ', x1, x2
 !print*, 'f1, f2 = ', f1, f2

 ! define initial constraints
 if(x1 < x2)then
  tempMin = x1
  tempMax = x2
 else
  tempMin = x2
  tempMax = x1
 endif
 !print*, 'tempMin, tempMax = ', tempMin, tempMax

 ! get starting trial
 xInc = huge(1._dp)
 xTry = 0.5_dp*(x1 + x2)
 fTry = resNrgFunc(xTry,scalarCanopyTemp,scalarBulkVolHeatCapVeg,snowfrz_scale)
 fDer = resNrgDer(xTry,scalarBulkVolHeatCapVeg,snowfrz_scale)
 !print*, 'xTry = ', xTry
 !print*, 'fTry = ', fTry

 ! check the functions at the limits (should be of opposing sign)
 !f1 = resNrgFunc(tempMax,scalarCanopyTemp,scalarBulkVolHeatCapVeg,snowfrz_scale)
 !f2 = resNrgFunc(tempMin,scalarCanopyTemp,scalarBulkVolHeatCapVeg,snowfrz_scale)
 !print*, 'f1, f2 = ', f1, f2

 ! -----------------------------------------------------------------------------------------------------------------------------------------------------
 ! iterate
 do iter=1,maxiter

  ! bisect if out of range
  if(xTry <= tempMin .or. xTry >= tempMax)then
   xTry = 0.5_dp*(tempMin + tempMax)  ! new value
   fBis = .true.

  ! value in range; use the newton step
  else
   xInc = fTry/fDer
   xTry = xTry + xInc
   fBis = .false.

  endif  ! (switch between bi-section and newton)

  ! compute new function and derivative
  fTry = resNrgFunc(xTry,scalarCanopyTemp,scalarBulkVolHeatCapVeg,snowfrz_scale)
  fDer = resNrgDer(xTry,scalarBulkVolHeatCapVeg,snowfrz_scale)
  !print*, 'tempMin, tempMax = ', tempMin, tempMax

  ! update limits
  if(fTry < 0._dp)then
   tempMax = min(xTry,tempMax)
  else
   tempMin = max(tempMin,xTry)
  endif

  ! check the functions at the limits (should be of opposing sign)
  !f1 = resNrgFunc(tempMax,scalarCanopyTemp,scalarBulkVolHeatCapVeg,snowfrz_scale)
  !f2 = resNrgFunc(tempMin,scalarCanopyTemp,scalarBulkVolHeatCapVeg,snowfrz_scale)
  !print*, 'f1, f2 = ', f1, f2

  ! print progress
  !write(*,'(a,1x,i4,1x,l1,1x,e20.10,1x,4(f20.10,1x))') 'iter, fBis, fTry, xTry, xInc, tempMin, tempMax = ', iter, fBis, fTry, xTry, xInc, tempMin, tempMax

  ! check convergence
  if(abs(fTry) < resNrgToler) exit 

  ! check non-convergence
  if(iter==maxiter)then
   ! (print out a 1-d x-section)
   do iTry=1,maxiter
    xTry = 1.0_dp*real(iTry,kind(1._dp))/real(maxiter,kind(1._dp)) + 272.5_dp
    fTry = resNrgFunc(xTry,scalarCanopyTemp,scalarBulkVolHeatCapVeg,snowfrz_scale)
    write(*,'(a,1x,i4,1x,e20.10,1x,4(f20.10,1x))') 'iTry, fTry, xTry = ', iTry, fTry, xTry
   end do
   ! (return with error)
   message=trim(message)//'unable to converge'
   err=20; return
  endif

 end do  ! iterating
 ! -----------------------------------------------------------------------------------------------------------------------------------------------------

 ! update state variables
 scalarCanopyTemp = xTry 
 scalarCanopyIce  = (1._dp - fracliquid(xTry,snowfrz_scale))*scalarCanopyWat
 scalarCanopyLiq  = scalarCanopyWat - scalarCanopyIce

 ! end association to variables in the data structure
 end associate


 contains

  ! ** internal functions

  ! calculate function
  function resNrgFunc(xTemp,xTemp0,bulkVolHeatCapVeg,snowfrz_scale)
  ! calculate the residual in energy (J m-3)
  implicit none
  real(dp),intent(in) :: xTemp              ! temperature (K)
  real(dp),intent(in) :: xTemp0             ! initial temperature (K)
  real(dp),intent(in) :: bulkVolHeatCapVeg  ! volumetric heat capacity of veg (J m-3 K-1)
  real(dp),intent(in) :: snowfrz_scale      ! scaling factor in freezing curve (K-1)
  real(dp)            :: xIce               ! canopy ice content (kg m-2)
  real(dp)            :: resNrgFunc         ! residual in energy (J m-3)
  xIce       = (1._dp - fracliquid(xTemp,snowfrz_scale))*scalarCanopyWat
  resNrgFunc = -bulkVolHeatCapVeg*(xTemp - xTemp0) + LH_fus*(xIce - scalarCanopyIceOld)/canopyDepth + nrgMeltFreeze
  return
  end function resNrgFunc

  ! calculate derivative
  function resNrgDer(xTemp,bulkVolHeatCapVeg,snowfrz_scale)
  ! calculate the derivatve (J m-3 K-1)
  implicit none
  real(dp),intent(in) :: xTemp              ! temperature (K)
  real(dp),intent(in) :: bulkVolHeatCapVeg  ! volumetric heat capacity of veg (J m-3 K-1)
  real(dp),intent(in) :: snowfrz_scale      ! scaling factor in freezing curve (K-1)
  real(dp)            :: dW_dT              ! derivative in canopy ice content w.r.t. temperature (kg m-2 K-1)
  real(dp)            :: resNrgDer          ! derivative (J m-3 K-1)
  dW_dT     = -scalarCanopyWat*dFracLiq_dTk(xTemp,snowfrz_scale)
  resNrgDer = bulkVolHeatCapVeg - dW_dT*LH_fus/canopyDepth
  return
  end function resNrgDer

 end subroutine tempAdjust


end module tempAdjust_module

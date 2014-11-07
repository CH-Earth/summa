module groundwatr_module
! numerical recipes data types
USE nrtype
! model constants
USE multiconst,only:iden_water ! density of water (kg m-3)
! access the number of snow and soil layers
USE data_struc,only:&
                    nSnow,   & ! number of snow layers  
                    nSoil,   & ! number of soil layers  
                    nLayers    ! total number of layers
implicit none
! constant parameters
real(dp),parameter     :: valueMissing=-9999._dp    ! missing value parameter
real(dp),parameter     :: verySmall=epsilon(1.0_dp) ! a very small number (used to avoid divide by zero)
real(dp),parameter     :: dx=1.e-8_dp               ! finite difference increment
private
public::groundwatr
contains

 ! ************************************************************************************************
 ! new subroutine: compute the groundwater sink term in Richards' equation
 ! ************************************************************************************************
 subroutine groundwatr(&

                       ! input: model control
                       dt,                                     & ! intent(in): length of the model time step
                       mLayerHydCond,                          & ! intent(in): hydraulic conductivity in each soil layer (m s-1)
                       dHydCond_dMatric,                       & ! intent(in): derivative in hydraulic conductivity w.r.t. matric head in each layer (s-1)
                       mLayerdTheta_dPsi,                      & ! intent(in): derivative in the soil water characteristic w.r.t. matric head in each layer (m-1)
                       mLayerMatricHeadLiq,                    & ! intent(in): liquid water matric potential (m)
                       mLayerVolFracLiq,                       & ! intent(in): volumetric fraction of liquid water (-)
                       mLayerVolFracIce,                       & ! intent(in): volumetric fraction of ice (-)

                       ! input/output: data structures
                       attr_data,                              & ! intent(in):    spatial attributes
                       mpar_data,                              & ! intent(in):    model parameters
                       mvar_data,                              & ! intent(inout): model variables for a local HRU

                       ! output: baseflow
                       mLayerBaseflow,                         & ! intent(out): baseflow from each soil layer (m s-1)
                       dBaseflow_dMatric,                      & ! intent(out): derivative in baseflow w.r.t. matric head (s-1)

                       ! output: error control
                       err,message)                              ! intent(out): error control
 ! ---------------------------------------------------------------------------------------
 ! provide access to the derived types to define the data structures
 USE data_struc,only:&
                     var_d,            & ! data vector (dp)
                     var_dlength         ! data vector with variable length dimension (dp)
 ! provide access to named variables defining elements in the data structures
 USE var_lookup,only:iLookATTR,iLookPARAM,iLookMVAR              ! named variables for structure elements
 ! utility modules
 USE soil_utils_module,only:volFracLiq          ! compute volumetric fraction of liquid water as a function of matric head
 USE soil_utils_module,only:hydCond_psi         ! compute hydraulic conductivity as a function of matric head
 implicit none
 ! ---------------------------------------------------------------------------------------
 ! * dummy variables
 ! ---------------------------------------------------------------------------------------
 ! model control
 real(dp),intent(in)              :: dt                           ! length of model time step (s)
 real(dp),intent(in)              :: mLayerHydCond(:)             ! hydraulic conductivity in each soil layer (m s-1)
 real(dp),intent(in)              :: dHydCond_dMatric(:)          ! derivative in hydraulic conductivity w.r.t matric head (s-1)
 real(dp),intent(in)              :: mLayerdTheta_dPsi(:)         ! derivative in the soil water characteristic w.r.t. matric head in each layer (m-1)
 real(dp),intent(in)              :: mLayerMatricHeadLiq(:)       ! matric head in each layer at the current iteration (m)
 real(dp),intent(in)              :: mLayerVolFracLiq(:)          ! volumetric fraction of liquid water (-)
 real(dp),intent(in)              :: mLayerVolFracIce(:)          ! volumetric fraction of ice (-)
 ! input/output: data structures
 type(var_d),intent(in)           :: attr_data                    ! spatial attributes
 type(var_d),intent(in)           :: mpar_data                    ! model parameters
 type(var_dlength),intent(inout)  :: mvar_data                    ! model variables for a local HRU
 ! output: baseflow
 real(dp),intent(out)             :: mLayerBaseflow(:)            ! baseflow from each soil layer (m s-1)
 real(dp),intent(out)             :: dBaseflow_dMatric(:,:)       ! derivative in baseflow w.r.t. matric head (s-1)
 ! output: error control
 integer(i4b),intent(out)         :: err                          ! error code
 character(*),intent(out)         :: message                      ! error message
 ! ---------------------------------------------------------------------------------------
 ! * variables in the data structures
 ! ---------------------------------------------------------------------------------------
 ! input: coordinate variables
 real(dp),dimension(nSoil)        :: mLayerDepth                  ! intent(in): depth of each soil layer (m)
 ! input: diagnostic variables
 real(dp)                         :: mLayerSatHydCondMP           ! intent(in): saturated hydraulic conductivity at the surface (m s-1)
 real(dp),dimension(nSoil)        :: mLayerColumnInflow           ! intent(in): inflow into each soil layer (m3/s)
 ! input: local attributes
 real(dp)                         :: HRUarea                      ! intent(in): HRU area (m2)
 real(dp)                         :: tan_slope                    ! intent(in): tan water table slope, taken as tan local ground surface slope (-)
 real(dp)                         :: contourLength                ! intent(in): length of contour at downslope edge of HRU (m)
 ! input: baseflow parameters
 real(dp)                         :: zScale_TOPMODEL              ! intent(in): TOPMODEL exponent (-)
 real(dp)                         :: kAnisotropic                 ! intent(in): anisotropy factor for lateral hydraulic conductivity (-)
 real(dp)                         :: fieldCapacity                ! intent(in): field capacity (-)
 real(dp)                         :: theta_sat                    ! intent(in): soil porosity (-) 
 real(dp)                         :: theta_res                    ! intent(in): residual volumetric water content (-) 
 ! input: van Genuchten soil parameters
 real(dp)                         :: vGn_alpha,vGn_n,vGn_m         ! van Genuchten parameters
 ! output: diagnostic variables
 real(dp)                         :: scalarExfiltration           ! intent(out): exfiltration from the soil profile (m s-1)
 real(dp),dimension(nSoil)        :: mLayerColumnOutflow          ! intent(out): column outflow from each soil layer (m3 s-1)
 ! ---------------------------------------------------------------------------------------
 ! * local variables
 ! ---------------------------------------------------------------------------------------
 ! general local variables
 character(LEN=256)              :: cmessage                      ! error message of downwind routine
 integer(i4b)                    :: iter                          ! iteration index
 integer(i4b),parameter          :: maxiter=100                   ! maximum number of iterations
 real(dp),parameter              :: checkSum=1.e-10_dp            ! tolerance to check that sum of fractional allocation sums to one
 ! local variables for the lateral flux among soil columns
 integer(i4b)                    :: ixSaturation                  ! index of the lowest saturated layer
 real(dp)                        :: subSurfaceStorage             ! sub surface storage (m)
 real(dp)                        :: subSurfaceStorageTrial        ! sub surface storage (m)
 real(dp)                        :: maximumSoilWater              ! maximum storage (m)
 real(dp)                        :: maximumFlowRate               ! flow rate under saturated conditions (m/s)
 real(dp)                        :: xFreeIce                      ! mass of ice in the "free" pore space (m)
 real(dp)                        :: totalColumnInflow             ! total column inflow (m/s)
 real(dp)                        :: totalColumnOutflow            ! total outflow from the soil column (m s-1)
 real(dp)                        :: totalOutflowDeriv             ! derivative in total outflow w.r.t. storage (s-1)
 real(dp)                        :: storageRes                    ! storage residual (m)
 real(dp)                        :: storageInc                    ! storage iteration increment (m)
 real(dp),parameter              :: tolRes=1.e-8_dp               ! convergence tolerance for the residual
 real(dp),parameter              :: tolInc=1.e-10_dp              ! convergence tolerance for the iteration increment
 real(dp)                        :: volFracLiq0,volFracLiq1

 real(qp)                        :: sumDepthAvgCond               ! total transmissivity for the "active" portion of the soil profile (m2 s-1)
 real(dp),dimension(nSoil)       :: fracTotalOutflow              ! fraction of outflow apportioned to each layer

 real(qp)                        :: fPart1,fPart2                 ! part of a function
 real(qp)                        :: dPart1,dPart2                 ! derivatives for part of a function
 real(dp),dimension(nSoil)       :: dOutflow_dMatric              ! derivative in total column outflow w.r.t matric head in each layer (s-1)
 real(dp),dimension(nSoil,nSoil) :: dFracTotalOutflow_dMatric     ! derivative in fraction of outflow apportioned to each layer w.r.t matric head in each layer (m-1)

 real(qp)                     :: hydCond1
 real(qp)                     :: sumDepthAvgCond1
 real(qp),dimension(nSoil)    :: mLayerHydCondCopy
 real(dp),pointer             :: scalarSatHydCond

 real(dp)                        :: xPore,xTemp,xTran,xFlow

 real(qp)                        :: f0,f1
 integer(i4b)                    :: iLayer,jLayer
 real(dp),dimension(nSoil)       :: mLayerVolFracLiqCopy  
 real(dp)                        :: subSurfaceStorage1
 real(dp)                        :: totalColumnOutflow1
 ! ***************************************************************************************
 ! ***************************************************************************************
 ! initialize error control
 err=0; message='groundwatr/'
 ! ---------------------------------------------------------------------------------------
 ! ---------------------------------------------------------------------------------------
 ! associate variables in data structures
 associate(&

 ! input: coordinate variables
 mLayerDepth             => mvar_data%var(iLookMVAR%mLayerDepth)%dat(1:nSoil),      & ! intent(in): [dp(:)] depth of each soil layer (m)
 
 ! input: diagnostic variables
 mLayerSatHydCondMP      => mvar_data%var(iLookMVAR%mLayerSatHydCondMP)%dat(1),     & ! intent(in): [dp]    saturated hydraulic conductivity at the surface (m s-1)
 mLayerColumnInflow      => mvar_data%var(iLookMVAR%mLayerColumnInflow)%dat,        & ! intent(in): [dp(:)] inflow into each soil layer (m3/s)

 ! input: local attributes
 HRUarea                 => attr_data%var(iLookATTR%HRUarea),                       & ! intent(in): [dp] HRU area (m2)
 tan_slope               => attr_data%var(iLookATTR%tan_slope),                     & ! intent(in): [dp] tan water table slope, taken as tan local ground surface slope (-)
 contourLength           => attr_data%var(iLookATTR%contourLength),                 & ! intent(in): [dp] length of contour at downslope edge of HRU (m)

 ! input: baseflow parameters
 zScale_TOPMODEL         => mpar_data%var(iLookPARAM%zScale_TOPMODEL),              & ! intent(in): [dp] TOPMODEL exponent (-)
 kAnisotropic            => mpar_data%var(iLookPARAM%kAnisotropic),                 & ! intent(in): [dp] anisotropy factor for lateral hydraulic conductivity (-
 fieldCapacity           => mpar_data%var(iLookPARAM%fieldCapacity),                & ! intent(in): [dp] field capacity (-)
 theta_sat               => mpar_data%var(iLookPARAM%theta_sat),                    & ! intent(in): [dp] soil porosity (-)
 theta_res               => mpar_data%var(iLookPARAM%theta_res),                    & ! intent(in): [dp] residual volumetric water content (-)

 ! input: van Genuchten soil parametrers
 vGn_alpha               => mpar_data%var(iLookPARAM%vGn_alpha),                    & ! intent(in): [dp] van Genutchen "alpha" parameter (m-1)
 vGn_n                   => mpar_data%var(iLookPARAM%vGn_n),                        & ! intent(in): [dp] van Genutchen "n" parameter (-)
 vGn_m                   => mvar_data%var(iLookMVAR%scalarVGn_m)%dat(1),            & ! intent(in): [dp] van Genutchen "m" parameter (-)

 ! output: diagnostic variables
 scalarExfiltration      => mvar_data%var(iLookMVAR%scalarExfiltration)%dat(1),     & ! intent(out):[dp]    exfiltration from the soil profile (m s-1)
 mLayerColumnOutflow     => mvar_data%var(iLookMVAR%mLayerColumnOutflow)%dat        & ! intent(out):[dp(:)] column outflow from each soil layer (m3 s-1)

 )  ! end association to variables in data structures

 ! ************************************************************************************************
 ! (1) estimate total column outflow and its derivative w.r.t. matric head in each layer
 ! ************************************************************************************************

 ! compute the total column inflow (m/s)
 totalColumnInflow = sum(mLayerColumnInflow(:))/HRUarea

 ! get index of the lowest saturated layer
 ixSaturation = nSoil+1  ! unsaturated profile when ixSaturation>nSoil
 do iLayer=nSoil,1,-1  ! start at the lowest soil layer and work upwards to the top layer
  if(mLayerVolFracLiq(iLayer) > fieldCapacity)then; ixSaturation = iLayer  ! index of saturated layer -- keeps getting over-written as move upwards
  else; exit; endif                                                        ! (only consider saturated layer at the bottom of the soil profile) 
 end do  ! (looping through soil layers)

 !write(*,'(a,i4,1x,10(f20.10,1x))') 'ixSaturation, mLayerVolFracLiq(:) = ', ixSaturation, mLayerVolFracLiq(:)



 !print*, 'HRUarea            = ', HRUarea
 !print*, 'tan_slope          = ', tan_slope
 !print*, 'mLayerSatHydCondMP = ', mLayerSatHydCondMP
 !print*, 'kAnisotropic       = ', kAnisotropic
 !print*, 'contourLength      = ', contourLength

 ! if the profile is active compute storage (m) and the depth-weighted hydraulic conductivity, i.e, the transmissivity (m2 s-1)
 if(ixSaturation <= nSoil)then
  subSurfaceStorage = sum((mLayerVolFracLiq(ixSaturation:nSoil) - fieldCapacity)*mLayerDepth(ixSaturation:nSoil))
  sumDepthAvgCond   = sum(mLayerHydCond(ixSaturation:nSoil)*mLayerDepth(ixSaturation:nSoil))
  !write(*,'(a,1x,10(e20.10,1x))') 'sumDepthAvgCond, verySmall = ', sumDepthAvgCond, verySmall
 else
  subSurfaceStorage = 0._dp
  sumDepthAvgCond   = 0._dp
 endif

 ! check for an early return
 if(ixSaturation > nSoil .or. sumDepthAvgCond < verySmall)then
  scalarExfiltration     = 0._dp   ! exfiltration from the soil profile (m s-1)
  mLayerColumnOutflow(:) = 0._dp   ! column outflow from each soil layer (m3 s-1)
  mLayerBaseflow(:)      = 0._dp   ! baseflow from each soil layer (m s-1)
  dBaseflow_dMatric(:,:) = 0._dp   ! derivative in baseflow w.r.t. matric head (s-1)
  return
 endif  ! if some layers are saturated

 ! compute the mass of ice in the "free" pore space
 ! NOTE: currently disabled to simplify the numerical solution
 xFreeIce = 0._dp  ! mass of ice in the "free" pore space (m)
 !forall(iLayer=1:nSoil, mLayerVolFracIce(iLayer) > fieldCapacity) xFreeIce = xFreeIce + (mLayerVolFracIce(iLayer) - fieldCapacity)*mLayerDepth(iLayer)

 ! compute maximum possible sub-surface free storage (m)
 ! NOTE: only need to calculate at the start of a sub-step
 maximumSoilWater = sum(mLayerDepth(1:nSoil))*(theta_sat - fieldCapacity) - xFreeIce    ! maximum aquifer storage (m)

 ! compute the flow rate under saturated conditions (m s-1)
 ! NOTE: only need to calculate at the start of a sub-step
 maximumFlowRate  = (1._dp/HRUarea)*tan_slope*mLayerSatHydCondMP*kAnisotropic*maximumSoilWater*contourLength &
                       / ((theta_sat - fieldCapacity)*zScale_TOPMODEL)   ! effective hydraulic conductivity (m/s)

 ! compute maximum transmissivity (m2 s-1)
 !xPore = theta_sat - fieldCapacity
 !xTemp = mLayerSatHydCondMP*(maximumSoilWater/xPore)/zScale_TOPMODEL

 ! compute transmissivity (m2 s-1)
 !xTran = xTemp * (subSurfaceStorage/maximumSoilWater)**zScale_TOPMODEL

 ! compute outflow (m3 s-1)
 !xFlow = tan_slope*contourLength*xTran

 !print*, 'xPore = ', xPore
 !print*, 'xTemp = ', xTemp
 !print*, 'xTran = ', xTran
 !write(*,'(a,1x,f20.10)') 'xFlow         = ', xFlow
 !write(*,'(a,1x,f20.10)') 'xFlow/HRUarea = ', xFlow/HRUarea


 !print*, 'maximumFlowRate = ', maximumFlowRate

 ! compute total column outflow
 totalColumnOutflow = maximumFlowRate*(subSurfaceStorage/maximumSoilWater)**zScale_TOPMODEL  ! m s-1
 totalOutflowDeriv  = (maximumFlowRate/maximumSoilWater)*zScale_TOPMODEL*(subSurfaceStorage/maximumSoilWater)**(zScale_TOPMODEL - 1._dp)  ! s-1
 !write(*,'(a,1x,10(f20.10,1x))') 'maximumSoilWater, subSurfaceStorage, totalColumnOutflow = ', &
 !                                 maximumSoilWater, subSurfaceStorage, totalColumnOutflow

 ! compute derivative in total column outflow w.r.t. matric head in each layer (s-1)
 dOutflow_dMatric(ixSaturation:nSoil) = totalOutflowDeriv*mLayerdTheta_dPsi(ixSaturation:nSoil)*mLayerDepth(ixSaturation:nSoil)
 if(ixSaturation>1) dOutflow_dMatric(1:ixSaturation-1) = 0._dp
 !print*, 'dOutflow_dMatric(:) = ', dOutflow_dMatric(:)

 ! **************************************************************************************************************
 ! (2) compute the fraction of total outflow apportioned to each layer, and its derivative w.r.t. matric head
 ! **************************************************************************************************************

 ! check that the depth average hydraulic conductivity is positive when baseflow occurs
 if(sumDepthAvgCond < tiny(theta_sat))then
  write(*,'(a,1x,10(e20.10,1x))') 'mLayerHydCond(ixSaturation:nSoil)       = ', mLayerHydCond(ixSaturation:nSoil)
  write(*,'(a,1x,10(e20.10,1x))') 'mLayerDepth(ixSaturation:nSoil)         = ', mLayerDepth(ixSaturation:nSoil)
  message=trim(message)//'zero depth-weighted conductivity when baseflow occurs'
  err=20; return
 endif

 ! loop through active soil layers
 do iLayer=ixSaturation,nSoil

  ! compute the function and derivative for the numerator
  fPart1 = mLayerHydCond(iLayer)*mLayerDepth(iLayer)        ! numerator = transmissivity of layer (m2 s-1)
  dPart1 = dHydCond_dMatric(iLayer)*mLayerDepth(iLayer)     ! derivative of the numerator (m s-1) -- needed for the diagonal terms

  ! compute fraction of total outflow apportioned to each layer
  fracTotalOutflow(iLayer) = fPart1/sumDepthAvgCond

  ! loop through active soil layers
  do jLayer=ixSaturation,nSoil

   ! compute the derivative for the denominator
   dPart2 = -mLayerDepth(jLayer)*dHydCond_dMatric(jLayer)/sumDepthAvgCond**2._qp ! derivative for the denominator (s m-3)

   ! diagonal terms
   if(iLayer==jLayer)then
    dFracTotalOutflow_dMatric(iLayer,jLayer) = fPart1*dPart2 + dPart1/sumDepthAvgCond  ! m-1

   ! x-derivative terms
   else
    dFracTotalOutflow_dMatric(iLayer,jLayer) = fPart1*dPart2

   end if ! (switch between diagonal and x-derivative terms)

  end do  ! looping through active soil layers
 end do  ! looping through active soil layers


 ! check derivatives
 !iLayer = 1
 !scalarSatHydCond => mvar_data%var(iLookMVAR%mLayerSatHydCond)%dat(iLayer)
 !hydCond1 = hydcond_psi(mLayerMatricHeadLiq(iLayer)+dx,scalarSatHydCond,vGn_alpha,vGn_n,vGn_m) 
 !f0 = mLayerHydCond(iLayer)
 !f1 = hydCond1
 !print*, 'f0 = ', f0
 !print*, 'f1 = ', f1
 !write(*,'(a,1x,10(e20.10,1x))') '(f1 - f0)/dx             = ', (f1 - f0)/dx
 !write(*,'(a,1x,10(e20.10,1x))') 'dHydCond_dMatric(iLayer) = ', dHydCond_dMatric(iLayer)   

 !mLayerHydCondCopy(:) = mLayerHydCond(:) 
 !mLayerHydCondCopy(iLayer) = hydCond1
 !sumDepthAvgCond1 = sum(mLayerHydCondCopy(ixSaturation:nSoil)*mLayerDepth(ixSaturation:nSoil))
 !f0 = mLayerHydCond(iLayer)*mLayerDepth(iLayer)/sumDepthAvgCond
 !f1 = hydCond1*mLayerDepth(iLayer)/sumDepthAvgCond1
 !write(*,'(a,1x,10(e20.10,1x))') '(f1 - f0)/dx                             = ', (f1 - f0)/dx
 !write(*,'(a,1x,10(e20.10,1x))') 'dFracTotalOutflow_dMatric(iLayer,iLayer) = ', dFracTotalOutflow_dMatric(iLayer,iLayer)
 !pause

 ! set unsaturated elements to zero
 if(ixSaturation > 1)then
  fracTotalOutflow(1:ixSaturation-1) = 0._dp
  dFracTotalOutflow_dMatric(1:ixSaturation-1,1:ixSaturation-1) = 0._dp
 endif

 ! check that fraction of baseflow apportioned to the different soil layers sums to 1
 if(abs(1._dp - sum(fracTotalOutflow)) > checkSum)then
  write(*,'(a,1x,10(f30.20,1x))') 'fracTotalOutflow      = ', fracTotalOutflow
  write(*,'(a,1x,10(f30.20,1x))') 'sum(fracTotalOutflow) = ', sum(fracTotalOutflow)
  message=trim(message)//'fraction of baseflow does not sum to 1'
  err=20; return
 endif

 ! check
 !write(*,'(a,1x,10(e20.10,1x))') 'fracTotalOutflow(ixSaturation:nSoil) = ', fracTotalOutflow(ixSaturation:nSoil)
 

 ! **************************************************************************************************************
 ! (3) compute baseflow, and its derivative w.r.t. matric head
 ! **************************************************************************************************************

 ! compute the outflow from each soil layer (m3 s-1)
 mLayerColumnOutflow(1:nSoil) = fracTotalOutflow(1:nSoil)*totalColumnOutflow*HRUarea

 ! compute the net baseflow from each soil layer (m s-1)
 mLayerBaseflow(1:nSoil) = (mLayerColumnOutflow(1:nSoil) - mLayerColumnInflow(1:nSoil))/HRUarea

 ! compute the derivative in baseflow w.r.t matric head (s-1)
 do iLayer=ixSaturation,nSoil
  do jLayer=ixSaturation,nSoil

   ! diagonal elements
   if(iLayer==jLayer)then
    dBaseflow_dMatric(iLayer,jLayer) = fracTotalOutflow(iLayer)*dOutflow_dMatric(iLayer) + totalColumnOutflow*dFracTotalOutflow_dMatric(iLayer,jLayer)

   ! off-diagonal elements
   else
    dBaseflow_dMatric(iLayer,jLayer) = totalColumnOutflow*dFracTotalOutflow_dMatric(iLayer,jLayer)

   endif   ! (switch between diagonal and off-diagonal elements

  end do  ! (looping through "active" soil layers)
 end do  ! (looping through "active" soil layers)


 ! start with zero exfiltration
 scalarExfiltration = 0._dp

 ! add exfiltration to the baseflow flux at the top layer
 mLayerBaseflow(1)      = mLayerBaseflow(1) + scalarExfiltration
 mLayerColumnOutflow(1) = mLayerColumnOutflow(1) + scalarExfiltration*HRUarea


!
!
!
!
!
! ! ************************************************************************************************
! ! (3) disaggregate total column outflow to estimate baseflow from each soil layer...
! ! ************************************************************************************************
! call disaggFlow(&
!                 ! input: model control
!                 dt,                         & ! intent(in): time step (s) -- used to calculate maximum possible inflow rate
!                 ixSaturation,               & ! intent(inout): index of the lowest saturated layer
!                 ! input: total column inflow and outflow
!                 totalColumnInflow,          & ! intent(in): total column inflow (m s-1)
!                 totalColumnOutflow,         & ! intent(in): total outflow from the soil column (m s-1)
!                 ! input: hydraulic conductivity in each soil layer
!                 mLayerVolFracLiq,           & ! intent(in): volumetric fraction of liquid water in each soil layer (-)
!                 mLayerVolFracIce,           & ! intent(in): volumetric fraction of ice in each soil layer (-)
!                 mLayerHydCond,              & ! intent(in): hydraulic conductivity in each soil layer (m s-1)
!                 dHydCond_dMatric,           & ! intent(in): derivative in hydraulic conductivity w.r.t. matric head in each layer (s-1)
!                 dOutflow_dMatric,           & ! intent(in): derivative in total column outflow w.r.t matric head in each layer (s-1)
!                 mLayerMatricHeadLiq,        & ! intent(in): liquid water matric potential (m)
!                 ! input: attributes and parameters
!                 HRUarea,                    & ! intent(in): HRU area (m2)
!                 theta_sat,                  & ! intent(in): soil porosity (-)
!                 fieldCapacity,              & ! intent(in): field capacity (-)
!                 ! input: storage and transmission properties in each layer
!                 mLayerDepth,                & ! intent(in): depth of each soil layer (m)
!                 mLayerColumnInflow,         & ! intent(in): inflow into each layer (m3 s-1)
!                 ! output: exfiltration and column outflow
!                 scalarExfiltration,         & ! intent(out): exfiltration (m s-1)
!                 mLayerColumnOutflow,        & ! intent(out): column outflow from each soil layer (m3 s-1)
!                 ! output: baseflow sink
!                 mLayerBaseflow,             & ! intent(out): baseflow from each soil layer (m s-1)
!                 dBaseflow_dMatric,          & ! intent(out): derivative in baseflow w.r.t. matric head (s-1)
!                 ! output: error control
!                 err,cmessage)                  ! intent(out): error control
! if(err/=0)then; message=trim(message)//trim(cmessage); return; endif

 ! end association to variables in data structures
 end associate

 !print*, 'in groundWatr: mLayerBaseflow = ', mLayerBaseflow
 !if(any(mLayerBaseflow < -verySmall)) pause 'negative baseflow'


 end subroutine groundwatr





 ! ************************************************************************************************
 ! ************************************************************************************************
 ! ************************************************************************************************
 ! ************************************************************************************************
 ! ***** PRIVATE SUBROUTINES **********************************************************************
 ! ************************************************************************************************
 ! ************************************************************************************************
 ! ************************************************************************************************
 ! ************************************************************************************************

 ! ************************************************************************************************
 ! private subroutine: compute saturated storage from individual layer state variables
 ! ************************************************************************************************
 subroutine satStorage(&
                       ! input
                       mLayerVolFracLiq,           & ! intent(in): volumetric liquid water content in each layer (-)
                       mLayerVolFracIce,           & ! intent(in): volumetric ice content in each layer (-)
                       mLayerDepth,                & ! intent(in): depth of each layer (m)
                       fieldCapacity,              & ! intent(in): field capacity (-)
                       theta_sat,                  & ! intent(in): soil porosity (-)
                       ! output
                       ixSaturation,               & ! intent(out): index of the lowest saturated layer
                       subSurfaceStorage,          & ! intent(out): sub surface storage (m)
                       maximumSoilWater,           & ! intent(out): maximum storage (m)
                       err,message)                  ! intent(out): error control
 ! ----------------------------------------------------------------------------------------------------------
 implicit none
 ! input
 real(dp),intent(in)          :: mLayerVolFracLiq(:)        ! volumetric liquid water content (-)
 real(dp),intent(in)          :: mLayerVolFracIce(:)        ! volumetric ice content (-)
 real(dp),intent(in)          :: mLayerDepth(:)             ! depth of each layer (m)
 real(dp),intent(in)          :: fieldCapacity              ! field capacity (-)
 real(dp),intent(in)          :: theta_sat                  ! soil porosity (-) 
 ! output
 integer(i4b),intent(out)     :: ixSaturation               ! index of the lowest saturated layer
 real(dp),intent(out)         :: subSurfaceStorage          ! sub surface storage (m)
 real(dp),intent(out)         :: maximumSoilWater           ! maximum storage (m)
 integer(i4b),intent(out)     :: err                        ! error code
 character(*),intent(out)     :: message                    ! error message
 ! local
 integer(i4b)                 :: iLayer                     ! index of model layer
 ! ----------------------------------------------------------------------------------------------------------
 ! initialize error control
 err=0; message='satStorage/'

 ! get index of the lowest saturated layer
 ixSaturation = nSoil+1  ! unsaturated profile when ixSaturation>nSoil
 do iLayer=nSoil,1,-1  ! start at the lowest soil layer and work upwards to the top layer
  if(mLayerVolFracLiq(iLayer) > fieldCapacity)then; ixSaturation = iLayer  ! index of saturated layer -- keeps getting over-written as move upwards
  else; exit; endif                                                        ! (only consider saturated layer at the bottom of the soil profile) 
 end do  ! (looping through soil layers)

 ! compute storage (m)
 if(ixSaturation < nSoil)then; subSurfaceStorage = sum((mLayerVolFracLiq(ixSaturation:nSoil) - fieldCapacity)*mLayerDepth(ixSaturation:nSoil))
 else;                         subSurfaceStorage = 0._dp
 endif  ! if some layers are saturated

 ! compute maximum possible sub-surface free storage (m)
 maximumSoilWater = sum(mLayerDepth(1:nSoil))*(theta_sat - fieldCapacity) &   ! maximum aquifer storage (m)
                         - sum(mLayerVolFracIce(1:nSoil)*mLayerDepth(1:nSoil))
 if(maximumSoilWater < 0._dp)then
  message=trim(message)//'ice content is greater than the sub-surface free storage -- need additional code to handle this situation'
  err=20; return
 endif

 end subroutine satStorage


 ! ************************************************************************************************
 ! private subroutine: disaggregate total inflow and total outflow to the baseflow sink term
 ! ************************************************************************************************
 subroutine disaggFlow(&
                       ! input: model control
                       dt,                         & ! intent(in): time step (s) -- used to calculate maximum possible inflow rate
                       ixSaturation,               & ! intent(inout): index of the lowest saturated layer
                       ! input: total column inflow and outflow
                       totalColumnInflow,          & ! intent(in): total column inflow (m s-1)
                       totalColumnOutflow,         & ! intent(in): total outflow from the soil column (m s-1)
                       ! input: hydraulic conductivity in each soil layer
                       mLayerVolFracLiq,           & ! intent(in): volumetric fraction of liquid water in each soil layer (-)
                       mLayerVolFracIce,           & ! intent(in): volumetric fraction of ice in each soil layer (-)
                       mLayerHydCond,              & ! intent(in): hydraulic conductivity in each soil layer (m s-1)
                       dHydCond_dMatric,           & ! intent(in): derivative in hydraulic conductivity w.r.t. matric head in each layer (s-1)
                       dOutflow_dMatric,           & ! intent(in): derivative in total column outflow w.r.t matric head in each layer (s-1)
                       mLayerMatricHeadLiq,        & ! intent(in): liquid water matric potential (m)
                       ! input: attributes and parameters
                       HRUarea,                    & ! intent(in): HRU area (m2)
                       theta_sat,                  & ! intent(in): soil porosity (-)
                       fieldCapacity,              & ! intent(in): field capacity (-)
                       ! input: storage and transmission properties in each layer
                       mLayerDepth,                & ! intent(in): depth of each soil layer (m)
                       mLayerColumnInflow,         & ! intent(in): inflow into each layer (m3 s-1)
                       ! output: exfiltration and column outflow
                       exfiltration,               & ! intent(out): exfiltration (m s-1)
                       mLayerColumnOutflow,        & ! intent(out): column outflow from each soil layer (m3 s-1)
                       ! output: baseflow sink
                       mLayerBaseflow,             & ! intent(out): baseflow from each soil layer (m s-1)
                       dBaseflow_dMatric,          & ! intent(out): derivative in baseflow w.r.t. matric head (s-1)
                       ! output: error control
                       err,message)                  ! intent(out): error control
 ! ----------------------------------------------------------------------------------------------------------
 ! utility modules
 USE soil_utils_module,only:hydCond_psi         ! compute hydraulic conductivity as a function of matric head
 ! data structures
 USE data_struc,only:mpar_data                  ! data structures
 USE data_struc,only:mvar_data                  ! data structures
 USE var_lookup,only:iLookPARAM,iLookMVAR       ! named variables for structure elements
 ! ----------------------------------------------------------------------------------------------------------
 implicit none
 ! input: model control
 real(dp),intent(in)          :: dt                         ! time step (s) -- used to calculate maximum possible inflow rate
 integer(i4b),intent(inout)   :: ixSaturation               ! index of the lowest saturated layer
 ! input: total column inflow and outflow
 real(dp),intent(in)          :: totalColumnInflow          ! total column inflow (m s-1)
 real(dp),intent(in)          :: totalColumnOutflow         ! total outflow from the soil column (m s-1)
 ! input: storage and transmission properties  in each soil layer
 real(dp),intent(in)          :: mLayerVolFracLiq(:)        ! volumetric fraction of liquid water in each soil layer (-)
 real(dp),intent(in)          :: mLayerVolFracIce(:)        ! volumetric fraction of ice in each soil layer (-)
 real(dp),intent(in)          :: mLayerHydCond(:)           ! hydraulic conductivity in each layer (m s-1)
 real(dp),intent(in)          :: dHydCond_dMatric(:)        ! derivative in hydraulic conductivity w.r.t matric head in each layer (s-1)
 real(dp),intent(in)          :: dOutflow_dMatric(:)        ! derivative in total column outflow w.r.t matric head in each layer (s-1)
 real(dp),intent(in)          :: mLayerMatricHeadLiq(:)     ! matric head in each layer at the current iteration (m)
 ! input: attributes and parameters
 real(dp),intent(in)          :: HRUarea                    ! HRU area (m2)
 real(dp),intent(in)          :: theta_sat                  ! soil porosity (-)
 real(dp),intent(in)          :: fieldCapacity              ! field capacity (-)
 ! input: storage and transmission properties in each layer
 real(dp),intent(in)          :: mLayerDepth(:)             ! depth of each layer (m)
 real(dp),intent(in)          :: mLayerColumnInflow(:)      ! inflow into each layer (m3 s-1)
 ! output: exfiltration and column outflow
 real(dp),intent(out)         :: exfiltration               ! exfiltration (m s-1)
 real(dp),intent(out)         :: mLayerColumnOutflow(:)     ! column outflow in each layer (m3 s-1)
 ! output: baseflow sink
 real(dp),intent(out)         :: mLayerBaseflow(:)          ! baseflow in each layer (m s-1)
 real(dp),intent(out)         :: dBaseflow_dMatric(:,:)     ! derivative in baseflow w.r.t. matric head (s-1)
 ! output: error control
 integer(i4b),intent(out)     :: err                        ! error code
 character(*),intent(out)     :: message                    ! error message
 ! ----------------------------------------------------------------------------------------------------------
 ! local variables
 real(dp)                     :: f0, f1
 real(dp)                     :: hydCond1
 real(dp)                     :: sumDepthAvgCond1
 real(dp),dimension(nSoil)    :: mLayerHydCondCopy
 real(dp),pointer             :: scalarSatHydCond,vGn_alpha,vGn_n,vGn_m
 real(dp)                     :: fPart1,fPart2              ! different parts of a function: used to compute analytical derivatives
 real(dp)                     :: dPart1,dPart2              ! derivatives for different parts of a function: used to compute analytical derivatives
 real(dp)                     :: sumDepthAvgCond            ! sum of depth-weighted hydraulic conductivity (m2 s-1)
 real(dp),dimension(nSoil)    :: fracTotalOutflow           ! fraction of outflow apportioned to each layer (-)
 real(dp),dimension(nSoil)    :: mLayerColumnInflowAdjusted ! adjusted column inflow to ensure no inflow into ice layers (m3 s-1)
 real(dp)                     :: totalInflowUnallocated     ! unallocated inflow (m3 s-1)
 real(dp)                     :: volTotalWater              ! volumetric fraction of total water (liquid + ice)
 real(dp)                     :: qMax                       ! max inflow rate (m s-1)
 real(dp)                     :: sink                       ! net source/sink (m s-1)
 integer(i4b)                 :: iLayer,jLayer              ! index of model layer
 ! ----------------------------------------------------------------------------------------------------------
 ! initialize error control
 err=0; message='disaggFlow/'

 ! point to variables in the data structures
 vGn_alpha               => mpar_data%var(iLookPARAM%vGn_alpha)                    ! intent(in): [dp] van Genutchen "alpha" parameter (m-1)
 vGn_n                   => mpar_data%var(iLookPARAM%vGn_n)                        ! intent(in): [dp] van Genutchen "n" parameter (-)
 vGn_m                   => mvar_data%var(iLookMVAR%scalarVGn_m)%dat(1)            ! intent(in): [dp] van Genutchen "m" parameter (-)

 ! initialize adjusted inflow to each layer (m3/s)
 mLayerColumnInflowAdjusted(1:nSoil) = 0._dp

 ! simple case with no flow
 if(totalColumnOutflow < tiny(theta_sat))then
  mLayerColumnOutflow(1:nSoil) = HRUarea*totalColumnOutflow/real(nSoil, kind(dp))
  mLayerBaseflow(1:nSoil)      = (mLayerColumnOutflow(1:nSoil) - mLayerColumnInflow(1:nSoil))/HRUarea
  return
 endif

 ! possible to start iterations with zero storage, and gain storage through inflow -- reset ixSaturation
 !if(ixSaturation > nSoil) ixSaturation=nSoil

 ! ****
 ! compute derivatives


 ! compute the depth-weighted hydraulic conductivity, i.e, the transmissivity (m2 s-1)
 sumDepthAvgCond = sum(mLayerHydCond(ixSaturation:nSoil)*mLayerDepth(ixSaturation:nSoil))

 ! check that the depth average hydraulic conductivity is positive when baseflow occurs
 if(sumDepthAvgCond < tiny(theta_sat))then
  write(*,'(a,1x,10(e20.10,1x))') 'mLayerHydCond(ixSaturation:nSoil)       = ', mLayerHydCond(ixSaturation:nSoil)
  write(*,'(a,1x,10(e20.10,1x))') 'mLayerDepth(ixSaturation:nSoil)         = ', mLayerDepth(ixSaturation:nSoil)
  message=trim(message)//'zero depth-weighted conductivity when baseflow occurs'
  err=20; return
 endif

 ! **************************************************************************************************************
 ! * compute the fraction of total outflow apportioned to each layer, and its derivative w.r.t. matric head
 ! **************************************************************************************************************

 ! loop through active soil layers
 !do iLayer=ixSaturation,nSoil

 ! ! compute the function and derivative for the numerator
 ! fPart1 = mLayerHydCond(iLayer)*mLayerDepth(iLayer)        ! numerator = transmissivity of layer (m2 s-1)
 ! dPart1 = dHydCond_dMatric(iLayer)*mLayerDepth(iLayer)     ! derivative of the numerator (m s-1) -- needed for the diagonal terms
!
!  ! compute fraction of total outflow apportioned to each layer
!  fracTotalOutflow(iLayer) = fPart1/sumDepthAvgCond
!
!  ! loop through active soil layers
!  do jLayer=ixSaturation,nSoil
!
!   ! compute the derivative for the denominator
!   dPart2 = -mLayerDepth(jLayer)*dHydCond_dMatric(jLayer)/sumDepthAvgCond**2._dp ! derivative for the denominator (s m-3)
!
!   ! diagonal terms
!   if(iLayer==jLayer)then
!    dFracTotalOutflow_dMatric(iLayer,jLayer) = fPart1*dPart2 + dPart1/sumDepthAvgCond  ! m-1
!
!   ! x-derivative terms
!   else
!    dFracTotalOutflow_dMatric(iLayer,jLayer) = fPart1*dPart2
!
!   end if ! (switch between diagonal and x-derivative terms)
!
!  end do  ! looping through active soil layers
! end do  ! looping through active soil layers
!
! ! set unsaturated elements to zero
! if(ixSaturation > 1)then
!  fracTotalOutflow(1:ixSaturation-1) = 0._dp
!  dFracTotalOutflow_dMatric(1:ixSaturation-1,1:ixSaturation-1) = 0._dp
! endif
!
! ! check that fraction of baseflow apportioned to the different soil layers sums to 1
! if(abs(1._dp - sum(fracTotalOutflow)) > verySmall)then
!  message=trim(message)//'fraction of baseflow does not sum to 1'
!  err=20; return
! endif
!
 ! **************************************************************************************************************

 ! compute the outflow from each soil layer (m3 s-1)
 mLayerColumnOutflow(1:nSoil) = fracTotalOutflow(1:nSoil)*totalColumnOutflow*HRUarea

 ! compute the exfiltration (m s-1)
 !exfiltration = totalInflowUnallocated/HRUarea

 ! compute the net baseflow from each soil layer (m s-1)
 mLayerBaseflow(1:nSoil) = fracTotalOutflow(1:nSoil) 
 !mLayerBaseflow(1:nSoil) = (mLayerColumnOutflow(1:nSoil) - mLayerColumnInflowAdjusted(1:nSoil))/HRUarea

 ! add exfiltration to the baseflow flux at the top layer
 mLayerBaseflow(1)      = mLayerBaseflow(1) + exfiltration
 mLayerColumnOutflow(1) = mLayerColumnOutflow(1) + exfiltration*HRUarea

 ! initialize total unallocated inflow (m3 s-1)
 totalInflowUnallocated = totalColumnInflow*HRUarea
 !print*, 'totalInflowUnallocated = ', totalInflowUnallocated

 ! adjust the inflow into each layer -- start at the bottom and work upwards
 do iLayer=nSoil,1,-1
  ! (first check if inflow exceeds outflow)
  if(totalInflowUnallocated < mLayerColumnOutflow(iLayer))then
   mLayerColumnInflowAdjusted(iLayer) = totalInflowUnallocated
   exit  ! allocated all inflow, so exit do loop because all other layers set to zero
  endif
  ! (next, try to allocate as much inflow as possible to the current layer)
  volTotalWater = mLayerVolFracLiq(iLayer) + mLayerVolFracIce(iLayer)
  qMax = HRUarea*mLayerDepth(iLayer)*(theta_sat - max(volTotalWater,fieldCapacity))/dt   ! maximum allowable net flux (m3 s-1)
  !sink = (totalInflowUnallocated + mLayerColumnOutflow(iLayer))                          ! maximum source/sink (m3 s-1)
  mLayerColumnInflowAdjusted(iLayer) = min(totalInflowUnallocated,qMax)                  ! inflow into the current layer (m3 s-1)
  ! (save remaining water for the next layer)
  totalInflowUnallocated = totalInflowUnallocated - mLayerColumnInflowAdjusted(iLayer)  ! m3 s-1
  ! (add as much water as possible to the current layer)
  !write(*,'(a,1x,i4,1x,10(e20.10,1x))') 'iLayer, qMax, mLayerColumnOutflow(iLayer), mLayerColumnInflowAdjusted(iLayer), totalInflowUnallocated = ', &
  !                                       iLayer, qMax, mLayerColumnOutflow(iLayer), mLayerColumnInflowAdjusted(iLayer), totalInflowUnallocated
 end do  ! (looping through soil layers)

 ! print inflow
 !print*, 'totalColumnInflow = ', totalColumnInflow
 !print*, 'mLayerColumnInflowAdjusted = ', mLayerColumnInflowAdjusted
 !pause

 ! compute the exfiltration (m s-1)
 exfiltration = totalInflowUnallocated/HRUarea

 ! compute the net baseflow from each soil layer (m s-1)
 mLayerBaseflow(1:nSoil) = fracTotalOutflow(1:nSoil) 
 !mLayerBaseflow(1:nSoil) = (mLayerColumnOutflow(1:nSoil) - mLayerColumnInflowAdjusted(1:nSoil))/HRUarea

 ! add exfiltration to the baseflow flux at the top layer
 mLayerBaseflow(1)      = mLayerBaseflow(1) + exfiltration
 mLayerColumnOutflow(1) = mLayerColumnOutflow(1) + exfiltration*HRUarea

 end subroutine disaggFlow


end module groundwatr_module

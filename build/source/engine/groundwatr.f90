module groundwatr_module
! numerical recipes data types
USE nrtype
! model constants
USE multiconst,only:iden_water ! density of water (kg m-3)
! access the number of snow and soil layers
USE data_struc,only:&
                    nSnow,     & ! number of snow layers  
                    nSoil,     & ! number of soil layers  
                    nLayers      ! total number of layers
! provide access to the derived types to define the data structures
USE data_struc,only:&
                    var_d,     & ! data vector (dp)
                    var_dlength  ! data vector with variable length dimension (dp)
! provide access to named variables defining elements in the data structures
USE var_lookup,only:iLookATTR,iLookPARAM,iLookMVAR     
 ! utility modules
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
 !
 ! Method
 ! ------
 !
 ! Here we assume that water avaialble for shallow groundwater flow includes is all water above
 ! "field capacity" below the depth zCrit, where zCrit is defined as the lowest point in the soil
 ! profile where the volumetric liquid water content is less than field capacity.
 !
 ! We further assume that transmssivity (m2 s-1) for each layer is defined asuming that the water
 ! available for saturated flow is located at the bottom of the soil profile. Specifically:
 !  trTotal(iLayer) = tran0*(zActive(iLayer)/soilDepth)**zScale_TOPMODEL
 !  trSoil(iLayer)  = trTotal(iLayer) - trTotal(iLayer+1)
 ! where zActive(iLayer) is the effective water table thickness for all layers up to and including
 ! the current layer (working from the bottom to the top).
 !
 ! The outflow from each layer is then (m3 s-1)
 !  mLayerOutflow(iLayer) = trSoil(iLayer)*tan_slope*contourLength
 ! where contourLength is the width of a hillslope (m) parallel to a stream
 !
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
 real(dp)                         :: soilDepth                    ! intent(in): total soil depth (m)
 real(dp),dimension(nSoil)        :: mLayerDepth                  ! intent(in): depth of each soil layer (m)
 ! input: diagnostic variables
 real(dp)                         :: surfaceHydCond               ! intent(in): saturated hydraulic conductivity at the surface (m s-1)
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
 real(dp)                         :: vGn_alpha,vGn_n,vGn_m        ! intent(in): van Genuchten parameters
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
 real(dp)                        :: activePorosity                ! "active" porosity associated with storage above a threshold (-)
 real(dp)                        :: tran0                         ! maximum transmissivity (m2 s-1)
 real(dp),dimension(nSoil)       :: zActive                       ! water table thickness associated with storage below and including the given layer (m)
 real(dp),dimension(nSoil)       :: trTotal                       ! total transmissivity associated with total water table depth zActive (m2 s-1)
 real(dp),dimension(nSoil)       :: trSoil                        ! transmissivity of water in a given layer (m2 s-1)
 ! local variables for the derivatives
 real(dp)                        :: length2area                   ! ratio of hillslope width to hillslope area (m m-2)
 real(dp),dimension(nSoil)       :: depth2capacity                ! ratio of layer depth to total subsurface storage capacity (-)
 real(dp),dimension(nSoil)       :: dXdS                          ! change in dimensionless flux w.r.t. change in dimensionless storage (-)
 ! local variables to compute the numerical Jacobian
 logical(lgt),parameter          :: doNumericalJacobian=.true.    ! flag to compute the numerical Jacobian
 real(dp),dimension(nSoil)       :: mLayerVolFracLiqPerturbed     ! perturbed volumetric fraction of liquid water (-)
 real(dp),dimension(nSoil)       :: mLayerBaseflowPerturbed       ! perturbed baseflow (m s-1)
 real(dp),dimension(nSoil,nSoil) :: nJac                          ! numerical Jacobian (s-1)






 real(dp)                        :: mLayerSatHydCondMP
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
 real(qp)                        :: dPart1,dPart2,dPart3          ! derivatives for part of a function
 real(dp),dimension(nSoil)       :: dOutflow_dMatric              ! derivative in total column outflow w.r.t matric head in each layer (s-1)
 real(dp),dimension(nSoil,nSoil) :: dFracTotalOutflow_dMatric     ! derivative in fraction of outflow apportioned to each layer w.r.t matric head in each layer (m-1)

 real(qp)                     :: hydCond1
 real(qp)                     :: sumDepthAvgCond1
 real(qp),dimension(nSoil)    :: mLayerHydCondCopy
 real(dp),pointer             :: scalarSatHydCond

 real(dp)                        :: xDepth,xPore,xTemp,xTran,xFlow

 real(qp)                        :: z0,z1
 real(qp)                        :: t0,t1,tOld
 real(qp)                        :: b0,b1


 real(qp)                        :: f0,f1
 integer(i4b)                    :: iLayer,jLayer, kLayer
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
 soilDepth               => mvar_data%var(iLookMVAR%iLayerHeight)%dat(nSoil),       & ! intent(in): [dp]    total soil depth (m)
 mLayerDepth             => mvar_data%var(iLookMVAR%mLayerDepth)%dat(1:nSoil),      & ! intent(in): [dp(:)] depth of each soil layer (m)
 
 ! input: diagnostic variables
 surfaceHydCond          => mvar_data%var(iLookMVAR%mLayerSatHydCondMP)%dat(1),     & ! intent(in): [dp]    saturated hydraulic conductivity at the surface (m s-1)
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
 ! (1) compute the "active" portion of the soil profile
 ! ************************************************************************************************

 ! get index of the lowest saturated layer
 ixSaturation = nSoil+1  ! unsaturated profile when ixSaturation>nSoil
 do iLayer=nSoil,1,-1  ! start at the lowest soil layer and work upwards to the top layer
  if(mLayerVolFracLiq(iLayer) > fieldCapacity)then; ixSaturation = iLayer  ! index of saturated layer -- keeps getting over-written as move upwards
  else; exit; endif                                                        ! (only consider saturated layer at the bottom of the soil profile) 
 end do  ! (looping through soil layers)

 ! check for an early return (no layers are "active")
 if(ixSaturation > nSoil)then
  scalarExfiltration     = 0._dp   ! exfiltration from the soil profile (m s-1)
  mLayerColumnOutflow(:) = 0._dp   ! column outflow from each soil layer (m3 s-1)
  mLayerBaseflow(:)      = 0._dp   ! baseflow from each soil layer (m s-1)
  dBaseflow_dMatric(:,:) = 0._dp   ! derivative in baseflow w.r.t. matric head (s-1)
  return
 endif  ! if some layers are saturated

 ! ************************************************************************************************
 ! (2) compute the baseflow flux and its derivative w.r.t volumetric liquid water content
 ! ************************************************************************************************

 ! use private subroutine to compute baseflow (for multiple calls for numerical Jacobian)
 call computeBaseflow(&
                      ! input: control and state variables
                      .true.,                  & ! intent(in): .true. if derivatives are desired
                      ixSaturation,            & ! intent(in): index of upper-most "saturated" layer
                      mLayerVolFracLiq,        & ! intent(in): volumetric fraction of liquid water in each soil layer (-)
                      ! input/output: data structures
                      attr_data,               & ! intent(in):    spatial attributes
                      mpar_data,               & ! intent(in):    model parameters
                      mvar_data,               & ! intent(inout): model variables for a local HRU
                      ! output: fluxes and derivatives
                      mLayerBaseflow,          & ! intent(out): baseflow flux in each soil layer (m s-1)
                      dBaseflow_dVolLiq)         ! intent(out): derivative in baseflow w.r.t. volumetric liquid water content (s-1)


 ! ************************************************************************************************
 ! (3) compute the derivatives in baseflow w.r.t volumetric liquid water content in each layer
 ! ************************************************************************************************
 
 ! initialize the derivative matrix
 dBaseflow_dMatric(:,:) = 0._dp

 ! compute ratio of hillslope width to hillslope area (m m-2)
 length2area = tan_slope*contourLength/HRUarea  

 ! compute the ratio of layer depth to maximum water holding capacity (-)
 depth2capacity(1:nSoil) = mLayerDepth(1:nSoil)/(activePorosity*soilDepth)

 ! compute the change in dimensionless flux w.r.t. change in dimensionless storage (-)
 dXdS(1:nSoil) = zScale_TOPMODEL*(zActive(1:nSoil)/SoilDepth)**(zScale_TOPMODEL - 1._dp) 

 ! loop through soil layers
 do iLayer=1,nSoil
  ! compute diagonal terms (s-1)
  dBaseflow_dMatric(iLayer,iLayer) = mLayerdTheta_dPsi(iLayer)*tran0*dXdS(iLayer)*depth2capacity(iLayer)*length2area
  ! compute off-diagonal terms
  do jLayer=iLayer+1,nSoil  ! (only dependent on layers below)
   dBaseflow_dMatric(iLayer,jLayer) = mLayerdTheta_dPsi(jLayer)*(tran0/soilDepth)*(dXdS(iLayer) - dXdS(iLayer+1))
  end do  ! looping through soil layers
 end do  ! looping through soil layers


 ! ************************************************************************************************
 ! (4) compute the numerical Jacobian
 ! ************************************************************************************************
 if(doNumericalJacobian)then

  ! first, print the analytical Jacobian
  print*, '** analytical Jacobian:'
  write(*,'(a4,1x,100(i12,1x))') 'xCol', (iLayer, iLayer=1,nSoil)
  do iLayer=1,nSoil; write(*,'(i4,1x,100(e12.5,1x))') iLayer, dBaseflow_dMatric(1:nSoil,iLayer); end do

  ! get a copy of the state vector to perturb
  mLayerVolFracLiqPerturbed(:) = mLayerVolFracLiq(:)

  ! loop through the state variables
  do iLayer=1,nSoil

   ! perturb state vector
   mLayerVolFracLiqPerturbed(iLayer) = mLayerVolFracLiq(iLayer) + dx       ! perturb state vector

   ! compute baseflow flux
   call computeBaseflow(&
                        ! input: control and state variables
                        .false.,                   & ! intent(in): .true. if derivatives are desired
                        ixSaturation,              & ! intent(in): index of upper-most "saturated" layer
                        mLayerVolFracLiqiPerturbed,& ! intent(in): volumetric fraction of liquid water in each soil layer (-)
                        ! input/output: data structures
                        attr_data,                 & ! intent(in):    spatial attributes
                        mpar_data,                 & ! intent(in):    model parameters
                        mvar_data,                 & ! intent(inout): model variables for a local HRU
                        ! output: fluxes and derivatives
                        mLayerBaseflowPerturbed,   & ! intent(out): baseflow flux in each soil layer (m s-1)
                        dBaseflow_dVolLiq)           ! intent(out): derivative in baseflow w.r.t. volumetric liquid water content (s-1)
   
   ! compute the numerical Jacobian
   nJac(:,iLayer) = (mLayerBaseflowPerturbed(:) - mLayerBaseflow(:))/dx    ! compute the Jacobian

   ! set the state back to the input value
   mLayerVolFracLiqPerturbed(iLayer) = mLayerVolFracLiq(iLayer)

  end do  ! looping through state variables

  ! print the numerical Jacobian
  print*, '** numerical Jacobian:'
  write(*,'(a4,1x,100(i12,1x))') 'xCol', (iLayer, iLayer=1,nSoil)
  do iLayer=1,nSoil; write(*,'(i4,1x,100(e12.5,1x))') iLayer, nJac(1:nSoil,iLayer); end do
  pause 'testing Jacobian'

 endif  ! if desire to compute the Jacobian

 pause 'checking derivatives'



 !iLayer = 1
 !jLayer = 4

 ! compute analytical derivatives for baseflow w.r.t. volumetric liqyud water content (m s-1)
 !dPart1 = mLayerDepth(iLayer)/(activePorosity*soilDepth)
 !dPart2 = tran0*zScale_TOPMODEL*(zActive(iLayer)/SoilDepth)**(zScale_TOPMODEL - 1._dp)
 !write(*,'(a,1x,e20.10,1x)') 'anal deriv   = ', dPart1*dPart2*tan_slope*contourLength/HRUarea 

 ! check x-derivative terms....
 !f0 = zActive(iLayer)
 !f1 = zActive(jLayer+1) + mLayerDepth(jLayer)*((mLayerVolFracLiq(jLayer)+dx) - fieldCapacity)/activePorosity
 !do kLayer=jLayer-1,iLayer,-1
 ! print*, 'kLayer = ', kLayer
 ! f1 = f1 + mLayerDepth(kLayer)*(mLayerVolFracLiq(kLayer) - fieldCapacity)/activePorosity
 ! if(kLayer==iLayer+1) tOld = tran0*(f1/soilDepth)**zScale_TOPMODEL
 !end do
 !write(*,'(a,1x,e20.10,1x)') '(f1 - f0)/dx = ', (f1 - f0)/dx

 ! check total transmissivity
 !t0 = tran0*(f0/soilDepth)**zScale_TOPMODEL - trTotal(iLayer+1)
 !t1 = tran0*(f1/soilDepth)**zScale_TOPMODEL - tOld
 !write(*,'(a,1x,e20.10,1x)') '(t1 - t0)/dx = ', (t1 - t0)/dx

 !dPart1 = tran0*zScale_TOPMODEL/SoilDepth
 !dPart2 = (zActive(iLayer)/SoilDepth)**(zScale_TOPMODEL - 1._dp)
 !dPart3 = (zActive(iLayer+1)/SoilDepth)**(zScale_TOPMODEL - 1._dp)

 !write(*,'(a,1x,e20.10,1x)') 'anal deriv   = ', dPart1*(dPart2 - dPart3)
 !pause ' check x-deriv'

 ! check depth
 !f0 = zActive(iLayer)
 !f1 = zActive(iLayer+1) + mLayerDepth(iLayer)*((mLayerVolFracLiq(iLayer)+dx) - fieldCapacity)/activePorosity
 !write(*,'(a,1x,e20.10,1x)') '(f1 - f0)/dx = ', (f1 - f0)/dx
 !write(*,'(a,1x,e20.10,1x)') 'anal deriv   = ', mLayerDepth(iLayer)/activePorosity


 ! check total transmissivity
 !t0 = tran0*(f0/soilDepth)**zScale_TOPMODEL - trTotal(iLayer+1)
 !t1 = tran0*(f1/soilDepth)**zScale_TOPMODEL - trTotal(iLayer+1)

 ! compute the analytical derivative
 !dPart1 = mLayerDepth(iLayer)/(activePorosity*soilDepth)
 !dPart2 = zScale_TOPMODEL*(f0/soilDepth)**(zScale_TOPMODEL - 1._dp)
 !write(*,'(a,1x,e20.10,1x)') '(t1 - t0)/dx = ', (t1 - t0)/dx
 !write(*,'(a,1x,e20.10,1x)') 'anal deriv   = ', tran0*dPart2*dPart1

 ! compute baseflow
 !b0 = (t0*tan_slope*contourLength - mLayerColumnInflow(iLayer))/HRUarea
 !b1 = (t1*tan_slope*contourLength - mLayerColumnInflow(iLayer))/HRUarea
 !write(*,'(a,1x,e20.10,1x)') '(b1 - b0)/dx = ', (b1 - b0)/dx
 !write(*,'(a,1x,e20.10,1x)') 'anal deriv   = ', tran0*dPart2*dPart1*tan_slope*contourLength/HRUarea

 !pause

 ! end association to variables in data structures
 end associate

 end subroutine groundwatr


 ! ***********************************************************************************************************************
 ! ***********************************************************************************************************************
 ! ***********************************************************************************************************************
 ! ***********************************************************************************************************************
 ! **** PRIVATE SUBROUTINES **********************************************************************************************
 ! ***********************************************************************************************************************
 ! ***********************************************************************************************************************
 ! ***********************************************************************************************************************
 ! ***********************************************************************************************************************

 ! ***********************************************************************************************************************
 ! * compute baseflow: private subroutine so can be used to test the numerical jacobian
 ! ***********************************************************************************************************************
 subroutine computeBaseflow(&
                            ! input: control and state variables
                            derivDesired,                  & ! intent(in): .true. if derivatives are desired
                            ixSaturation,                  & ! intent(in): index of upper-most "saturated" layer
                            mLayerVolFracLiq,              & ! intent(in): volumetric fraction of liquid water in each soil layer (-)
                            ! input/output: data structures
                            attr_data,                     & ! intent(in):    spatial attributes
                            mpar_data,                     & ! intent(in):    model parameters
                            mvar_data,                     & ! intent(inout): model variables for a local HRU
                            ! output: fluxes and derivatives
                            mLayerBaseflow,                & ! intent(out): baseflow flux in each soil layer (m s-1)
                            dBaseflow_dVolLiq)               ! intent(out): derivative in baseflow w.r.t. volumetric liquid water content (s-1)
 implicit none
 ! ---------------------------------------------------------------------------------------
 ! * dummy variables
 ! ---------------------------------------------------------------------------------------
 ! input: control and state variables
 logical(lgt),intent(in)          :: derivDesired            ! .true. if derivatives are desired
 integer(i4b),intent(in)          :: ixSaturation            ! index of upper-most "saturated" layer
 real(dp),intent(in)              :: mLayerVolFracLiq(:)     ! volumetric fraction of liquid water (-)
 real(dp),intent(out)             :: mLayerBaseflow(:)       ! baseflow from each soil layer (m s-1)
 ! input/output: data structures
 type(var_d),intent(in)           :: attr_data               ! spatial attributes
 type(var_d),intent(in)           :: mpar_data               ! model parameters
 type(var_dlength),intent(inout)  :: mvar_data               ! model variables for a local HRU
 ! output: baseflow
 real(dp),intent(out)             :: mLayerBaseflow(:)       ! baseflow from each soil layer (m s-1)
 real(dp),intent(out)             :: dBaseflow_dVolLiq(:,:)  ! derivative in baseflow w.r.t. matric head (s-1)
 ! ---------------------------------------------------------------------------------------
 ! * variables in the data structures
 ! ---------------------------------------------------------------------------------------
 ! input: coordinate variables
 real(dp)                         :: soilDepth               ! intent(in): total soil depth (m)
 real(dp),dimension(nSoil)        :: mLayerDepth             ! intent(in): depth of each soil layer (m)
 ! input: diagnostic variables
 real(dp)                         :: surfaceHydCond          ! intent(in): saturated hydraulic conductivity at the surface (m s-1)
 real(dp),dimension(nSoil)        :: mLayerColumnInflow      ! intent(in): inflow into each soil layer (m3/s)
 ! input: local attributes
 real(dp)                         :: HRUarea                 ! intent(in): HRU area (m2)
 real(dp)                         :: tan_slope               ! intent(in): tan water table slope, taken as tan local ground surface slope (-)
 real(dp)                         :: contourLength           ! intent(in): length of contour at downslope edge of HRU (m)
 ! input: baseflow parameters
 real(dp)                         :: zScale_TOPMODEL         ! intent(in): TOPMODEL exponent (-)
 real(dp)                         :: kAnisotropic            ! intent(in): anisotropy factor for lateral hydraulic conductivity (-)
 real(dp)                         :: fieldCapacity           ! intent(in): field capacity (-)
 real(dp)                         :: theta_sat               ! intent(in): soil porosity (-) 
 ! output: diagnostic variables
 real(dp)                         :: scalarExfiltration      ! intent(out): exfiltration from the soil profile (m s-1)
 real(dp),dimension(nSoil)        :: mLayerColumnOutflow     ! intent(out): column outflow from each soil layer (m3 s-1)
 ! ---------------------------------------------------------------------------------------
 ! * local variables
 ! ---------------------------------------------------------------------------------------
 ! general local variables
 integer(i4b)                    :: iLayer,jLayer            ! index of model layer
 ! local variables for the lateral flux among soil columns
 integer(i4b)                    :: ixSaturation             ! index of the lowest saturated layer
 real(dp)                        :: activePorosity           ! "active" porosity associated with storage above a threshold (-)
 real(dp)                        :: tran0                    ! maximum transmissivity (m2 s-1)
 real(dp),dimension(nSoil)       :: zActive                  ! water table thickness associated with storage below and including the given layer (m)
 real(dp),dimension(nSoil)       :: trTotal                  ! total transmissivity associated with total water table depth zActive (m2 s-1)
 real(dp),dimension(nSoil)       :: trSoil                   ! transmissivity of water in a given layer (m2 s-1)
 ! local variables for the derivatives
 real(dp)                        :: length2area              ! ratio of hillslope width to hillslope area (m m-2)
 real(dp),dimension(nSoil)       :: depth2capacity           ! ratio of layer depth to total subsurface storage capacity (-)
 real(dp),dimension(nSoil)       :: dXdS                     ! change in dimensionless flux w.r.t. change in dimensionless storage (-)
 ! local variables for testing (debugging)
 logical(lgt),parameter          :: printFlag=.false.        ! flag for printing (debugging)
 real(dp)                        :: xDepth,xTran,xFlow       ! temporary variables (depth, transmissivity, flow)
 ! ---------------------------------------------------------------------------------------
 ! * association to data in structures
 ! ---------------------------------------------------------------------------------------
 associate(&

 ! input: coordinate variables
 soilDepth               => mvar_data%var(iLookMVAR%iLayerHeight)%dat(nSoil),       & ! intent(in): [dp]    total soil depth (m)
 mLayerDepth             => mvar_data%var(iLookMVAR%mLayerDepth)%dat(1:nSoil),      & ! intent(in): [dp(:)] depth of each soil layer (m)

 ! input: diagnostic variables
 surfaceHydCond          => mvar_data%var(iLookMVAR%mLayerSatHydCondMP)%dat(1),     & ! intent(in): [dp]    saturated hydraulic conductivity at the surface (m s-1)
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

 ! output: diagnostic variables
 scalarExfiltration      => mvar_data%var(iLookMVAR%scalarExfiltration)%dat(1),     & ! intent(out):[dp]    exfiltration from the soil profile (m s-1)
 mLayerColumnOutflow     => mvar_data%var(iLookMVAR%mLayerColumnOutflow)%dat        & ! intent(out):[dp(:)] column outflow from each soil layer (m3 s-1)

 )  ! end association to variables in data structures
 ! ***********************************************************************************************************************
 ! ***********************************************************************************************************************
 ! start routine here

 ! ***********************************************************************************************************************
 ! (1) compute the baseflow flux in each soil layer
 ! ***********************************************************************************************************************

 ! compute the porosity and the maximum transmissivity
 ! NOTE: this can be done as a pre-processing step
 activePorosity = theta_sat - fieldCapacity                               ! "active" porosity (-)
 tran0          = kAnisotropic*surfaceHydCond*soilDepth/zScale_TOPMODEL   ! maximum transmissivity (m2 s-1)

 ! compute the water table thickness (m) and transmissivity in each layer (m2 s-1)
 do iLayer=nSoil,ixSaturation,-1  ! loop through "active" soil layers, from lowest to highest
  if(iLayer==nSoil)then
   zActive(iLayer) = mLayerDepth(iLayer)*(mLayerVolFracLiq(iLayer) - fieldCapacity)/activePorosity  ! water table thickness associated with storage in a given layer (m)
   trTotal(iLayer) = tran0*(zActive(iLayer)/soilDepth)**zScale_TOPMODEL   ! total transmissivity for total depth zActive (m2 s-1)
   trSoil(iLayer)  = trTotal(iLayer)                                      ! transmissivity of water in a given layer (m2 s-1)
  else
   zActive(iLayer) = zActive(iLayer+1) + mLayerDepth(iLayer)*(mLayerVolFracLiq(iLayer) - fieldCapacity)/activePorosity
   trTotal(iLayer) = tran0*(zActive(iLayer)/soilDepth)**zScale_TOPMODEL
   trSoil(iLayer)  = trTotal(iLayer) - trTotal(iLayer+1)
  endif
  !write(*,'(a,1x,i4,1x,5(f20.10,1x))') 'iLayer, zActive(iLayer), trTotal(iLayer), trSoil(iLayer) = ', iLayer, zActive(iLayer), trTotal(iLayer), trSoil(iLayer)
 end do  ! looping through soil layers

 ! set un-used portions of the vectors to zero
 if(ixSaturation>1)then
  zActive(1:ixSaturation-1) = 0._dp
  trTotal(1:ixSaturation-1) = 0._dp
  trSoil(1:ixSaturation-1)  = 0._dp
 endif

 ! compute the outflow from each layer (m3 s-1)
 mLayerColumnOutflow(1:nSoil) = trSoil(1:nSoil)*tan_slope*contourLength

 ! compute the baseflow in each layer (m s-1)
 mLayerBaseflow(1:nSoil) = (mLayerColumnOutflow(1:nSoil) - mLayerColumnInflow(1:nSoil))/HRUarea

 ! test
 if(printFlag)then
  xDepth = sum(mLayerDepth(ixSaturation:nSoil)*(mLayerVolFracLiq(ixSaturation:nSoil) - fieldCapacity))/activePorosity  ! "effective" water table thickness (m)
  xTran  = tran0*(xDepth/soilDepth)**zScale_TOPMODEL  ! transmissivity for the entire aquifer (m2 s-1)
  xFlow  = xTran*tan_slope*contourLength/HRUarea   ! total column outflow (m s-1)
  write(*,'(a,1x,5(e20.10,1x))') 'xDepth, zActive(ixSaturation) = ', xDepth, zActive(ixSaturation)
  write(*,'(a,1x,5(e20.10,1x))') 'xTran, trTotal(ixSaturation)  = ', xTran, trTotal(ixSaturation)
  write(*,'(a,1x,5(e20.10,1x))') 'xFlow, totalColumnOutflow     = ', xFlow, sum(mLayerColumnOutflow(:))/HRUarea
 endif

 ! ***********************************************************************************************************************
 ! (2) compute the derivative in the baseflow flux w.r.t. volumetric liquid water content (m s-1)
 ! ***********************************************************************************************************************

 ! initialize the derivative matrix
 dBaseflow_dMatric(:,:) = 0._dp

 ! check if derivatives are actually required
 if(.not.derivDesired) return

 ! compute ratio of hillslope width to hillslope area (m m-2)
 length2area = tan_slope*contourLength/HRUarea

 ! compute the ratio of layer depth to maximum water holding capacity (-)
 depth2capacity(1:nSoil) = mLayerDepth(1:nSoil)/(activePorosity*soilDepth)

 ! compute the change in dimensionless flux w.r.t. change in dimensionless storage (-)
 dXdS(1:nSoil) = zScale_TOPMODEL*(zActive(1:nSoil)/SoilDepth)**(zScale_TOPMODEL - 1._dp)

 ! loop through soil layers
 do iLayer=1,nSoil
  ! compute diagonal terms (s-1)
  dBaseflow_dMatric(iLayer,iLayer) = tran0*dXdS(iLayer)*depth2capacity(iLayer)*length2area
  ! compute off-diagonal terms
  do jLayer=iLayer+1,nSoil  ! (only dependent on layers below)
   dBaseflow_dVolLiq(iLayer,jLayer) = (tran0/soilDepth)*(dXdS(iLayer) - dXdS(iLayer+1))
  end do  ! looping through soil layers
 end do  ! looping through soil layers

 ! end association to data in structures
 end associate

 end subroutine computeBaseflow

 ! -- end of private subroutines

end module groundwatr_module

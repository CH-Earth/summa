module groundwatr_module
USE nrtype
implicit none
private
public::groundwatr
public::waterTableHeight
contains

 ! ************************************************************************************************
 ! new subroutine: compute water table depth and baseflow
 ! ************************************************************************************************
 subroutine groundwatr(&
                       ! input (state variables)
                       scalarWaterTableDepthIter,                       & ! input: trial value of water table depth (m)
                       mLayerVolFracLiqIter,                            & ! input: trial value of volumetric fraction of liquid water in each soil layer (-)
                       ! input (depth-weighted hydraulic conductivity)
                       mLayerHydCond,                                   & ! input: hydraulic conductivity in each soil layer (m s-1)
                       mLayerDepth,                                     & ! input: depth of each soil layer (m)
                       ! input (diagnostic variables and parameters)
                       theta_sat,                                       & ! input: porosity (-)
                       k_surf,                                          & ! input: saturated hydraulic conductivity at the surface (m s-1)
                       kAnisotropic,                                    & ! input: anisotropy factor for lateral hydraulic conductivity (-)
                       specificYield,                                   & ! input: specific yield (-)
                       zScale_TOPMODEL,                                 & ! input: scale factor for TOPMODEL-ish baseflow parameterization (m)
                       ! output (states)
                       scalarWaterTableDepthNew,                        & ! output: water table depth at the end of the time step (m)
                       mLayerBaseflow,                                  & ! output: baseflow from each soil layer (m s-1)
                       scalarBaseflow,                                  & ! output: total baseflow (m s-1)
                       ! output: error control
                       err,message)                                       ! output: error control
 implicit none
 ! -----------------------------------------------------------------------------------------------------------------------------------------
 ! input (state variables)
 real(dp),intent(in)                 :: scalarWaterTableDepthIter ! trial value of water table depth (m)
 real(dp),intent(in)                 :: mLayerVolFracLiqIter      ! trial value of volumetric fraction of liquid water in each soil layer (-)
 ! input (depth-weighted hydraulic conductivity)
 real(dp),intent(in)                 :: mLayerHydCond(:)          ! hydraulic conductivity in each soil layer (m s-1)
 real(dp),intent(in)                 :: mLayerDepth(:)            ! depth of each soil layer (m)
 ! input (diagnostic variables and parameters)
 real(dp),intent(in)                 :: theta_sat                 ! porosity (-)
 real(dp),intent(in)                 :: k_surf                    ! saturated hydraulic conductivity at the surface (m s-1) 
 real(dp),intent(in)                 :: kAnisotropic              ! anisotropy factor for lateral hydraulic conductivity (-)
 real(dp),intent(in)                 :: specificYield             ! fraction of water volume drained by gravity in an unconfined aquifer (-)
 real(dp),intent(in)                 :: zScale_TOPMODEL           ! scale factor for TOPMODEL-ish baseflow parameterization (m)
 ! output (states)
 real(dp),intent(out)                :: scalarWaterTableDepthNew  ! water table depth at the end of the time step (m)
 real(dp),intent(out)                :: mLayerBaseflow            ! baseflow from each soil layer (m s-1)
 real(dp),intent(out)                :: scalarBaseflow            ! total baseflow (m s-1)
 ! output: error control
 integer(i4b),intent(out)            :: err                       ! error code
 character(*),intent(out)            :: message                   ! error message
 ! -----------------------------------------------------------------------------------------------------------------------------------------
 ! algorithmic control parameters
 integer(i4b),parameter              :: maxiter=10                ! maximum number of iterations
 real(dp),parameter                  :: incTol= 1.e-6_dp          ! convergence tolerance for the iteration increment (m)
 real(dp),parameter                  :: resTol= 1.e-6_dp          ! convergence tolerance for the residual (m)
 ! define local variables
 character(len=256)                  :: cmessage                  ! error message for downwind routine
 integer(i4b)                        :: iter                      ! iteration index
 integer(i4b)                        :: iSoil                     ! index of each soil layer
 integer(i4b)                        :: nSoil                     ! number of soil layers
 integer(i4b)                        :: nUnsat                    ! number of unsaturated layers
 integer(i4b)                        :: nUnsat_init               ! number of unsaturated layers at the start of the iterations
 real(dp)                            :: baseflowMax               ! maximum baseflow rate from the aquifer (m s-1)
 real(dp),dimension(0:size(mLayerVolFracLiqIter)) :: waterReq2FillPore  ! water required to fill pore space (m)
 real(dp)                            :: drainablePorosity         ! drainable porosity (-)
 real(dp)                            :: aquiferStorageTrial       ! trial value of aquifer storage (m) -- local variable to preserve input
 real(dp)                            :: zWaterTrial               ! depth of the water table (m)
 real(dp)                            :: dBaseflow_dStorage        ! derivative in the baseflow term w.r.t. aquifer storage (s-1)
 real(dp)                            :: netFlux                   ! recharge minus baseflow (m s-1)
 real(dp)                            :: residual                  ! residual in aquifer storage (m)
 real(dp)                            :: dStorage                  ! iteration increment for aquifer storage (m)
 integer(i4b)                        :: idiff                     ! trial integers used to find the bracket
 integer(i4b),parameter              :: maxdiff=100               ! maximum trial values to find the bracket
 real(dp),parameter                  :: xdiff=0.00001_dp          ! increments used to find the bracket
 real(dp)                            :: xbeg,xend                 ! lower and upper points in the bracket
 real(dp)                            :: rbeg,rend                 ! residuals for the lower and upper points in the bracket
 integer(i4b)                        :: itest                     ! trial integers used to test along a line
 integer(i4b),parameter              :: ntest=100000              ! number of values to test along a line
 real(dp)                            :: xinc                      ! increment in storage to test along a line
 real(dp)                            :: xtry,rtry                 ! trial value and residual
 real(dp)                            :: zWaterTrial_init          ! initial value of the water table depth


 

 ! initialize error control
 err=0; message="groundwatr/"

 ! identify the number of soil layers
 nSoil = size(mLayerDepth)

 ! compute the total water deficit (m)
 waterDeficit = sum( (theta_sat - mLayerVolFracLiqIter(:)) * mLayerDepth(:) )






 end subroutine groundwatr


end module groundwatr_module

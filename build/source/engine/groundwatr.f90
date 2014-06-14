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
! look-up values for method used to compute derivative
USE mDecisions_module,only:  &
 numerical,                  & ! numerical solution
 analytical                    ! analytical solution
! look-up values for the choice of groundwater parameterization
USE mDecisions_module,only:  &
 qbaseTopmodel,              & ! TOPMODEL-ish baseflow parameterization
 bigBucket,                  & ! a big bucket (lumped aquifer model)
 noExplicit                    ! no explicit groundwater parameterization
! look-up values for the choice of groundwater representation (local-column, or single-basin)
USE mDecisions_module,only:  &
 localColumn,                & ! separate groundwater representation in each local soil column
 singleBasin                   ! single groundwater store over the entire basin
implicit none
! constant parameters
real(dp),parameter     :: valueMissing=-9999._dp    ! missing value parameter
real(dp),parameter     :: verySmall=epsilon(1.0_dp) ! a very small number (used to avoid divide by zero)
real(dp),parameter     :: dx=1.e-8_dp               ! finite difference increment
private
public::satStorage
public::soilBsFlow
public::disaggFlow
public::aquifrFlux
public::groundwatr
contains

 ! ************************************************************************************************
 ! new subroutine: compute the "saturated" storage at the bottom of the soil profile
 ! ************************************************************************************************
 subroutine satStorage(&
                       ! input
                       mLayerVolFracLiqTrial,      & ! intent(in): volumetric liquid water content in each layer (-)
                       mLayerVolFracIceTrial,      & ! intent(in): volumetric ice content in each layer (-)
                       ! output
                       ixSaturation,               & ! intent(out): index of the lowest saturated layer
                       subSurfaceStorage,          & ! intent(out): sub surface storage (m)
                       maximumSoilWater,           & ! intent(out): maximum storage (m)
                       err,message)                  ! intent(out): error control
 ! ----------------------------------------------------------------------------------------------------------
 ! model variables, parameters, forcing data, etc.
 USE data_struc,only:attr_data,type_data,mpar_data,forc_data,mvar_data,indx_data    ! data structures
 USE var_lookup,only:iLookATTR,iLookTYPE,iLookPARAM,iLookFORCE,iLookMVAR,iLookINDEX ! named variables for structure elements
 implicit none
 ! input
 real(dp),intent(in)          :: mLayerVolFracLiqTrial(:)   ! volumetric liquid water content (-)
 real(dp),intent(in)          :: mLayerVolFracIceTrial(:)   ! volumetric ice content (-)
 ! output
 integer(i4b),intent(out)     :: ixSaturation               ! index of the lowest saturated layer
 real(dp),intent(out)         :: subSurfaceStorage          ! sub surface storage (m)
 real(dp),intent(out)         :: maximumSoilWater           ! maximum storage (m)
 integer(i4b),intent(out)     :: err                        ! error code
 character(*),intent(out)     :: message                    ! error message
 ! local
 character(LEN=256)           :: cmessage                   ! error message of downwind routine
 ! ----------------------------------------------------------------------------------------------------------
 ! initialize error control
 err=0; message='satStorage/'
 call satStorage_muster(&
                        ! input
                        mLayerVolFracLiqTrial,                     & ! intent(in): volumetric liquid water content in each layer (-)
                        mLayerVolFracIceTrial,                     & ! intent(in): volumetric ice content in each layer (-)
                        mvar_data%var(iLookMVAR%mLayerDepth)%dat,  & ! intent(in): depth of each layer (m)
                        mpar_data%var(iLookPARAM%fieldCapacity),   & ! intent(in): field capacity (-)
                        mpar_data%var(iLookPARAM%theta_sat),       & ! intent(in): soil porosity (-)
                        ! output
                        ixSaturation,                              & ! intent(out): index of the lowest saturated layer
                        subSurfaceStorage,                         & ! intent(out): sub surface storage (m)
                        maximumSoilWater,                          & ! intent(out): maximum storage (m)
                        err,cmessage)                                ! intent(out): error control
 if(err/=0)then; message=trim(message)//trim(cmessage); return; endif
 end subroutine satStorage


 ! ************************************************************************************************
 ! new subroutine: compute the flux and derivative
 ! ************************************************************************************************
 subroutine soilBsFlow(&
                       ! input: model control
                       iter,                      & ! intent(in): iteration index
                       ! input: storage
                       subSurfaceStorage,         & ! intent(in): sub surface storage (m)
                       maximumSoilWater,          & ! intent(in): maximum possible sub surface storage (m)
                       ! input/output: diagnostic variables and fluxes constant over iterations
                       maximumFlowRate,           & ! intent(inout): flow rate under saturated conditions (m/s)
                       totalColumnInflow,         & ! intent(inout): total column inflow (m/s)
                       ! output: outflow and its derivative w.r.t. storage
                       totalColumnOutflow,        & ! intent(out): total outflow from the soil column (m s-1)
                       totalOutflowDeriv,         & ! intent(out): derivative in total outflow w.r.t. storage (s-1)
                       ! output: error control
                       err,message)                 ! intent(out): error control
 ! model variables, parameters, forcing data, etc.
 USE data_struc,only:attr_data,type_data,mpar_data,forc_data,mvar_data,indx_data    ! data structures
 USE var_lookup,only:iLookATTR,iLookTYPE,iLookPARAM,iLookFORCE,iLookMVAR,iLookINDEX ! named variables for structure elements
 ! ----------------------------------------------------------------------------------------------------------
 implicit none
 ! input: model control
 integer(i4b),intent(in)      :: iter               ! iteration index
 ! input: storage
 real(dp),intent(in)          :: subSurfaceStorage  ! sub surface storage (m)
 real(dp),intent(in)          :: maximumSoilWater   ! maximum possible sub surface storage (m)
 ! input/output: diagnostic variables and fluxes constant over iterations
 real(dp),intent(inout)       :: maximumFlowRate    ! flow rate under saturated conditions (m/s)
 real(dp),intent(inout)       :: totalColumnInflow  ! total column inflow (m/s)
 ! output: outflow and its derivative w.r.t. storage
 real(dp),intent(out)         :: totalColumnOutflow ! total outflow from the soil column (m s-1)
 real(dp),intent(out)         :: totalOutflowDeriv  ! derivative in total outflow w.r.t. storage (s-1)
 ! output: error control
 integer(i4b),intent(out)     :: err                ! error code
 character(*),intent(out)     :: message            ! error message
 ! local
 character(LEN=256)           :: cmessage           ! error message of downwind routine
 ! ----------------------------------------------------------------------------------------------------------
 ! initialize error control
 err=0; message='soilBsFlow/'
 call soilBsFlow_muster(&
                        ! input: model control
                        iter,                                              & ! intent(in): iteration index
                        ! input: storage
                        subSurfaceStorage,                                 & ! intent(in): sub surface storage (m)
                        maximumSoilWater,                                  & ! intent(in): maximum possible sub surface storage (m)
                        ! input: inflow
                        mvar_data%var(iLookMVAR%mLayerColumnInflow)%dat,   & ! intent(in): inflow into each soil layer (m3/s)
                        ! input/output: diagnostic variables and fluxes constant over iterations
                        maximumFlowRate,                                   & ! intent(inout): flow rate under saturated conditions (m/s)
                        totalColumnInflow,                                 & ! intent(inout): total column inflow (m/s)
                        ! input: storage and transmission properties
                        mvar_data%var(iLookMVAR%mLayerSatHydCondMP)%dat(0),& ! intent(in): saturated hydraulic conductivity at the surface (m s-1)
                        mpar_data%var(iLookPARAM%kAnisotropic),            & ! intent(in): anisotropy factor for lateral hydraulic conductivity (-)
                        mpar_data%var(iLookPARAM%fieldCapacity),           & ! intent(in): field capacity (-)
                        mpar_data%var(iLookPARAM%theta_sat),               & ! intent(in): soil porosity (-)
                        mpar_data%var(iLookPARAM%zScale_TOPMODEL),         & ! intent(in): TOPMODEL exponent (-)
                        ! input: local attributes
                        attr_data%var(iLookATTR%HRUarea),                  & ! intent(in): HRU area (m2)
                        attr_data%var(iLookATTR%tan_slope),                & ! intent(in): tan water table slope, taken as tan local ground surface slope (-)
                        attr_data%var(iLookATTR%contourLength),            & ! intent(in): length of contour at downslope edge of HRU (m)
                        ! output: outflow and its derivative w.r.t. storage
                        totalColumnOutflow,                                & ! intent(out): total outflow from the soil column (m s-1)
                        totalOutflowDeriv,                                 & ! intent(out): derivative in total outflow w.r.t. storage (s-1)
                        ! output: error control
                        err,cmessage)                                        ! intent(out): error control
 if(err/=0)then; message=trim(message)//trim(cmessage); return; endif
 end subroutine soilBsFlow


 ! ************************************************************************************************
 ! new subroutine: disaggregate total inflow and total outflow to the baseflow sink term
 ! ************************************************************************************************
 subroutine disaggFlow(&
                       ! input
                       dt,                         & ! intent(in): time step (s) -- used to calculate maximum possible inflow rate
                       ixSaturation,               & ! intent(inout): index of the lowest saturated layer
                       totalColumnInflow,          & ! intent(in): total column inflow (m s-1)
                       totalColumnOutflow,         & ! intent(in): total outflow from the soil column (m s-1)
                       mLayerVolFracLiq,           & ! intent(in): volumetric fraction of liquid water in each soil layer (-)
                       mLayerVolFracIce,           & ! intent(in): volumetric fraction of ice in each soil layer (-)
                       mLayerHydCond,              & ! intent(in): hydraulic conductivity in each soil layer (m s-1)
                       ! output
                       mLayerBaseflow,             & ! intent(out): baseflow from each soil layer (m s-1)
                       err,message)                  ! intent(out): error control
 ! model variables, parameters, forcing data, etc.
 USE data_struc,only:attr_data,type_data,mpar_data,forc_data,mvar_data,indx_data    ! data structures
 USE var_lookup,only:iLookATTR,iLookTYPE,iLookPARAM,iLookFORCE,iLookMVAR,iLookINDEX ! named variables for structure elements
 ! ----------------------------------------------------------------------------------------------------------
 implicit none
 ! input
 real(dp),intent(in)          :: dt                         ! time step (s) -- used to calculate maximum possible inflow rate
 integer(i4b),intent(inout)   :: ixSaturation               ! index of the lowest saturated layer
 real(dp),intent(in)          :: totalColumnInflow          ! total column inflow (m s-1)
 real(dp),intent(in)          :: totalColumnOutflow         ! total outflow from the soil column (m s-1)
 real(dp),intent(in)          :: mLayerVolFracLiq(:)        ! volumetric fraction of liquid water in each soil layer (-)
 real(dp),intent(in)          :: mLayerVolFracIce(:)        ! volumetric fraction of ice in each soil layer (-)
 real(dp),intent(in)          :: mLayerHydCond(:)           ! hydraulic conductivity in each soil layer (m s-1)
 ! output
 real(dp),intent(out)         :: mLayerBaseflow(:)          ! baseflow in each layer (m s-1)
 integer(i4b),intent(out)     :: err                        ! error code
 character(*),intent(out)     :: message                    ! error message
 ! local
 character(LEN=256)           :: cmessage           ! error message of downwind routine
 ! ----------------------------------------------------------------------------------------------------------
 ! initialize error control
 err=0; message='disaggFlow/'
 call disaggFlow_muster(&
                        ! input: model control
                        dt,                                                          & ! intent(in): time step (s) -- used to calculate maximum possible inflow rate
                        ixSaturation,                                                & ! intent(inout): index of the lowest saturated layer
                        ! input: total column inflow and outflow
                        totalColumnInflow,                                           & ! intent(in): total column inflow (m s-1)
                        totalColumnOutflow,                                          & ! intent(in): total outflow from the soil column (m s-1)
                        ! input: hydraulic conductivity in each soil layer
                        mLayerVolFracLiq,                                            & ! intent(in): volumetric fraction of liquid water in each soil layer (-)
                        mLayerVolFracIce,                                            & ! intent(in): volumetric fraction of ice in each soil layer (-)
                        mLayerHydCond,                                               & ! intent(in): hydraulic conductivity in each soil layer (m s-1)
                        ! input: attributes and parameters
                        attr_data%var(iLookATTR%HRUarea),                            & ! intent(in): HRU area (m2)
                        mpar_data%var(iLookPARAM%theta_sat),                         & ! intent(in): soil porosity (-)
                        mpar_data%var(iLookPARAM%fieldCapacity),                     & ! intent(in): field capacity (-)
                        ! input: storage and transmission properties in each layer
                        mvar_data%var(iLookMVAR%mLayerDepth)%dat,                    & ! intent(in): depth of each layer (m)
                        mvar_data%var(iLookMVAR%mLayerColumnInflow)%dat(1:nSoil),    & ! total inflow to each layer in the soil column (m3 s-1)
                        ! output: exfiltration and column outflow
                        mvar_data%var(iLookMVAR%scalarExfiltration)%dat(1),          & ! intent(out): (scalar) exfiltration at the end-of-step (m s-1)
                        mvar_data%var(iLookMVAR%mLayerColumnOutflow)%dat(1:nSoil),   & ! intent(out): total outflow from each layer of the soil column (m3 s-1)
                        ! output: baseflow sink
                        mLayerBaseflow,                                              & ! intent(out): baseflow from each soil layer (m s-1)
                        ! output: error control
                        err,cmessage)                                                  ! intent(out): error control
 if(err/=0)then; message=trim(message)//trim(cmessage); return; endif
 end subroutine disaggFlow



 ! ************************************************************************************************
 ! new subroutine: compute fluxes for the aquifer
 ! ************************************************************************************************
 subroutine aquifrFlux(&
                       ! input
                       ixGroundwater,                             & ! intent(in): index defining the choice of groundwater parameterization
                       scalarAquiferStorageTrial,                 & ! intent(in): aquifer storage (m)
                       scalarCanopyTranspiration,                 & ! intent(in): canopy transpiration (kg m-2 s-1)
                       scalarSoilDrainage,                        & ! intent(in): drainage from the bottom of the soil profile
                       ! output
                       scalarAquiferTranspire,                    & ! intent(out): transpiration loss from the aquifer (m s-1)
                       scalarAquiferRecharge,                     & ! intent(out): recharge to the aquifer (m s-1)
                       scalarAquiferBaseflow,                     & ! intent(out): total baseflow from the aquifer (m s-1)
                       scalarAquiferBaseflowDeriv,                & ! intent(out): derivative in baseflow from the aquifer w.r.t. aquifer storage (s-1)
                       err,message)                                 ! intent(out): error control
 ! model decisions
 USE data_struc,only:model_decisions                                ! model decision structure
 USE var_lookup,only:iLookDECISIONS                                 ! named variables for elements of the decision structure
 ! model variables, parameters, forcing data, etc.
 USE data_struc,only:attr_data,type_data,mpar_data,forc_data,mvar_data,indx_data    ! data structures
 USE var_lookup,only:iLookATTR,iLookTYPE,iLookPARAM,iLookFORCE,iLookMVAR,iLookINDEX ! named variables for structure elements
 ! ----------------------------------------------------------------------------------------------------------
 implicit none
 ! input
 integer(i4b),intent(in)          :: ixGroundwater                  ! index defining the choice of groundwater parameterization
 real(dp),intent(in)              :: scalarAquiferStorageTrial      ! aquifer storage (m)
 real(dp),intent(in)              :: scalarCanopyTranspiration      ! canopy transpiration (kg m-2 s-1)
 real(dp),intent(in)              :: scalarSoilDrainage             ! drainage from the bottom of the soil profile (m s-1)
 ! output
 real(dp),intent(out)             :: scalarAquiferTranspire         ! transpiration loss from the aquifer (m s-1)
 real(dp),intent(out)             :: scalarAquiferRecharge          ! recharge to the aquifer (m s-1)
 real(dp),intent(out)             :: scalarAquiferBaseflow          ! total baseflow from the aquifer (m s-1)
 real(dp),intent(out)             :: scalarAquiferBaseflowDeriv     ! derivative in baseflow from the aquifer w.r.t. aquifer storage (s-1)
 integer(i4b),intent(out)         :: err                            ! error code
 character(*),intent(out)         :: message                        ! error message
 ! local
 character(LEN=256)               :: cmessage                       ! error message of downwind routine
 ! ----------------------------------------------------------------------------------------------------------
 ! initialize error control
 err=0; message='aquifrFlux/'
 call aquifrFlux_muster(&
                        ! input: model decisions
                        ixGroundwater,                                          & ! intent(in): index defining the choice of groundwater parameterization
                        ! input: model state variable
                        scalarAquiferStorageTrial,                              & ! intent(in): aquifer storage (m)
                        scalarCanopyTranspiration,                              & ! intent(in): canopy transpiration (kg m-2 s-1)
                        scalarSoilDrainage,                                     & ! intent(in): drainage from the bottom of the soil profile
                        ! input: factors limiting transpiration (from vegFlux routine)
                        mvar_data%var(iLookMVAR%scalarTranspireLim)%dat(1),     & ! intent(in): weighted average of the transpiration limiting factor (-)
                        mvar_data%var(iLookMVAR%scalarAquiferRootFrac)%dat(1),  & ! intent(in): fraction of roots below the lowest soil layer (-)
                        mvar_data%var(iLookMVAR%scalarTranspireLimAqfr)%dat(1), & ! intent(in): transpiration limiting factor for the aquifer (-)
                        ! input: factors controlling baseflow
                        mpar_data%var(iLookPARAM%kAnisotropic),                 & ! intent(in): anisotropy factor for lateral hydraulic conductivity (-)
                        mvar_data%var(iLookMVAR%iLayerSatHydCond)%dat(0),       & ! intent(in): saturated hydraulic conductivity at the surface (m s-1)
                        mpar_data%var(iLookPARAM%aquiferScaleFactor),           & ! intent(in): scaling factor for aquifer storage in the big bucket (m)
                        mpar_data%var(iLookPARAM%aquiferBaseflowExp),           & ! intent(in): baseflow exponent for the big bucket (-)
                        ! output: model fluxes and derivatives
                        scalarAquiferTranspire,                                 & ! intent(out): transpiration loss from the aquifer (m s-1)
                        scalarAquiferRecharge,                                  & ! intent(out): recharge to the aquifer (m s-1)
                        scalarAquiferBaseflow,                                  & ! intent(out): total baseflow from the aquifer (m s-1)
                        scalarAquiferBaseflowDeriv,                             & ! intent(out): derivative in baseflow from the aquifer w.r.t. aquifer storage (s-1)
                        ! output: error control
                        err,cmessage)                                             ! intent(out): error control
 if(err/=0)then; message=trim(message)//trim(cmessage); return; endif
 end subroutine aquifrFlux




 ! ************************************************************************************************
 ! new subroutine: compute water balance for the basin-wide aquifer (split from other routines)
 ! ************************************************************************************************
 subroutine groundwatr(&
                       ! input: model control
                       dt,                   & ! intent(in): time step (s)
                       ixDerivMethod,        & ! intent(in): method used to calculate derivatives
                       ! input: effective parameters
                       aquiferHydCond,       & ! intent(in): effective hydraulic conductivity (m s-1)
                       aquiferScaleFactor,   & ! intent(in): scaling factor for aquifer storage (m)
                       aquiferBaseflowExp,   & ! intent(in): exponent in bucket baseflow parameterization (-)
                       ! input: aquifer fluxes
                       aquiferRecharge,      & ! intent(in): aquifer recharge (m s-1)
                       aquiferTranspire,     & ! intent(in): aquifer transpiration (m s-1)
                       ! input-output
                       aquiferStorage,       & ! intent(inout): aquifer storage (m)
                       ! output
                       aquiferBaseflow,      & ! intent(out): aquifer baseflow (m s-1)
                       err,message)            ! intent(out): error control
 implicit none
 ! -----------------------------------------------------------------------------------------------------------------------------------------
 ! input: model control
 real(dp),intent(in)              :: dt                      ! time step (s)
 integer(i4b),intent(in)          :: ixDerivMethod           ! method used to calculate derivatives
 ! input: effective parameters
 real(dp),intent(in)              :: aquiferHydCond          ! effective hydraulic conductivity (m s-1)
 real(dp),intent(in)              :: aquiferScaleFactor      ! scaling factor for aquifer storage (m)
 real(dp),intent(in)              :: aquiferBaseflowExp      ! exponent in bucket baseflow parameterization (-)
 ! input: aquifer fluxes
 real(dp),intent(in)              :: aquiferRecharge         ! aquifer recharge (m s-1)
 real(dp),intent(in)              :: aquiferTranspire        ! aquifer transpiration (m s-1)
 ! input-output
 real(dp),intent(inout)           :: aquiferStorage          ! aquifer storage (m)
 ! output
 real(dp),intent(out)             :: aquiferBaseflow         ! aquifer baseflow (m s-1)
 integer(i4b),intent(out)         :: err                     ! error code
 character(*),intent(out)         :: message                 ! error message
 ! -----------------------------------------------------------------------------------------------------------------------------------------
 ! local variables
 character(LEN=256)               :: cmessage                ! error message of downwind routine
 integer(i4b)                     :: iter                    ! iteration index
 integer(i4b),parameter           :: maxiter=20              ! maximum number of iterations
 real(dp)                         :: aquiferStorageTrial     ! trial value of aquifer storage
 real(dp)                         :: scalarBaseflow          ! baseflow (m s-1)
 real(dp)                         :: scalarBaseflowDeriv     ! derivative in baseflow w.r.t. aquifer storage (s-1)
 real(dp)                         :: res                     ! residual in water balance (m)
 real(dp)                         :: aquiferIncr             ! iteration increment (m)
 real(dp),parameter               :: tolRes=1.e-8_dp         ! convergence tolerance for the residual
 real(dp),parameter               :: tolInc=1.e-10_dp        ! convergence tolerance for the iteration increment
 ! -----------------------------------------------------------------------------------------------------------------------------------------
 ! initialize error control
 err=0; message='groundwatr/'

 ! initialize state variables
 aquiferStorageTrial   = aquiferStorage ! initialize as start-of-step value

 ! iterate
 do iter=1,maxiter

  ! compute baseflow
  call q_baseflow(&
                  ! input: model decisions
                  .true.,               & ! intent(in): flag indicating if derivatives are desired
                  ixDerivMethod,        & ! intent(in): method used to calculate derivatives
                  ! input: effective parameters
                  aquiferHydCond,       & ! intent(in): effective hydraulic conductivity (m s-1)
                  aquiferScaleFactor,   & ! intent(in): scaling factor for aquifer storage (m)
                  aquiferBaseflowExp,   & ! intent(in): exponent in bucket baseflow parameterization (-)
                  ! input: state variables
                  aquiferStorageTrial,  & ! intent(in): choice of groundwater parameterization
                  ! output
                  scalarBaseflow,       & ! intent(out): total baseflow (m s-1)
                  scalarBaseflowDeriv,  & ! intent(out): derivative in baseflow flux w.r.t. water table depth (m s-1)
                  err,message)            ! intent(out): error control
  if(err/=0)then; message=trim(message)//trim(cmessage); return; endif

  ! compute the residual
  res = (aquiferRecharge + aquiferTranspire - scalarBaseflow)*dt - (aquiferStorageTrial - aquiferStorage)

  ! compute the iteration increment
  aquiferIncr = res/(1._dp + scalarBaseflowDeriv*dt)

  ! print progress
  !print*, 'scalarBaseflowDeriv = ', scalarBaseflowDeriv
  !write(*,'(a,i4,1x,2(f20.10,1x),2(e20.10,1x))') 'iter, aquiferStorageTrial, scalarBaseflow, res, aquiferIncr = ', &
  !                                                iter, aquiferStorageTrial, scalarBaseflow, res, aquiferIncr

  ! update the aquifer
  aquiferStorageTrial = aquiferStorageTrial + aquiferIncr
  !if(iter > 10) pause

  ! check convergence
  if(res < tolRes .or. aquiferIncr < tolInc) exit

  ! check that we converged
  if(iter==maxiter)then; err=20; message=trim(message)//'failed to converge'; return; endif

 end do ! iterating

 ! return aquifer storaage and baseflow
 aquiferStorage  = aquiferStorageTrial
 aquiferBaseflow = scalarBaseflow

 end subroutine groundwatr


 ! ************************************************************************************************
 ! ************************************************************************************************
 ! ***** PRIVATE SUBROUTINES **********************************************************************
 ! ************************************************************************************************
 ! ************************************************************************************************

 ! ************************************************************************************************
 ! private subroutine: compute fluxes for the aquifer
 ! ************************************************************************************************
 subroutine aquifrFlux_muster(&
                              ! input: model decisions
                              ixGroundwater,                      & ! intent(in): index defining the choice of groundwater parameterization
                              ! input: model state variable
                              scalarAquiferStorageTrial,          & ! intent(in): aquifer storage (m)
                              scalarCanopyTranspiration,          & ! intent(in): canopy transpiration (kg m-2 s-1)
                              scalarSoilDrainage,                 & ! intent(in): drainage from the bottom of the soil profile
                              ! input: factors limiting transpiration (from vegFlux routine)
                              scalarTranspireLim,                 & ! intent(in): weighted average of the transpiration limiting factor (-)
                              scalarAquiferRootFrac,              & ! intent(in): fraction of roots below the lowest soil layer (-)
                              scalarTranspireLimAqfr,             & ! intent(in): transpiration limiting factor for the aquifer (-)
                              ! input: factors controlling baseflow
                              kAnisotropic,                       & ! intent(in): anisotropy factor for lateral hydraulic conductivity (-)
                              surfaceHydCond,                     & ! intent(in): hydraulic conductivity at the surface (m s-1)
                              aquiferScaleFactor,                 & ! intent(in): scaling factor for aquifer storage in the big bucket (m)
                              aquiferBaseflowExp,                 & ! intent(in): baseflow exponent for the big bucket (-)
                              ! output: model fluxes and derivatives
                              scalarAquiferTranspire,             & ! intent(out): transpiration loss from the aquifer (m s-1)
                              scalarAquiferRecharge,              & ! intent(out): recharge to the aquifer (m s-1)
                              scalarAquiferBaseflow,              & ! intent(out): total baseflow from the aquifer (m s-1)
                              scalarAquiferBaseflowDeriv,         & ! intent(out): derivative in baseflow from the aquifer w.r.t. aquifer storage (s-1)
                              ! output: error control
                              err,message)                          ! intent(out): error control
 implicit none
 ! input: model decisions
 integer(i4b),intent(in)          :: ixGroundwater                  ! index defining the choice of groundwater parameterization
 ! input: model state variable
 real(dp),intent(in)              :: scalarAquiferStorageTrial      ! aquifer storage (m)
 real(dp),intent(in)              :: scalarCanopyTranspiration      ! canopy transpiration (kg m-2 s-1)
 real(dp),intent(in)              :: scalarSoilDrainage             ! drainage from the bottom of the soil profile
 ! input: factors limiting transpiration (from vegFlux routine)
 real(dp),intent(in)              :: scalarTranspireLim             ! weighted average of the transpiration limiting factor (-)
 real(dp),intent(in)              :: scalarAquiferRootFrac          ! fraction of roots below the lowest soil layer (-)
 real(dp),intent(in)              :: scalarTranspireLimAqfr         ! transpiration limiting factor for the aquifer (-)
 ! input: factors controlling baseflow
 real(dp),intent(in)              :: kAnisotropic                   ! anisotropy factor for lateral hydraulic conductivity (-)
 real(dp),intent(in)              :: surfaceHydCond                 ! hydraulic conductivity at the surface (m s-1)
 real(dp),intent(in)              :: aquiferScaleFactor             ! scaling factor for aquifer storage in the big bucket (m)
 real(dp),intent(in)              :: aquiferBaseflowExp             ! baseflow exponent for the big bucket (-)
 ! output: model fluxes and derivatives
 real(dp),intent(out)             :: scalarAquiferTranspire         ! transpiration loss from the aquifer (m s-1)
 real(dp),intent(out)             :: scalarAquiferRecharge          ! recharge to the aquifer (m s-1)
 real(dp),intent(out)             :: scalarAquiferBaseflow          ! total baseflow from the aquifer (m s-1)
 real(dp),intent(out)             :: scalarAquiferBaseflowDeriv     ! derivative in baseflow from the aquifer w.r.t. aquifer storage (s-1)
 ! output: error control
 integer(i4b),intent(out)         :: err                            ! error code
 character(*),intent(out)         :: message                        ! error message
 ! -------------------------------------------------------------------------------------------------------------------------------
 ! local
 character(LEN=256)               :: cmessage                       ! error message of downwind routine
 real(dp)                         :: aquiferHydCond                 ! hydraulic conductivity of the aquifer (m s-1)
 real(dp)                         :: aquiferTranspireFrac           ! fraction of transpiration from the aquifer (-)
 ! -------------------------------------------------------------------------------------------------------------------------------
 ! initialize error control
 err=0; message='aquifrFlux/'

 ! -------------------------------------------------------------------------------------------------------------------------------------------------
 ! * compute the aquifer transpiration
 ! -------------------------------------------------------------------------------------------------------------------------------------------------

 ! compute the fraction of transpiration loss from the aquifer
 if(scalarTranspireLim > epsilon(scalarTranspireLim))then ! (transpiration may be non-zero even if the soil moisture limiting factor is zero)
  aquiferTranspireFrac = scalarAquiferRootFrac*scalarTranspireLimAqfr/scalarTranspireLim
 else
  aquiferTranspireFrac = scalarAquiferRootFrac
 endif

 ! compute transpiration loss from the aquifer (kg m-2 s-1 --> m s-1)
 scalarAquiferTranspire = aquiferTranspireFrac*scalarCanopyTranspiration/iden_water

 ! -------------------------------------------------------------------------------------------------------------------------------------------------
 ! * compute the recharge to the aquifer
 ! -------------------------------------------------------------------------------------------------------------------------------------------------
 select case(ixGroundwater)

  ! using the big bucket...
  case(bigBucket)
   scalarAquiferRecharge = scalarSoilDrainage  ! recharge = drainage flux from the bottom of the soil profile (m s-1)

  ! no explicit aquifer...
  case(qbaseTopmodel,noExplicit)
   scalarAquiferRecharge = 0._dp

  ! error checking...
  case default; err=20; message=trim(message)//'unknown groundwater parameterization'; return

 end select  ; ! (choice of groundwater parameterization)

 ! -------------------------------------------------------------------------------------------------------------------------------------------------
 ! * compute the baseflow flux and its derivative
 ! -------------------------------------------------------------------------------------------------------------------------------------------------
 select case(ixGroundwater)

  ! the big bucket
  case(bigBucket)

   ! compute effective hydraulic conductivity
   aquiferHydCond = kAnisotropic*scalarSoilDrainage

   ! compute baseflow
   call q_baseflow(&
                   ! input: model decisions
                   .true.,                      & ! intent(in): flag indicating if derivatives are desired
                   analytical,                  & ! intent(in): choice of method used to compute derivative
                   ! input: effective parameters
                   aquiferHydCond,              & ! intent(in): effective hydraulic conductivity (m s-1)
                   aquiferScaleFactor,          & ! intent(in): scaling factor for aquifer storage (m)
                   aquiferBaseflowExp,          & ! intent(in): exponent in bucket baseflow parameterization (-)
                   ! input: state variables
                   scalarAquiferStorageTrial,   & ! intent(in): aquifer storage (m)
                   ! output
                   scalarAquiferBaseflow,       & ! intent(out): total baseflow (m s-1)
                   scalarAquiferBaseflowDeriv,  & ! intent(out): derivative in baseflow flux w.r.t. water table depth (m s-1)
                   err,cmessage)                  ! intent(out): error control
   if(err/=0)then; message=trim(message)//trim(cmessage); return; endif

  ! no representation of the aquifer
  case(noExplicit,qbaseTopmodel)

   scalarAquiferBaseflow      = 0._dp
   scalarAquiferBaseflowDeriv = 0._dp

  ! error check
  case default
   message=trim(message)//'unable to identify baseflow parameterization'
   err=20; return

 end select

 end subroutine aquifrFlux_muster



 ! ************************************************************************************************
 ! private subroutine: compute baseflow
 ! ************************************************************************************************
 subroutine q_baseflow(&
                       ! input: model decisions
                       deriv_desired,        & ! intent(in): flag indicating if derivatives are desired
                       ixDerivMethod,        & ! intent(in): method used to calculate derivatives
                       ! input: effective parameters
                       aquiferHydCond,       & ! intent(in): effective hydraulic conductivity (m s-1)
                       aquiferScaleFactor,   & ! intent(in): scaling factor for aquifer storage (m)
                       aquiferBaseflowExp,   & ! intent(in): exponent in bucket baseflow parameterization (-)
                       ! input: state variables
                       scalarAquiferStorage, & ! intent(in): aquifer storage (m)
                       ! output
                       scalarBaseflow,       & ! intent(out): total baseflow (m s-1)
                       scalarBaseflowDeriv,  & ! intent(out): derivative in baseflow flux w.r.t. water table depth (m s-1)
                       err,message)            ! intent(out): error control
 implicit none
 ! -----------------------------------------------------------------------------------------------------------------------------------------
 ! input: model decisions
 logical(lgt),intent(in)          :: deriv_desired           ! flag indicating if derivatives are desired
 integer(i4b),intent(in)          :: ixDerivMethod           ! method used to calculate derivatives
 ! input: effective parameters
 real(dp),intent(in)              :: aquiferHydCond          ! effective hydraulic conductivity (m s-1)
 real(dp),intent(in)              :: aquiferScaleFactor      ! scaling factor for aquifer storage (m)
 real(dp),intent(in)              :: aquiferBaseflowExp      ! exponent in bucket baseflow parameterization (-)
 ! input: state variables
 real(dp),intent(in)              :: scalarAquiferStorage    ! trial value of aquifer strorage (m)
 ! output
 real(dp),intent(out)             :: scalarBaseflow          ! baseflow (m s-1)
 real(dp),intent(out)             :: scalarBaseflowDeriv     ! derivative in baseflow flux w.r.t. water table depth (m s-1)
 integer(i4b),intent(out)         :: err                     ! error code
 character(*),intent(out)         :: message                 ! error message
 ! -----------------------------------------------------------------------------------------------------------------------------------------
 ! local variables
 real(dp)                         :: aquiferStorageTrial     ! trial value of aquifer strorage (m)
 real(dp)                         :: scaledStorage           ! scaled storage (-)
 integer(i4b)                     :: itry                    ! index of different flux calculations
 integer(i4b)                     :: nFlux                   ! number of flux calculations required (>1 = numerical derivatives)
 integer(i4b),parameter           :: unperturbed=0           ! named variable to identify the case of unperturbed state variables
 integer(i4b),parameter           :: perturbState=1          ! named variable to identify the case where we perturb the state in the current layer
 integer(i4b),parameter           :: perturbStateAbove=2     ! named variable to identify the case where we perturb the state layer above
 integer(i4b),parameter           :: perturbStateBelow=3     ! named variable to identify the case where we perturb the state layer below
 real(dp)                         :: scalarFlux              ! baseflow flux (m s-1)
 real(dp)                         :: scalarFlux_dStateAbove  ! baseflow flux with perturbation to the state above (m s-1)
 ! -----------------------------------------------------------------------------------------------------------------------------------------
 ! initialize error control
 err=0; message='q_baseflow/'
 
 ! identify number of ADDITIONAL flux evaluations (used when computing numerical derivatives)
 if(ixDerivMethod==numerical .and. deriv_desired)then
  nFlux=3   ! NOTE: we cycle through undesired perturbations, so only actually do one additional flux eval
 else
  nFlux=0
 endif

 ! *** loop to compute numerical derivatives
 ! either one or multiple flux calls, depending on if using analytical or numerical derivatives
 do itry=nFlux,0,-1  ! (work backwards to ensure all computed fluxes come from the un-perturbed case)

  ! skip undesired perturbations
  if(itry==perturbState .or. itry==perturbStateBelow) cycle

  ! define the trial value for aquifer storage
  select case(itry)
   case(unperturbed);       aquiferStorageTrial  = scalarAquiferStorage
   case(perturbStateAbove); aquiferStorageTrial  = scalarAquiferStorage + dx
    case(perturbStateBelow,perturbState); err=10; message=trim(message)//'only perturb aquifer storage when computing baseflow flux -- should not get here'; return
   case default; err=10; message=trim(message)//"unknown perturbation"; return
  end select ! (type of perturbation)

  ! compute scaled storage (-)
  scaledStorage  = aquiferStorageTrial/aquiferScaleFactor

  ! compute baseflow
  scalarBaseflow = aquiferHydCond*(scaledStorage**aquiferBaseflowExp)

  ! compute derivative in baseflow
  if(ixDerivMethod==analytical .and. deriv_desired)then
   scalarBaseflowDeriv = (aquiferHydCond/aquiferScaleFactor)*aquiferBaseflowExp*scaledStorage**(aquiferBaseflowExp - 1._dp)
  else
   scalarBaseflowDeriv = valueMissing
  endif

  ! get copies of baseflow flux to compute derivatives
  if(deriv_desired .and. ixDerivMethod==numerical)then
   select case(itry)
    case(unperturbed);       scalarFlux             = scalarBaseflow
    case(perturbStateAbove); scalarFlux_dStateAbove = scalarBaseflow
    case(perturbStateBelow,perturbState); err=10; message=trim(message)//'only perturb aquifer storage when computing baseflow flux -- should not get here'; return
    case default; err=10; message=trim(message)//'unknown perturbation'; return
   end select
  endif

 end do  ! (multiple flux calls for computing  numerical derivatives

 ! * compute derivatives
 ! NOTE: baseflow derivatives w.r.t. state below are *actually* w.r.t. water table depth, so need to be corrected for aquifer storage
 if(deriv_desired)then
  if(ixDerivMethod==numerical) scalarBaseflowDeriv = (scalarFlux_dStateAbove - scalarFlux)/dx
 else
  scalarBaseflowDeriv = valueMissing
 endif

 end subroutine q_baseflow


 ! ************************************************************************************************
 ! private subroutine: compute baseflow fluxes
 ! ************************************************************************************************
 subroutine satStorage_muster(&
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
 err=0; message='satStorage_muster/'

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

 end subroutine satStorage_muster




 ! ************************************************************************************************
 ! private subroutine: compute baseflow fluxes
 ! ************************************************************************************************
 subroutine soilBsFlow_muster(&
                              ! input: model control
                              iter,                     & ! intent(in): iteration index
                              ! input: storage
                              subSurfaceStorage,        & ! intent(in): sub surface storage (m)
                              maximumSoilWater,         & ! intent(in): maximum possible sub surface storage (m)
                              ! input: inflow
                              mLayerColumnInflow,       & ! intent(in): inflow into each soil layer (m3/s)
                              ! input/output: diagnostic variables and fluxes constant over iterations
                              maximumFlowRate,          & ! intent(inout): flow rate under saturated conditions (m/s)
                              totalColumnInflow,        & ! intent(inout): total column inflow (m/s)
                              ! input: storage and transmission properties
                              srfSatHydCond,            & ! intent(in): saturated hydraulic conductivity at the surface (m s-1) 
                              kAnisotropic,             & ! intent(in): anisotropy factor for lateral hydraulic conductivity (-)
                              fieldCapacity,            & ! intent(in): field capacity (-)
                              theta_sat,                & ! intent(in): soil porosity (-)
                              TOPMODELexp,              & ! intent(in): TOPMODEL exponent (-)
                              ! input: local attributes
                              HRUarea,                  & ! intent(in): HRU area (m2)
                              tan_slope,                & ! intent(in): tan water table slope, taken as tan local ground surface slope (-)
                              contourLength,            & ! intent(in): length of contour at downslope edge of HRU (m)
                              ! output: outflow and its derivative w.r.t. storage
                              totalColumnOutflow,       & ! intent(out): total outflow from the soil column (m s-1)
                              totalOutflowDeriv,        & ! intent(out): derivative in total outflow w.r.t. storage (s-1)
                              ! output: error control
                              err,message)                ! intent(out): error control
 ! ----------------------------------------------------------------------------------------------------------
 implicit none
 ! input: model control
 integer(i4b),intent(in)      :: iter                     ! iteration index
 ! input: storage
 real(dp),intent(in)          :: subSurfaceStorage        ! sub surface storage (m)
 real(dp),intent(in)          :: maximumSoilWater         ! maximum possible sub surface storage (m)
 ! input: inflow
 real(dp),intent(in)          :: mLayerColumnInflow(:)    ! inflow into each soil layer (m3/s)
 ! input/output: diagnostic variables and fluxes constant over iterations
 real(dp),intent(inout)       :: maximumFlowRate          ! flow rate under saturated conditions (m/s)
 real(dp),intent(inout)       :: totalColumnInflow        ! total column inflow (m/s)
 ! input: storage and transmission properties
 real(dp),intent(in)          :: srfSatHydCond            ! saturated hydraulic conductivity at the surface (m s-1) 
 real(dp),intent(in)          :: kAnisotropic             ! anisotropy factor for lateral hydraulic conductivity (-)
 real(dp),intent(in)          :: fieldCapacity            ! field capacity (-)
 real(dp),intent(in)          :: theta_sat                ! soil porosity (-)
 real(dp),intent(in)          :: TOPMODELexp              ! TOPMODEL exponent (-)
 ! input: local attributes
 real(dp),intent(in)          :: HRUarea                  ! HRU area (m2)
 real(dp),intent(in)          :: tan_slope                ! tan water table slope, taken as tan local ground surface slope (-)
 real(dp),intent(in)          :: contourLength            ! length of contour at downslope edge of HRU (m)
 ! output: outflow and its derivative w.r.t. storage
 real(dp),intent(out)         :: totalColumnOutflow       ! total outflow from the soil column (m s-1)
 real(dp),intent(out)         :: totalOutflowDeriv        ! derivative in total outflow w.r.t. storage (s-1)
 ! output: error control
 integer(i4b),intent(out)     :: err                      ! error code
 character(*),intent(out)     :: message                  ! error message
 ! ----------------------------------------------------------------------------------------------------------
 ! initialize error control
 err=0; message='soilBsFlow_muster/'

 ! *** preliminaries (only need to do once; intent = inout)...
 if(iter==1)then

  ! compute maximum flow rate -- i.e., the flow rate under saturated conditions (m/s)
  maximumFlowRate  = (1._dp/HRUarea)*tan_slope*srfSatHydCond*kAnisotropic*maximumSoilWater*contourLength &
                       / ((theta_sat - fieldCapacity)*TOPMODELexp)   ! effective hydraulic conductivity (m/s)

  ! compute the total column inflow (m/s)
  totalColumnInflow = sum(mLayerColumnInflow(:))/HRUarea

 endif  ! (if the first iteration)

 ! compute the total column outflow (m/s) and its derivative w.r.t. storage (s-1)
 if(subSurfaceStorage > tiny(theta_sat))then
  totalColumnOutflow = (1._dp/HRUarea)*maximumFlowRate*(subSurfaceStorage/maximumSoilWater)**TOPMODELexp  ! m s-1
  totalOutflowDeriv  = (1._dp/HRUarea)*(maximumFlowRate/maximumSoilWater)*TOPMODELexp*(subSurfaceStorage/maximumSoilWater)**(TOPMODELexp - 1._dp)  ! s-1
 else
  totalColumnOutflow = 0._dp
  totalOutflowDeriv  = 0._dp
 endif

 end subroutine soilBsFlow_muster





 ! ************************************************************************************************
 ! private subroutine: disaggregate total inflow and total outflow to the baseflow sink term
 ! ************************************************************************************************
 subroutine disaggFlow_muster(&
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
                              ! output: error control
                              err,message)                  ! intent(out): error control
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
 ! output: error control
 integer(i4b),intent(out)     :: err                        ! error code
 character(*),intent(out)     :: message                    ! error message
 ! ----------------------------------------------------------------------------------------------------------
 ! local variables
 real(dp)                     :: sumDepthAvgCond            ! sum of depth-weighted hydraulic conductivity (m2 s-1)
 real(dp),dimension(nSoil)    :: fracTotalOutflow           ! fraction of outflow apportioned to each layer (-)
 real(dp),dimension(nSoil)    :: mLayerColumnInflowAdjusted ! adjusted column inflow to ensure no inflow into ice layers (m3 s-1)
 real(dp)                     :: totalInflowUnallocated     ! unallocated inflow (m3 s-1)
 real(dp)                     :: volTotalWater              ! volumetric fraction of total water (liquid + ice)
 real(dp)                     :: qMax                       ! max inflow rate (m s-1)
 real(dp)                     :: sink                       ! net source/sink (m s-1)
 integer(i4b)                 :: iLayer                     ! index of model layer
 ! ----------------------------------------------------------------------------------------------------------
 ! initialize error control
 err=0; message='disaggFlow_muster/'

 ! initialize adjusted inflow to each layer (m3/s)
 mLayerColumnInflowAdjusted(1:nSoil) = 0._dp

 ! simple case with no flow
 if(totalColumnOutflow < tiny(theta_sat))then
  mLayerColumnOutflow(1:nSoil) = HRUarea*totalColumnOutflow/real(nSoil, kind(dp))
  mLayerBaseflow(1:nSoil)      = (mLayerColumnOutflow(1:nSoil) - mLayerColumnInflow(1:nSoil))/HRUarea
  return
 endif

 ! possible to start iterations with zero storage, and gain storage through inflow -- reset ixSaturation
 if(ixSaturation > nSoil) ixSaturation=nSoil

 ! compute the depth-weighted hydraulic conductivity (m2 s-1)
 sumDepthAvgCond = sum(mLayerHydCond(ixSaturation:nSoil)*mLayerDepth(ixSaturation:nSoil))
 if(sumDepthAvgCond < tiny(theta_sat))then
  write(*,'(a,1x,10(e20.10,1x))') 'mLayerHydCond(ixSaturation:nSoil)       = ', mLayerHydCond(ixSaturation:nSoil)
  write(*,'(a,1x,10(e20.10,1x))') 'mLayerDepth(ixSaturation:nSoil)         = ', mLayerDepth(ixSaturation:nSoil)
  message=trim(message)//'zero depth-weighted conductivity when baseflow occurs'
  err=20; return
 endif

 ! compute the fraction of outflow apportioned to each layer (-)
 fracTotalOutflow(ixSaturation:nSoil) = (mLayerHydCond(ixSaturation:nSoil)*mLayerDepth(ixSaturation:nSoil)) / sumDepthAvgCond
 if(ixSaturation > 1) fracTotalOutflow(1:ixSaturation-1) = 0._dp
 if(abs(1._dp - sum(fracTotalOutflow)) > verySmall)then
  message=trim(message)//'fraction of baseflow does not sum to 1'
  err=20; return
 endif

 ! compute the outflow from each soil layer (m3 s-1)
 mLayerColumnOutflow(1:nSoil) = fracTotalOutflow(1:nSoil)*totalColumnOutflow*HRUarea

 ! initialize total unallocated inflow (m3 s-1)
 totalInflowUnallocated = totalColumnInflow

 ! adjust the inflow into each layer -- start at the bottom and work upwards
 do iLayer=nSoil,1,-1
  ! (define volumetric fraction of total water)
  volTotalWater = mLayerVolFracLiq(iLayer) + mLayerVolFracIce(iLayer)
  ! (define maximum flow into the layer)
  qMax = mLayerDepth(iLayer)*(theta_sat - max(volTotalWater,fieldCapacity))/dt   ! max inflow rate (m s-1)
  sink = (totalInflowUnallocated - mLayerColumnOutflow(iLayer))/HRUarea          ! maximum source/sink (m s-1)
  ! (add as much water as possible to the current layer)
  mLayerColumnInflowAdjusted(iLayer) = min(qMax,sink)*HRUarea                    ! inflow into the current layer (m3 s-1)
  ! (save remaining water for the next layer)
  totalInflowUnallocated = totalInflowUnallocated - mLayerColumnInflowAdjusted(iLayer)  ! m3 s-1
  if(totalInflowUnallocated < tiny(dt)) exit ! exit do loop as initialized at zero earlier
 end do  ! (looping through soil layers)

 ! compute the exfiltration (m s-1)
 exfiltration = totalInflowUnallocated/HRUarea

 ! compute the net baseflow from each soil layer (m s-1)
 mLayerBaseflow(1:nSoil) = (mLayerColumnOutflow(1:nSoil) - mLayerColumnInflowAdjusted(1:nSoil))/HRUarea

 ! add exfiltration to the baseflow flux at the top layer
 mLayerBaseflow(1)      = mLayerBaseflow(1) + exfiltration
 mLayerColumnOutflow(1) = mLayerColumnOutflow(1) + exfiltration*HRUarea

 end subroutine disaggFlow_muster


end module groundwatr_module

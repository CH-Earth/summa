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

module run_oneGRU_module

! numerical recipes data types
USE nrtype

! access missing values
USE globalData,only:realMissing            ! real missing value
USE globalData,only:integerMissing         ! integer missing value

! access vegetation data
USE globalData,only:greenVegFrac_monthly   ! fraction of green vegetation in each month (0-1)

! provide access to Noah-MP constants
USE module_sf_noahmplsm,only:isWater       ! parameter for water land cover type

! define data types
USE data_types,only:&
                    ! mapping
                    gru2hru_map,         & ! mapping between the GRUs and HRUs
                    ! no spatial dimension
                    var_i,               & ! x%var(:)            (i4b)
                    var_d,               & ! x%var(:)            (dp)
                    var_ilength,         & ! x%var(:)%dat        (i4b)
                    var_dlength,         & ! x%var(:)%dat        (dp)
                    ! no variable dimension
                    hru_i,               & ! x%hru(:)            (i4b)
                    hru_d,               & ! x%hru(:)            (dp)
                    ! gru dimension
                    gru_int,             & ! x%gru(:)%var(:)     (i4b)
                    gru_double,          & ! x%gru(:)%var(:)     (dp)
                    gru_intVec,          & ! x%gru(:)%var(:)%dat (i4b)
                    gru_doubleVec,       & ! x%gru(:)%var(:)%dat (dp)
                    ! gru+hru dimension
                    gru_hru_int,         & ! x%gru(:)%hru(:)%var(:)     (i4b)
                    gru_hru_double,      & ! x%gru(:)%hru(:)%var(:)     (dp)
                    gru_hru_intVec,      & ! x%gru(:)%hru(:)%var(:)%dat (i4b)
                    gru_hru_doubleVec      ! x%gru(:)%hru(:)%var(:)%dat (dp)

! provide access to the named variables that describe elements of parameter structures
USE var_lookup,only:iLookTYPE              ! look-up values for classification of veg, soils etc.
USE var_lookup,only:iLookATTR              ! look-up values for local attributes
USE var_lookup,only:iLookFORCE             ! look-up values for forcing data structures
USE var_lookup,only:iLookPARAM             ! look-up values for local column model parameters
USE var_lookup,only:iLookBPAR              ! look-up values for basin-average model parameters

! provide access to the named variables that describe elements of variable structures
USE var_lookup,only:iLookINDEX             ! look-up values for local column index variables
USE var_lookup,only:iLookPROG              ! look-up values for local column model prognostic (state) variables
USE var_lookup,only:iLookDIAG              ! look-up values for local column model diagnostic variables
USE var_lookup,only:iLookFLUX              ! look-up values for local column model fluxes
USE var_lookup,only:iLookBVAR              ! look-up values for basin-average model variables

! provide access to model decisions
USE globalData,only:model_decisions        ! model decision structure
USE var_lookup,only:iLookDECISIONS         ! look-up values for model decisions

! provide access to the named variables that describe model decisions
USE mDecisions_module,only:&               ! look-up values for LAI decisions
 monthlyTable,& ! LAI/SAI taken directly from a monthly table for different vegetation classes
 specified      ! LAI/SAI computed from green vegetation fraction and winterSAI and summerLAI parameters
USE mDecisions_module,only:&               ! look-up values for the choice of method for the spatial representation of groundwater
 localColumn, & ! separate groundwater representation in each local soil column
 singleBasin    ! single groundwater store over the entire basin

! -----------------------------------------------------------------------------------------------------------------------------------
! -----------------------------------------------------------------------------------------------------------------------------------
! ----- global variables that are modified ------------------------------------------------------------------------------------------
! -----------------------------------------------------------------------------------------------------------------------------------
! -----------------------------------------------------------------------------------------------------------------------------------

! Noah-MP parameters
USE NOAHMP_VEG_PARAMETERS,only:SAIM,LAIM   ! 2-d tables for stem area index and leaf area index (vegType,month)
USE NOAHMP_VEG_PARAMETERS,only:HVT,HVB     ! height at the top and bottom of vegetation (vegType)
USE noahmp_globals,only:RSMIN              ! minimum stomatal resistance (vegType)

! urban vegetation category (could be local)
USE globalData,only:urbanVegCategory       ! vegetation category for urban areas

! -----------------------------------------------------------------------------------------------------------------------------------
! -----------------------------------------------------------------------------------------------------------------------------------
! -----------------------------------------------------------------------------------------------------------------------------------
implicit none
private
public::run_oneGRU

contains

 ! ************************************************************************************************
 ! public subroutine run_oneGRU: simulation for a single GRU
 ! ************************************************************************************************

 ! simulation for a single GRU
 subroutine run_oneGRU(&
                       ! model control
                       iGRU,               & ! intent(in):    GRU index
                       dt_init,            & ! intent(inout): used to initialize the length of the sub-step for each HRU
                       gru_struc,          & ! intent(inout): mapping between the GRUs and the HRUs
                       ixComputeVegFlux,   & ! intent(inout): flag to indicate if we are computing fluxes over vegetation (false=no, true=yes)
                       ! data structures (input)
                       timeStruct,         & ! intent(in):    model time data
                       typeStruct,         & ! intent(in):    local classification of soil veg etc. for each HRU
                       attrStruct,         & ! intent(in):    local attributes for each HRU
                       mparStruct,         & ! intent(in):    local model parameters
                       ! data structures (input-output)
                       indxStruct,         & ! intent(inout): model indices
                       forcStruct,         & ! intent(inout): model forcing data
                       progStruct,         & ! intent(inout): prognostic variables for a local HRU
                       diagStruct,         & ! intent(inout): diagnostic variables for a local HRU
                       fluxStruct,         & ! intent(inout): model fluxes for a local HRU
                       bvarStruct,         & ! intent(inout): basin-average variables
                       ! error control
                       err,message)         ! intent(out):   error control

 ! ----- define downstream subroutines -----------------------------------------------------------------------------------

 USE module_sf_noahmplsm,only:redprm                               ! module to assign more Noah-MP parameters
 USE derivforce_module,only:derivforce                             ! module to compute derived forcing data
 USE coupled_em_module,only:coupled_em                             ! module to run the coupled energy and mass model
 USE qTimeDelay_module,only:qOverland                              ! module to route water through an "unresolved" river network
 
 ! ----- define dummy variables ------------------------------------------------------------------------------------------
 
 implicit none

 ! model control
 integer(i4b)            , intent(in)    :: iGRU                   ! GRU index
 type(hru_d)             , intent(inout) :: dt_init(:)             ! used to initialize the length of the sub-step for each HRU
 type(gru2hru_map)       , intent(inout) :: gru_struc(:)           ! mapping between the GRUs and the HRUs
 type(hru_i)             , intent(inout) :: ixComputeVegFlux(:)    ! flag to indicate if we are computing fluxes over vegetation (false=no, true=yes)

 ! data structures (input)
 type(var_i)             , intent(in)    :: timeStruct             ! x%var(:)                   -- model time data
 type(gru_hru_int)       , intent(in)    :: typeStruct             ! x%gru(:)%hru(:)%var(:)     -- local classification of soil veg etc. for each HRU
 type(gru_hru_double)    , intent(in)    :: attrStruct             ! x%gru(:)%hru(:)%var(:)     -- local attributes for each HRU
 type(gru_hru_doubleVec) , intent(in)    :: mparStruct             ! x%gru(:)%hru(:)%var(:)%dat -- local (HRU) model parameters
 ! data structures (input-output)
 type(gru_hru_intVec)    , intent(inout) :: indxStruct             ! x%gru(:)%hru(:)%var(:)%dat -- model indices
 type(gru_hru_double)    , intent(inout) :: forcStruct             ! x%gru(:)%hru(:)%var(:)     -- model forcing data
 type(gru_hru_doubleVec) , intent(inout) :: progStruct             ! x%gru(:)%hru(:)%var(:)%dat -- model prognostic (state) variables
 type(gru_hru_doubleVec) , intent(inout) :: diagStruct             ! x%gru(:)%hru(:)%var(:)%dat -- model diagnostic variables
 type(gru_hru_doubleVec) , intent(inout) :: fluxStruct             ! x%gru(:)%hru(:)%var(:)%dat -- model fluxes
 type(gru_doubleVec)     , intent(inout) :: bvarStruct             ! x%gru(:)%var(:)%dat        -- basin-average variables
 ! error control
 integer(i4b)            , intent(out)   :: err                    ! error code
 character(*)            , intent(out)   :: message                ! error message

 ! ----- define local variables ------------------------------------------------------------------------------------------

 ! general local variables
 character(len=256)                      :: cmessage               ! error message
 integer(i4b)                            :: iHRU                   ! HRU index
 ! HRU initialization
 real(dp)                                :: fracHRU                ! fractional area of a given HRU (-)
 logical(lgt)            , parameter     :: overwriteRSMIN=.false. ! flag to overwrite RSMIN
 integer(i4b)            , parameter     :: maxSoilLayers=10       ! Maximum Number of Soil Layers
 real(dp)                , allocatable   :: zSoilReverseSign(:)    ! height at bottom of each soil layer, negative downwards (m)
 ! run the model
 integer(i4b),parameter                  :: no=0                   ! .false.
 integer(i4b),parameter                  :: yes=1                  ! .true.
 logical(lgt)                            :: computeVegFluxFlag     ! flag to indicate if we are computing fluxes over vegetation (.false. means veg is buried with snow)
 ! compute fluxes across HRUs
 integer(i4b)                            :: jHRU,kHRU              ! index of the hydrologic response unit 

 ! initialize error control
 err=0; message='run_oneGRU/'
 
 ! ----- basin initialization --------------------------------------------------------------------------------------------

 ! initialize runoff variables
 bvarStruct%gru(iGRU)%var(iLookBVAR%basin__SurfaceRunoff)%dat(1)    = 0._dp  ! surface runoff (m s-1)
 bvarStruct%gru(iGRU)%var(iLookBVAR%basin__ColumnOutflow)%dat(1)    = 0._dp  ! outflow from all "outlet" HRUs (those with no downstream HRU)

 ! initialize baseflow variables
 bvarStruct%gru(iGRU)%var(iLookBVAR%basin__AquiferRecharge)%dat(1)  = 0._dp ! recharge to the aquifer (m s-1)
 bvarStruct%gru(iGRU)%var(iLookBVAR%basin__AquiferBaseflow)%dat(1)  = 0._dp ! baseflow from the aquifer (m s-1)
 bvarStruct%gru(iGRU)%var(iLookBVAR%basin__AquiferTranspire)%dat(1) = 0._dp ! transpiration loss from the aquifer (m s-1)

 ! initialize total inflow for each layer in a soil column
 do iHRU=1,gru_struc(iGRU)%hruCount
  fluxStruct%gru(iGRU)%hru(iHRU)%var(iLookFLUX%mLayerColumnInflow)%dat(:) = 0._dp
 end do

 ! ***********************************************************************************************************************
 ! ********** LOOP THROUGH HRUS ******************************************************************************************
 ! ***********************************************************************************************************************

 ! loop through HRUs
 do iHRU=1,gru_struc(iGRU)%hruCount

  ! ----- hru initialization ---------------------------------------------------------------------------------------------

  ! make associations to variables in data structures
  associate(&
    nSnow   =>indxStruct%gru(iGRU)%hru(iHRU)%var(iLookINDEX%nSnow)%dat(1)    , &  ! number of snow layers
    nSoil   =>indxStruct%gru(iGRU)%hru(iHRU)%var(iLookINDEX%nSoil)%dat(1)    , &  ! number of soil layers
    nLayers =>indxStruct%gru(iGRU)%hru(iHRU)%var(iLookINDEX%nLayers)%dat(1)    &  ! total number of layers
  ) ! associations

  ! cycle water pixel
  if (typeStruct%gru(iGRU)%hru(iHRU)%var(iLookTYPE%vegTypeIndex) == isWater) cycle

  ! identify the area covered by the current HRU
  fracHRU =  attrStruct%gru(iGRU)%hru(iHRU)%var(iLookATTR%HRUarea) / bvarStruct%gru(iGRU)%var(iLookBVAR%basin__totalArea)%dat(1)

  ! assign model layers
  ! NOTE: layer structure is different for each HRU
  gru_struc(iGRU)%hruInfo(iHRU)%nSnow = nSnow
  gru_struc(iGRU)%hruInfo(iHRU)%nSoil = nSoil

  ! get height at bottom of each soil layer, negative downwards (used in Noah MP)
  allocate(zSoilReverseSign(gru_struc(iGRU)%hruInfo(iHRU)%nSoil),stat=err)
  if(err/=0)then
   message=trim(message)//'problem allocating space for zSoilReverseSign'
   err=20; return
  endif
  zSoilReverseSign(:) = -progStruct%gru(iGRU)%hru(iHRU)%var(iLookPROG%iLayerHeight)%dat(gru_struc(iGRU)%hruInfo(iHRU)%nSnow+1:nLayers)

  ! get NOAH-MP parameters
  ! Passing a maxSoilLayer in order to pass the check for NROOT, that is done to avoid making any changes to Noah-MP code.
  !  --> NROOT from Noah-MP veg tables (as read here) is not used in SUMMA
  call REDPRM(typeStruct%gru(iGRU)%hru(iHRU)%var(iLookTYPE%vegTypeIndex),      & ! vegetation type index
              typeStruct%gru(iGRU)%hru(iHRU)%var(iLookTYPE%soilTypeIndex),     & ! soil type
              typeStruct%gru(iGRU)%hru(iHRU)%var(iLookTYPE%slopeTypeIndex),    & ! slope type index
              zSoilReverseSign,                                                & ! * not used: height at bottom of each layer [NOTE: negative] (m)
              maxSoilLayers,                                                   & ! number of soil layers
              urbanVegCategory)                                                  ! vegetation category for urban areas

  ! deallocate height at bottom of each soil layer(used in Noah MP)
  deallocate(zSoilReverseSign,stat=err)
  if(err/=0)then
   message=trim(message)//'problem deallocating space for zSoilReverseSign'
   err=20; return
  endif

  ! overwrite the minimum resistance
  if(overwriteRSMIN) RSMIN = mparStruct%gru(iGRU)%hru(iHRU)%var(iLookPARAM%minStomatalResistance)%dat(1)

  ! overwrite the vegetation height
  HVT(typeStruct%gru(iGRU)%hru(iHRU)%var(iLookTYPE%vegTypeIndex)) = mparStruct%gru(iGRU)%hru(iHRU)%var(iLookPARAM%heightCanopyTop)%dat(1)
  HVB(typeStruct%gru(iGRU)%hru(iHRU)%var(iLookTYPE%vegTypeIndex)) = mparStruct%gru(iGRU)%hru(iHRU)%var(iLookPARAM%heightCanopyBottom)%dat(1)

  ! overwrite the tables for LAI and SAI
  if(model_decisions(iLookDECISIONS%LAI_method)%iDecision == specified)then
   SAIM(typeStruct%gru(iGRU)%hru(iHRU)%var(iLookTYPE%vegTypeIndex),:) = mparStruct%gru(iGRU)%hru(iHRU)%var(iLookPARAM%winterSAI)%dat(1)
   LAIM(typeStruct%gru(iGRU)%hru(iHRU)%var(iLookTYPE%vegTypeIndex),:) = mparStruct%gru(iGRU)%hru(iHRU)%var(iLookPARAM%summerLAI)%dat(1)*greenVegFrac_monthly
  end if

  ! compute derived forcing variables
  call derivforce(timeStruct%var,                    & ! vector of time information
                  forcStruct%gru(iGRU)%hru(iHRU)%var,& ! vector of model forcing data
                  attrStruct%gru(iGRU)%hru(iHRU)%var,& ! vector of model attributes
                  mparStruct%gru(iGRU)%hru(iHRU),    & ! vector of model parameters
                  progStruct%gru(iGRU)%hru(iHRU),    & ! data structure of model prognostic variables
                  diagStruct%gru(iGRU)%hru(iHRU),    & ! data structure of model diagnostic variables
                  fluxStruct%gru(iGRU)%hru(iHRU),    & ! data structure of model fluxes
                  err,cmessage)                        ! error control
  if(err/=0)then; err=20; message=trim(message)//trim(cmessage); return; endif

  ! ----- run the model --------------------------------------------------------------------------------------------------

  ! set the flag to compute the vegetation flux
  computeVegFluxFlag = (ixComputeVegFlux(iGRU)%hru(iHRU) == yes)

  ! initialize the number of flux calls
  diagStruct%gru(iGRU)%hru(iHRU)%var(iLookDIAG%numFluxCalls)%dat(1) = 0._dp

  ! run the model for a single parameter set and time step
  call coupled_em(&
                  ! model control
                  gru_struc(iGRU)%hruInfo(iHRU)%hru_id,    & ! intent(in):    hruId
                  dt_init(iGRU)%hru(iHRU),                 & ! intent(inout): initial time step
                  computeVegFluxFlag,                      & ! intent(inout): flag to indicate if we are computing fluxes over vegetation
                  ! data structures (input)
                  typeStruct%gru(iGRU)%hru(iHRU),          & ! intent(in):    local classification of soil veg etc. for each HRU
                  attrStruct%gru(iGRU)%hru(iHRU),          & ! intent(in):    local attributes for each HRU
                  forcStruct%gru(iGRU)%hru(iHRU),          & ! intent(in):    model forcing data
                  mparStruct%gru(iGRU)%hru(iHRU),          & ! intent(in):    model parameters
                  bvarStruct%gru(iGRU),                    & ! intent(in):    basin-average model variables
                  ! data structures (input-output)
                  indxStruct%gru(iGRU)%hru(iHRU),          & ! intent(inout): model indices
                  progStruct%gru(iGRU)%hru(iHRU),          & ! intent(inout): model prognostic variables for a local HRU
                  diagStruct%gru(iGRU)%hru(iHRU),          & ! intent(inout): model diagnostic variables for a local HRU
                  fluxStruct%gru(iGRU)%hru(iHRU),          & ! intent(inout): model fluxes for a local HRU
                  ! error control
                  err,message)            ! intent(out): error control
  if(err/=0)then; err=20; message=trim(message)//trim(cmessage); return; endif

  ! update layer numbers that could be changed in coupled_em()
  gru_struc(iGRU)%hruInfo(iHRU)%nSnow = indxStruct%gru(iGRU)%hru(iHRU)%var(iLookINDEX%nSnow)%dat(1)
  gru_struc(iGRU)%hruInfo(iHRU)%nSoil = indxStruct%gru(iGRU)%hru(iHRU)%var(iLookINDEX%nSoil)%dat(1)

  ! save the flag for computing the vegetation fluxes
  if(computeVegFluxFlag)      ixComputeVegFlux(iGRU)%hru(iHRU) = yes
  if(.not.computeVegFluxFlag) ixComputeVegFlux(iGRU)%hru(iHRU) = no

  ! ----- compute fluxes across HRUs --------------------------------------------------------------------------------------------------

  kHRU = 0
  ! identify the downslope HRU
  dsHRU: do jHRU=1,gru_struc(iGRU)%hruCount
   if(typeStruct%gru(iGRU)%hru(iHRU)%var(iLookTYPE%downHRUindex) == typeStruct%gru(iGRU)%hru(jHRU)%var(iLookTYPE%hruId))then
    if(kHRU==0)then  ! check there is a unique match
     kHRU=jHRU
     exit dsHRU
    end if  ! (check there is a unique match)
   end if  ! (if identified a downslope HRU)
  end do dsHRU

  ! add inflow to the downslope HRU
  if(kHRU > 0)then  ! if there is a downslope HRU
   fluxStruct%gru(iGRU)%hru(kHRU)%var(iLookFLUX%mLayerColumnInflow)%dat(:) = fluxStruct%gru(iGRU)%hru(kHRU)%var(iLookFLUX%mLayerColumnInflow)%dat(:)  + fluxStruct%gru(iGRU)%hru(iHRU)%var(iLookFLUX%mLayerColumnOutflow)%dat(:)

  ! increment basin column outflow (m3 s-1)
  else
   bvarStruct%gru(iGRU)%var(iLookBVAR%basin__ColumnOutflow)%dat(1)   = bvarStruct%gru(iGRU)%var(iLookBVAR%basin__ColumnOutflow)%dat(1) + sum(fluxStruct%gru(iGRU)%hru(iHRU)%var(iLookFLUX%mLayerColumnOutflow)%dat(:))
  end if

  ! increment basin surface runoff (m s-1)
  bvarStruct%gru(iGRU)%var(iLookBVAR%basin__SurfaceRunoff)%dat(1)    = bvarStruct%gru(iGRU)%var(iLookBVAR%basin__SurfaceRunoff)%dat(1)     + fluxStruct%gru(iGRU)%hru(iHRU)%var(iLookFLUX%scalarSurfaceRunoff)%dat(1)    * fracHRU

  ! increment basin-average baseflow input variables (m s-1)
  bvarStruct%gru(iGRU)%var(iLookBVAR%basin__AquiferRecharge)%dat(1)  = bvarStruct%gru(iGRU)%var(iLookBVAR%basin__AquiferRecharge)%dat(1)   + fluxStruct%gru(iGRU)%hru(iHRU)%var(iLookFLUX%scalarSoilDrainage)%dat(1)     * fracHRU
  bvarStruct%gru(iGRU)%var(iLookBVAR%basin__AquiferTranspire)%dat(1) = bvarStruct%gru(iGRU)%var(iLookBVAR%basin__AquiferTranspire)%dat(1)  + fluxStruct%gru(iGRU)%hru(iHRU)%var(iLookFLUX%scalarAquiferTranspire)%dat(1) * fracHRU

  ! increment aquifer baseflow -- ONLY if baseflow is computed individually for each HRU
  ! NOTE: groundwater computed later for singleBasin
  if(model_decisions(iLookDECISIONS%spatial_gw)%iDecision == localColumn)then
   bvarStruct%gru(iGRU)%var(iLookBVAR%basin__AquiferBaseflow)%dat(1)  =  bvarStruct%gru(iGRU)%var(iLookBVAR%basin__AquiferBaseflow)%dat(1)  &
           +  fluxStruct%gru(iGRU)%hru(iHRU)%var(iLookFLUX%scalarAquiferBaseflow)%dat(1) * fracHRU  &
           +  fluxStruct%gru(iGRU)%hru(iHRU)%var(iLookFLUX%scalarSoilDrainage)%dat(1)    * fracHRU
  end if

  ! end HRU associations
  end associate

 end do  ! (looping through HRUs)

 ! ***********************************************************************************************************************
 ! ********** END LOOP THROUGH HRUS **************************************************************************************
 ! ***********************************************************************************************************************

 ! compute water balance for the basin aquifer
 if(model_decisions(iLookDECISIONS%spatial_gw)%iDecision == singleBasin)then
  message=trim(message)//'multi_driver/bigBucket groundwater code not transferred from old code base yet'
  err=20; return
 end if

 ! perform the routing
 associate(totalArea => bvarStruct%gru(iGRU)%var(iLookBVAR%basin__totalArea)%dat(1) )
 call qOverland(&
                ! input
                model_decisions(iLookDECISIONS%subRouting)%iDecision,                      &  ! intent(in): index for routing method
                bvarStruct%gru(iGRU)%var(iLookBVAR%basin__SurfaceRunoff)%dat(1),           &  ! intent(in): surface runoff (m s-1)
                bvarStruct%gru(iGRU)%var(iLookBVAR%basin__ColumnOutflow)%dat(1)/totalArea, &  ! intent(in): outflow from all "outlet" HRUs (those with no downstream HRU)
                bvarStruct%gru(iGRU)%var(iLookBVAR%basin__AquiferBaseflow)%dat(1),         &  ! intent(in): baseflow from the aquifer (m s-1)
                bvarStruct%gru(iGRU)%var(iLookBVAR%routingFractionFuture)%dat,             &  ! intent(in): fraction of runoff in future time steps (m s-1)
                bvarStruct%gru(iGRU)%var(iLookBVAR%routingRunoffFuture)%dat,               &  ! intent(in): runoff in future time steps (m s-1)
                ! output
                bvarStruct%gru(iGRU)%var(iLookBVAR%averageInstantRunoff)%dat(1),           &  ! intent(out): instantaneous runoff (m s-1)
                bvarStruct%gru(iGRU)%var(iLookBVAR%averageRoutedRunoff)%dat(1),            &  ! intent(out): routed runoff (m s-1)
                err,message)                                                                  ! intent(out): error control
 if(err/=0)then; err=20; message=trim(message)//trim(cmessage); return; endif
 end associate

 end subroutine run_oneGRU

end module run_oneGRU_module

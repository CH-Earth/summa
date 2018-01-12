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

! access integers to define "yes" and "no"
USE globalData,only:yes,no                 ! .true. and .false.

! access the mapping betweeen GRUs and HRUs
USE globalData,only:gru_struc              ! gru-hru mapping structures

! define data types
USE data_types,only:&
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
USE var_lookup,only:iLookINDEX             ! look-up values for local column index variables
USE var_lookup,only:iLookFLUX              ! look-up values for local column model fluxes
USE var_lookup,only:iLookBVAR              ! look-up values for basin-average model variables

! provide access to model decisions
USE globalData,only:model_decisions        ! model decision structure
USE var_lookup,only:iLookDECISIONS         ! look-up values for model decisions

! provide access to the named variables that describe model decisions
USE mDecisions_module,only:&               ! look-up values for the choice of method for the spatial representation of groundwater
 localColumn, & ! separate groundwater representation in each local soil column
 singleBasin    ! single groundwater store over the entire basin

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

 USE run_oneHRU_module,only:run_oneHRU                             ! module to run for one HRU
 USE qTimeDelay_module,only:qOverland                              ! module to route water through an "unresolved" river network
 implicit none

 ! ----- define dummy variables ------------------------------------------------------------------------------------------
 
 ! model control
 integer(i4b)            , intent(in)    :: iGRU                   ! GRU index
 type(hru_d)             , intent(inout) :: dt_init(:)             ! used to initialize the length of the sub-step for each HRU
 type(hru_i)             , intent(inout) :: ixComputeVegFlux(:)    ! flag to indicate if we are computing fluxes over vegetation (false=no, true=yes)
 ! data structures (input)
 type(var_i)             , intent(in)    :: timeStruct             ! x%var(:)                   -- model time data
 type(gru_hru_int)       , intent(in)    :: typeStruct             ! x%gru(:)%hru(:)%var(:)     -- local classification of soil veg etc. for each HRU
 type(gru_hru_double)    , intent(in)    :: attrStruct             ! x%gru(:)%hru(:)%var(:)     -- local attributes for each HRU
 ! data structures (input-output)
 type(gru_hru_doubleVec) , intent(inout) :: mparStruct             ! x%gru(:)%hru(:)%var(:)%dat -- local (HRU) model parameters
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
 integer(i4b)                            :: jHRU,kHRU              ! index of the hydrologic response unit 
 real(dp)                                :: fracHRU                ! fractional area of a given HRU (-)
 logical(lgt)                            :: computeVegFluxFlag     ! flag to indicate if we are computing fluxes over vegetation (.false. means veg is buried with snow)

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
 ! ********** RUN FOR ONE HRU ********************************************************************************************
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

  ! set the flag to compute the vegetation flux
  computeVegFluxFlag = (ixComputeVegFlux(iGRU)%hru(iHRU) == yes)

  ! ----- run the model --------------------------------------------------------------------------------------------------

  ! simulation for a single HRU
  call run_oneHRU(&
                  ! model control
                  gru_struc(iGRU)%hruInfo(iHRU)%hru_id,    & ! intent(in):    hruId
                  dt_init(iGRU)%hru(iHRU),                 & ! intent(inout): initial time step
                  computeVegFluxFlag,                      & ! intent(inout): flag to indicate if we are computing fluxes over vegetation (false=no, true=yes)
                  nSnow,nSoil,nLayers,                     & ! intent(inout): number of snow and soil layers
                  ! data structures (input)
                  timeStruct%var,                          & ! intent(in):    model time data
                  typeStruct%gru(iGRU)%hru(iHRU),          & ! intent(in):    local classification of soil veg etc. for each HRU
                  attrStruct%gru(iGRU)%hru(iHRU),          & ! intent(in):    local attributes for each HRU
                  bvarStruct%gru(iGRU),                    & ! intent(in):    basin-average model variables
                  ! data structures (input-output)
                  mparStruct%gru(iGRU)%hru(iHRU),          & ! intent(inout): model parameters
                  indxStruct%gru(iGRU)%hru(iHRU),          & ! intent(inout): model indices
                  forcStruct%gru(iGRU)%hru(iHRU),          & ! intent(inout): model forcing data
                  progStruct%gru(iGRU)%hru(iHRU),          & ! intent(inout): model prognostic variables for a local HRU
                  diagStruct%gru(iGRU)%hru(iHRU),          & ! intent(inout): model diagnostic variables for a local HRU
                  fluxStruct%gru(iGRU)%hru(iHRU),          & ! intent(inout): model fluxes for a local HRU
                  ! error control
                  err,message)                               ! intent(out):   error control
  if(err/=0)then; err=20; message=trim(message)//trim(cmessage); return; endif

  ! update layer numbers that could be changed in run_oneHRU -- needed for model output
  gru_struc(iGRU)%hruInfo(iHRU)%nSnow = nSnow
  gru_struc(iGRU)%hruInfo(iHRU)%nSoil = nSoil

  ! save the flag for computing the vegetation fluxes
  if(computeVegFluxFlag)      ixComputeVegFlux(iGRU)%hru(iHRU) = yes
  if(.not.computeVegFluxFlag) ixComputeVegFlux(iGRU)%hru(iHRU) = no

  ! ----- compute fluxes across HRUs --------------------------------------------------------------------------------------------------

  ! identify the area covered by the current HRU
  fracHRU =  attrStruct%gru(iGRU)%hru(iHRU)%var(iLookATTR%HRUarea) / bvarStruct%gru(iGRU)%var(iLookBVAR%basin__totalArea)%dat(1)

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

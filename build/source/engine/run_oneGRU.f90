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

! access missing values
USE globalData,only:realMissing       ! real missing value
USE globalData,only:integerMissing    ! integer missing value

! access the mapping betweeen GRUs and HRUs
USE globalData,only:gru_struc         ! gru-hru mapping structures

! access the minimum and maximum HRUs in the file
USE globalData,only:ixHRUfile_min,ixHRUfile_max

! define data types
USE data_types,only:gru_hru_double    ! x%gru(:)%hru(:)%var(:)     (dp)

implicit none
private
public::run_oneGRU

contains

 ! ************************************************************************************************
 ! public subroutine run_oneGRU: simulation for a single GRU
 ! ************************************************************************************************

 subroutine run_oneGRU(&
                       ! model control
                       iGRU,               & ! intent(in):    GRU index
                       ! data structures (input)
                       typeStruct,         & ! intent(in):    local classification of soil veg etc. for each HRU
                       attrStruct,         & ! intent(in):    local attributes for each HRU
                       forcStruct,         & ! intent(in):    model forcing data
                       mparStruct,         & ! intent(in):    model parameters
                       bvarStruct,         & ! intent(in):    basin-average variables
                       ! data structures (input-output)
                       indxStruct,         & ! intent(inout): model indices
                       progStruct,         & ! intent(inout): prognostic variables for a local HRU
                       diagStruct,         & ! intent(inout): diagnostic variables for a local HRU
                       fluxStruct,         & ! intent(inout): model fluxes for a local HRU
                       ! error control
                       err,message)         ! intent(out):   error control


! define the primary data structures (scalars)
type(var_i)                      :: timeStruct                 ! x%var(:)                   -- model time data
type(gru_hru_double)             :: forcStruct                 ! x%gru(:)%hru(:)%var(:)     -- model forcing data
type(gru_hru_double)             :: attrStruct                 ! x%gru(:)%hru(:)%var(:)     -- local attributes for each HRU
type(gru_hru_int)                :: typeStruct                 ! x%gru(:)%hru(:)%var(:)     -- local classification of soil veg etc. for each HRU
! define the primary data structures (variable length vectors)
type(gru_hru_intVec)             :: indxStruct                 ! x%gru(:)%hru(:)%var(:)%dat -- model indices
type(gru_hru_doubleVec)          :: mparStruct                 ! x%gru(:)%hru(:)%var(:)%dat -- model parameters
type(gru_hru_doubleVec)          :: progStruct                 ! x%gru(:)%hru(:)%var(:)%dat -- model prognostic (state) variables
type(gru_hru_doubleVec)          :: diagStruct                 ! x%gru(:)%hru(:)%var(:)%dat -- model diagnostic variables
type(gru_hru_doubleVec)          :: fluxStruct                 ! x%gru(:)%hru(:)%var(:)%dat -- model fluxes
! define the basin-average structures
type(gru_double)                 :: bparStruct                 ! x%gru(:)%var(:)            -- basin-average parameters
type(gru_doubleVec)              :: bvarStruct                 ! x%gru(:)%var(:)%dat        -- basin-average variables




 ! initialize error control
 err=0; message='run_oneGRU/'


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

  ! loop through HRUs
  do iHRU=1,gru_struc(iGRU)%hruCount

   ! identify the area covered by the current HRU
   fracHRU =  attrStruct%gru(iGRU)%hru(iHRU)%var(iLookATTR%HRUarea) / bvarStruct%gru(iGRU)%var(iLookBVAR%basin__totalArea)%dat(1)

   ! assign model layers
   ! NOTE: layer structure is different for each HRU
   gru_struc(iGRU)%hruInfo(iHRU)%nSnow = indxStruct%gru(iGRU)%hru(iHRU)%var(iLookINDEX%nSnow)%dat(1)
   gru_struc(iGRU)%hruInfo(iHRU)%nSoil = indxStruct%gru(iGRU)%hru(iHRU)%var(iLookINDEX%nSoil)%dat(1)
   nLayers                                 = indxStruct%gru(iGRU)%hru(iHRU)%var(iLookINDEX%nLayers)%dat(1)

   ! get height at bottom of each soil layer, negative downwards (used in Noah MP)
   allocate(zSoilReverseSign(gru_struc(iGRU)%hruInfo(iHRU)%nSoil),stat=err); call handle_err(err,'problem allocating space for zSoilReverseSign')
   zSoilReverseSign(:) = -progStruct%gru(iGRU)%hru(iHRU)%var(iLookPROG%iLayerHeight)%dat(gru_struc(iGRU)%hruInfo(iHRU)%nSnow+1:nLayers)

   ! get NOAH-MP parameters
   ! Passing a maxSoilLayer in order to pass the check for NROOT, that is done to avoid making any changes to Noah-MP code. NROOT from Noah-MP veg tables (as read here) is not used in SUMMA
   call REDPRM(typeStruct%gru(iGRU)%hru(iHRU)%var(iLookTYPE%vegTypeIndex),      & ! vegetation type index
               typeStruct%gru(iGRU)%hru(iHRU)%var(iLookTYPE%soilTypeIndex),     & ! soil type
               typeStruct%gru(iGRU)%hru(iHRU)%var(iLookTYPE%slopeTypeIndex),    & ! slope type index
               zSoilReverseSign,                                                & ! * not used: height at bottom of each layer [NOTE: negative] (m)
               maxSoilLayers,                                                   & ! number of soil layers
               urbanVegCategory)                                                  ! vegetation category for urban areas

   ! deallocate height at bottom of each soil layer(used in Noah MP)
   deallocate(zSoilReverseSign,stat=err); call handle_err(err,'problem deallocating space for zSoilReverseSign')

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

   ! cycle water pixel
   if (typeStruct%gru(iGRU)%hru(iHRU)%var(iLookTYPE%vegTypeIndex) == isWater) cycle

   ! compute derived forcing variables
   call derivforce(timeStruct%var,                    & ! vector of time information
                   forcStruct%gru(iGRU)%hru(iHRU)%var,& ! vector of model forcing data
                   attrStruct%gru(iGRU)%hru(iHRU)%var,& ! vector of model attributes
                   mparStruct%gru(iGRU)%hru(iHRU),    & ! vector of model parameters
                   progStruct%gru(iGRU)%hru(iHRU),    & ! data structure of model prognostic variables
                   diagStruct%gru(iGRU)%hru(iHRU),    & ! data structure of model diagnostic variables
                   fluxStruct%gru(iGRU)%hru(iHRU),    & ! data structure of model fluxes
                   err,message)                         ! error control
   call handle_err(err,message)

   ! ****************************************************************************
   ! *** run the model
   ! ****************************************************************************
   ! set the flag to compute the vegetation flux
   computeVegFluxFlag = (computeVegFlux(iGRU)%hru(iHRU) == yes)

   !print*, 'iHRU = ', iHRU

   ! initialize the number of flux calls
   diagStruct%gru(iGRU)%hru(iHRU)%var(iLookDIAG%numFluxCalls)%dat(1) = 0._dp

   ! run the model for a single parameter set and time step
   call coupled_em(&
                   ! model control
                   gru_struc(iGRU)%hruInfo(iHRU)%hru_id,    & ! intent(in):    hruId
                   dt_init(iGRU)%hru(iHRU),                 & ! intent(inout): initial time step
                   computeVegFluxFlag,                          & ! intent(inout): flag to indicate if we are computing fluxes over vegetation (.false. means veg is buried with snow)
                   ! data structures (input)
                   typeStruct%gru(iGRU)%hru(iHRU),          & ! intent(in):    local classification of soil veg etc. for each HRU
                   attrStruct%gru(iGRU)%hru(iHRU),          & ! intent(in):    local attributes for each HRU
                   forcStruct%gru(iGRU)%hru(iHRU),          & ! intent(in):    model forcing data
                   mparStruct%gru(iGRU)%hru(iHRU),          & ! intent(in):    model parameters
                   bvarStruct%gru(iGRU),                        & ! intent(in):    basin-average model variables
                   ! data structures (input-output)
                   indxStruct%gru(iGRU)%hru(iHRU),          & ! intent(inout): model indices
                   progStruct%gru(iGRU)%hru(iHRU),          & ! intent(inout): model prognostic variables for a local HRU
                   diagStruct%gru(iGRU)%hru(iHRU),          & ! intent(inout): model diagnostic variables for a local HRU
                   fluxStruct%gru(iGRU)%hru(iHRU),          & ! intent(inout): model fluxes for a local HRU
                   ! error control
                   err,message)            ! intent(out): error control
   call handle_err(err,message)

   ! update layer numbers that could be changed in coupled_em()
   gru_struc(iGRU)%hruInfo(iHRU)%nSnow = indxStruct%gru(iGRU)%hru(iHRU)%var(iLookINDEX%nSnow)%dat(1)
   gru_struc(iGRU)%hruInfo(iHRU)%nSoil = indxStruct%gru(iGRU)%hru(iHRU)%var(iLookINDEX%nSoil)%dat(1)

!   ! check feasibiility of certain states
!   call check_icond(nGRU,nHRU,                     & ! number of response units
!                    progStruct,                    & ! model prognostic (state) variables
!                    mparStruct,                    & ! model parameters
!                    indxStruct,                    & ! layer indexes
!                    err,message)                     ! error control
!   call handle_err(err,message)

   ! save the flag for computing the vegetation fluxes
   if(computeVegFluxFlag)      computeVegFlux(iGRU)%hru(iHRU) = yes
   if(.not.computeVegFluxFlag) computeVegFlux(iGRU)%hru(iHRU) = no

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

   ! print progress
   !print*, 'resetStats     = ', resetStats
   !print*, 'finalizeStats  = ', finalizeStats
   !print*, 'statCounter    = ', statCounter
   !print*, 'outputTimeStep = ', outputTimeStep 

   ! calculate output Statistics
   call calcStats(forcStat%gru(iGRU)%hru(iHRU)%var,forcStruct%gru(iGRU)%hru(iHRU)%var,statForc_meta,resetStats,finalizeStats,statCounter,err,message); call handle_err(err,message)
   call calcStats(progStat%gru(iGRU)%hru(iHRU)%var,progStruct%gru(iGRU)%hru(iHRU)%var,statProg_meta,resetStats,finalizeStats,statCounter,err,message); call handle_err(err,message)
   call calcStats(diagStat%gru(iGRU)%hru(iHRU)%var,diagStruct%gru(iGRU)%hru(iHRU)%var,statDiag_meta,resetStats,finalizeStats,statCounter,err,message); call handle_err(err,message)
   call calcStats(fluxStat%gru(iGRU)%hru(iHRU)%var,fluxStruct%gru(iGRU)%hru(iHRU)%var,statFlux_meta,resetStats,finalizeStats,statCounter,err,message); call handle_err(err,message)
   call calcStats(indxStat%gru(iGRU)%hru(iHRU)%var,indxStruct%gru(iGRU)%hru(iHRU)%var,statIndx_meta,resetStats,finalizeStats,statCounter,err,message); call handle_err(err,message)

  end do  ! (looping through HRUs)

  ! compute water balance for the basin aquifer
  if(model_decisions(iLookDECISIONS%spatial_gw)%iDecision == singleBasin)then
   call handle_err(20,'multi_driver/bigBucket groundwater code not transferred from old code base yet')
  end if

  ! perform the routing
  associate(totalArea => bvarStruct%gru(iGRU)%var(iLookBVAR%basin__totalArea)%dat(1) )
  call qOverland(&
                 ! input
                 model_decisions(iLookDECISIONS%subRouting)%iDecision,            &  ! intent(in): index for routing method
                 bvarStruct%gru(iGRU)%var(iLookBVAR%basin__SurfaceRunoff)%dat(1),           &  ! intent(in): surface runoff (m s-1)
                 bvarStruct%gru(iGRU)%var(iLookBVAR%basin__ColumnOutflow)%dat(1)/totalArea, &  ! intent(in): outflow from all "outlet" HRUs (those with no downstream HRU)
                 bvarStruct%gru(iGRU)%var(iLookBVAR%basin__AquiferBaseflow)%dat(1),         &  ! intent(in): baseflow from the aquifer (m s-1)
                 bvarStruct%gru(iGRU)%var(iLookBVAR%routingFractionFuture)%dat,             &  ! intent(in): fraction of runoff in future time steps (m s-1)
                 bvarStruct%gru(iGRU)%var(iLookBVAR%routingRunoffFuture)%dat,               &  ! intent(in): runoff in future time steps (m s-1)
                 ! output
                 bvarStruct%gru(iGRU)%var(iLookBVAR%averageInstantRunoff)%dat(1),           &  ! intent(out): instantaneous runoff (m s-1)
                 bvarStruct%gru(iGRU)%var(iLookBVAR%averageRoutedRunoff)%dat(1),            &  ! intent(out): routed runoff (m s-1)
                 err,message)                                                        ! intent(out): error control
  call handle_err(err,message)
  end associate

  ! calc basin stats
  call calcStats(bvarStat%gru(iGRU)%var(:),bvarStruct%gru(iGRU)%var(:),statBvar_meta,resetStats,finalizeStats,statCounter,err,message); call handle_err(err,message)

  ! write basin-average variables
  call writeBasin(iGRU,finalizeStats,outputTimeStep,bvar_meta,bvarStat%gru(iGRU)%var,bvarStruct%gru(iGRU)%var,bvarChild_map,err,message); call handle_err(err,message)

 end do  ! (looping through GRUs)

 call WriteTime(finalizeStats,outputTimeStep,time_meta,timeStruct%var,err,message)

 ! write the model output to the NetCDF file
 ! Passes the full metadata structure rather than the stats metadata structure because
 !  we have the option to write out data of types other than statistics.
 !  Thus, we must also pass the stats parent->child maps from childStruct.
 call writeData(finalizeStats,outputTimeStep,nHRUrun,maxLayers,forc_meta,forcStat,forcStruct,forcChild_map,indxStruct,err,message); call handle_err(err,message)
 call writeData(finalizeStats,outputTimeStep,nHRUrun,maxLayers,prog_meta,progStat,progStruct,progChild_map,indxStruct,err,message); call handle_err(err,message)
 call writeData(finalizeStats,outputTimeStep,nHRUrun,maxLayers,diag_meta,diagStat,diagStruct,diagChild_map,indxStruct,err,message); call handle_err(err,message)
 call writeData(finalizeStats,outputTimeStep,nHRUrun,maxLayers,flux_meta,fluxStat,fluxStruct,fluxChild_map,indxStruct,err,message); call handle_err(err,message)
 call writeData(finalizeStats,outputTimeStep,nHRUrun,maxLayers,indx_meta,indxStat,indxStruct,indxChild_map,indxStruct,err,message); call handle_err(err,message)

 ! increment output file timestep
 do iFreq = 1,maxvarFreq
  statCounter(iFreq) = statCounter(iFreq)+1
  if(finalizeStats(iFreq)) outputTimeStep(iFreq) = outputTimeStep(iFreq) + 1
 end do

 ! increment forcingStep
 forcingStep=forcingStep+1

 ! if finalized stats, then reset stats on the next time step
 resetStats(:) = finalizeStats(:)

 ! save time vector
 oldTimeVec(:) = timeStruct%var

 !print*, 'PAUSE: in driver: testing differences'; read(*,*)
 !stop 'end of time step'

 ! *****************************************************************************
 ! *** create a new NetCDF output file, and write parameters and forcing data
 ! *****************************************************************************

 ! define the need to create a new output file
 select case(newOutputFile)
  ! (don't ever create a new output file)
  case(noNewFiles); defNewOutputFile=.false.
  ! (check for the start of the USA water year)
  case(newFileEveryOct1)
   defNewOutputFile = (timeStruct%var(iLookTIME%im)  ==10 .and. &   ! month = October
                       timeStruct%var(iLookTIME%id)  ==1  .and. &   ! day = 1
                       timeStruct%var(iLookTIME%ih)  ==0  .and. &   ! hour = 1
                       timeStruct%var(iLookTIME%imin)==0)           ! minute = 0
  ! (check that we found the option)
  case default; call handle_err(20,'unable to identify the option to define new output files')
 end select

 ! create hte new output file
 if(defNewOutputFile)then

  ! close any output files that are already open
  do iFreq = 1,maxvarFreq
   if (ncid(iFreq)/=integerMissing) then
    call nc_file_close(ncid(iFreq),err,message)
    call handle_err(err,message)
   end if
  end do

  ! define the filename
  write(fileout,'(a,i0,a,i0,a)') trim(OUTPUT_PATH)//trim(OUTPUT_PREFIX),&
                                 timeStruct%var(iLookTIME%iyyy),'-',timeStruct%var(iLookTIME%iyyy)+1,&
                                 trim(output_fileSuffix)

  ! define the file
  call def_output(summaVersion,buildTime,gitBranch,gitHash,nGRU,nHRU,gru_struc(1)%hruInfo(1)%nSoil,fileout,err,message)
  call handle_err(err,message)

  ! write parameters for each HRU, and re-set indices
  do iGRU=1,nGRU
   do iHRU=1,gru_struc(iGRU)%hruCount
    call writeParm(gru_struc(iGRU)%hruInfo(iHRU)%hru_ix,attrStruct%gru(iGRU)%hru(iHRU),attr_meta,err,message); call handle_err(err,'[attr]/'//message)
    call writeParm(gru_struc(iGRU)%hruInfo(iHRU)%hru_ix,typeStruct%gru(iGRU)%hru(iHRU),type_meta,err,message); call handle_err(err,'[type]/'//message)
    call writeParm(gru_struc(iGRU)%hruInfo(iHRU)%hru_ix,mparStruct%gru(iGRU)%hru(iHRU),mpar_meta,err,message); call handle_err(err,'[mpar]'//message)
    ! re-initalize the indices for model writing
    outputTimeStep(:)=1
   end do  ! (looping through HRUs)
   call writeParm(integerMissing,bparStruct%gru(iGRU),bpar_meta,err,message); call handle_err(err,message)
  end do  ! (looping through GRUs)

 end if  ! if defining a new file

 ! *****************************************************************************
 ! *** write restart file
 ! *****************************************************************************

 ! query whether this timestep requires a re-start file
 select case(ixRestart)
  case(ixRestart_iy);    printRestart = (timeStruct%var(iLookTIME%im) == 1 .and. timeStruct%var(iLookTIME%id) == 1 .and. timeStruct%var(iLookTIME%ih) == 0  .and. timeStruct%var(iLookTIME%imin) == 0)
  case(ixRestart_im);    printRestart = (timeStruct%var(iLookTIME%id) == 1 .and. timeStruct%var(iLookTIME%ih) == 0 .and. timeStruct%var(iLookTIME%imin) == 0)
  case(ixRestart_id);    printRestart = (timeStruct%var(iLookTIME%ih) == 0 .and. timeStruct%var(iLookTIME%imin) == 0)
  case(ixRestart_end);   printRestart = (timeStruct%var(iLookTIME%im) == finshTime%var(2) .and. timeStruct%var(iLookTIME%id) == finshTime%var(3) .and. timeStruct%var(iLookTIME%ih) == finshTime%var(4)  .and. timeStruct%var(iLookTIME%imin) == finshTime%var(5))
  case(ixRestart_never); printRestart = .false.
  case default; call handle_err(20,'unable to identify option for the restart file')
 end select

 ! print a restart file if requested
 if(printRestart)then
  write(timeString,'(a,i4,3(a,i2.2))') '_',timeStruct%var(iLookTIME%iyyy),'-',timeStruct%var(iLookTIME%im),'-',timeStruct%var(iLookTIME%id),'-',timeStruct%var(iLookTIME%ih)
  restartFile=trim(OUTPUT_PATH)//trim(OUTPUT_PREFIX)//'_'//trim('summaRestart')//trim(timeString)//trim(output_fileSuffix)//'.nc'
  call writeRestart(restartFile,nGRU,nHRU,prog_meta,progStruct,maxLayers,maxSnowLayers,indx_meta,indxStruct,err,message)
  call handle_err(err,message)
 end if

end do  ! (looping through time)

! close any remaining output files
do iFreq = 1,maxvarFreq
 if (ncid(iFreq).ne.integerMissing) then
  call nc_file_close(ncid(iFreq),err,message)
  call handle_err(err,message)
 end if
end do

! deallocate space for dt_init and upArea
deallocate(dt_init,upArea,stat=err); call handle_err(err,'unable to deallocate space for dt_init and upArea')

call stop_program('finished simulation successfully.')

contains

 ! **************************************************************************************************
 ! internal function to obtain the command line arguments
 ! **************************************************************************************************
 subroutine getCommandArguments()
 implicit none
 integer(i4b)                     :: iArgument                  ! index of command line argument
 integer(i4b)                     :: nArgument                  ! number of command line arguments
 character(len=256),allocatable   :: argString(:)               ! string to store command line arguments
 integer(i4b)                     :: nLocalArgument             ! number of command line arguments to read for a switch
 character(len=70), parameter     :: spaces = ''
 nArgument = command_argument_count()
 ! check numbers of command-line arguments and obtain all arguments
 if (nArgument < 1) then
  call printCommandHelp()
 end if

 allocate(argString(nArgument))
 do iArgument = 1,nArgument
  call get_command_argument(iArgument,argString(iArgument))
  ! print versions if needed
  if (trim(argString(iArgument)) == '-v' .or. trim(argString(iArgument)) == '--version') then
   ! print version numbers

   print "(A)", '----------------------------------------------------------------------'
   print "(A)", '     SUMMA - Structure for Unifying Multiple Modeling Alternatives    '
   print "(A)", spaces(1:int((70 - len_trim(summaVersion) - 9) / 2))//'Version: '   //trim(summaVersion)
   print "(A)", spaces(1:int((70 - len_trim(buildTime) - 12) / 2))  //'Build Time: '//trim(buildTime)
   print "(A)", spaces(1:int((70 - len_trim(gitBranch) - 12) / 2))  //'Git Branch: '//trim(gitBranch)
   print "(A)", spaces(1:int((70 - len_trim(gitHash) - 10) / 2))    //'Git Hash: '  //trim(gitHash)
   print "(A)", '----------------------------------------------------------------------'
   if (nArgument == 1) stop
  end if
 end do

 ! initialize command line argument variables
 startGRU = integerMissing; checkHRU = integerMissing
 nGRU = integerMissing; nHRU = integerMissing
 newOutputFile = noNewFiles
 iRunMode = iRunModeFull

 ! loop through all command arguments
 nLocalArgument = 0
 do iArgument = 1,nArgument
  if (nLocalArgument>0) then; nLocalArgument = nLocalArgument -1; cycle; end if ! skip the arguments have been read
  select case (trim(argString(iArgument)))

   case ('-m', '--master')
    ! update arguments
    nLocalArgument = 1
    if (iArgument+nLocalArgument>nArgument) call handle_err(1,"missing argument file_suffix; type 'summa.exe --help' for correct usage")
    ! get name of master control file
    summaFileManagerFile=trim(argString(iArgument+1))
    print "(A)", "file_master is '"//trim(summaFileManagerFile)//"'."

   ! define the formation of new output files
   case ('-n', '--newFile')
    ! check that the number of command line arguments is correct
    nLocalArgument = 1  ! expect just one argument for new output files
    if (iArgument+nLocalArgument>nArgument) call handle_err(1,"missing argument file_suffix; type 'summa.exe --help' for correct usage")
    ! get the decision for the formation of new output files
    select case( trim(argString(iArgument+1)) )
     case('noNewFiles');       newOutputFile = noNewFiles
     case('newFileEveryOct1'); newOutputFile = newFileEveryOct1
     case default;             call handle_err(1,'unknown option for new output file: expect "noNewFiles" or "newFileEveryOct1"')
    end select

   case ('-s', '--suffix')
    ! define file suffix
    nLocalArgument = 1
    ! check if the number of command line arguments is correct
    if (iArgument+nLocalArgument>nArgument) call handle_err(1,"missing argument file_suffix; type 'summa.exe --help' for correct usage")
    output_fileSuffix=trim(argString(iArgument+1))
    print "(A)", "file_suffix is '"//trim(output_fileSuffix)//"'."

   case ('-h', '--hru')
    ! define a single HRU run
    if (iRunMode == iRunModeGRU) call handle_err(1,"single-HRU run and GRU-parallelization run cannot be both selected.")
    iRunMode=iRunModeHRU
    nLocalArgument = 1
    ! check if the number of command line arguments is correct
    if (iArgument+nLocalArgument>nArgument) call handle_err(1,"missing argument checkHRU; type 'summa.exe --help' for correct usage")
    read(argString(iArgument+1),*) checkHRU ! read the index of the HRU for a single HRU run
    nHRU=1; nGRU=1                          ! nHRU and nGRU are both one in this case
    ! examines the checkHRU is correct
    if (checkHRU<1) then
     call handle_err(1,"illegal iHRU specification; type 'summa.exe --help' for correct usage")
    else
     print '(A)',' Single-HRU run activated. HRU '//trim(argString(iArgument+1))//' is selected for simulation.'
    end if

   case ('-g','--gru')
    ! define a GRU parallelization run; get the starting GRU and countGRU
    if (iRunMode == iRunModeHRU) call handle_err(1,"single-HRU run and GRU-parallelization run cannot be both selected.")
    iRunMode=iRunModeGRU
    nLocalArgument = 2
    ! check if the number of command line arguments is correct
    if (iArgument+nLocalArgument>nArgument) call handle_err(1,"missing argument startGRU or countGRU; type 'summa.exe --help' for correct usage")
    read(argString(iArgument+1),*) startGRU ! read the argument of startGRU
    read(argString(iArgument+2),*) nGRU     ! read the argument of countGRU
    if (startGRU<1 .or. nGRU<1) then
     call handle_err(1,'startGRU and countGRU must be larger than 1.')
    else
     print '(A)', ' GRU-Parallelization run activated. '//trim(argString(iArgument+2))//' GRUs are selected for simulation.'
    end if

   case ('-p', '--progress')
    ! define the frequency to print progress
    nLocalArgument = 1
    ! check if the number of command line arguments is correct
    if (iArgument+nLocalArgument>nArgument) call handle_err(1, "missing argument freqProgress; type 'summa.exe --help' for correct usage")
    select case (trim(argString(iArgument+1)))
     case ('m' , 'month'); ixProgress = ixProgress_im
     case ('d' , 'day');   ixProgress = ixProgress_id
     case ('h' , 'hour');  ixProgress = ixProgress_ih
     case ('n' , 'never'); ixProgress = ixProgress_never
     case default;         call handle_err(1,'unknown frequency to print progress')
    end select

   case ('-r', '--restart')
    ! define the frequency to write restart files
    nLocalArgument = 1
    ! check if the number of command line arguments is correct
    if (iArgument+nLocalArgument>nArgument) call handle_err(1, "missing argument freqRestart; type 'summa.exe --help' for correct usage")
    select case (trim(argString(iArgument+1)))
     case ('y' , 'year');  ixRestart = ixRestart_iy
     case ('m' , 'month'); ixRestart = ixRestart_im
     case ('d' , 'day');   ixRestart = ixRestart_id
     case ('e' , 'end');   ixRestart = ixRestart_end
     case ('n' , 'never'); ixRestart = ixRestart_never
     case default;         call handle_err(1,'unknown frequency to write restart files')
    end select

   ! do nothing
   case ('-v','--version')

   ! print help message
   case ('--help')
    call printCommandHelp

   case default
    call printCommandHelp
    call handle_err(1, 'unknown command line option')

  end select
 end do  ! looping through command line arguments

 ! check if master_file has been received.
 if (len(trim(summaFileManagerFile))==0) call handle_err(1, "master_file is not received; type 'summa.exe --help' for correct usage")

 ! set startGRU for full run
 if (iRunMode==iRunModeFull) startGRU=1

 end subroutine getCommandArguments

 ! **************************************************************************************************
 ! internal subroutine to print the correct command line usage of SUMMA
 ! **************************************************************************************************
 subroutine printCommandHelp()
 implicit none
 ! command line usage
 print "(//A)",'Usage: summa.exe -m master_file [-s fileSuffix] [-g startGRU countGRU] [-h iHRU] [-r freqRestart] [-p freqProgress] [-c]'
 print "(A,/)",  ' summa.exe          summa executable'
 print "(A)",  'Running options:'
 print "(A)",  ' -m --master        Define path/name of master file (required)'
 print "(A)",  ' -n --newFile       Define frequency [noNewFiles,newFileEveryOct1] of new output files'
 print "(A)",  ' -s --suffix        Add fileSuffix to the output files'
 print "(A)",  ' -g --gru           Run a subset of countGRU GRUs starting from index startGRU'
 print "(A)",  ' -h --hru           Run a single HRU with index of iHRU'
 print "(A)",  ' -r --restart       Define frequency [y,m,d,e,never] to write restart files'
 print "(A)",  ' -p --progress      Define frequency [m,d,h,never] to print progress'
 print "(A)",  ' -v --version       Display version information of the current built'
 stop
 end subroutine printCommandHelp

 ! **************************************************************************************************
 ! internal subroutine handle_err: error handler
 ! **************************************************************************************************
 subroutine handle_err(err,message)
 ! used to handle error codes
 USE var_lookup,only:iLookPROG,iLookDIAG,iLookFLUX,iLookPARAM,iLookINDEX    ! named variables defining elements in data structure
 implicit none
 ! dummy variables
 integer(i4b),intent(in) :: err             ! error code
 character(*),intent(in) :: message         ! error message
 ! local variables
 integer(i4b)            :: nc_err          ! error code of nc_close
 character(len=256)      :: cmessage        ! error message of the downwind routine

 ! return if A-OK
 if(err==0) return
 ! process error messages
 if (err>0) then
  write(*,'(//a/)') 'FATAL ERROR: '//trim(message)
 else
  write(*,'(//a/)') 'WARNING: '//trim(message); print*,'(can keep going, but stopping anyway)'
 endif
 ! dump variables
 print*, 'error, variable dump:'
 if(allocated(timeStruct%var))then
  ! print time step
  print*, 'modelTimeStep = ', modelTimeStep
  ! print information for the HRUs
  if(iGRU<=nGRU)then
   if(iHRU<=gru_struc(iGRU)%hruCount)then
    print*, 'initial time step  = ', dt_init(iGRU)%hru(iHRU)
    print*, 'HRU index          = ', typeStruct%gru(iGRU)%hru(iHRU)%var(iLookTYPE%hruId)
    print*, 'pptrate            = ', forcStruct%gru(iGRU)%hru(iHRU)%var(iLookFORCE%pptrate)
    print*, 'airtemp            = ', forcStruct%gru(iGRU)%hru(iHRU)%var(iLookFORCE%airtemp)
    print*, 'theta_res          = ', mparStruct%gru(iGRU)%hru(iHRU)%var(iLookPARAM%theta_res)%dat(1)            ! soil residual volumetric water content (-)
    print*, 'theta_sat          = ', mparStruct%gru(iGRU)%hru(iHRU)%var(iLookPARAM%theta_sat)%dat(1)            ! soil porosity (-)
    print*, 'plantWiltPsi       = ', mparStruct%gru(iGRU)%hru(iHRU)%var(iLookPARAM%plantWiltPsi)%dat(1)         ! matric head at wilting point (m)
    print*, 'soilStressParam    = ', mparStruct%gru(iGRU)%hru(iHRU)%var(iLookPARAM%soilStressParam)%dat(1)      ! parameter in the exponential soil stress function (-)
    print*, 'critSoilWilting    = ', mparStruct%gru(iGRU)%hru(iHRU)%var(iLookPARAM%critSoilWilting)%dat(1)      ! critical vol. liq. water content when plants are wilting (-)
    print*, 'critSoilTranspire  = ', mparStruct%gru(iGRU)%hru(iHRU)%var(iLookPARAM%critSoilTranspire)%dat(1)    ! critical vol. liq. water content when transpiration is limited (-)
    print*, 'scalarSWE          = ', progStruct%gru(iGRU)%hru(iHRU)%var(iLookPROG%scalarSWE)%dat(1)
    print*, 'scalarSnowDepth    = ', progStruct%gru(iGRU)%hru(iHRU)%var(iLookPROG%scalarSnowDepth)%dat(1)
    print*, 'scalarCanopyTemp   = ', progStruct%gru(iGRU)%hru(iHRU)%var(iLookPROG%scalarCanopyTemp)%dat(1)
    print*, 'scalarRainPlusMelt = ', fluxStruct%gru(iGRU)%hru(iHRU)%var(iLookFLUX%scalarRainPlusMelt)%dat(1)
    write(*,'(a,100(i4,1x))'   ) 'layerType          = ', indxStruct%gru(iGRU)%hru(iHRU)%var(iLookINDEX%layerType)%dat
    write(*,'(a,100(f11.5,1x))') 'mLayerDepth        = ', progStruct%gru(iGRU)%hru(iHRU)%var(iLookPROG%mLayerDepth)%dat
    write(*,'(a,100(f11.5,1x))') 'mLayerTemp         = ', progStruct%gru(iGRU)%hru(iHRU)%var(iLookPROG%mLayerTemp)%dat
    write(*,'(a,100(f11.5,1x))') 'mLayerVolFracIce   = ', progStruct%gru(iGRU)%hru(iHRU)%var(iLookPROG%mLayerVolFracIce)%dat
    write(*,'(a,100(f11.5,1x))') 'mLayerVolFracLiq   = ', progStruct%gru(iGRU)%hru(iHRU)%var(iLookPROG%mLayerVolFracLiq)%dat
    print*, 'mLayerMatricHead   = ', progStruct%gru(iGRU)%hru(iHRU)%var(iLookPROG%mLayerMatricHead)%dat
    print*, 'column inflow      = ', fluxStruct%gru(iGRU)%hru(iHRU)%var(iLookFLUX%mLayerColumnInflow)%dat
   endif  ! if HRU is valid
  endif  ! if GRU is valid
 endif  ! if the time structure is allocated
 print*,'error code = ', err
 if(allocated(timeStruct%var)) print*, timeStruct%var
 !write(*,'(a)') trim(message)

 ! close any remaining output files
 do iFreq = 1,maxvarFreq
  if (ncid(iFreq).ne.integerMissing) then
   call nc_file_close(ncid(iFreq),nc_err,cmessage)
   if(nc_err/=0) print*, trim(cmessage)
  end if
 end do

 stop 1
 end subroutine handle_err

 ! **************************************************************************************************
 ! private subroutine stop_program: stop program execution
 ! **************************************************************************************************
 subroutine stop_program(message)
 ! used to stop program execution
 implicit none
 ! define dummy variables
 character(*),intent(in)::message
 ! define the local variables
 integer(i4b),parameter :: outunit=6               ! write to screen
 integer(i4b)           :: ctime2(8)               ! final time
 real(dp)               :: elpSec                  ! elapsed seconds

 ! close any remaining output files
 ! NOTE: use the direct NetCDF call with no error checking since the file may already be closed
 do iFreq = 1,maxvarFreq
  if (ncid(iFreq).ne.integerMissing) then
   err = nf90_close(ncid(iFreq))
  end if
 end do

 ! get the final date and time
 call date_and_time(values=ctime2)

 elpSec = elapsedSec(ctime1,ctime2)

 ! print initial and final date and time
 write(outunit,"(A,I4,'-',I2.2,'-',I2.2,2x,I2,':',I2.2,':',I2.2,'.',I3.3)") 'initial date/time = ',ctime1(1:3),ctime1(5:8)
 write(outunit,"(A,I4,'-',I2.2,'-',I2.2,2x,I2,':',I2.2,':',I2.2,'.',I3.3)") '  final date/time = ',ctime2(1:3),ctime2(5:8)
 ! print elapsed time
 write(outunit,"(/,A,1PG15.7,A)")                                           '     elapsed time = ', elpSec,          ' s'
 write(outunit,"(A,1PG15.7,A)")                                             '       or           ', elpSec/60_dp,    ' m'
 write(outunit,"(A,1PG15.7,A)")                                             '       or           ', elpSec/3600_dp,  ' h'
 write(outunit,"(A,1PG15.7,A/)")                                            '       or           ', elpSec/86400_dp, ' d'
 ! stop with message
 print*,'FORTRAN STOP: '//trim(message)
 stop
 end subroutine

end program multi_driver


 ! **************************************************************************************************
 ! private subroutine SOIL_VEG_GEN_PARM: Read soil, vegetation and other model parameters (from NOAH)
 ! **************************************************************************************************
!-----------------------------------------------------------------
SUBROUTINE SOIL_VEG_GEN_PARM(FILENAME_VEGTABLE, FILENAME_SOILTABLE, FILENAME_GENERAL, MMINLU, MMINSL)
!-----------------------------------------------------------------
  use module_sf_noahlsm, only : shdtbl, nrotbl, rstbl, rgltbl, &
       &                        hstbl, snuptbl, maxalb, laimintbl, &
       &                        bb, drysmc, f11, maxsmc, laimaxtbl, &
       &                        emissmintbl, emissmaxtbl, albedomintbl, &
       &                        albedomaxtbl, wltsmc, qtz, refsmc, &
       &                        z0mintbl, z0maxtbl, &
       &                        satpsi, satdk, satdw, &
       &                        theta_res, theta_sat, vGn_alpha, vGn_n, k_soil, &  ! MPC add van Genutchen parameters
       &                        fxexp_data, lvcoef_data, &
       &                        lutype, maxalb, &
       &                        slope_data, frzk_data, bare, cmcmax_data, &
       &                        cfactr_data, csoil_data, czil_data, &
       &                        refkdt_data, natural, refdk_data, &
       &                        rsmax_data, salp_data, sbeta_data, &
       &                        zbot_data, smhigh_data, smlow_data, &
       &                        lucats, topt_data, slcats, slpcats, sltype

  IMPLICIT NONE

  CHARACTER(LEN=*), INTENT(IN) :: FILENAME_VEGTABLE, FILENAME_SOILTABLE, FILENAME_GENERAL
  CHARACTER(LEN=*), INTENT(IN) :: MMINLU, MMINSL
  integer :: LUMATCH, IINDEX, LC, NUM_SLOPE
  integer :: ierr
  INTEGER , PARAMETER :: OPEN_OK = 0

  character*128 :: mess , message

!-----SPECIFY VEGETATION RELATED CHARACTERISTICS :
!             ALBBCK: SFC albedo (in percentage)
!                 Z0: Roughness length (m)
!             SHDFAC: Green vegetation fraction (in percentage)
!  Note: The ALBEDO, Z0, and SHDFAC values read from the following table
!          ALBEDO, amd Z0 are specified in LAND-USE TABLE; and SHDFAC is
!          the monthly green vegetation data
!             CMXTBL: MAX CNPY Capacity (m)
!             NROTBL: Rooting depth (layer)
!              RSMIN: Mimimum stomatal resistance (s m-1)
!              RSMAX: Max. stomatal resistance (s m-1)
!                RGL: Parameters used in radiation stress function
!                 HS: Parameter used in vapor pressure deficit functio
!               TOPT: Optimum transpiration air temperature. (K)
!             CMCMAX: Maximum canopy water capacity
!             CFACTR: Parameter used in the canopy inteception calculati
!               SNUP: Threshold snow depth (in water equivalent m) that
!                     implies 100% snow cover
!                LAI: Leaf area index (dimensionless)
!             MAXALB: Upper bound on maximum albedo over deep snow
!
!-----READ IN VEGETAION PROPERTIES FROM VEGPARM.TBL
!

  OPEN(19, FILE=trim(FILENAME_VEGTABLE),FORM='FORMATTED',STATUS='OLD',IOSTAT=ierr)
  IF(ierr .NE. OPEN_OK ) THEN
     WRITE(message,FMT='(A)') &
          'module_sf_noahlsm.F: soil_veg_gen_parm: failure opening VEGPARM.TBL'
     CALL wrf_error_fatal ( message )
  END IF

  LUMATCH=0

  FIND_LUTYPE : DO WHILE (LUMATCH == 0)
     READ (19,*,END=2002)
     READ (19,*,END=2002)LUTYPE
     READ (19,*)LUCATS,IINDEX

     IF(LUTYPE.EQ.MMINLU)THEN
        WRITE( mess , * ) 'LANDUSE TYPE = ' // TRIM ( LUTYPE ) // ' FOUND', LUCATS,' CATEGORIES'
        ! CALL wrf_message( mess )
        LUMATCH=1
     ELSE
        call wrf_message ( "Skipping over LUTYPE = " // TRIM ( LUTYPE ) )
        DO LC = 1, LUCATS+12
           read(19,*)
        ENDDO
     ENDIF
  ENDDO FIND_LUTYPE
! prevent possible array overwrite, Bill Bovermann, IBM, May 6, 2008
  IF ( SIZE(SHDTBL)       < LUCATS .OR. &
       SIZE(NROTBL)       < LUCATS .OR. &
       SIZE(RSTBL)        < LUCATS .OR. &
       SIZE(RGLTBL)       < LUCATS .OR. &
       SIZE(HSTBL)        < LUCATS .OR. &
       SIZE(SNUPTBL)      < LUCATS .OR. &
       SIZE(MAXALB)       < LUCATS .OR. &
       SIZE(LAIMINTBL)    < LUCATS .OR. &
       SIZE(LAIMAXTBL)    < LUCATS .OR. &
       SIZE(Z0MINTBL)     < LUCATS .OR. &
       SIZE(Z0MAXTBL)     < LUCATS .OR. &
       SIZE(ALBEDOMINTBL) < LUCATS .OR. &
       SIZE(ALBEDOMAXTBL) < LUCATS .OR. &
       SIZE(EMISSMINTBL ) < LUCATS .OR. &
       SIZE(EMISSMAXTBL ) < LUCATS ) THEN
     CALL wrf_error_fatal('Table sizes too small for value of LUCATS in module_sf_noahdrv.F')
  ENDIF

  IF(LUTYPE.EQ.MMINLU)THEN
     DO LC=1,LUCATS
        READ (19,*)IINDEX,SHDTBL(LC),                        &
             NROTBL(LC),RSTBL(LC),RGLTBL(LC),HSTBL(LC), &
             SNUPTBL(LC),MAXALB(LC), LAIMINTBL(LC),     &
             LAIMAXTBL(LC),EMISSMINTBL(LC),             &
             EMISSMAXTBL(LC), ALBEDOMINTBL(LC),         &
             ALBEDOMAXTBL(LC), Z0MINTBL(LC), Z0MAXTBL(LC)
     ENDDO

     READ (19,*)
     READ (19,*)TOPT_DATA
     READ (19,*)
     READ (19,*)CMCMAX_DATA
     READ (19,*)
     READ (19,*)CFACTR_DATA
     READ (19,*)
     READ (19,*)RSMAX_DATA
     READ (19,*)
     READ (19,*)BARE
     READ (19,*)
     READ (19,*)NATURAL
  ENDIF

2002 CONTINUE

  CLOSE (19)
  IF (LUMATCH == 0) then
     CALL wrf_error_fatal ("Land Use Dataset '"//MMINLU//"' not found in VEGPARM.TBL.")
  ENDIF

!
!-----READ IN SOIL PROPERTIES FROM SOILPARM.TBL
!
  OPEN(19, FILE=trim(FILENAME_SOILTABLE),FORM='FORMATTED',STATUS='OLD',IOSTAT=ierr)
  IF(ierr .NE. OPEN_OK ) THEN
     WRITE(message,FMT='(A)') &
          'module_sf_noahlsm.F: soil_veg_gen_parm: failure opening SOILPARM.TBL'
     CALL wrf_error_fatal ( message )
  END IF

  WRITE(mess,*) 'INPUT SOIL TEXTURE CLASSIFICATION = ', TRIM ( MMINSL )
  ! CALL wrf_message( mess )

  LUMATCH=0

  ! MPC add a new soil table
  FIND_soilTYPE : DO WHILE (LUMATCH == 0)
   READ (19,*)
   READ (19,*,END=2003)SLTYPE
   READ (19,*)SLCATS,IINDEX
   IF(SLTYPE.EQ.MMINSL)THEN
     WRITE( mess , * ) 'SOIL TEXTURE CLASSIFICATION = ', TRIM ( SLTYPE ) , ' FOUND', &
          SLCATS,' CATEGORIES'
     ! CALL wrf_message ( mess )
     LUMATCH=1
   ELSE
    call wrf_message ( "Skipping over SLTYPE = " // TRIM ( SLTYPE ) )
    DO LC = 1, SLCATS
     read(19,*)
    ENDDO
   ENDIF
  ENDDO FIND_soilTYPE
  ! prevent possible array overwrite, Bill Bovermann, IBM, May 6, 2008
  IF ( SIZE(BB    ) < SLCATS .OR. &
       SIZE(DRYSMC) < SLCATS .OR. &
       SIZE(F11   ) < SLCATS .OR. &
       SIZE(MAXSMC) < SLCATS .OR. &
       SIZE(REFSMC) < SLCATS .OR. &
       SIZE(SATPSI) < SLCATS .OR. &
       SIZE(SATDK ) < SLCATS .OR. &
       SIZE(SATDW ) < SLCATS .OR. &
       SIZE(WLTSMC) < SLCATS .OR. &
       SIZE(QTZ   ) < SLCATS  ) THEN
     CALL wrf_error_fatal('Table sizes too small for value of SLCATS in module_sf_noahdrv.F')
  ENDIF

  ! MPC add new soil table
  select case(trim(SLTYPE))
   case('STAS','STAS-RUC')  ! original soil tables
     DO LC=1,SLCATS
        READ (19,*) IINDEX,BB(LC),DRYSMC(LC),F11(LC),MAXSMC(LC),&
             REFSMC(LC),SATPSI(LC),SATDK(LC), SATDW(LC),   &
             WLTSMC(LC), QTZ(LC)
     ENDDO
   case('ROSETTA')          ! new soil table
     DO LC=1,SLCATS
        READ (19,*) IINDEX,&
             ! new soil parameters (from Rosetta)
             theta_res(LC), theta_sat(LC),        &
             vGn_alpha(LC), vGn_n(LC), k_soil(LC), &
             ! original soil parameters
             BB(LC),DRYSMC(LC),F11(LC),MAXSMC(LC),&
             REFSMC(LC),SATPSI(LC),SATDK(LC), SATDW(LC),   &
             WLTSMC(LC), QTZ(LC)
     ENDDO
   case default
     CALL wrf_message( 'SOIL TEXTURE IN INPUT FILE DOES NOT ' )
     CALL wrf_message( 'MATCH SOILPARM TABLE'                 )
     CALL wrf_error_fatal ( 'INCONSISTENT OR MISSING SOILPARM FILE' )
  end select

2003 CONTINUE

  CLOSE (19)

  IF(LUMATCH.EQ.0)THEN
     CALL wrf_message( 'SOIL TEXTURE IN INPUT FILE DOES NOT ' )
     CALL wrf_message( 'MATCH SOILPARM TABLE'                 )
     CALL wrf_error_fatal ( 'INCONSISTENT OR MISSING SOILPARM FILE' )
  ENDIF

!
!-----READ IN GENERAL PARAMETERS FROM GENPARM.TBL
!
  OPEN(19, FILE=trim(FILENAME_GENERAL),FORM='FORMATTED',STATUS='OLD',IOSTAT=ierr)
  IF(ierr .NE. OPEN_OK ) THEN
     WRITE(message,FMT='(A)') &
          'module_sf_noahlsm.F: soil_veg_gen_parm: failure opening GENPARM.TBL'
     CALL wrf_error_fatal ( message )
  END IF

  READ (19,*)
  READ (19,*)
  READ (19,*) NUM_SLOPE

  SLPCATS=NUM_SLOPE
! prevent possible array overwrite, Bill Bovermann, IBM, May 6, 2008
  IF ( SIZE(slope_data) < NUM_SLOPE ) THEN
     CALL wrf_error_fatal('NUM_SLOPE too large for slope_data array in module_sf_noahdrv')
  ENDIF

  DO LC=1,SLPCATS
     READ (19,*)SLOPE_DATA(LC)
  ENDDO

  READ (19,*)
  READ (19,*)SBETA_DATA
  READ (19,*)
  READ (19,*)FXEXP_DATA
  READ (19,*)
  READ (19,*)CSOIL_DATA
  READ (19,*)
  READ (19,*)SALP_DATA
  READ (19,*)
  READ (19,*)REFDK_DATA
  READ (19,*)
  READ (19,*)REFKDT_DATA
  READ (19,*)
  READ (19,*)FRZK_DATA
  READ (19,*)
  READ (19,*)ZBOT_DATA
  READ (19,*)
  READ (19,*)CZIL_DATA
  READ (19,*)
  READ (19,*)SMLOW_DATA
  READ (19,*)
  READ (19,*)SMHIGH_DATA
  READ (19,*)
  READ (19,*)LVCOEF_DATA
  CLOSE (19)

!-----------------------------------------------------------------
END SUBROUTINE SOIL_VEG_GEN_PARM
!-----------------------------------------------------------------

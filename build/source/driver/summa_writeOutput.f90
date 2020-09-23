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

module summa_writeOutput
! used to define/write output files

! named variables to define new output files
USE globalData, only: noNewFiles              ! no new output files
USE globalData, only: newFileEveryOct1        ! create a new file on Oct 1 every year (start of the USA water year)

! metadata
USE globalData,only:time_meta                 ! metadata on the model time
USE globalData,only:forc_meta                 ! metadata on the model forcing data
USE globalData,only:diag_meta                 ! metadata on the model diagnostic variables
USE globalData,only:prog_meta                 ! metadata on the model prognostic variables
USE globalData,only:flux_meta                 ! metadata on the model fluxes
USE globalData,only:indx_meta                 ! metadata on the model index variables
USE globalData,only:bvar_meta                 ! metadata on basin-average variables

! child metadata for stats
USE globalData,only:statForc_meta             ! child metadata for stats
USE globalData,only:statProg_meta             ! child metadata for stats
USE globalData,only:statDiag_meta             ! child metadata for stats
USE globalData,only:statFlux_meta             ! child metadata for stats
USE globalData,only:statIndx_meta             ! child metadata for stats
USE globalData,only:statBvar_meta             ! child metadata for stats

! index of the child data structure
USE globalData,only:forcChild_map             ! index of the child data structure: stats forc
USE globalData,only:progChild_map             ! index of the child data structure: stats prog
USE globalData,only:diagChild_map             ! index of the child data structure: stats diag
USE globalData,only:fluxChild_map             ! index of the child data structure: stats flux
USE globalData,only:indxChild_map             ! index of the child data structure: stats indx
USE globalData,only:bvarChild_map             ! index of the child data structure: stats bvar

! named variables
USE var_lookup,only:maxvarFreq                ! maximum number of output files
USE var_lookup,only:iLookTIME                 ! named variables for time data structure
USE var_lookup,only:iLookDIAG                 ! named variables for local column model diagnostic variables
USE var_lookup,only:iLookPROG                 ! named variables for local column model prognostic variables
USE var_lookup,only:iLookINDEX                ! named variables for local column index variables
USE var_lookup,only:iLookFreq                 ! named variables for the frequency structure

! safety: set private unless specified otherwise
implicit none
private
public::summa_writeOutputFiles
contains

 ! used to define/write output files
 subroutine summa_writeOutputFiles(modelTimeStep, summa1_struc, err, message)
 ! ---------------------------------------------------------------------------------------
 ! * desired modules
 ! ---------------------------------------------------------------------------------------
 ! data types
 USE nrtype                                                  ! variable types, etc.
 USE summa_type,only:summa1_type_dec                         ! master summa data type
 ! subroutines and functions
 USE time_utils_module,only:elapsedSec                       ! calculate the elapsed time
 USE summa_alarms,only:summa_setWriteAlarms                  ! set alarms to control model output
 USE summa_defineOutput,only:summa_defineOutputFiles         ! define summa output files
 USE modelwrite_module,only:writeRestart                     ! module to write model Restart
 USE modelwrite_module,only:writeData,writeBasin             ! module to write model output
 USE modelwrite_module,only:writeTime                        ! module to write model time
 USE output_stats,only:calcStats                             ! module for compiling output statistics
 ! global data: general
 USE globalData,only:forcingStep                             ! index of current time step in current forcing file
 USE globalData,only:gru_struc                               ! gru-hru mapping structures
 USE globalData,only:structInfo                              ! information on the data structures
 ! global data: time structures
 USE globalData,only:oldTime                                 ! time from the previous time step
 USE globalData,only:finshTime                               ! end time of simulation
 ! global data: decisions for model alarms
 USE globalData,only:ixProgress                              ! define frequency to write progress
 USE globalData,only:ixRestart                               ! define frequency to write restart files
 USE globalData,only:newOutputFile                           ! define option for new output files
 ! controls on statistics output
 USE globalData,only:statCounter                             ! time counter for stats
 USE globalData,only:resetStats                              ! flags to reset statistics
 USE globalData,only:finalizeStats                           ! flags to finalize statistics
 USE globalData,only:outputTimeStep                          ! timestep in output files
 ! output constraints
 USE globalData,only:maxLayers                               ! maximum number of layers
 USE globalData,only:maxSnowLayers                           ! maximum number of snow layers
 ! timing variables
 USE globalData,only:startWrite,endWrite                     ! date/time for the start and end of the model writing
 USE globalData,only:elapsedWrite                            ! elapsed time to write data
 ! file information
 USE summaFileManager,only:OUTPUT_PATH,OUTPUT_PREFIX         ! define output file
 USE summaFileManager,only:STATE_PATH                        ! optional path to state output files (defaults to OUTPUT_PATH)
 USE globalData,only:output_fileSuffix                       ! suffix for the output & state files (optional summa argument)
 USE globalData,only:nHRUrun                                 ! number of HRU in the run
 USE globalData,only:nGRUrun                                 ! number of GRU in the run
 ! ---------------------------------------------------------------------------------------
 ! * variables
 ! ---------------------------------------------------------------------------------------
 implicit none
 ! dummy variables
 integer(i4b),intent(in)               :: modelTimeStep      ! time step index
 type(summa1_type_dec),intent(inout)   :: summa1_struc       ! master summa data structure
 integer(i4b),intent(out)              :: err                ! error code
 character(*),intent(out)              :: message            ! error message
 ! local variables
 character(LEN=256)                    :: cmessage           ! error message of downwind routine
 character(len=256)                    :: timeString                 ! portion of restart file name that contains the write-out time
 character(len=256)                    :: restartFile                ! restart file name
 logical(lgt)                          :: printRestart=.false.       ! flag to print a re-start file
 logical(lgt)                          :: printProgress=.false.      ! flag to print simulation progress
 logical(lgt)                          :: defNewOutputFile=.false.   ! flag to define new output files
 integer(i4b)                          :: iGRU,iHRU          ! indices of GRUs and HRUs
 integer(i4b)                          :: iStruct            ! index of model structure
 integer(i4b)                          :: iFreq              ! index of the output frequency
 ! ---------------------------------------------------------------------------------------
 ! associate to elements in the data structure
 summaVars: associate(&

  ! statistics structures
  forcStat             => summa1_struc%forcStat    , & ! x%gru(:)%hru(:)%var(:)%dat -- model forcing data
  progStat             => summa1_struc%progStat    , & ! x%gru(:)%hru(:)%var(:)%dat -- model prognostic (state) variables
  diagStat             => summa1_struc%diagStat    , & ! x%gru(:)%hru(:)%var(:)%dat -- model diagnostic variables
  fluxStat             => summa1_struc%fluxStat    , & ! x%gru(:)%hru(:)%var(:)%dat -- model fluxes
  indxStat             => summa1_struc%indxStat    , & ! x%gru(:)%hru(:)%var(:)%dat -- model indices
  bvarStat             => summa1_struc%bvarStat    , & ! x%gru(:)%var(:)%dat        -- basin-average variabl

  ! primary data structures
  timeStruct           => summa1_struc%timeStruct  , & ! x%var(:)                   -- model time data
  forcStruct           => summa1_struc%forcStruct  , & ! x%gru(:)%hru(:)%var(:)     -- model forcing data
  indxStruct           => summa1_struc%indxStruct  , & ! x%gru(:)%hru(:)%var(:)%dat -- model indices
  progStruct           => summa1_struc%progStruct  , & ! x%gru(:)%hru(:)%var(:)%dat -- model prognostic (state) variables
  diagStruct           => summa1_struc%diagStruct  , & ! x%gru(:)%hru(:)%var(:)%dat -- model diagnostic variables
  fluxStruct           => summa1_struc%fluxStruct  , & ! x%gru(:)%hru(:)%var(:)%dat -- model fluxes
  bvarStruct           => summa1_struc%bvarStruct  , & ! x%gru(:)%var(:)%dat        -- basin-average variables

  ! miscellaneous variables
  nGRU                 => summa1_struc%nGRU        , & ! number of grouped response units
  nHRU                 => summa1_struc%nHRU          & ! number of global hydrologic response units

 ) ! assignment to variables in the data structures
 ! ---------------------------------------------------------------------------------------
 ! initialize error control
 err=0; message='summa_manageOutputFiles/'

 ! identify the start of the writing
 call date_and_time(values=startWrite)

 ! *****************************************************************************
 ! *** initialize statistics
 ! *****************************************************************************

 ! initialize the statistics flags
 if(modelTimeStep==1)then

  ! initialize time step index
  statCounter(1:maxVarFreq) = 1
  outputTimeStep(1:maxVarFreq) = 1

  ! initialize flags to reset/finalize statistics
  resetStats(:)    = .true.   ! start by resetting statistics
  finalizeStats(:) = .false.  ! do not finalize stats on the first time step

  ! set stats flag for the timestep-level output
  finalizeStats(iLookFreq%timestep)=.true.

  ! initialize number of hru and gru in global data
  nGRUrun = nGRU
  nHRUrun = nHRU
 endif  ! if the first time step

 ! *****************************************************************************
 ! *** set alarms for writing data
 ! *****************************************************************************

 ! set alarms to control model output
 call summa_setWriteAlarms(oldTime%var, timeStruct%var, finshTime%var, &   ! time vectors
                           newOutputFile, defNewOutputFile,            &   ! flag to define new output file
                           ixRestart,     printRestart,                &   ! flag to print the restart file
                           ixProgress,    printProgress,               &   ! flag to print simulation progress
                           resetStats,    finalizeStats,               &   ! flags to reset and finalize stats
                           statCounter,                                &   ! statistics counter
                           err, cmessage)                                  ! error control
 if(err/=0)then; message=trim(message)//trim(cmessage); return; endif

 ! print progress
 if(printProgress) write(*,'(i4,1x,5(i2,1x))') timeStruct%var(1:5)

 ! *****************************************************************************
 ! *** define summa output files
 ! *****************************************************************************

 ! check the need to create a new output file
 if(defNewOutputFile .or. modelTimeStep==1)then

  ! define summa output files
  call summa_defineOutputFiles(modelTimeStep, summa1_struc, err, cmessage)
  if(err/=0)then; message=trim(message)//trim(cmessage); return; endif

  ! re-initalize the indices for model writing
  outputTimeStep(:)=1

 end if  ! if defining a new file

 ! ****************************************************************************
 ! *** calculate output statistics
 ! ****************************************************************************

 ! loop through GRUs and HRUs
 do iGRU=1,nGRU
  do iHRU=1,gru_struc(iGRU)%hruCount

   ! calculate output Statistics
   do iStruct=1,size(structInfo)
    select case(trim(structInfo(iStruct)%structName))
     case('forc'); call calcStats(forcStat%gru(iGRU)%hru(iHRU)%var,forcStruct%gru(iGRU)%hru(iHRU)%var,statForc_meta,resetStats,finalizeStats,statCounter,err,cmessage)
     case('prog'); call calcStats(progStat%gru(iGRU)%hru(iHRU)%var,progStruct%gru(iGRU)%hru(iHRU)%var,statProg_meta,resetStats,finalizeStats,statCounter,err,cmessage)
     case('diag'); call calcStats(diagStat%gru(iGRU)%hru(iHRU)%var,diagStruct%gru(iGRU)%hru(iHRU)%var,statDiag_meta,resetStats,finalizeStats,statCounter,err,cmessage)
     case('flux'); call calcStats(fluxStat%gru(iGRU)%hru(iHRU)%var,fluxStruct%gru(iGRU)%hru(iHRU)%var,statFlux_meta,resetStats,finalizeStats,statCounter,err,cmessage)
     case('indx'); call calcStats(indxStat%gru(iGRU)%hru(iHRU)%var,indxStruct%gru(iGRU)%hru(iHRU)%var,statIndx_meta,resetStats,finalizeStats,statCounter,err,cmessage)
    end select
    if(err/=0)then; message=trim(message)//trim(cmessage)//'['//trim(structInfo(iStruct)%structName)//']'; return; endif
   end do  ! (looping through structures)

  end do  ! (looping through HRUs)

  ! calc basin stats
  call calcStats(bvarStat%gru(iGRU)%var(:),bvarStruct%gru(iGRU)%var(:),statBvar_meta,resetStats,finalizeStats,statCounter,err,cmessage)
  if(err/=0)then; message=trim(message)//trim(cmessage)//'[bvar stats]'; return; endif

  ! write basin-average variables
  call writeBasin(iGRU,finalizeStats,outputTimeStep,bvar_meta,bvarStat%gru(iGRU)%var,bvarStruct%gru(iGRU)%var,bvarChild_map,err,cmessage)
  if(err/=0)then; message=trim(message)//trim(cmessage)//'[bvar]'; return; endif

 end do  ! (looping through GRUs)

 ! ****************************************************************************
 ! *** write data
 ! ****************************************************************************

 ! get the number of HRUs in the run domain
 nHRUrun = sum(gru_struc%hruCount)

 ! write time information
 call writeTime(finalizeStats,outputTimeStep,time_meta,timeStruct%var,err,message)

 ! write the model output to the NetCDF file
 ! Passes the full metadata structure rather than the stats metadata structure because
 !  we have the option to write out data of types other than statistics.
 !  Thus, we must also pass the stats parent->child maps from childStruct.
 do iStruct=1,size(structInfo)
  select case(trim(structInfo(iStruct)%structName))
   case('forc'); call writeData(finalizeStats,outputTimeStep,nHRUrun,maxLayers,forc_meta,forcStat,forcStruct,forcChild_map,indxStruct,err,cmessage)
   case('prog'); call writeData(finalizeStats,outputTimeStep,nHRUrun,maxLayers,prog_meta,progStat,progStruct,progChild_map,indxStruct,err,cmessage)
   case('diag'); call writeData(finalizeStats,outputTimeStep,nHRUrun,maxLayers,diag_meta,diagStat,diagStruct,diagChild_map,indxStruct,err,cmessage)
   case('flux'); call writeData(finalizeStats,outputTimeStep,nHRUrun,maxLayers,flux_meta,fluxStat,fluxStruct,fluxChild_map,indxStruct,err,cmessage)
   case('indx'); call writeData(finalizeStats,outputTimeStep,nHRUrun,maxLayers,indx_meta,indxStat,indxStruct,indxChild_map,indxStruct,err,cmessage)
  end select
  if(err/=0)then; message=trim(message)//trim(cmessage)//'['//trim(structInfo(iStruct)%structName)//']'; return; endif
 end do  ! (looping through structures)

 ! *****************************************************************************
 ! *** write restart file
 ! *****************************************************************************

 ! print a restart file if requested
 if(printRestart)then
  write(timeString,'(i4,3(i2.2))') timeStruct%var(iLookTIME%iyyy),timeStruct%var(iLookTIME%im),timeStruct%var(iLookTIME%id),timeStruct%var(iLookTIME%ih)
  
  if(STATE_PATH == '') then
    restartFile=trim(OUTPUT_PATH)//trim(OUTPUT_PREFIX)//'_restart_'//trim(timeString)//trim(output_fileSuffix)//'.nc'
  else
    restartFile=trim(STATE_PATH)//trim(OUTPUT_PREFIX)//'_restart_'//trim(timeString)//trim(output_fileSuffix)//'.nc'
  endif

  call writeRestart(restartFile,nGRU,nHRU,prog_meta,progStruct,bvar_meta,bvarStruct,maxLayers,maxSnowLayers,indx_meta,indxStruct,err,cmessage)  
  if(err/=0)then; message=trim(message)//trim(cmessage); return; endif
 end if

 ! *****************************************************************************
 ! *** update counters
 ! *****************************************************************************

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
 oldTime%var(:) = timeStruct%var(:)

 ! *****************************************************************************
 ! *** finalize
 ! *****************************************************************************

 ! identify the end of the writing
 call date_and_time(values=endWrite)

 ! aggregate the elapsed time for model writing
 elapsedWrite = elapsedWrite + elapsedSec(startWrite, endWrite)

 ! end associate statements
 end associate summaVars

 end subroutine summa_writeOutputFiles
end module summa_writeOutput



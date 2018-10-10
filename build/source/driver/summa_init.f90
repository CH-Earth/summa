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

module summa_init
! used to declare and allocate summa data structures and initialize model state to known values

! access missing values
USE globalData,only:integerMissing   ! missing integer
USE globalData,only:realMissing      ! missing double precision number

! named variables for run time options
USE globalData,only:iRunModeFull,iRunModeGRU,iRunModeHRU

! size of data structures
USE var_lookup,only:maxvarTime                              ! size of variable vectors
USE var_lookup,only:maxvarForc,maxvarProg,maxvarDiag        ! size of variable vectors
USE var_lookup,only:maxvarFlux,maxvarIndx,maxvarBvar        ! size of variable vectors

! metadata structures
USE globalData,only:time_meta,forc_meta,attr_meta,type_meta ! metadata structures
USE globalData,only:prog_meta,diag_meta,flux_meta           ! metadata structures
USE globalData,only:mpar_meta,indx_meta                     ! metadata structures
USE globalData,only:bpar_meta,bvar_meta                     ! metadata structures
USE globalData,only:averageFlux_meta                        ! metadata for time-step average fluxes

! statistics metadata structures
USE globalData,only:statForc_meta                           ! child metadata for stats
USE globalData,only:statProg_meta                           ! child metadata for stats
USE globalData,only:statDiag_meta                           ! child metadata for stats
USE globalData,only:statFlux_meta                           ! child metadata for stats
USE globalData,only:statIndx_meta                           ! child metadata for stats
USE globalData,only:statBvar_meta                           ! child metadata for stats

! mapping from original to child structures
USE globalData,only:forcChild_map                           ! index of the child data structure: stats forc
USE globalData,only:progChild_map                           ! index of the child data structure: stats prog
USE globalData,only:diagChild_map                           ! index of the child data structure: stats diag
USE globalData,only:fluxChild_map                           ! index of the child data structure: stats flux
USE globalData,only:indxChild_map                           ! index of the child data structure: stats indx
USE globalData,only:bvarChild_map                           ! index of the child data structure: stats bvar

! provide access to file paths
USE summaFileManager,only:SETNGS_PATH                       ! define path to settings files (e.g., Noah vegetation tables)
USE summaFileManager,only:MODEL_INITCOND                    ! name of model initial conditions file
USE summaFileManager,only:LOCAL_ATTRIBUTES                  ! name of model initial attributes file

! safety: set private unless specified otherwise
implicit none
private
public::summa_initialize
contains

 subroutine summa_initialize(summa1_struc, err, message)
 ! ---------------------------------------------------------------------------------------
 ! * desired modules
 ! ---------------------------------------------------------------------------------------
 ! data types
 USE nrtype                                                  ! variable types, etc.
 USE summa_type, only:summa1_type_dec                        ! master summa data type
 ! subroutines and functions: initial priming
 USE summa_util, only:getCommandArguments                    ! process command line arguments
 USE summaFileManager,only:summa_SetDirsUndPhiles            ! sets directories and filenames
 USE time_utils_module,only:elapsedSec                       ! calculate the elapsed time
 USE,intrinsic :: ieee_arithmetic                            ! IEEE arithmetic (obviously)
 ! subroutines and functions: define metadata structures
 USE popMetadat_module,only:popMetadat                       ! module to populate metadata structures
 USE flxMapping_module,only:flxMapping                       ! module to map fluxes to states
 USE checkStruc_module,only:checkStruc                       ! module to check metadata structures
 USE childStruc_module,only:childStruc                       ! module to create a child data structure
 ! subroutines and functions: read dimensions (NOTE: NetCDF)
 USE read_attrb_module,only:read_dimension                   ! module to read dimensions of GRU and HRU
 USE read_icond_module,only:read_icond_nlayers               ! module to read initial condition dimensions
 ! subroutines and functions: allocate space
 USE allocspace_module,only:allocGlobal                      ! module to allocate space for global data structures
 USE allocspace_module,only:allocLocal                       ! module to allocate space for local data structures
 ! timing variables
 USE globalData,only:startInit,endInit                       ! date/time for the start and end of the initialization
 USE globalData,only:elapsedInit                             ! elapsed time for the initialization
 ! model time structures
 USE globalData,only:refTime                                 ! reference time
 USE globalData,only:startTime                               ! start time
 USE globalData,only:finshTime                               ! end time
 ! run time options
 USE globalData,only:startGRU                                ! index of the starting GRU for parallelization run
 USE globalData,only:checkHRU                                ! index of the HRU for a single HRU run
 USE globalData,only:iRunMode                                ! define the current running mode
 ! miscellaneous global data
 USE globalData,only:dNaN                                    ! double precision NaN
 USE globalData,only:ncid                                    ! file id of netcdf output file
 USE globalData,only:gru_struc                               ! gru-hru mapping structures
 USE globalData,only:doJacobian                              ! flag to compute the Jacobian
 USE globalData,only:structInfo                              ! information on the data structures
 ! named variables that describe elements of child  model structures
 USE var_lookup,only:iLookVarType                            ! look-up values for variable type structure
 USE var_lookup,only:childFLUX_MEAN                          ! look-up values for timestep-average model fluxes
 ! ---------------------------------------------------------------------------------------
 ! * variables
 ! ---------------------------------------------------------------------------------------
 implicit none
 ! dummy variables
 type(summa1_type_dec),intent(inout)   :: summa1_struc       ! master summa data structure
 integer(i4b),intent(out)              :: err                ! error code
 character(*),intent(out)              :: message            ! error message
 ! local variables
 character(LEN=256)                    :: cmessage           ! error message of downwind routine
 logical(lgt), dimension(maxvarFlux)   :: flux_mask          ! mask defining desired flux variables
 logical(lgt), dimension(maxvarForc)   :: statForc_mask      ! mask defining forc stats
 logical(lgt), dimension(maxvarProg)   :: statProg_mask      ! mask defining prog stats
 logical(lgt), dimension(maxvarDiag)   :: statDiag_mask      ! mask defining diag stats
 logical(lgt), dimension(maxvarFlux)   :: statFlux_mask      ! mask defining flux stats
 logical(lgt), dimension(maxvarIndx)   :: statIndx_mask      ! mask defining indx stats
 logical(lgt), dimension(maxvarBvar)   :: statBvar_mask      ! mask defining bvar stats
 character(len=256)                    :: restartFile        ! restart file name
 character(len=256)                    :: attrFile           ! attributes file name
 integer(i4b)                          :: iStruct,iGRU       ! looping variables
 integer(i4b)                          :: fileGRU            ! [used for filenames] number of GRUs in the input file
 integer(i4b)                          :: fileHRU            ! [used for filenames] number of HRUs in the input file
 ! ---------------------------------------------------------------------------------------
 ! associate to elements in the data structure
 summaVars: associate(&

  ! statistics structures
  forcStat             => summa1_struc%forcStat            , & ! x%gru(:)%hru(:)%var(:)%dat -- model forcing data
  progStat             => summa1_struc%progStat            , & ! x%gru(:)%hru(:)%var(:)%dat -- model prognostic (state) variables
  diagStat             => summa1_struc%diagStat            , & ! x%gru(:)%hru(:)%var(:)%dat -- model diagnostic variables
  fluxStat             => summa1_struc%fluxStat            , & ! x%gru(:)%hru(:)%var(:)%dat -- model fluxes
  indxStat             => summa1_struc%indxStat            , & ! x%gru(:)%hru(:)%var(:)%dat -- model indices
  bvarStat             => summa1_struc%bvarStat            , & ! x%gru(:)%var(:)%dat        -- basin-average variabl

  ! primary data structures (scalars)
  timeStruct           => summa1_struc%timeStruct          , & ! x%var(:)                   -- model time data
  forcStruct           => summa1_struc%forcStruct          , & ! x%gru(:)%hru(:)%var(:)     -- model forcing data
  attrStruct           => summa1_struc%attrStruct          , & ! x%gru(:)%hru(:)%var(:)     -- local attributes for each HRU
  typeStruct           => summa1_struc%typeStruct          , & ! x%gru(:)%hru(:)%var(:)     -- local classification of soil veg etc. for each HRU

  ! primary data structures (variable length vectors)
  indxStruct           => summa1_struc%indxStruct          , & ! x%gru(:)%hru(:)%var(:)%dat -- model indices
  mparStruct           => summa1_struc%mparStruct          , & ! x%gru(:)%hru(:)%var(:)%dat -- model parameters
  progStruct           => summa1_struc%progStruct          , & ! x%gru(:)%hru(:)%var(:)%dat -- model prognostic (state) variables
  diagStruct           => summa1_struc%diagStruct          , & ! x%gru(:)%hru(:)%var(:)%dat -- model diagnostic variables
  fluxStruct           => summa1_struc%fluxStruct          , & ! x%gru(:)%hru(:)%var(:)%dat -- model fluxes

  ! basin-average structures
  bparStruct           => summa1_struc%bparStruct          , & ! x%gru(:)%var(:)            -- basin-average parameters
  bvarStruct           => summa1_struc%bvarStruct          , & ! x%gru(:)%var(:)%dat        -- basin-average variables

  ! ancillary data structures
  dparStruct           => summa1_struc%dparStruct          , & ! x%gru(:)%hru(:)%var(:)     -- default model parameters

  ! run time variables
  computeVegFlux       => summa1_struc%computeVegFlux      , & ! flag to indicate if we are computing fluxes over vegetation (.false. means veg is buried with snow)
  dt_init              => summa1_struc%dt_init             , & ! used to initialize the length of the sub-step for each HRU
  upArea               => summa1_struc%upArea              , & ! area upslope of each HRU

  ! miscellaneous variables
  summa1open           => summa1_struc%summa1open          , & ! flag to define if the summa file is open??
  numout               => summa1_struc%numout              , & ! number of output variables??
  ts                   => summa1_struc%ts                  , & ! model time step ??
  nGRU                 => summa1_struc%nGRU                , & ! number of grouped response units
  nHRU                 => summa1_struc%nHRU                , & ! number of global hydrologic response units
  hruCount             => summa1_struc%hruCount            , & ! number of local hydrologic response units
  greenVegFrac_monthly => summa1_struc%greenVegFrac_monthly, & ! fraction of green vegetation in each month (0-1)
  summaFileManagerFile => summa1_struc%summaFileManagerFile  & ! path/name of file defining directories and files

 ) ! assignment to variables in the data structures
 ! ---------------------------------------------------------------------------------------
 ! initialize error control
 err=0; message='summa_initialize/'

 ! initialize the start of the initialization
 call date_and_time(values=startInit)

 ! *****************************************************************************
 ! *** inital priming -- get command line arguments, identify files, etc.
 ! *****************************************************************************

 ! initialize the Jacobian flag
 doJacobian=.false.        ! initialize the Jacobian flag
 ncid(:) = integerMissing  ! initialize netcdf file id

 ! define double precision NaNs (shared in globalData)
 dNaN = ieee_value(1._dp, ieee_quiet_nan)

 ! get the command line arguments
 call getCommandArguments(summa1_struc,err,cmessage)
 if(err/=0)then; message=trim(message)//trim(cmessage); return; endif

 ! set directories and files -- summaFileManager used as command-line argument
 call summa_SetDirsUndPhiles(summaFileManagerFile,err,cmessage)
 if(err/=0)then; message=trim(message)//trim(cmessage); return; endif

 ! *****************************************************************************
 ! *** populate/check metadata structures
 ! *****************************************************************************

 ! populate metadata for all model variables
 call popMetadat(err,cmessage)
 if(err/=0)then; message=trim(message)//trim(cmessage); return; endif

 ! define mapping between fluxes and states
 call flxMapping(err,cmessage)
 if(err/=0)then; message=trim(message)//trim(cmessage); return; endif

 ! check data structures
 call checkStruc(err,cmessage)
 if(err/=0)then; message=trim(message)//trim(cmessage); return; endif

 ! define the mask to identify the subset of variables in the "child" data structure (just scalar variables)
 flux_mask = (flux_meta(:)%vartype==iLookVarType%scalarv)

 ! create the averageFlux metadata structure
 call childStruc(flux_meta, flux_mask, averageFlux_meta, childFLUX_MEAN, err, cmessage)
 if(err/=0)then; message=trim(message)//trim(cmessage); return; endif

 ! *****************************************************************************
 ! *** read the number of GRUs and HRUs
 ! *****************************************************************************
 ! obtain the HRU and GRU dimensions in the LocalAttribute file
 attrFile = trim(SETNGS_PATH)//trim(LOCAL_ATTRIBUTES)
 select case (iRunMode)
  case(iRunModeFull); call read_dimension(trim(attrFile),fileGRU,fileHRU,nGRU,nHRU,err,cmessage)
  case(iRunModeGRU ); call read_dimension(trim(attrFile),fileGRU,fileHRU,nGRU,nHRU,err,cmessage,startGRU=startGRU)
  case(iRunModeHRU ); call read_dimension(trim(attrFile),fileGRU,fileHRU,nGRU,nHRU,err,cmessage,checkHRU=checkHRU)
 end select
 if(err/=0)then; message=trim(message)//trim(cmessage); return; endif

 ! *****************************************************************************
 ! *** read the number of snow and soil layers
 ! *****************************************************************************
 ! obtain the number of snow and soil layers from the initial conditions file
 restartFile = trim(SETNGS_PATH)//trim(MODEL_INITCOND)
 call read_icond_nlayers(trim(restartFile),nGRU,indx_meta,err,cmessage)
 if(err/=0)then; message=trim(message)//trim(cmessage); return; endif

 ! *****************************************************************************
 ! *** allocate space for data structures
 ! *****************************************************************************

 ! allocate time structures
 do iStruct=1,3
  select case(iStruct)
   case(1); call allocLocal(time_meta, refTime,   err=err, message=cmessage)  ! reference time for the model simulation
   case(2); call allocLocal(time_meta, startTime, err=err, message=cmessage)  ! start time for the model simulation
   case(3); call allocLocal(time_meta, finshTime, err=err, message=cmessage)  ! end time for the model simulation
  end select
  if(err/=0)then; message=trim(message)//trim(cmessage); return; endif
 end do  ! looping through time structures

 ! allocate other data structures
 do iStruct=1,size(structInfo)
  ! allocate space
  select case(trim(structInfo(iStruct)%structName))
   case('time'); call allocGlobal(time_meta,  timeStruct,  err, cmessage)   ! model forcing data
   case('forc'); call allocGlobal(forc_meta,  forcStruct,  err, cmessage)   ! model forcing data
   case('attr'); call allocGlobal(attr_meta,  attrStruct,  err, cmessage)   ! local attributes for each HRU
   case('type'); call allocGlobal(type_meta,  typeStruct,  err, cmessage)   ! local classification of soil veg etc. for each HRU
   case('mpar'); call allocGlobal(mpar_meta,  mparStruct,  err, cmessage)   ! model parameters
   case('indx'); call allocGlobal(indx_meta,  indxStruct,  err, cmessage)   ! model variables
   case('prog'); call allocGlobal(prog_meta,  progStruct,  err, cmessage)   ! model prognostic (state) variables
   case('diag'); call allocGlobal(diag_meta,  diagStruct,  err, cmessage)   ! model diagnostic variables
   case('flux'); call allocGlobal(flux_meta,  fluxStruct,  err, cmessage)   ! model fluxes
   case('bpar'); call allocGlobal(bpar_meta,  bparStruct,  err, cmessage)   ! basin-average parameters
   case('bvar'); call allocGlobal(bvar_meta,  bvarStruct,  err, cmessage)   ! basin-average variables
   case('deriv'); cycle
   case default; err=20; message='unable to find structure name: '//trim(structInfo(iStruct)%structName)
  end select
  ! check errors
  if(err/=0)then
   message=trim(cmessage)//'[structure =  '//trim(structInfo(iStruct)%structName)//']'
   return
  endif
 end do  ! looping through data structures

 ! allocate space for default model parameters
 ! NOTE: This is done here, rather than in the loop above, because dpar is not one of the "standard" data structures
 call allocGlobal(mpar_meta,dparStruct,err,cmessage)   ! default model parameters
 if(err/=0)then
  message=trim(cmessage)//' [problem allocating dparStruct]'
  return
 endif

 ! allocate space for the time step and computeVegFlux flags (recycled for each GRU for subsequent model calls)
 allocate(dt_init%gru(nGRU),upArea%gru(nGRU),computeVegFlux%gru(nGRU),stat=err)
 if(err/=0)then
  message='problem allocating space for dt_init, upArea, or computeVegFlux [GRU]'
  return
 endif

 ! allocate space for the HRUs
 do iGRU=1,nGRU
  hruCount = gru_struc(iGRU)%hruCount  ! gru_struc populated in "read_dimension"
  allocate(dt_init%gru(iGRU)%hru(hruCount),upArea%gru(iGRU)%hru(hruCount),computeVegFlux%gru(iGRU)%hru(hruCount),stat=err)
  if(err/=0)then
   message='problem allocating space for dt_init, upArea, or computeVegFlux [HRU]'
   return
  endif
 end do

 ! *****************************************************************************
 ! *** allocate space for output statistics data structures
 ! *****************************************************************************

 ! child metadata structures - so that we do not carry full stats structures around everywhere
 ! only carry stats for variables with output frequency > model time step
 statForc_mask = (forc_meta(:)%vartype==iLookVarType%scalarv.and.forc_meta(:)%varDesire)
 statProg_mask = (prog_meta(:)%vartype==iLookVarType%scalarv.and.prog_meta(:)%varDesire)
 statDiag_mask = (diag_meta(:)%vartype==iLookVarType%scalarv.and.diag_meta(:)%varDesire)
 statFlux_mask = (flux_meta(:)%vartype==iLookVarType%scalarv.and.flux_meta(:)%varDesire)
 statIndx_mask = (indx_meta(:)%vartype==iLookVarType%scalarv.and.indx_meta(:)%varDesire)
 statBvar_mask = (bvar_meta(:)%vartype==iLookVarType%scalarv.and.bvar_meta(:)%varDesire)

 ! create the stats metadata structures
 do iStruct=1,size(structInfo)
  select case (trim(structInfo(iStruct)%structName))
   case('forc'); call childStruc(forc_meta,statForc_mask,statForc_meta,forcChild_map,err,cmessage)
   case('prog'); call childStruc(prog_meta,statProg_mask,statProg_meta,progChild_map,err,cmessage)
   case('diag'); call childStruc(diag_meta,statDiag_mask,statDiag_meta,diagChild_map,err,cmessage)
   case('flux'); call childStruc(flux_meta,statFlux_mask,statFlux_meta,fluxChild_map,err,cmessage)
   case('indx'); call childStruc(indx_meta,statIndx_mask,statIndx_meta,indxChild_map,err,cmessage)
   case('bvar'); call childStruc(bvar_meta,statBvar_mask,statBvar_meta,bvarChild_map,err,cmessage)
  end select
  ! check errors
  if(err/=0)then
   message=trim(cmessage)//'[statistics for =  '//trim(structInfo(iStruct)%structName)//']'
   return
  endif
 end do ! iStruct

 ! set all stats metadata to correct var types
 statForc_meta(:)%vartype = iLookVarType%outstat
 statProg_meta(:)%vartype = iLookVarType%outstat
 statDiag_meta(:)%vartype = iLookVarType%outstat
 statFlux_meta(:)%vartype = iLookVarType%outstat
 statIndx_meta(:)%vartype = iLookVarType%outstat
 statBvar_meta(:)%vartype = iLookVarType%outstat

 ! loop through data structures
 do iStruct=1,size(structInfo)

  ! allocate space
  select case(trim(structInfo(iStruct)%structName))
   case('forc'); call allocGlobal(statForc_meta(:)%var_info,forcStat,err,cmessage)   ! model forcing data
   case('prog'); call allocGlobal(statProg_meta(:)%var_info,progStat,err,cmessage)   ! model prognostic (state) variables
   case('diag'); call allocGlobal(statDiag_meta(:)%var_info,diagStat,err,cmessage)   ! model diagnostic variables
   case('flux'); call allocGlobal(statFlux_meta(:)%var_info,fluxStat,err,cmessage)   ! model fluxes
   case('indx'); call allocGlobal(statIndx_meta(:)%var_info,indxStat,err,cmessage)   ! index vars
   case('bvar'); call allocGlobal(statBvar_meta(:)%var_info,bvarStat,err,cmessage)   ! basin-average variables
   case default; cycle
  end select

  ! check errors
  if(err/=0)then
   message=trim(cmessage)//'[statistics for =  '//trim(structInfo(iStruct)%structName)//']'
   return
  endif

 end do ! iStruct

 ! identify the end of the initialization
 call date_and_time(values=endInit)

 ! aggregate the elapsed time for the initialization
 elapsedInit = elapsedSec(startInit, endInit)

 ! end associate statements
 end associate summaVars

 end subroutine summa_initialize
end module summa_init

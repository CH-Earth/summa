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

module summa_init
! used to declare and allocate summa data structures and initialize model state to known values

! check if mizuroute is active
use build_options, only: mizuroute_active
use build_options, only: ngen_forcing_active

#ifdef MIZUROUTE_ACTIVE
USE mizuroute_coupling, only: init_mizuroute_from_summa
#endif

! access missing values
USE globalData,only:integerMissing   ! missing integer
USE globalData,only:realMissing      ! missing real number

! global data 
USE globalData,only:data_step        ! length of the data step (s)
USE globalData, only: iulog          ! I/O unit for logging messages 

! output constraints
USE globalData,only:maxLayers        ! maximum number of layers
USE globalData,only:maxSoilLayers    ! maximum number of soil layers
USE globalData,only:maxSnowLayers    ! maximum number of snow layers

! named variables for run time options
USE globalData,only:iRunModeFull,iRunModeGRU,iRunModeHRU

! metadata structures
USE globalData,only:time_meta,forc_meta,attr_meta,type_meta ! metadata structures
USE globalData,only:prog_meta,diag_meta,flux_meta,id_meta   ! metadata structures
USE globalData,only:mpar_meta,indx_meta                     ! metadata structures
USE globalData,only:bpar_meta,bvar_meta                     ! metadata structures
USE globalData,only:lookup_meta

! statistics metadata structures
USE globalData,only:statForc_meta                           ! child metadata for stats
USE globalData,only:statProg_meta                           ! child metadata for stats
USE globalData,only:statDiag_meta                           ! child metadata for stats
USE globalData,only:statFlux_meta                           ! child metadata for stats
USE globalData,only:statIndx_meta                           ! child metadata for stats
USE globalData,only:statBvar_meta                           ! child metadata for stats

! provide access to file paths
USE summaFileManager,only:SETTINGS_PATH                     ! define path to settings files (e.g., parameters, soil and veg. tables)
USE summaFileManager,only:STATE_PATH                        ! optional path to state/init. condition files (defaults to SETTINGS_PATH)
USE summaFileManager,only:MODEL_INITCOND                    ! name of model initial conditions file
USE summaFileManager,only:LOCAL_ATTRIBUTES                  ! name of model initial attributes file

! model decisions
USE globalData,only:model_decisions                         ! model decision structure
USE var_lookup,only:iLookDECISIONS                          ! look-up values for model decisions

! named variables to define the decisions for snow layers
USE mDecisions_module,only:&
  sameRulesAllLayers,&                  ! SNTHERM option: same combination/sub-dividion rules applied to all layers
  rulesDependLayerIndex                 ! CLM option: combination/sub-dividion rules depend on layer index

! named variables for the output buffer
USE mDecisions_module,only:&
 writePerStep,   &                      ! write data per time step (default)
 writeFullSeries                        ! write all data for a given output file

! safety: set private unless specified otherwise
implicit none
private
public::summa_initialize
contains

! used to declare and allocate summa data structures and initialize model state to known values
subroutine summa_initialize(summa1_struc, err, message)
  ! ---------------------------------------------------------------------------------------
  ! * desired modules
  ! ---------------------------------------------------------------------------------------
  ! data types
  USE nr_type                                                 ! variable types, etc.
  USE summa_type, only:summa1_type_dec                        ! master summa data type
  ! subroutines and functions: initial priming
  USE summa_util, only:getCommandArguments                    ! process command line arguments
  USE summaFileManager,only:summa_SetTimesDirsAndFiles        ! sets directories and filenames
  USE summa_globalData,only:summa_defineGlobalData            ! used to define global summa data structures
  USE summa_config,only:load_summa_config                     ! load TOML configuration settings (for parsing later)
  USE time_utils_module,only:elapsedSec                       ! calculate the elapsed time
  ! subroutines and functions: parallelization
  USE summa_work_balance,only:balance_even                    ! module to identify start/end indices for a given rank
  ! subroutines and functions: read dimensions (NOTE: NetCDF)
  USE read_attrb_module,only:read_dimension                   ! module to read dimensions of GRU and HRU
  USE read_attrb_module,only:read_mapping_vectors             ! module to define mapping between GRU and HRU
  USE read_icond_module,only:read_icond_nlayers               ! module to read initial condition dimensions
  ! subroutines and functions: allocate space
  USE allocspace_module,only:alloc_driver_work                ! module to allocate space for work structures
  USE allocspace_module,only:allocGlobal                      ! module to allocate space for global data structures
  USE allocspace_module,only:allocLocal                       ! module to allocate space for local data structures
  ! subroutines and functions: model decisions and forcing file information
  USE mDecisions_module,only:mDecisions                       ! module to read model decisions
  USE ffile_info_module,only:ffile_info                       ! module to read information on forcing datafile
  ! timing variables
  USE globalData,only:startInit,endInit                       ! date/time for the start and end of the initialization
  USE globalData,only:elapsedInit                             ! elapsed time for the initialization
  USE globalData,only:elapsedRead                             ! elapsed time for the data read
  USE globalData,only:elapsedWrite                            ! elapsed time for the stats/write
  USE globalData,only:elapsedPhysics                          ! elapsed time for the physics
  ! model time structures
  USE globalData,only:startTime                               ! start time
  USE globalData,only:finshTime                               ! end time
  USE globalData,only:refTime                                 ! reference time
  USE globalData,only:oldTime                                 ! time from previous step
  ! buffered write
  USE globalData,only:numtim                                  ! number of time steps
  ! run time options
  USE globalData,only:startGRU_user => startGRU               ! index of the starting GRU defined using the -g runtime option
  USE globalData,only:checkHRU                                ! index of the HRU for a single HRU run using the -h runtime option
  USE globalData,only:iRunMode                                ! define the current running mode
  ! miscellaneous global data
  USE globalData,only:ncid                                    ! file id of netcdf output file
  USE globalData,only:gru_struc                               ! gru-hru mapping structures (constructed in read_mapping_vectors)
  USE globalData,only:structInfo                              ! information on the data structures
  USE globalData,only:output_fileSuffix                       ! suffix for the output file
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
  character(len=256)                    :: restartFile        ! restart file name
  character(len=256)                    :: attrFile           ! attributes file name
  character(len=128)                    :: fmtGruOutput       ! a format string used to write start and end GRU in output file names
  integer(i4b)                          :: iStruct,iGRU,iHRU  ! looping variables
  integer(i4b)                          :: nGRU_file          ! number of GRUs in the complete input file
  integer(i4b)                          :: nHRU_file          ! number of HRUs in the complete input file
  integer(i4b)                          :: startGRU_local     ! file index of first GRU assigned to the current rank
  integer(i4b)                          :: startGRU_domain    ! file index of first GRU in the run domain
  integer(i4b)                          :: nGRU_domain        ! number of GRUs in the run domain
  ! ---------------------------------------------------------------------------------------
  ! associate to elements in the data structure
  summaVars: associate(&
    
    ! parallel execution context
    parallel             => summa1_struc%parallel            , & ! x%comm, x%rank, x%size -- parallel execution context

    ! lookup table data structure
    lookupStruct         => summa1_struc%lookupStruct        , & ! x%gru(:)%hru(:)%z(:)%var(:)%lookup(:) -- lookup tables
    
    ! statistics structures
    forcStat             => summa1_struc%forcStat            , & ! x%gru(:)%hru(:)%var(:)%dat -- model forcing data
    progStat             => summa1_struc%progStat            , & ! x%gru(:)%hru(:)%var(:)%dat -- model prognostic (state) variables
    diagStat             => summa1_struc%diagStat            , & ! x%gru(:)%hru(:)%var(:)%dat -- model diagnostic variables
    fluxStat             => summa1_struc%fluxStat            , & ! x%gru(:)%hru(:)%var(:)%dat -- model fluxes
    indxStat             => summa1_struc%indxStat            , & ! x%gru(:)%hru(:)%var(:)%dat -- model indices
    bvarStat             => summa1_struc%bvarStat            , & ! x%gru(:)%var(:)%dat        -- basin-average variables

    ! primary data structures (scalars)
    timeStruct           => summa1_struc%timeStruct          , & ! x%var(:)                   -- model time data
    forcStruct           => summa1_struc%forcStruct          , & ! x%gru(:)%hru(:)%var(:)     -- model forcing data
    attrStruct           => summa1_struc%attrStruct          , & ! x%gru(:)%hru(:)%var(:)     -- local attributes for each HRU
    typeStruct           => summa1_struc%typeStruct          , & ! x%gru(:)%hru(:)%var(:)     -- local classification of soil veg etc. for each HRU
    idStruct             => summa1_struc%idStruct            , & ! x%gru(:)%hru(:)%var(:)     --

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

    ! GRU and HRU dimensions
    nGRU_user            => summa1_struc%nGRU_user           , & ! number of GRUs assigned to the current rank 
    nGRU_local           => summa1_struc%nGRU_local          , & ! number of GRUs assigned to the current rank 
    nHRU_local           => summa1_struc%nHRU_local          , & ! number of HRUs assigned to the current rank

    ! manager file
    summaFileManagerFile => summa1_struc%summaFileManagerFile, & ! path/name of file defining directories and files
    summaConfigFile      => summa1_struc%summaConfigFile       & ! path/name of summa configuration file

    ) ! assignment to variables in the data structures
    ! ---------------------------------------------------------------------------------------
    ! initialize error control
    err=0; message='summa_initialize/'

    ! initialize the start of the initialization
    call date_and_time(values=startInit)

    ! *****************************************************************************
    ! *** inital priming -- get command line arguments, identify files, etc.
    ! *****************************************************************************

    ! initialize the netcdf file id
    ncid(:) = integerMissing

    ! initialize the elapsed time for cumulative quantities
    elapsedRead=0._rkind
    elapsedWrite=0._rkind
    elapsedPhysics=0._rkind

    ! get the command line arguments
    call getCommandArguments(summa1_struc,err,cmessage)
    if(err/=0)then; message=trim(message)//trim(cmessage); return; endif

    ! load configuration settings from TOML file
    call load_summa_config(trim(summaConfigFile), summa1_struc, err, cmessage)
    if(err/=0)then; message=trim(message)//trim(cmessage); return; endif

    ! set directories and files -- summaFileManager used as command-line argument
    call summa_SetTimesDirsAndFiles(summaFileManagerFile,err,cmessage)
    if(err/=0)then; message=trim(message)//trim(cmessage); return; endif

    ! define global data (parameters, metadata)
    call summa_defineGlobalData(err, cmessage)
    if(err/=0)then; message=trim(message)//trim(cmessage); return; endif

    ! *****************************************************************************
    ! *** define spatial indexing and the run domain
    ! *****************************************************************************

    ! Spatial indexing uses three distinct reference domains:
    !
    !   file   : the complete spatial domain in the input files
    !            nGRU_file, nHRU_file
    !
    !   domain : the subset of GRUs selected for the model run
    !            startGRU_domain, nGRU_domain
    !
    !   local  : the subset of the run domain assigned to the current rank
    !            startGRU_local, nGRU_local, nHRU_local
    !
    ! startGRU_domain and startGRU_local are indices into the input file.

    ! Read the number of GRUs and HRUs in the LocalAttributes file and identify
    ! the GRU range selected for the model run. startGRU_domain is the file index
    ! of the first selected GRU and nGRU_domain is the number of selected GRUs.
    ! No parallel decomposition is performed here.

    ! obtain the HRU and GRU dimensions in the LocalAttributes file
    attrFile = trim(SETTINGS_PATH)//trim(LOCAL_ATTRIBUTES)
    
    select case (iRunMode)
    
      case (iRunModeFull)
        call read_dimension(trim(attrFile), nGRU_file, nHRU_file, &
                            startGRU_domain, nGRU_domain, err, cmessage)
    
      case (iRunModeGRU)

        nGRU_domain = nGRU_user

        call read_dimension(trim(attrFile), nGRU_file, nHRU_file, &
                            startGRU_domain, nGRU_domain, err, cmessage, &
                            startGRU_user=startGRU_user)
    
      case (iRunModeHRU)
        call read_dimension(trim(attrFile), nGRU_file, nHRU_file, &
                            startGRU_domain, nGRU_domain, err, cmessage, &
                            checkHRU=checkHRU)
    
    end select

    if(err/=0)then; message=trim(message)//trim(cmessage); return; endif

    ! *****************************************************************************
    ! *** define the local spatial domain
    ! *****************************************************************************

    ! Partition the GRUs in the run domain evenly across the available ranks.
    ! startGRU_local is the file index of the first GRU assigned to this rank
    ! and nGRU_local is the number of GRUs assigned to this rank. For a serial
    ! run, the local GRU range is identical to the run domain.

    write(iulog,*) 'Parallel context: rank =', parallel%rank, ' size =', parallel%size, &
                                    ' comm =', parallel%comm

    ! define start and count indices for each local rank
    if(iRunMode /= iRunModeHRU)then
    
      call balance_even(startGRU_domain, nGRU_domain, &
                        parallel%rank, parallel%size, &
                        startGRU_local, nGRU_local, &
                        err, cmessage)
      if(err/=0)then; message=trim(message)//trim(cmessage); return; endif

    else
      startGRU_local = integerMissing  ! assigned in read_mapping_vectors
      nGRU_local     = nGRU_domain     ! =1
    endif

    write(iulog,*) 'Run mode =', iRunMode, iRunModeFull
    write(iulog,*) 'File dimensions:  nGRU_file =', nGRU_file, '  nHRU_file =', nHRU_file
    write(iulog,*) 'Run domain:  startGRU =', startGRU_domain, ' nGRU =', nGRU_domain
    write(iulog,*) 'Local rank:  startGRU =', startGRU_local,  ' nGRU =', nGRU_local

    ! *****************************************************************************
    ! *** construct the local GRU-HRU mapping
    ! *****************************************************************************

    ! Read the GRU and HRU identifiers needed for this rank and construct the
    ! local GRU-HRU and HRU-GRU mapping structures. nHRU_local is determined
    ! from the HRUs belonging to the GRUs assigned to this rank.

    call read_mapping_vectors(attrFile,                                                 &
                              nGRU_file, nHRU_file,                                     &
                              startGRU_local, nGRU_local, nHRU_local,                   &
                              merge(checkHRU, integerMissing, iRunMode == iRunModeHRU), &
                              summa1_struc%gru_struc, summa1_struc%index_map,           & 
                              err, cmessage)
    if(err/=0)then; message=trim(message)//trim(cmessage); return; endif

    ! *****************************************************************************
    ! *** read the number of snow and soil layers
    ! *****************************************************************************
    ! set restart filename and read the number of snow and soil layers from the initial conditions (restart) file
    if(STATE_PATH == '') then
      restartFile = trim(SETTINGS_PATH)//trim(MODEL_INITCOND)
    else
      restartFile = trim(STATE_PATH)//trim(MODEL_INITCOND)
    endif
    call read_icond_nlayers(trim(restartFile),nGRU_local,indx_meta,err,cmessage)
    if(err/=0)then; message=trim(message)//trim(cmessage); return; endif

    ! *****************************************************************************
    ! *** allocate space for data structures
    ! *****************************************************************************

    ! Allocate the non-spatial time structures and the spatial model structures.

    ! Spatial structures are allocated for the GRUs and HRUs assigned to this rank,
    ! as defined by the local gru_struc mapping. For a serial run, the local spatial
    ! domain is identical to the complete run domain.

    ! allocate time structures
    do iStruct=1,4
      select case(iStruct)
        case(1); call allocLocal(time_meta, startTime, err=err, message=cmessage)  ! start time for the model simulation
        case(2); call allocLocal(time_meta, finshTime, err=err, message=cmessage)  ! end time for the model simulation
        case(3); call allocLocal(time_meta, refTime,   err=err, message=cmessage)  ! reference time for the model simulation
        case(4); call allocLocal(time_meta, oldTime,   err=err, message=cmessage)  ! time from the previous step
      end select
      if(err/=0)then; message=trim(message)//trim(cmessage); return; endif
    end do  ! looping through time structures

    ! allocate other data structures
    do iStruct=1,size(structInfo)
      ! allocate space
      select case(trim(structInfo(iStruct)%structName))
        case('time'  ); call allocGlobal(time_meta,    timeStruct,    err, cmessage)   ! model time data
        case('forc'  ); call allocGlobal(forc_meta,    forcStruct,    err, cmessage)   ! model forcing data
        case('attr'  ); call allocGlobal(attr_meta,    attrStruct,    err, cmessage)   ! local attributes for each HRU
        case('type'  ); call allocGlobal(type_meta,    typeStruct,    err, cmessage)   ! local classification of soil veg etc. for each HRU
        case('id'    ); call allocGlobal(id_meta,      idStruct,      err, cmessage)   ! local values of hru and gru IDs
        case('mpar'  ); call allocGlobal(mpar_meta,    mparStruct,    err, cmessage)   ! model parameters
        case('indx'  ); call allocGlobal(indx_meta,    indxStruct,    err, cmessage)   ! model variables
        case('prog'  ); call allocGlobal(prog_meta,    progStruct,    err, cmessage)   ! model prognostic (state) variables
        case('diag'  ); call allocGlobal(diag_meta,    diagStruct,    err, cmessage)   ! model diagnostic variables
        case('flux'  ); call allocGlobal(flux_meta,    fluxStruct,    err, cmessage)   ! model fluxes
        case('bpar'  ); call allocGlobal(bpar_meta,    bparStruct,    err, cmessage)   ! basin-average parameters
        case('bvar'  ); call allocGlobal(bvar_meta,    bvarStruct,    err, cmessage)   ! basin-average variables
        case('lookup'); call allocGlobal(lookup_meta,  lookupStruct,  err, cmessage)   ! lookup tables
        case('deriv' ); cycle ! derivatives are not stored in the data structure, but are instead computed on the fly and stored in local variables
        case default; err=20; message='unable to find structure name: '//trim(structInfo(iStruct)%structName)
      end select
      ! check errors
      if(err/=0)then
        message=trim(message)//trim(cmessage)//'[structure =  '//trim(structInfo(iStruct)%structName)//']'
        return
      endif
    end do  ! looping through data structures

    ! allocate space for default model parameters
    ! NOTE: This is done here, rather than in the loop above, because dpar is not one of the "standard" data structures
    call allocGlobal(mpar_meta,dparStruct,err,cmessage)   ! default model parameters
    if(err/=0)then
      message=trim(message)//trim(cmessage)//' [problem allocating dparStruct]'
      return
    endif

    ! allocate driver work structures
    call alloc_driver_work(nGRU_local, dt_init, upArea, computeVegFlux, err, cmessage)
    if(err/=0)then; message=trim(message)//trim(cmessage); return; endif

    ! *****************************************************************************
    ! *** allocate space for output statistics structures assigned to this rank
    ! *****************************************************************************

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
        message=trim(message)//trim(cmessage)//'[statistics for =  '//trim(structInfo(iStruct)%structName)//']'
        return
      endif

    end do ! iStruct

    ! *****************************************************************************
    ! if using NGEN forcing only need to set the hourly data_step (fixed)
    ! *****************************************************************************
    if (ngen_forcing_active) then
      data_step = 3600._rkind
    
    ! *****************************************************************************
    ! *** read description of model forcing datafile used in each HRU
    ! *****************************************************************************
    else
      call ffile_info(nGRU_local,err,cmessage)
      if(err/=0)then; message=trim(message)//trim(cmessage); return; endif
    endif
    
    ! *****************************************************************************
    ! *** read model decisions
    ! *****************************************************************************
    ! NOTE: Must be after ffile_info because mDecisions uses the data_step
    call mDecisions(err,cmessage)
    if(err/=0)then; message=trim(message)//trim(cmessage); return; endif
    
    ! get the maximum number of snow layers
    select case(model_decisions(iLookDECISIONS%snowLayers)%iDecision)
     case(sameRulesAllLayers);    maxSnowLayers = 100
     case(rulesDependLayerIndex); maxSnowLayers = 5
     case default; err=20; message=trim(message)//'unable to identify option to combine/sub-divide snow layers'; return
    end select ! (option to combine/sub-divide snow layers)
    
    ! get the maximum number of layers
    maxLayers     = 0
    maxSoilLayers = 0
    do iGRU=1,nGRU_local
     do iHRU=1,gru_struc(iGRU)%hruCount
      maxSoilLayers = max(maxSoilLayers, gru_struc(iGRU)%hruInfo(iHRU)%nSoil)
      maxLayers = max(maxLayers, maxSnowLayers+gru_struc(iGRU)%hruInfo(iHRU)%nSoil)
     end do
    end do
   
    ! get the number of time steps in the output buffer
    select case(model_decisions(iLookDECISIONS%write_buff)%iDecision)
     case (writePerStep);    summa1_struc%n_write = 1
     case (writeFullSeries); summa1_struc%n_write = numtim
     case default; err=20; message=trim(message)//'unable to identify option for output buffer'; return
    end select

    ! save the length of the data window
    summa1_struc%data_step = data_step

    ! *****************************************************************************
    ! *** initialize mizuRoute (if mizuRoute is active)
    ! *****************************************************************************
    
    if (mizuroute_active) then

      ! allocate data structure for mizuroute coupling
      allocate(summa1_struc%coupling(nGRU_local), stat=err)
      if(err/=0)then
        message=trim(message)//' [problem allocating mizuroute coupling structure]'
        return
      endif

      ! populate mizuroute coupling IDs
      summa1_struc%coupling(:)%id = summa1_struc%gru_struc(:)%gru_id

      call init_mizuroute_from_summa(summa1_struc, err, cmessage) 
      if(err/=0)then; message=trim(message)//trim(cmessage); return; endif

    endif

    ! *****************************************************************************
    ! *** define the suffix for the model output file
    ! *****************************************************************************
    
    if (output_fileSuffix(1:1) /= '_') output_fileSuffix='_'//trim(output_fileSuffix)   ! separate output_fileSuffix from others by underscores
    if (output_fileSuffix(len_trim(output_fileSuffix):len_trim(output_fileSuffix)) == '_') output_fileSuffix(len_trim(output_fileSuffix):len_trim(output_fileSuffix)) = ' '
    
    select case (iRunMode)
    
      case (iRunModeGRU, iRunModeFull)

        ! add GRU range text string for GRU-subset runs or parallel full-domain runs
        if (iRunMode == iRunModeGRU .or. parallel%size > 1)then

          ! left zero padding for startGRU and endGRU
          write(fmtGruOutput,"(i0)") ceiling(log10(real(nGRU_file)+0.1))                    ! maximum width of startGRU and endGRU
          fmtGruOutput = "i"//trim(fmtGruOutput)//"."//trim(fmtGruOutput)                   ! construct the format string for startGRU and endGRU
          fmtGruOutput = "('_G',"//trim(fmtGruOutput)//",'-',"//trim(fmtGruOutput)//")"
          
          write(output_fileSuffix((len_trim(output_fileSuffix)+1):len(output_fileSuffix)),fmtGruOutput) &
                                   startGRU_local, startGRU_local+nGRU_local-1
  
        endif

      case (iRunModeHRU)
        write(output_fileSuffix((len_trim(output_fileSuffix)+1):len(output_fileSuffix)),"('_H',i0)") checkHRU

    end select

    ! identify the end of the initialization
    call date_and_time(values=endInit)

    ! aggregate the elapsed time for the initialization
    elapsedInit = elapsedSec(startInit, endInit)

  ! end associate statements
  end associate summaVars

  !stop 'end of summa_initialize'

end subroutine summa_initialize

end module summa_init

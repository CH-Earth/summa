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

program multi_driver
! used to evaluate different methods for simulating snow processes
! *****************************************************************************
! use desired modules
! *****************************************************************************
USE nrtype                                                  ! variable types, etc.
! provide access to subroutines and functions
USE summaFileManager,only:summa_SetDirsUndPhiles            ! sets directories and filenames
USE module_sf_noahmplsm,only:read_mp_veg_parameters         ! module to read NOAH vegetation tables
USE module_sf_noahmplsm,only:redprm                         ! module to assign more Noah-MP parameters
USE module_sf_noahmplsm,only:isWater                        ! parameter for water land cover type
USE nr_utility_module,only:arth                             ! get a sequence of numbers
USE ascii_util_module,only:file_open                        ! open ascii file
USE ascii_util_module,only:get_vlines                       ! read a vector of non-comment lines from an ASCII file
USE ascii_util_module,only:split_line                       ! extract the list of variable names from the character string
use time_utils_module,only:elapsedSec                       ! calculate the elapsed time
USE allocspace_module,only:allocate_gru_struc               ! read/allocate space for gru-hru mask structures
USE allocspace_module,only:allocGlobal                      ! module to allocate space for global data structures
USE allocspace_module,only:allocLocal                       ! module to allocate space for local data structures
USE childStruc_module,only:childStruc                       ! module to create a child data structure
USE mDecisions_module,only:mDecisions                       ! module to read model decisions
USE popMetadat_module,only:popMetadat                       ! module to populate metadata structures
USE checkStruc_module,only:checkStruc                       ! module to check metadata structures
USE def_output_module,only:def_output                       ! module to define model output
USE ffile_info_module,only:ffile_info                       ! module to read information on forcing datafile
USE read_attrb_module,only:read_attrb                       ! module to read local attributes
USE read_pinit_module,only:read_pinit                       ! module to read initial model parameter values
USE paramCheck_module,only:paramCheck                       ! module to check consistency of model parameters
USE read_icond_module,only:read_icond                       ! module to read initial conditions
USE pOverwrite_module,only:pOverwrite                       ! module to overwrite default parameter values with info from the Noah tables
USE read_param_module,only:read_param                       ! module to read model parameter sets
USE ConvE2Temp_module,only:E2T_lookup                       ! module to calculate a look-up table for the temperature-enthalpy conversion
USE var_derive_module,only:calcHeight                       ! module to calculate height at layer interfaces and layer mid-point
USE var_derive_module,only:v_shortcut                       ! module to calculate "short-cut" variables
USE var_derive_module,only:rootDensty                       ! module to calculate the vertical distribution of roots
USE var_derive_module,only:satHydCond                       ! module to calculate the saturated hydraulic conductivity in each soil layer
USE var_derive_module,only:fracFuture                       ! module to calculate the fraction of runoff in future time steps (time delay histogram)
USE read_force_module,only:read_force                       ! module to read model forcing data
USE derivforce_module,only:derivforce                       ! module to compute derived forcing data
USE modelwrite_module,only:writeParm,writeTime              ! module to write model attributes and parameters
USE modelwrite_module,only:writeData,writeBasin             ! module to write model output
USE vegPhenlgy_module,only:vegPhenlgy                       ! module to compute vegetation phenology
USE coupled_em_module,only:coupled_em                       ! module to run the coupled energy and mass model
USE groundwatr_module,only:groundwatr                       ! module to simulate regional groundwater balance
USE qTimeDelay_module,only:qOverland                        ! module to route water through an "unresolved" river network
USE netcdf_util_module,only:nc_file_close                   ! module to handle netcdf stuff for inputs and outputs
! provide access to file paths
USE summaFileManager,only:SETNGS_PATH                       ! define path to settings files (e.g., Noah vegetation tables)
USE summaFileManager,only:MODEL_INITCOND                    ! name of model initial conditions file
USE summaFileManager,only:OUTPUT_PATH,OUTPUT_PREFIX         ! define output file
USE summaFileManager,only:LOCALPARAM_INFO,BASINPARAM_INFO   ! files defining the default values and constraints for model parameters
! provide access to the derived types to define the data structures
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
USE data_types,only:extended_info          ! extended metadata structure
! provide access to metadata structures
USE globalData,only:time_meta,forc_meta,attr_meta,type_meta ! metadata structures
USE globalData,only:prog_meta,diag_meta,flux_meta           ! metadata structures
USE globalData,only:mpar_meta,indx_meta                     ! metadata structures
USE globalData,only:bpar_meta,bvar_meta                     ! metadata structures
USE globalData,only:averageFlux_meta                        ! metadata for time-step average fluxes
USE globalData,only:model_decisions                         ! model decision structure
! provide access to global data
USE globalData,only:refTime                                 ! reference time
USE globalData,only:startTime                               ! start time
USE globalData,only:finshTime                               ! end time
USE globalData,only:doJacobian                              ! flag to compute the Jacobian
USE globalData,only:gru_struc                               ! gru-hru mapping structures
USE globalData,only:index_map                               ! hru-hru mapping structures
USE globalData,only:localParFallback                        ! local column default parameters
USE globalData,only:basinParFallback                        ! basin-average default parameters
USE globalData,only:structInfo                              ! information on the data structures                  
USE globalData,only:numtim                                  ! number of time steps
USE globalData,only:urbanVegCategory                        ! vegetation category for urban areas
USE globalData,only:globalPrintFlag                         ! global print flag
USE multiconst,only:integerMissing                          ! missing integer value
! provide access to Noah-MP parameters
USE NOAHMP_VEG_PARAMETERS,only:SAIM,LAIM                    ! 2-d tables for stem area index and leaf area index (vegType,month)
USE NOAHMP_VEG_PARAMETERS,only:HVT,HVB                      ! height at the top and bottom of vegetation (vegType)
USE noahmp_globals,only:RSMIN                               ! minimum stomatal resistance (vegType)
! provide access to the named variables that describe elements of parent model structures
USE var_lookup,only:iLookTIME,iLookFORCE                    ! look-up values for time and forcing data structures
USE var_lookup,only:iLookTYPE                               ! look-up values for classification of veg, soils etc.
USE var_lookup,only:iLookATTR                               ! look-up values for local attributes
USE var_lookup,only:iLookPARAM                              ! look-up values for local column model parameters
USE var_lookup,only:iLookINDEX                              ! look-up values for local column index variables
USE var_lookup,only:iLookPROG                               ! look-up values for local column model prognostic (state) variables
USE var_lookup,only:iLookDIAG                               ! look-up values for local column model diagnostic variables 
USE var_lookup,only:iLookFLUX                               ! look-up values for local column model fluxes 
USE var_lookup,only:iLookBVAR                               ! look-up values for basin-average model variables
USE var_lookup,only:iLookBPAR                               ! look-up values for basin-average model parameters
USE var_lookup,only:iLookDECISIONS                          ! look-up values for model decisions
USE var_lookup,only:iLookVarType                            ! look-up values for variable type structure
! provide access to the named variables that describe elements of child  model structures
USE var_lookup,only:childFLUX_MEAN                          ! look-up values for timestep-average model fluxes
! provide access to the named variables that describe model decisions
USE mDecisions_module,only:  &                              ! look-up values for method used to compute derivative
 numerical,   & ! numerical solution
 analytical     ! analytical solution
USE mDecisions_module,only:&                                ! look-up values for LAI decisions
 monthlyTable,& ! LAI/SAI taken directly from a monthly table for different vegetation classes
 specified      ! LAI/SAI computed from green vegetation fraction and winterSAI and summerLAI parameters
USE mDecisions_module,only:&                                ! look-up values for the choice of method for the spatial representation of groundwater
 localColumn, & ! separate groundwater representation in each local soil column
 singleBasin    ! single groundwater store over the entire basin
USE output_stats,only:allocStat,calcStats                   ! module for compiling output statistics
USE globalData,only:nFreq,outFreq                           ! model output files
USE globalData,only:ncid                                    ! file id of netcdf output file
USE var_lookup,only:maxFreq                                 ! maximum # of output files 
implicit none

! *****************************************************************************
! (0) variable definitions
! *****************************************************************************
type(gru_hru_doubleVec)          :: forcStat                   ! x%gru(:)%hru(:)%var(:)%dat -- model forcing data
type(gru_hru_doubleVec)          :: progStat                   ! x%gru(:)%hru(:)%var(:)%dat -- model prognostic (state) variables
type(gru_hru_doubleVec)          :: diagStat                   ! x%gru(:)%hru(:)%var(:)%dat -- model diagnostic variables
type(gru_hru_doubleVec)          :: fluxStat                   ! x%gru(:)%hru(:)%var(:)%dat -- model fluxes
type(gru_hru_doubleVec)          :: indxStat                   ! x%gru(:)%hru(:)%var(:)%dat -- model indices
type(gru_doubleVec)              :: bvarStat                   ! x%gru(:)%var(:)%dat        -- basin-average variabl
! define the primary data structures (scalars)
type(var_i)                      :: timeStruct                 ! x%var(:)                   -- model time data
type(gru_hru_double)             :: forcStruct                 ! x%gru(:)%hru(:)%var(:)     -- model forcing data
type(gru_hru_double)             :: attrStruct                 ! x%gru(:)%hru(:)%var(:)     -- local attributes for each HRU
type(gru_hru_int)                :: typeStruct                 ! x%gru(:)%hru(:)%var(:)     -- local classification of soil veg etc. for each HRU
type(gru_hru_double)             :: mparStruct                 ! x%gru(:)%hru(:)%var(:)     -- model parameters
! define the primary data structures (variable length vectors)
type(gru_hru_intVec)             :: indxStruct                 ! x%gru(:)%hru(:)%var(:)%dat -- model indices
type(gru_hru_doubleVec)          :: progStruct                 ! x%gru(:)%hru(:)%var(:)%dat -- model prognostic (state) variables
type(gru_hru_doubleVec)          :: diagStruct                 ! x%gru(:)%hru(:)%var(:)%dat -- model diagnostic variables
type(gru_hru_doubleVec)          :: fluxStruct                 ! x%gru(:)%hru(:)%var(:)%dat -- model fluxes
! define the basin-average structures
type(gru_double)                 :: bparStruct                 ! x%gru(:)%var(:)            -- basin-average parameters
type(gru_doubleVec)              :: bvarStruct                 ! x%gru(:)%var(:)%dat        -- basin-average variables
! define the ancillary data structures
type(gru_hru_double)             :: dparStruct                 ! x%gru(:)%hru(:)%var(:)     -- default model parameters
! define indices
integer(i4b)                     :: iVar                       ! index of a model variable 
integer(i4b)                     :: iStruct                    ! loop through data structures
integer(i4b)                     :: iGRU
integer(i4b)                     :: localHRU,jHRU,kHRU         ! index of the hydrologic response unit
integer(i4b)                     :: nGRU                       ! number of grouped response units
integer(i4b)                     :: nHRU                       ! number of global hydrologic response units
integer(i4b)                     :: hruCount                   ! number of local hydrologic response units
integer(i4b)                     :: modelTimeStep=0            ! index of model time step
integer(i4b)                     :: waterYearTimeStep=0        ! index of water year
integer(i4b),dimension(maxFreq)  :: outputTimeStep=0           ! timestep in output files
integer(i4b)                     :: ix_gru                     ! index of GRU that corresponds to the global HRU
integer(i4b)                     :: ix_hru                     ! index of local HRU that corresponds to the global HRU
! define the time output
logical(lgt)                     :: printProgress              ! flag to print progress
integer(i4b),parameter           :: ixProgress_im=1000         ! named variable to print progress once per month
integer(i4b),parameter           :: ixProgress_id=1001         ! named variable to print progress once per day
integer(i4b),parameter           :: ixProgress_ih=1002         ! named variable to print progress once per hour
integer(i4b),parameter           :: ixProgress_never=1003      ! named variable to print progress never
integer(i4b)                     :: ixProgress=ixProgress_ih   ! define frequency to write progress
! define the re-start file
logical(lgt)                     :: printRestart               ! flag to print a re-start file
integer(i4b),parameter           :: ixRestart_iy=1000          ! named variable to print a re-start file once per year
integer(i4b),parameter           :: ixRestart_im=1001          ! named variable to print a re-start file once per month
integer(i4b),parameter           :: ixRestart_id=1002          ! named variable to print a re-start file once per day
integer(i4b),parameter           :: ixRestart_never=1003       ! named variable to print a re-start file never
integer(i4b)                     :: ixRestart=ixRestart_im     ! define frequency to write restart files
! define output file
integer(i4b)                     :: ctime1(8)                  ! initial time
character(len=64)                :: output_fileSuffix=''       ! suffix for the output file
character(len=256)               :: summaFileManagerFile=''    ! path/name of file defining directories and files
character(len=256)               :: fileout=''                 ! output filename
! define model control structures
integer(i4b)                     :: nLayers                    ! total number of layers
integer(i4b),parameter           :: no=0                       ! .false.
integer(i4b),parameter           :: yes=1                      ! .true.
logical(lgt)                     :: computeVegFluxFlag         ! flag to indicate if we are computing fluxes over vegetation (.false. means veg is buried with snow) 
type(hru_i),allocatable          :: computeVegFlux(:)          ! flag to indicate if we are computing fluxes over vegetation (.false. means veg is buried with snow) 
type(hru_d),allocatable          :: dt_init(:)                 ! used to initialize the length of the sub-step for each HRU
type(hru_d),allocatable          :: upArea(:)                  ! area upslope of each HRU 
! general local variables        
real(dp)                         :: fracHRU                    ! fractional area of a given HRU (-)
integer(i4b)                     :: fileUnit                   ! file unit (output from file_open; a unit not currently used)
logical(lgt),parameter           :: printTime=.true.           ! flag to print the time information
character(LEN=256),allocatable   :: dataLines(:)               ! vector of character strings from non-comment lines
character(LEN=256),allocatable   :: chardata(:)                ! vector of character data
integer(i4b)                     :: iWord                      ! loop through words in a string
logical(lgt)                     :: flux_mask(size(flux_meta)) ! mask defining desired flux variables
integer(i4b)                     :: forcNcid=integerMissing    ! netcdf id for current netcdf forcing file
integer(i4b)                     :: iFile=1                    ! index of current forcing file from forcing file list
integer(i4b)                     :: forcingStep=-999           ! index of current time step in current forcing file
real(dp),allocatable             :: zSoilReverseSign(:)        ! height at bottom of each soil layer, negative downwards (m)
real(dp),dimension(12)           :: greenVegFrac_monthly       ! fraction of green vegetation in each month (0-1)
real(dp),parameter               :: doubleMissing=-9999._dp    ! missing value
logical(lgt),parameter           :: overwriteRSMIN=.false.     ! flag to overwrite RSMIN
real(dp)                         :: notUsed_canopyDepth        ! NOT USED: canopy depth (m)
real(dp)                         :: notUsed_exposedVAI         ! NOT USED: exposed vegetation area index (m2 m-2)
! error control
integer(i4b)                     :: err=0                      ! error code
character(len=1024)              :: message=''                 ! error message
! output control 
integer(i4b)                     :: iFreq                      ! index for looping through output files

! parallelize the model run
integer(i4b)                     :: startGRU                   ! index of the starting GRU for parallelization run
integer(i4b)                     :: checkHRU                   ! index of the HRU for a single HRU run
character(len=256),allocatable   :: argString(:)               ! string to store command line arguments
integer(i4b)                     :: nArgument                  ! number of command line arguments
integer(i4b)                     :: skipArgument               ! number of command line arguments to skip
integer(i4b)                     :: iArgument                  ! index of command line argument
! *****************************************************************************
! (1) inital priming -- get command line arguments, identify files, etc.
! *****************************************************************************
! get the initial time
call date_and_time(values=ctime1)
print "(A,I2.2,':',I2.2,':',I2.2)", 'start at ',ctime1(5:7)

! check numbers of command-line arguments and obtain all arguments 
nArgument = getCommandArguments(argString)

if (nArgument < 1) then 
 call printCommandHelp()
end if

! initialize command line argument variables
startGRU=integerMissing; nGRU=integerMissing; checkHRU=integerMissing; nHRU=integerMissing
summaFileManagerFile = ''; output_fileSuffix = ''

! loop through all command arguments
skipArgument = 0
do iArgument= 1, nArgument
 if (skipArgument>0) then; skipArgument = skipArgument -1; cycle; end if ! skip the arguments have been read 
 select case (trim(argString(iArgument)))
  case ('-s', '--suffix')
   ! define file suffix
   skipArgument = 1
   if (iArgument+skipArgument>nArgument) call handle_err(1,"missing argument file_suffix; type 'summa.exe --help' for correct usage")   ! check if the number of command line arguments is correct
   output_fileSuffix=trim(argString(iArgument+1))
   print "(A)", "file_suffix is '"//trim(output_fileSuffix)//"'."
  case ('-h', '--hru')
   ! define a single HRU run
   skipArgument = 1
   if (iArgument+skipArgument>nArgument) call handle_err(1,"missing argument checkHRU; type 'summa.exe --help' for correct usage")   ! check if the number of command line arguments is correct
   read(argString(iArgument+1),*) checkHRU ! read the index of the HRU for a single HRU run 
   nHRU=1; nGRU=1                          ! nHRU and nGRU are both one in this case
   ! examines the checkHRU is correct 
   if (checkHRU<1) then
    call handle_err(1,"illegal checkHRU specification; type 'summa.exe --help' for correct usage") 
   else
    print '(A)',' Single-HRU run activated. HRU '//trim(argString(iArgument+1))//' is selected for simulation.'
   end if
  case ('-g','--gru')
   ! define a GRU parallelization run; get the starting GRU and countGRU
   skipArgument = 2
   if (iArgument+skipArgument>nArgument) call handle_err(1,"missing argument startGRU or countGRU; type 'summa.exe --help' for correct usage")   ! check if the number of command line arguments is correct
   read(argString(iArgument+1),*) startGRU ! read the argument of startGRU
   read(argString(iArgument+2),*) nGRU     ! read the argument of countGRU 
   if (startGRU<1 .or. nGRU<1) then
    call handle_err(1,'startGRU and countGRU must be larger than 1.') 
   else
    print '(A)', ' GRU-Parallelization run activated. '//trim(argString(iArgument+2))//' GRUs are selected for simulation.'
   end if
  case ('--help')
   call printCommandHelp
  case default
   ! this is the master file argument, the file will be opened in subroutine summa_SetDirsUndPhiles
   if (len(trim(summaFileManagerFile))>0) then
    ! check if summaFileManagerFile has been defined
    print "(A)", "master_file has been set to '"//trim(summaFileManagerFile)//"'. Argument '"//trim(argString(iArgument))//"' is ignored."
   else
    summaFileManagerFile=trim(argString(iArgument))
    print "(A)", "master_file is '"//trim(summaFileManagerFile)//"'."
   end if           
 end select
end do 
! check if master_file has been received.
if (len(trim(summaFileManagerFile))==0) then
 call handle_err(1,"master_file is not received; type 'summa.exe --help' for correct usage")
end if
if (startGRU/=integerMissing .and. checkHRU/=integerMissing) then 
 call handle_err(1,"single-HRU run and GRU-parallelization run cannot be both selected.")
end if

! set directories and files -- summaFileManager used as command-line argument
call summa_SetDirsUndPhiles(summaFileManagerFile,err,message); call handle_err(err,message)
! initialize the Jacobian flag
doJacobian=.false.

! allocate time structures
call allocLocal(time_meta, refTime,   err=err, message=message); call handle_err(err,message)  ! reference time for the model simulation
call allocLocal(time_meta, startTime, err=err, message=message); call handle_err(err,message)  ! start time for the model simulation 
call allocLocal(time_meta, finshTime, err=err, message=message); call handle_err(err,message)  ! end time for the model simulation

! *****************************************************************************
! (2) populate/check metadata structures
! *****************************************************************************
! populate metadata for all model variables
call popMetadat(err,message); call handle_err(err,message)
! check d2ta structures
call checkStruc(err,message); call handle_err(err,message)

! define the mask to identify the subset of variables in the "child" data structure
! NOTE: The mask is true if (1) the variable is a scalar; *OR* (2) the variable is a flux at the layer interfaces.
!       The interface variables are included because there is a need to calculate the mean flux of surface infiltration (at interface=0)
!        and soil drainage (at interface=nSoil)
flux_mask = ((flux_meta(:)%vartype==iLookVarType%scalarv).or.(flux_meta(:)%vartype==iLookVarType%ifcSoil))

! create the averageFlux metadata structure
call childStruc(flux_meta, flux_mask, averageFlux_meta, childFLUX_MEAN, err, message)
call handle_err(err,message)

! *****************************************************************************
! (3a) read the number of GRUs and HRUs, and allocate the gru-hru mapping structures
! *****************************************************************************

! read and allocate for gru-hru mapping structures and get global variables, e.g, nGRU, nHRU,
! needed in the consequent allocations
! nGRU-is the total number of GRUs of the simulation domain
! nHRU-is the total number of HRUs of the simulation domain
! hruCount-is a local variable for the total number of HRUs in a GRU
call allocate_gru_struc(nGRU,startGRU,nHRU,err,message); call handle_err(err,message)

! *****************************************************************************
! (3b) read the number of snow and soil layers
! *****************************************************************************

! *** TEMPORARY CODE ***
! code will be replaced once merge with the NetCDF branch

! NOTE: currently the same initial conditions for all HRUs and GRUs; will change when shift to NetCDF

! loop through GRUs
do iGRU=1,nGRU

 ! check the GRU-HRU mapping structure is allocated 
 if(.not.allocated(gru_struc(iGRU)%hruInfo)) call handle_err(err,'gru_struc(iGRU)%hruInfo is not allocated')

 ! loop through HRUs
 do localHRU=1,gru_struc(iGRU)%hruCount  ! loop through HRUs within a given GRU

  ! get a vector of non-commented lines
  call file_open(trim(SETNGS_PATH)//trim(MODEL_INITCOND),fileUnit,err,message); call handle_err(err,message)
  call get_vlines(fileUnit,dataLines,err,message); call handle_err(err,message)
  close(fileUnit)
  
  ! get the number of snow and soil layers for each HRU
  gru_struc(iGRU)%hruInfo(localHRU)%nSnow = 0 ! initialize the number of snow layers
  gru_struc(iGRU)%hruInfo(localHRU)%nSoil = 0 ! initialize the number of soil layers
  do iVar=1,size(dataLines)
   ! split the line into an array of words
   call split_line(dataLines(iVar),chardata,err,message); call handle_err(err,message)
   ! check if the line contains initial conditions data (contains the word "snow" or "soil")
   do iword=1,size(chardata)
    if(chardata(iword)=='snow') gru_struc(iGRU)%hruInfo(localHRU)%nSnow = gru_struc(iGRU)%hruInfo(localHRU)%nSnow+1
    if(chardata(iword)=='soil') gru_struc(iGRU)%hruInfo(localHRU)%nSoil = gru_struc(iGRU)%hruInfo(localHRU)%nSoil+1
    if(chardata(iword)=='snow' .or. chardata(iword)=='soil') exit ! exit once read the layer type
   end do
   deallocate(chardata)
  end do
  deallocate(dataLines)

 end do  ! looping through HRUs
end do  ! looping through HRUs

! **** END OF TEMPORARY CODE ***

! *****************************************************************************
! (3b) allocate space for data structures
! *****************************************************************************

! loop through data structures
do iStruct=1,size(structInfo)
 ! allocate space
 select case(trim(structInfo(iStruct)%structName))
  case('time'); call allocGlobal(time_meta,  timeStruct,  err, message)   ! model forcing data
  case('forc'); call allocGlobal(forc_meta,  forcStruct,  err, message)   ! model forcing data
  case('attr'); call allocGlobal(attr_meta,  attrStruct,  err, message)   ! local attributes for each HRU
  case('type'); call allocGlobal(type_meta,  typeStruct,  err, message)   ! local classification of soil veg etc. for each HRU
  case('mpar'); call allocGlobal(mpar_meta,  mparStruct,  err, message)   ! model parameters
  case('indx'); call allocGlobal(indx_meta,  indxStruct,  err, message)   ! model variables
  case('prog'); call allocGlobal(prog_meta,  progStruct,  err, message)   ! model prognostic (state) variables
  case('diag'); call allocGlobal(diag_meta,  diagStruct,  err, message)   ! model diagnostic variables
  case('flux'); call allocGlobal(flux_meta,  fluxStruct,  err, message)   ! model fluxes
  case('bpar'); call allocGlobal(bpar_meta,  bparStruct,  err, message)   ! basin-average parameters
  case('bvar'); call allocGlobal(bvar_meta,  bvarStruct,  err, message)   ! basin-average variables
  case('deriv');cycle   ! model derivatives -- no need to have the GRU dimension
  case default; call handle_err(20,'unable to identify lookup structure')
 end select
 ! check errors
 call handle_err(err,trim(message)//'[structure =  '//trim(structInfo(iStruct)%structName)//']')
enddo  ! looping through data structures

! allocate space for default model parameters
! NOTE: This is done here, rather than in the loop above, because dpar is not one of the "standard" data structures
call allocGlobal(mpar_meta,  dparStruct,  err, message)   ! default model parameters
call handle_err(err,trim(message)//' [problem allocating dparStruct]')

! allocate space for the time step and computeVegFlux flags (recycled for each GRU for subsequent calls to coupled_em)
allocate(dt_init(nGRU),upArea(nGRU),computeVegFlux(nGRU),stat=err)
call handle_err(err,'problem allocating space for dt_init, upArea, or computeVegFlux [GRU]')

! allocate space for the HRUs
do iGRU=1,nGRU
 hruCount = gru_struc(iGRU)%hruCount
 allocate(dt_init(iGRU)%hru(hruCount),upArea(iGRU)%hru(hruCount),computeVegFlux(iGRU)%hru(hruCount),stat=err)
 call handle_err(err,'problem allocating space for dt_init, upArea, or computeVegFlux [HRU]')
end do

! read local attributes for each HRU
call read_attrb(nGRU,startGRU,checkHRU,attrStruct,typeStruct,err,message); call handle_err(err,message)

! *****************************************************************************
! (4a) read description of model forcing datafile used in each HRU
! *****************************************************************************
call ffile_info(nGRU,err,message); call handle_err(err,message)

! *****************************************************************************
! (4b) read model decisions
! *****************************************************************************
call mDecisions(err,message); call handle_err(err,message)

! *****************************************************************************
! (3c) allocate space for output statistics data structures
! *****************************************************************************
! loop through data structures
do iStruct=1,size(structInfo)
 ! allocate space
 select case(trim(structInfo(iStruct)%structName))
  case('forc'); call allocStat(forc_meta,forcStat,err,message)   ! model forcing data
  case('prog'); call allocStat(prog_meta,progStat,err,message)   ! model prognostic (state) variables
  case('diag'); call allocStat(diag_meta,diagStat,err,message)   ! model diagnostic variables
  case('flux'); call allocStat(flux_meta,fluxStat,err,message)   ! model fluxes
  case('indx'); call allocStat(flux_meta,indxStat,err,message)   ! index vars
  case('bvar'); call allocStat(bvar_meta,bvarStat,err,message)   ! basin-average variables
  case default; cycle;
 endselect  
 ! check errors
 call handle_err(err,trim(message)//'[statistics for =  '//trim(structInfo(iStruct)%structName)//']')
enddo ! iStruct

! *****************************************************************************
! (5a) read default model parameters
! *****************************************************************************
! read default values and constraints for model parameters (local column, and basin-average)
call read_pinit(LOCALPARAM_INFO,.TRUE., mpar_meta,localParFallback,err,message); call handle_err(err,message)
call read_pinit(BASINPARAM_INFO,.FALSE.,bpar_meta,basinParFallback,err,message); call handle_err(err,message)

! *****************************************************************************
! (5b) read Noah vegetation and soil tables
! *****************************************************************************

! define monthly fraction of green vegetation
!                           J        F        M        A        M        J        J        A        S        O        N        D
greenVegFrac_monthly = (/0.01_dp, 0.02_dp, 0.03_dp, 0.07_dp, 0.50_dp, 0.90_dp, 0.95_dp, 0.96_dp, 0.65_dp, 0.24_dp, 0.11_dp, 0.02_dp/)

! read Noah soil and vegetation tables
call soil_veg_gen_parm(trim(SETNGS_PATH)//'VEGPARM.TBL',                              & ! filename for vegetation table
                       trim(SETNGS_PATH)//'SOILPARM.TBL',                             & ! filename for soils table
                       trim(SETNGS_PATH)//'GENPARM.TBL',                              & ! filename for general table
                       trim(model_decisions(iLookDECISIONS%vegeParTbl)%cDecision),    & ! classification system used for vegetation
                       trim(model_decisions(iLookDECISIONS%soilCatTbl)%cDecision))      ! classification system used for soils

! read Noah-MP vegetation tables
call read_mp_veg_parameters(trim(SETNGS_PATH)//'MPTABLE.TBL',                         & ! filename for Noah-MP table
                            trim(model_decisions(iLookDECISIONS%vegeParTbl)%cDecision)) ! classification system used for vegetation

! define urban vegetation category
select case(trim(model_decisions(iLookDECISIONS%vegeParTbl)%cDecision))
 case('USGS');                     urbanVegCategory =    1
 case('MODIFIED_IGBP_MODIS_NOAH'); urbanVegCategory =   13
 case('plumberCABLE');             urbanVegCategory = -999
 case('plumberCHTESSEL');          urbanVegCategory = -999
 case('plumberSUMMA');             urbanVegCategory = -999
 case default; call handle_err(30,'unable to identify vegetation category')
end select

! set default model parameters
do iGRU=1,nGRU
 do localHRU=1,gru_struc(iGRU)%hruCount
  ! set parmameters to their default value
  dparStruct%gru(iGRU)%hru(localHRU)%var(:) = localParFallback(:)%default_val         ! x%hru(:)%var(:)
  ! overwrite default model parameters with information from the Noah-MP tables
  call pOverwrite(typeStruct%gru(iGRU)%hru(localHRU)%var(iLookTYPE%vegTypeIndex),  &  ! vegetation category
                  typeStruct%gru(iGRU)%hru(localHRU)%var(iLookTYPE%soilTypeIndex), &  ! soil category
                  dparStruct%gru(iGRU)%hru(localHRU)%var,                          &  ! default model parameters
                  err,message); call handle_err(err,message)            ! error control
  ! copy over to the parameter structure
  mparStruct%gru(iGRU)%hru(localHRU)%var(:) = dparStruct%gru(iGRU)%hru(localHRU)%var(:)
 end do  ! looping through HRUs
 ! set default for basin-average parameters
 bparStruct%gru(iGRU)%var(:) = basinParFallback(:)%default_val
end do  ! looping through GRUs

! *****************************************************************************
! (5c) read trial model parameter values for each HRU, and populate initial data structures
! *****************************************************************************
call read_param(nHRU,typeStruct,mparStruct,err,message); call handle_err(err,message)

! *****************************************************************************
! (5d) compute derived model variables that are pretty much constant for the basin as a whole
! *****************************************************************************

! loop through GRUs
do iGRU=1,nGRU
 ! calculate the fraction of runoff in future time steps
 call fracFuture(bparStruct%gru(iGRU)%var,    &  ! vector of basin-average model parameters
                 bvarStruct%gru(iGRU),        &  ! data structure of basin-average variables
                 err,message)                    ! error control
 call handle_err(err,message)

 ! loop through local HRUs
 do localHRU=1,gru_struc(iGRU)%hruCount
  kHRU=0
  ! check the network topology (only expect there to be one downslope HRU)
  do jHRU=1,gru_struc(iGRU)%hruCount
   if(typeStruct%gru(iGRU)%hru(localHRU)%var(iLookTYPE%downHRUindex) == typeStruct%gru(iGRU)%hru(jHRU)%var(iLookTYPE%hruIndex))then
    if(kHRU==0)then  ! check there is a unique match
     kHRU=jHRU
    else
     call handle_err(20,'multi_driver: only expect there to be one downslope HRU')
    endif  ! (check there is a unique match)
   endif  ! (if identified a downslope HRU)
  end do
  
  ! check that the parameters are consistent
  call paramCheck(mparStruct%gru(iGRU)%hru(localHRU)%var,err,message); call handle_err(err,message)
  
  ! calculate a look-up table for the temperature-enthalpy conversion 
  call E2T_lookup(mparStruct%gru(iGRU)%hru(localHRU)%var,err,message); call handle_err(err,message)
 
  ! read description of model initial conditions -- also initializes model structure components
  ! NOTE: at this stage the same initial conditions are used for all HRUs -- need to modify
  call read_icond(gru_struc(iGRU)%hruInfo(localHRU)%nSnow, & ! number of snow layers
                  gru_struc(iGRU)%hruInfo(localHRU)%nSoil, & ! number of soil layers
                  mparStruct%gru(iGRU)%hru(localHRU)%var,  & ! vector of model parameters
                  indxStruct%gru(iGRU)%hru(localHRU),      & ! data structure of model indices
                  progStruct%gru(iGRU)%hru(localHRU),      & ! model prognostic (state) variables
                  err,message)                           ! error control
  call handle_err(err,message)
  
  ! re-calculate height of each layer
  call calcHeight(&
                  ! input/output: data structures
                  indxStruct%gru(iGRU)%hru(localHRU),   & ! intent(in): layer type
                  progStruct%gru(iGRU)%hru(localHRU),   & ! intent(inout): model prognostic (state) variables for a local HRU
                  ! output: error control
                  err,message); call handle_err(err,message)
  
  ! calculate vertical distribution of root density
  call rootDensty(mparStruct%gru(iGRU)%hru(localHRU)%var,& ! vector of model parameters
                  indxStruct%gru(iGRU)%hru(localHRU),    & ! data structure of model indices
                  progStruct%gru(iGRU)%hru(localHRU),    & ! data structure of model prognostic (state) variables
                  diagStruct%gru(iGRU)%hru(localHRU),    & ! data structure of model diagnostic variables
                  err,message)                         ! error control
  call handle_err(err,message)
  
  ! calculate saturated hydraulic conductivity in each soil layer
  call satHydCond(mparStruct%gru(iGRU)%hru(localHRU)%var,& ! vector of model parameters
                  indxStruct%gru(iGRU)%hru(localHRU),    & ! data structure of model indices
                  progStruct%gru(iGRU)%hru(localHRU),    & ! data structure of model prognostic (state) variables
                  fluxStruct%gru(iGRU)%hru(localHRU),    & ! data structure of model fluxes 
                  err,message)                         ! error control
  call handle_err(err,message)
  
  ! calculate "short-cut" variables such as volumetric heat capacity
  call v_shortcut(mparStruct%gru(iGRU)%hru(localHRU)%var,& ! vector of model parameters
                  diagStruct%gru(iGRU)%hru(localHRU),    & ! data structure of model diagnostic variables
                  err,message)                         ! error control
  call handle_err(err,message)
  
  ! overwrite the vegetation height
  HVT(typeStruct%gru(iGRU)%hru(localHRU)%var(iLookTYPE%vegTypeIndex)) = mparStruct%gru(iGRU)%hru(localHRU)%var(iLookPARAM%heightCanopyTop)
  HVB(typeStruct%gru(iGRU)%hru(localHRU)%var(iLookTYPE%vegTypeIndex)) = mparStruct%gru(iGRU)%hru(localHRU)%var(iLookPARAM%heightCanopyBottom)
  
  ! overwrite the tables for LAI and SAI
  if(model_decisions(iLookDECISIONS%LAI_method)%iDecision == specified)then
   SAIM(typeStruct%gru(iGRU)%hru(localHRU)%var(iLookTYPE%vegTypeIndex),:) = mparStruct%gru(iGRU)%hru(localHRU)%var(iLookPARAM%winterSAI)
   LAIM(typeStruct%gru(iGRU)%hru(localHRU)%var(iLookTYPE%vegTypeIndex),:) = mparStruct%gru(iGRU)%hru(localHRU)%var(iLookPARAM%summerLAI)*greenVegFrac_monthly
  endif
  
  ! initialize canopy drip
  ! NOTE: canopy drip from the previous time step is used to compute throughfall for the current time step
  fluxStruct%gru(iGRU)%hru(localHRU)%var(iLookFLUX%scalarCanopyLiqDrainage)%dat(1) = 0._dp  ! not used  
 end do  ! (looping through HRUs)

 ! compute total area of the upstream HRUS that flow into each HRU
 do localHRU=1,gru_struc(iGRU)%hruCount
  upArea(iGRU)%hru(localHRU) = 0._dp
  do jHRU=1,gru_struc(iGRU)%hruCount
   ! check if jHRU flows into localHRU; assume no exchange between GRUs
   if(typeStruct%gru(iGRU)%hru(jHRU)%var(iLookTYPE%downHRUindex)==typeStruct%gru(iGRU)%hru(localHRU)%var(iLookTYPE%hruIndex))then
    upArea(iGRU)%hru(localHRU) = upArea(iGRU)%hru(localHRU) + attrStruct%gru(iGRU)%hru(jHRU)%var(iLookATTR%HRUarea)
   endif   ! (if jHRU is an upstream HRU)
  end do  ! jHRU
 end do  ! localHRU

 ! identify the total basin area for a GRU (m2)
 associate(totalArea => bvarStruct%gru(iGRU)%var(iLookBVAR%basin__totalArea)%dat(1) )
 totalArea = 0._dp
 do localHRU=1,gru_struc(iGRU)%hruCount
  totalArea = totalArea + attrStruct%gru(iGRU)%hru(localHRU)%var(iLookATTR%HRUarea)
 end do
 end associate

 ! initialize aquifer storage
 ! NOTE: this is ugly: need to add capabilities to initialize basin-wide state variables
 ! There are two options for groundwater:
 !  (1) where groundwater is included in the local column (i.e., the HRUs); and
 !  (2) where groundwater is included for the single basin (i.e., the GRUS, where multiple HRUS drain into a GRU).
 ! For water balance calculations it is important to ensure that the local aquifer storage is zero if groundwater is treated as a basin-average state variable (singleBasin);
 !  and ensure that basin-average aquifer storage is zero when groundwater is included in the local columns (localColumn).
 select case(model_decisions(iLookDECISIONS%spatial_gw)%iDecision)
  ! the basin-average aquifer storage is not used if the groundwater is included in the local column
  case(localColumn)
   bvarStruct%gru(iGRU)%var(iLookBVAR%basin__AquiferStorage)%dat(1) = 0._dp ! set to zero to be clear that there is no basin-average aquifer storage in this configuration 
  ! NOTE: the local column aquifer storage is not used if the groundwater is basin-average
  ! (i.e., where multiple HRUs drain to a basin-average aquifer)
  case(singleBasin)
   bvarStruct%gru(iGRU)%var(iLookBVAR%basin__AquiferStorage)%dat(1) = 1._dp
   do localHRU=1,gru_struc(iGRU)%hruCount
    progStruct%gru(iGRU)%hru(localHRU)%var(iLookPROG%scalarAquiferStorage)%dat(1) = 0._dp  ! set to zero to be clear that there is no local aquifer storage in this configuration
   end do
  case default; call handle_err(20,'unable to identify decision for regional representation of groundwater')
 end select

 ! initialize time step length for each HRU
 do localHRU=1,gru_struc(iGRU)%hruCount
  dt_init(iGRU)%hru(localHRU) = progStruct%gru(iGRU)%hru(localHRU)%var(iLookPROG%dt_init)%dat(1) ! seconds
 end do

end do  ! (looping through GRUs)

! *****************************************************************************
! (5e) initialize first output sequence 
! *****************************************************************************
! define the file
! NOTE: currently assumes that nSoil is constant across the model domain
write(fileout,'(a,i0,a,i0,a,i0,a,i0)') trim(OUTPUT_PATH)//trim(OUTPUT_PREFIX)//'_spinup'//trim(output_fileSuffix),startGRU,'-',startGRU+nGRU-1
call def_output(nHRU,gru_struc(1)%hruInfo(1)%nSoil,fileout,err,message); call handle_err(err,message)

! write local model attributes and parameters to the model output file
do iGRU=1,nGRU
 do localHRU=1,gru_struc(iGRU)%hruCount
  call writeParm(gru_struc(iGRU)%hruInfo(localHRU)%hru_ix,attrStruct%gru(iGRU)%hru(localHRU)%var,attr_meta,err,message); call handle_err(err,message)
  call writeParm(gru_struc(iGRU)%hruInfo(localHRU)%hru_ix,typeStruct%gru(iGRU)%hru(localHRU)%var,type_meta,err,message); call handle_err(err,message)
  call writeParm(gru_struc(iGRU)%hruInfo(localHRU)%hru_ix,mparStruct%gru(iGRU)%hru(localHRU)%var,mpar_meta,err,message); call handle_err(err,message)
 enddo ! HRU
 call writeParm(integerMissing,bparStruct%gru(iGRU)%var,bpar_meta,err,message); call handle_err(err,message)
enddo ! GRU

! ****************************************************************************
! (6) loop through time
! ****************************************************************************
! initialize time step index
waterYearTimeStep = 1
outputTimeStep(1:nFreq) = 1

do modelTimeStep=1,numtim

 ! set print flag
 globalPrintFlag=.false.

 ! read forcing data 
 do iGRU=1,nGRU
  do localHRU=1,gru_struc(iGRU)%hruCount
   
   ! read forcing data
   call read_force(&
                   ! input
                   modelTimeStep,                              & ! intent(in):    time step index
                   iGRU,                                       & ! intent(in):    index of gru
                   localHRU,                                   & ! intent(in):    index of LOCAL hru
                   gru_struc(iGRU)%hruInfo(localHRU)%hru_ncdf, & ! intent(in):    index of hru in netcdf
                   ! input-output
                   iFile,                                      & ! intent(inout): index of current forcing file in forcing file list
                   forcingStep,                                & ! intent(inout): index of read position in time dimension in current netcdf file
                   forcNcid,                                   & ! intent(inout): netcdf file identifier for the current forcing file
                   ! output
                   timeStruct%var,                             & ! intent(out):   time data structure (integer)
                   forcStruct%gru(iGRU)%hru(localHRU)%var,     & ! intent(out):   forcing data structure (double precision)
                   err, message)                               ! intent(out):   error control
   call handle_err(err,message)
  end do 
 end do  ! (end looping through global HRUs)

 ! print progress
 select case(ixProgress)
  case(ixProgress_im);    printProgress = (timeStruct%var(iLookTIME%id)   == 1 .and. timeStruct%var(iLookTIME%ih)   == 0 .and. timeStruct%var(iLookTIME%imin) == 0)
  case(ixProgress_id);    printProgress = (timeStruct%var(iLookTIME%ih)   == 0 .and. timeStruct%var(iLookTIME%imin) == 0)
  case(ixProgress_ih);    printProgress = (timeStruct%var(iLookTIME%imin) == 0)
  case(ixProgress_never); printProgress = .false.
  case default; call handle_err(20,'unable to identify option for the restart file')
 end select
 if(printProgress) write(*,'(i4,1x,5(i2,1x))') timeStruct%var

 ! NOTE: this is done because of the check in coupled_em if computeVegFlux changes in subsequent time steps
 !  (if computeVegFlux changes, then the number of state variables changes, and we need to reoranize the data structures)
 ! compute the exposed LAI and SAI and whether veg is buried by snow
 if(modelTimeStep==1)then  
  do iGRU=1,nGRU
   do localHRU=1,gru_struc(iGRU)%hruCount

    ! get vegetation phenology
    call vegPhenlgy(&
                    ! input/output: data structures
                    model_decisions,                & ! intent(in):    model decisions
                    typeStruct%gru(iGRU)%hru(localHRU), & ! intent(in):    type of vegetation and soil
                    attrStruct%gru(iGRU)%hru(localHRU), & ! intent(in):    spatial attributes
                    mparStruct%gru(iGRU)%hru(localHRU), & ! intent(in):    model parameters
                    progStruct%gru(iGRU)%hru(localHRU), & ! intent(in):    model prognostic variables for a local HRU
                    diagStruct%gru(iGRU)%hru(localHRU), & ! intent(inout): model diagnostic variables for a local HRU
                    ! output
                    computeVegFluxFlag,             & ! intent(out): flag to indicate if we are computing fluxes over vegetation (.false. means veg is buried with snow)
                    notUsed_canopyDepth,            & ! intent(out): NOT USED: canopy depth (m)
                    notUsed_exposedVAI,             & ! intent(out): NOT USED: exposed vegetation area index (m2 m-2)
                    err,message)                      ! intent(out): error control
    call handle_err(err,message)

    ! save the flag for computing the vegetation fluxes
    if(computeVegFluxFlag)      computeVegFlux(iGRU)%hru(localHRU) = yes
    if(.not.computeVegFluxFlag) computeVegFlux(iGRU)%hru(localHRU) = no

    ! define the green vegetation fraction of the grid box (used to compute LAI)
    diagStruct%gru(iGRU)%hru(localHRU)%var(iLookDIAG%scalarGreenVegFraction)%dat(1) = greenVegFrac_monthly(timeStruct%var(iLookTIME%im))

   end do  ! looping through HRUs
  end do  ! looping through GRUs
 endif  ! if the first time step

 ! *****************************************************************************
 ! (7) create a new NetCDF output file, and write parameters and forcing data
 ! *****************************************************************************
 ! check the start of a new water year
 if(timeStruct%var(iLookTIME%im)  ==10 .and. &   ! month = October
    timeStruct%var(iLookTIME%id)  ==1  .and. &   ! day = 1
    timeStruct%var(iLookTIME%ih)  ==1  .and. &   ! hour = 1
    timeStruct%var(iLookTIME%imin)==0)then       ! minute = 0

  ! close any output files that are already open
  do iFreq = 1,nFreq
   if (ncid(iFreq).ne.integerMissing) then
    call nc_file_close(ncid(iFreq),err,message)
    call handle_err(err,message)
   endif
  enddo
 
  ! define the filename
  write(fileout,'(a,i0,a,i0,a,i0,a,i0)') &
                                 trim(OUTPUT_PATH)//trim(OUTPUT_PREFIX)//'_',&
                                 timeStruct%var(iLookTIME%iyyy),'-',timeStruct%var(iLookTIME%iyyy)+1,&
                                 trim(output_fileSuffix),startGRU,'-',startGRU+nGRU-1

  ! define the file
  call def_output(nHRU,gru_struc(1)%hruInfo(1)%nSoil,fileout,err,message); call handle_err(err,message)

  ! write parameters for each HRU, and re-set indices
  do iGRU=1,nGRU
   do localHRU=1,gru_struc(iGRU)%hruCount
    call writeParm(localHRU,attrStruct%gru(iGRU)%hru(localHRU)%var,attr_meta,err,message); call handle_err(err,message)
    call writeParm(localHRU,typeStruct%gru(iGRU)%hru(localHRU)%var,type_meta,err,message); call handle_err(err,message)
    call writeParm(localHRU,mparStruct%gru(iGRU)%hru(localHRU)%var,mpar_meta,err,message); call handle_err(err,message)
    ! re-initalize the indices for midSnow, midSoil, midToto, and ifcToto
    waterYearTimeStep=1
    outputTimeStep=1
    indxStruct%gru(iGRU)%hru(localHRU)%var(iLookINDEX%midSnowStartIndex)%dat(1) = 1
    indxStruct%gru(iGRU)%hru(localHRU)%var(iLookINDEX%midSoilStartIndex)%dat(1) = 1
    indxStruct%gru(iGRU)%hru(localHRU)%var(iLookINDEX%midTotoStartIndex)%dat(1) = 1
    indxStruct%gru(iGRU)%hru(localHRU)%var(iLookINDEX%ifcSnowStartIndex)%dat(1) = 1
    indxStruct%gru(iGRU)%hru(localHRU)%var(iLookINDEX%ifcSoilStartIndex)%dat(1) = 1
    indxStruct%gru(iGRU)%hru(localHRU)%var(iLookINDEX%ifcTotoStartIndex)%dat(1) = 1
   enddo  ! (looping through HRUs)
   call writeParm(integerMissing,bparStruct%gru(iGRU)%var,bpar_meta,err,message); call handle_err(err,message)
  enddo  ! (looping through GRUs)

 endif  ! if start of a new water year, and defining a new file

 ! ****************************************************************************
 ! (8) loop through HRUs and GRUs
 ! ****************************************************************************

 ! initialize variables
 do iGRU=1,nGRU

  ! initialize runoff variables
  bvarStruct%gru(iGRU)%var(iLookBVAR%basin__SurfaceRunoff)%dat(1)    = 0._dp  ! surface runoff (m s-1)
  bvarStruct%gru(iGRU)%var(iLookBVAR%basin__ColumnOutflow)%dat(1)    = 0._dp  ! outflow from all "outlet" HRUs (those with no downstream HRU)
  
  ! initialize baseflow variables
  bvarStruct%gru(iGRU)%var(iLookBVAR%basin__AquiferRecharge)%dat(1)  = 0._dp ! recharge to the aquifer (m s-1)
  bvarStruct%gru(iGRU)%var(iLookBVAR%basin__AquiferBaseflow)%dat(1)  = 0._dp ! baseflow from the aquifer (m s-1)
  bvarStruct%gru(iGRU)%var(iLookBVAR%basin__AquiferTranspire)%dat(1) = 0._dp ! transpiration loss from the aquifer (m s-1)
  
  ! initialize total inflow for each layer in a soil column
  do localHRU=1,gru_struc(iGRU)%hruCount
   fluxStruct%gru(iGRU)%hru(localHRU)%var(iLookFLUX%mLayerColumnInflow)%dat(:) = 0._dp
  end do

  ! loop through HRUs
  do localHRU=1,gru_struc(iGRU)%hruCount
   
   ! This line is used for debugging 
   !write(*,'(4(A5,I10))') 'STEP',istep,'GRU',gru_struc(iGRU)%gru_id,'HRU',gru_struc(iGRU)%hruInfo(localHRU)%hru_id,'TYPE',typeStruct%gru(iGRU)%hru(localHRU)%var(iLookTYPE%vegTypeIndex)

   ! identify the area covered by the current HRU
   fracHRU =  attrStruct%gru(iGRU)%hru(localHRU)%var(iLookATTR%HRUarea) / bvarStruct%gru(iGRU)%var(iLookBVAR%basin__totalArea)%dat(1)
  
   ! assign model layers
   ! NOTE: layer structure is different for each HRU
   gru_struc(iGRU)%hruInfo(localHRU)%nSnow = indxStruct%gru(iGRU)%hru(localHRU)%var(iLookINDEX%nSnow)%dat(1)
   gru_struc(iGRU)%hruInfo(localHRU)%nSoil = indxStruct%gru(iGRU)%hru(localHRU)%var(iLookINDEX%nSoil)%dat(1)
   nLayers                                 = indxStruct%gru(iGRU)%hru(localHRU)%var(iLookINDEX%nLayers)%dat(1)
  
   ! get height at bottom of each soil layer, negative downwards (used in Noah MP)
   allocate(zSoilReverseSign(gru_struc(iGRU)%hruInfo(localHRU)%nSoil),stat=err); call handle_err(err,'problem allocating space for zSoilReverseSign')
   zSoilReverseSign(:) = -progStruct%gru(iGRU)%hru(localHRU)%var(iLookPROG%iLayerHeight)%dat(gru_struc(iGRU)%hruInfo(localHRU)%nSnow+1:nLayers)
  
   ! get NOAH-MP parameters
   call REDPRM(typeStruct%gru(iGRU)%hru(localHRU)%var(iLookTYPE%vegTypeIndex),      & ! vegetation type index
               typeStruct%gru(iGRU)%hru(localHRU)%var(iLookTYPE%soilTypeIndex),     & ! soil type
               typeStruct%gru(iGRU)%hru(localHRU)%var(iLookTYPE%slopeTypeIndex),    & ! slope type index
               zSoilReverseSign,                                                & ! * not used: height at bottom of each layer [NOTE: negative] (m)
               gru_struc(iGRU)%hruInfo(localHRU)%nSoil,                             & ! number of soil layers
               urbanVegCategory)                                                  ! vegetation category for urban areas
  
   ! deallocate height at bottom of each soil layer(used in Noah MP)
   deallocate(zSoilReverseSign,stat=err); call handle_err(err,'problem deallocating space for zSoilReverseSign')
   
   ! overwrite the minimum resistance
   if(overwriteRSMIN) RSMIN = mparStruct%gru(iGRU)%hru(localHRU)%var(iLookPARAM%minStomatalResistance)
  
   ! overwrite the vegetation height
   HVT(typeStruct%gru(iGRU)%hru(localHRU)%var(iLookTYPE%vegTypeIndex)) = mparStruct%gru(iGRU)%hru(localHRU)%var(iLookPARAM%heightCanopyTop)
   HVB(typeStruct%gru(iGRU)%hru(localHRU)%var(iLookTYPE%vegTypeIndex)) = mparStruct%gru(iGRU)%hru(localHRU)%var(iLookPARAM%heightCanopyBottom)
  
   ! overwrite the tables for LAI and SAI
   if(model_decisions(iLookDECISIONS%LAI_method)%iDecision == specified)then
    SAIM(typeStruct%gru(iGRU)%hru(localHRU)%var(iLookTYPE%vegTypeIndex),:) = mparStruct%gru(iGRU)%hru(localHRU)%var(iLookPARAM%winterSAI)
    LAIM(typeStruct%gru(iGRU)%hru(localHRU)%var(iLookTYPE%vegTypeIndex),:) = mparStruct%gru(iGRU)%hru(localHRU)%var(iLookPARAM%summerLAI)*greenVegFrac_monthly
   endif

   ! cycle water pixel
   if (typeStruct%gru(iGRU)%hru(localHRU)%var(iLookTYPE%vegTypeIndex) == isWater) cycle
   
   ! compute derived forcing variables
   call derivforce(timeStruct%var,                    & ! vector of time information
                   forcStruct%gru(iGRU)%hru(localHRU)%var,& ! vector of model forcing data
                   attrStruct%gru(iGRU)%hru(localHRU)%var,& ! vector of model attributes
                   mparStruct%gru(iGRU)%hru(localHRU)%var,& ! vector of model parameters
                   diagStruct%gru(iGRU)%hru(localHRU),    & ! data structure of model diagnostic variables
                   fluxStruct%gru(iGRU)%hru(localHRU),    & ! data structure of model fluxes
                   err,message)                         ! error control
   call handle_err(err,message)
  
   ! ****************************************************************************
   ! (9) run the model
   ! ****************************************************************************
   ! define the need to calculate the re-start file
   select case(ixRestart)
    case(ixRestart_iy);    printRestart = (timeStruct%var(iLookTIME%im) == 1 .and. timeStruct%var(iLookTIME%id) == 1 .and. timeStruct%var(iLookTIME%ih) == 0  .and. timeStruct%var(iLookTIME%imin) == 0)
    case(ixRestart_im);    printRestart = (timeStruct%var(iLookTIME%id) == 1 .and. timeStruct%var(iLookTIME%ih) == 0 .and. timeStruct%var(iLookTIME%imin) == 0)
    case(ixRestart_id);    printRestart = (timeStruct%var(iLookTIME%ih) == 0 .and. timeStruct%var(iLookTIME%imin) == 0)
    case(ixRestart_never); printRestart = .false.
    case default; call handle_err(20,'unable to identify option for the restart file')
   end select
 
   ! set the flag to compute the vegetation flux
   computeVegFluxFlag = (computeVegFlux(iGRU)%hru(localHRU) == yes)
 
   ! run the model for a single parameter set and time step
   call coupled_em(&
                   ! model control
                   modelTimeStep,                               & ! intent(in):    time step index
                   gru_struc(iGRU)%hruInfo(localHRU)%hru_id,    & ! intent(in):    hruId
                   printRestart,                                & ! intent(in):    flag to print a re-start file
                   output_fileSuffix,                           & ! intent(in):    name of the experiment used in the restart file
                   dt_init(iGRU)%hru(localHRU),                 & ! intent(inout): initial time step
                   computeVegFluxFlag,                          & ! intent(inout): flag to indicate if we are computing fluxes over vegetation (.false. means veg is buried with snow)
                   ! data structures (input)
                   timeStruct,                                  & ! intent(in):    model time data
                   typeStruct%gru(iGRU)%hru(localHRU),          & ! intent(in):    local classification of soil veg etc. for each HRU
                   attrStruct%gru(iGRU)%hru(localHRU),          & ! intent(in):    local attributes for each HRU
                   forcStruct%gru(iGRU)%hru(localHRU),          & ! intent(in):    model forcing data
                   mparStruct%gru(iGRU)%hru(localHRU),          & ! intent(in):    model parameters
                   bvarStruct%gru(iGRU),                        & ! intent(in):    basin-average model variables
                   ! data structures (input-output)
                   indxStruct%gru(iGRU)%hru(localHRU),          & ! intent(inout): model indices
                   progStruct%gru(iGRU)%hru(localHRU),          & ! intent(inout): model prognostic variables for a local HRU
                   diagStruct%gru(iGRU)%hru(localHRU),          & ! intent(inout): model diagnostic variables for a local HRU
                   fluxStruct%gru(iGRU)%hru(localHRU),          & ! intent(inout): model fluxes for a local HRU
                   ! error control
                   err,message)            ! intent(out): error control
   call handle_err(err,message)
 
   ! save the flag for computing the vegetation fluxes
   if(computeVegFluxFlag)      computeVegFlux(iGRU)%hru(localHRU) = yes
   if(.not.computeVegFluxFlag) computeVegFlux(iGRU)%hru(localHRU) = no

   kHRU = 0
   ! identify the downslope HRU
   dsHRU: do jHRU=1,gru_struc(iGRU)%hruCount
    if(typeStruct%gru(iGRU)%hru(localHRU)%var(iLookTYPE%downHRUindex) == typeStruct%gru(iGRU)%hru(jHRU)%var(iLookTYPE%hruIndex))then
     if(kHRU==0)then  ! check there is a unique match
      kHRU=jHRU
      exit dsHRU
     endif  ! (check there is a unique match)
    endif  ! (if identified a downslope HRU)
   end do dsHRU
  
   ! add inflow to the downslope HRU
   if(kHRU > 0)then  ! if there is a downslope HRU
    fluxStruct%gru(iGRU)%hru(kHRU)%var(iLookFLUX%mLayerColumnInflow)%dat(:) = fluxStruct%gru(iGRU)%hru(kHRU)%var(iLookFLUX%mLayerColumnInflow)%dat(:)  + fluxStruct%gru(iGRU)%hru(localHRU)%var(iLookFLUX%mLayerColumnOutflow)%dat(:)
  
   ! increment basin column outflow (m3 s-1)
   else
    bvarStruct%gru(iGRU)%var(iLookBVAR%basin__ColumnOutflow)%dat(1)   = bvarStruct%gru(iGRU)%var(iLookBVAR%basin__ColumnOutflow)%dat(1) + sum(fluxStruct%gru(iGRU)%hru(localHRU)%var(iLookFLUX%mLayerColumnOutflow)%dat(:))
   endif
  
   ! increment basin surface runoff (m s-1)
   bvarStruct%gru(iGRU)%var(iLookBVAR%basin__SurfaceRunoff)%dat(1)    = bvarStruct%gru(iGRU)%var(iLookBVAR%basin__SurfaceRunoff)%dat(1)     + fluxStruct%gru(iGRU)%hru(localHRU)%var(iLookFLUX%scalarSurfaceRunoff)%dat(1)    * fracHRU
  
   ! increment basin-average baseflow input variables (m s-1)
   bvarStruct%gru(iGRU)%var(iLookBVAR%basin__AquiferRecharge)%dat(1)  = bvarStruct%gru(iGRU)%var(iLookBVAR%basin__AquiferRecharge)%dat(1)   + fluxStruct%gru(iGRU)%hru(localHRU)%var(iLookFLUX%scalarSoilDrainage)%dat(1)     * fracHRU
   bvarStruct%gru(iGRU)%var(iLookBVAR%basin__AquiferTranspire)%dat(1) = bvarStruct%gru(iGRU)%var(iLookBVAR%basin__AquiferTranspire)%dat(1)  + fluxStruct%gru(iGRU)%hru(localHRU)%var(iLookFLUX%scalarAquiferTranspire)%dat(1) * fracHRU
  
   ! increment aquifer baseflow -- ONLY if baseflow is computed individually for each HRU
   ! NOTE: groundwater computed later for singleBasin
   if(model_decisions(iLookDECISIONS%spatial_gw)%iDecision == localColumn)then
    bvarStruct%gru(iGRU)%var(iLookBVAR%basin__AquiferBaseflow)%dat(1)  = bvarStruct%gru(iGRU)%var(iLookBVAR%basin__AquiferBaseflow)%dat(1)  + fluxStruct%gru(iGRU)%hru(localHRU)%var(iLookFLUX%scalarAquiferBaseflow)%dat(1) * fracHRU
   endif


   ! calculate output Statistics
   call calcStats(forcStat%gru(iGRU)%hru(localHRU)%var,forcStruct%gru(iGRU)%hru(localHRU)%var,forc_meta,waterYearTimeStep,err,message);       call handle_err(err,message)
   call calcStats(progStat%gru(iGRU)%hru(localHRU)%var,progStruct%gru(iGRU)%hru(localHRU)%var,prog_meta,waterYearTimeStep,err,message);       call handle_err(err,message)
   call calcStats(diagStat%gru(iGRU)%hru(localHRU)%var,diagStruct%gru(iGRU)%hru(localHRU)%var,diag_meta,waterYearTimeStep,err,message);       call handle_err(err,message)
   call calcStats(fluxStat%gru(iGRU)%hru(localHRU)%var,fluxStruct%gru(iGRU)%hru(localHRU)%var,flux_meta,waterYearTimeStep,err,message);       call handle_err(err,message)
   call calcStats(indxStat%gru(iGRU)%hru(localHRU)%var,indxStruct%gru(iGRU)%hru(localHRU)%var,indx_meta,waterYearTimeStep,err,message);       call handle_err(err,message)

   ! write the model output to the NetCDF file
   call writeData(waterYearTimeStep,outputTimeStep,forc_meta,forcStat%gru(iGRU)%hru(localHRU)%var,forcStruct%gru(iGRU)%hru(localHRU)%var,indxStruct%gru(iGRU)%hru(localHRU)%var,gru_struc(iGRU)%hruInfo(localHRU)%hru_ix,err,message); call handle_err(err,message)
   call writeData(waterYearTimeStep,outputTimeStep,prog_meta,progStat%gru(iGRU)%hru(localHRU)%var,progStruct%gru(iGRU)%hru(localHRU)%var,indxStruct%gru(iGRU)%hru(localHRU)%var,gru_struc(iGRU)%hruInfo(localHRU)%hru_ix,err,message); call handle_err(err,message)
   call writeData(waterYearTimeStep,outputTimeStep,diag_meta,diagStat%gru(iGRU)%hru(localHRU)%var,diagStruct%gru(iGRU)%hru(localHRU)%var,indxStruct%gru(iGRU)%hru(localHRU)%var,gru_struc(iGRU)%hruInfo(localHRU)%hru_ix,err,message); call handle_err(err,message)
   call writeData(waterYearTimeStep,outputTimeStep,flux_meta,fluxStat%gru(iGRU)%hru(localHRU)%var,fluxStruct%gru(iGRU)%hru(localHRU)%var,indxStruct%gru(iGRU)%hru(localHRU)%var,gru_struc(iGRU)%hruInfo(localHRU)%hru_ix,err,message); call handle_err(err,message)
   call writeData(waterYearTimeStep,outputTimeStep,indx_meta,indxStat%gru(iGRU)%hru(localHRU)%var,indxStruct%gru(iGRU)%hru(localHRU)%var,indxStruct%gru(iGRU)%hru(localHRU)%var,gru_struc(iGRU)%hruInfo(localHRU)%hru_ix,err,message); call handle_err(err,message)
  
   ! increment the model indices
   nLayers = gru_struc(iGRU)%hruInfo(localHRU)%nSnow + gru_struc(iGRU)%hruInfo(localHRU)%nSoil
   indxStruct%gru(iGRU)%hru(localHRU)%var(iLookINDEX%midSnowStartIndex)%dat(1) = indxStruct%gru(iGRU)%hru(localHRU)%var(iLookINDEX%midSnowStartIndex)%dat(1) + gru_struc(iGRU)%hruInfo(localHRU)%nSnow
   indxStruct%gru(iGRU)%hru(localHRU)%var(iLookINDEX%midSoilStartIndex)%dat(1) = indxStruct%gru(iGRU)%hru(localHRU)%var(iLookINDEX%midSoilStartIndex)%dat(1) + gru_struc(iGRU)%hruInfo(localHRU)%nSoil
   indxStruct%gru(iGRU)%hru(localHRU)%var(iLookINDEX%midTotoStartIndex)%dat(1) = indxStruct%gru(iGRU)%hru(localHRU)%var(iLookINDEX%midTotoStartIndex)%dat(1) + nLayers
   indxStruct%gru(iGRU)%hru(localHRU)%var(iLookINDEX%ifcSnowStartIndex)%dat(1) = indxStruct%gru(iGRU)%hru(localHRU)%var(iLookINDEX%ifcSnowStartIndex)%dat(1) + gru_struc(iGRU)%hruInfo(localHRU)%nSnow+1
   indxStruct%gru(iGRU)%hru(localHRU)%var(iLookINDEX%ifcSoilStartIndex)%dat(1) = indxStruct%gru(iGRU)%hru(localHRU)%var(iLookINDEX%ifcSoilStartIndex)%dat(1) + gru_struc(iGRU)%hruInfo(localHRU)%nSoil+1
   indxStruct%gru(iGRU)%hru(localHRU)%var(iLookINDEX%ifcTotoStartIndex)%dat(1) = indxStruct%gru(iGRU)%hru(localHRU)%var(iLookINDEX%ifcTotoStartIndex)%dat(1) + nLayers+1
  

  end do  ! (looping through HRUs)

  ! compute water balance for the basin aquifer
  if(model_decisions(iLookDECISIONS%spatial_gw)%iDecision == singleBasin)then
   call handle_err(20,'multi_driver/bigBucket groundwater code not transferred from old code base yet')
  endif
  
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
  call calcStats(bvarStat%gru(iGRU)%var(:),bvarStruct%gru(iGRU)%var(:),bvar_meta,waterYearTimeStep,err,message); call handle_err(err,message)

  ! write basin-average variables
  call writeBasin(waterYearTimeStep,outputTimeStep,bvar_meta,bvarStat%gru(iGRU)%var,bvarStruct%gru(iGRU)%var,err,message); call handle_err(err,message)

 enddo  ! (looping through GRUs)

 ! write current time to all files
 call WriteTime(waterYearTimeStep,outputTimeStep,time_meta,timeStruct%var,err,message)
 
 ! increment output file timestep
 do iFreq = 1,nFreq
  if (mod(waterYearTimeStep,outFreq(iFreq))==0) then
   outputTimeStep(iFreq) = outputTimeStep(iFreq) + 1
  endif
 enddo
 
 ! increment forcingStep
 forcingStep=forcingStep+1

 ! increment the time index
 waterYearTimeStep = waterYearTimeStep+1

 !print*, 'PAUSE: in driver: testing differences'; read(*,*)
 !stop 'end of time step'

end do  ! (looping through time)

! close any remaining output files
do iFreq = 1,nFreq
 if (ncid(iFreq).ne.integerMissing) then
  call nc_file_close(ncid(iFreq),err,message)
  call handle_err(err,message)
 endif
enddo
 
! deallocate space for dt_init and upArea
deallocate(dt_init,upArea,stat=err); call handle_err(err,'unable to deallocate space for dt_init and upArea')

call stop_program('finished simulation successfully.')

contains

 ! **************************************************************************************************
 ! private subroutine handle_err: error handler
 ! **************************************************************************************************
 subroutine handle_err(err,message)
 ! used to handle error codes
 USE var_lookup,only:iLookPROG,iLookDIAG,iLookFLUX,iLookPARAM,iLookINDEX    ! named variables defining elements in data structure
 implicit none
 ! define dummy variables
 integer(i4b),intent(in)::err             ! error code
 character(*),intent(in)::message         ! error message
 ! return if A-OK
 if(err==0) return
 ! process error messages
 if (err>0) then
  write(*,'(a)') 'FATAL ERROR: '//trim(message)
 else
  write(*,'(a)') 'WARNING: '//trim(message); print*,'(can keep going, but stopping anyway)'
 endif
 if(iGRU>=1 .and. iGRU<=nGRU)then
  if(localHRU>=1 .and. localHRU<=gru_struc(iGRU)%hruCount)then
   ! dump variables
   print*, 'error, variable dump:'
   print*, 'modelTimeStep      = ', modelTimeStep
   print*, 'HRU index          = ', typeStruct%gru(iGRU)%hru(localHRU)%var(iLookTYPE%hruIndex)
   print*, 'pptrate            = ', forcStruct%gru(iGRU)%hru(localHRU)%var(iLookFORCE%pptrate)
   print*, 'airtemp            = ', forcStruct%gru(iGRU)%hru(localHRU)%var(iLookFORCE%airtemp)
   print*, 'theta_res          = ', mparStruct%gru(iGRU)%hru(localHRU)%var(iLookPARAM%theta_res)            ! soil residual volumetric water content (-)
   print*, 'theta_sat          = ', mparStruct%gru(iGRU)%hru(localHRU)%var(iLookPARAM%theta_sat)            ! soil porosity (-)
   print*, 'plantWiltPsi       = ', mparStruct%gru(iGRU)%hru(localHRU)%var(iLookPARAM%plantWiltPsi)         ! matric head at wilting point (m)
   print*, 'soilStressParam    = ', mparStruct%gru(iGRU)%hru(localHRU)%var(iLookPARAM%soilStressParam)      ! parameter in the exponential soil stress function (-)
   print*, 'critSoilWilting    = ', mparStruct%gru(iGRU)%hru(localHRU)%var(iLookPARAM%critSoilWilting)      ! critical vol. liq. water content when plants are wilting (-)
   print*, 'critSoilTranspire  = ', mparStruct%gru(iGRU)%hru(localHRU)%var(iLookPARAM%critSoilTranspire)    ! critical vol. liq. water content when transpiration is limited (-)
   print*, 'scalarSWE          = ', progStruct%gru(iGRU)%hru(localHRU)%var(iLookPROG%scalarSWE)%dat(1)
   print*, 'scalarSnowDepth    = ', progStruct%gru(iGRU)%hru(localHRU)%var(iLookPROG%scalarSnowDepth)%dat(1)
   print*, 'scalarCanopyTemp   = ', progStruct%gru(iGRU)%hru(localHRU)%var(iLookPROG%scalarCanopyTemp)%dat(1)
   print*, 'scalarRainPlusMelt = ', fluxStruct%gru(iGRU)%hru(localHRU)%var(iLookFLUX%scalarRainPlusMelt)%dat(1)
   write(*,'(a,100(i4,1x))'   ) 'layerType          = ', indxStruct%gru(iGRU)%hru(localHRU)%var(iLookINDEX%layerType)%dat
   write(*,'(a,100(f11.5,1x))') 'mLayerDepth        = ', progStruct%gru(iGRU)%hru(localHRU)%var(iLookPROG%mLayerDepth)%dat
   write(*,'(a,100(f11.5,1x))') 'mLayerTemp         = ', progStruct%gru(iGRU)%hru(localHRU)%var(iLookPROG%mLayerTemp)%dat
   write(*,'(a,100(f11.5,1x))') 'mLayerVolFracIce   = ', progStruct%gru(iGRU)%hru(localHRU)%var(iLookPROG%mLayerVolFracIce)%dat
   write(*,'(a,100(f11.5,1x))') 'mLayerVolFracLiq   = ', progStruct%gru(iGRU)%hru(localHRU)%var(iLookPROG%mLayerVolFracLiq)%dat
   print*, 'mLayerMatricHead   = ', progStruct%gru(iGRU)%hru(localHRU)%var(iLookPROG%mLayerMatricHead)%dat
   print*, 'column inflow      = ', fluxStruct%gru(iGRU)%hru(localHRU)%var(iLookFLUX%mLayerColumnInflow)%dat
  endif  ! if HRU is valid
 endif  ! if GRU is valid
 print*,'error code = ', err
 if(allocated(timeStruct%var)) print*, timeStruct%var
 !write(*,'(a)') trim(message)
 stop
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
 
 ! get the final date and time
 call date_and_time(values=ctime2)
 
 elpSec = elapsedSec(ctime1,ctime2)
 
 ! print initial and final date and time
 write(outunit,"(A,I4,'-',I2.2,'-',I2.2,2x,I2,':',I2.2,':',I2.2,'.',I3.3)"),'initial date/time = ',ctime1(1:3),ctime1(5:8)
 write(outunit,"(A,I4,'-',I2.2,'-',I2.2,2x,I2,':',I2.2,':',I2.2,'.',I3.3)"),'  final date/time = ',ctime2(1:3),ctime2(5:8)
 ! print elapsed time
 write(outunit,"(/,A,1PG15.7,A)"),                                          '     elapsed time = ', elpSec,          ' s'
 write(outunit,"(A,1PG15.7,A)"),                                            '       or           ', elpSec/60_sp,    ' m'
 write(outunit,"(A,1PG15.7,A)"),                                            '       or           ', elpSec/3600_sp,  ' h'
 write(outunit,"(A,1PG15.7,A/)"),                                           '       or           ', elpSec/86400_sp, ' d'
 ! stop with message
 print*,'FORTRAN STOP: '//trim(message)
 stop
 end subroutine
 
 ! **************************************************************************************************
 ! private function to obtain all the command line arguments
 ! **************************************************************************************************
 function getCommandArguments(argString)
  character(len=*), allocatable  :: argString(:)
  integer(i4b)                   :: getCommandArguments
  integer(i4b)                   :: iArgument
  
  getCommandArguments = command_argument_count()
  if (getCommandArguments > 0) then
   allocate(argString(getCommandArguments))
   do iArgument = 1, getCommandArguments
    call get_command_argument(iArgument,argString(iArgument))   
   end do  
  end if
  return
 end function getCommandArguments
 
 subroutine printCommandHelp()  
  print "(A)", 'Usage: summa.exe master_file [-s file_suffix] [-g startGRU countGRU] [-h checkHRU]'
  print "(A)", '  summa.exe   -- summa executable'
  print "(A)", '  file_suffix -- text string defining the output file suffix'
  print "(A)", '  master_file -- path/name of master file'
  print "(A)", '  startGRU    -- the index of the first GRU for a parallelization run'
  print "(A)", '  countGRU    -- the number of GRUs for a parallelization run'
  print "(A)", '  checkHRU    -- the index of the HRU for a single HRU run'
  stop 
 end subroutine printCommandHelp

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
!
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
!
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

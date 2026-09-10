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

! ----------------------------------------------------------------------------------------------------------------
! * part 1: parameters that are fixed across multiple instantiations
! ----------------------------------------------------------------------------------------------------------------

MODULE globalData
  ! data types
  USE nr_type
  USE netcdf
  USE data_types,only:gru2hru_map         ! mapping between the GRUs and HRUs
  USE data_types,only:hru2gru_map         ! mapping between the GRUs and HRUs
  USE data_types,only:model_options       ! the model decision structure
  USE data_types,only:file_info           ! metadata for model forcing datafile
  USE data_types,only:par_info            ! default parameter values and parameter bounds
  USE data_types,only:var_info            ! metadata for variables in each model structure
  USE data_types,only:flux2state          ! extended metadata to define flux-to-state mapping
  USE data_types,only:extended_info       ! extended metadata for variables in each model structure
  USE data_types,only:struct_info         ! summary information on all data structures
  USE data_types,only:var_i               ! vector of integers
  USE data_types,only:gru_hru_dom_int     ! x%gru(:)%hru(:)%dom(:)%var(:)  (i4b)
  USE data_types,only:gru_hru_dom_double  ! x%gru(:)%hru(:)%dom(:)%var(:)  (rkind)
  USE data_types,only:gru_hru_double      ! x%gru(:)%hru(:)%var(:)         (rkind)
  USE data_types,only:gru_double          ! x%gru(:)%var(:)                (rkind)
  ! number of variables in each data structure
  USE var_lookup,only:maxvarTime          ! time:                     maximum number variables
  USE var_lookup,only:maxvarForc          ! forcing data:             maximum number variables
  USE var_lookup,only:maxvarAttr          ! attributes:               maximum number variables
  USE var_lookup,only:maxvarType          ! type index:               maximum number variables
  USE var_lookup,only:maxvarId            ! IDs index:                maximum number variables
  USE var_lookup,only:maxvarProg          ! prognostic variables:     maximum number variables
  USE var_lookup,only:maxvarDiag          ! diagnostic variables:     maximum number variables
  USE var_lookup,only:maxvarFlux          ! model fluxes:             maximum number variables
  USE var_lookup,only:maxvarDeriv         ! model derivatives:        maximum number variables
  USE var_lookup,only:maxvarIndx          ! model indices:            maximum number variables
  USE var_lookup,only:maxvarMpar          ! model parameters:         maximum number variables
  USE var_lookup,only:maxvarBvar          ! basin-average variables:  maximum number variables
  USE var_lookup,only:maxvarBpar          ! basin-average parameters: maximum number variables
  USE var_lookup,only:maxvarGrid          ! glacier grid variables:   maximum number variables
  USE var_lookup,only:maxvarDecisions     ! maximum number of decisions
  USE var_lookup,only:maxvarFreq          ! maximum number of output files
  USE var_lookup,only:maxvarLookup        ! maximum number of variables in the lookup
  implicit none
  private

  ! ----------------------------------------------------------------------------------------------------------------
  ! * part 1: parameters that are fixed across multiple instantiations
  ! ----------------------------------------------------------------------------------------------------------------

  ! define missing values
  real(rkind),parameter,public                :: quadMissing    = nr_quadMissing    ! (from nr_type) missing quadruple precision number
  real(rkind),parameter,public                :: realMissing    = nr_realMissing    ! (from nr_type) missing real number
  integer(i4b),parameter,public               :: integerMissing = nr_integerMissing ! (from nr_type) missing integer
  integer(i4b),parameter,public               :: int8Missing    = nr_int8Missing    ! (from nr_type) missing 8-bit integer
  ! define run modes
  integer(i4b),parameter,public               :: iRunModeFull=1                     ! named variable defining running mode as full run (all GRUs)
  integer(i4b),parameter,public               :: iRunModeGRU=2                      ! named variable defining running mode as GRU-parallelization run (GRU subset)
  integer(i4b),parameter,public               :: iRunModeHRU=3                      ! named variable defining running mode as single-HRU run (ONE HRU)
  ! define progress modes
  integer(i4b),parameter,public               :: ixProgress_im=1000                 ! named variable to print progress once per month
  integer(i4b),parameter,public               :: ixProgress_id=1001                 ! named variable to print progress once per day
  integer(i4b),parameter,public               :: ixProgress_ih=1002                 ! named variable to print progress once per hour
  integer(i4b),parameter,public               :: ixProgress_never=1003              ! named variable to print progress never
  integer(i4b),parameter,public               :: ixProgress_it=1004                 ! named variable to print progress every timestep
  ! define restart frequency
  integer(i4b),parameter,public               :: ixRestart_iy=1000                  ! named variable to print a re-start file once per year
  integer(i4b),parameter,public               :: ixRestart_im=1001                  ! named variable to print a re-start file once per month
  integer(i4b),parameter,public               :: ixRestart_id=1002                  ! named variable to print a re-start file once per day
  integer(i4b),parameter,public               :: ixRestart_end=1003                 ! named variable to print a re-start file at the end of a run
  integer(i4b),parameter,public               :: ixRestart_never=1004               ! named variable to print a re-start file never
  ! define output file frequency
  integer(i4b),parameter,public               :: noNewFiles=1001                    ! no new output files
  integer(i4b),parameter,public               :: newFileEveryOct1=1002              ! create a new file on Oct 1 every year (start of the USA water year)
  ! define named variables for "yes" and "no"
  integer(i4b),parameter,public               :: no=0                               ! .false.
  integer(i4b),parameter,public               :: yes=1                              ! .true.
  ! define named variables to describe the domain type
  integer(i4b),parameter,public               :: iname_cas =1000                    ! named variable to denote a canopy air space state variable
  integer(i4b),parameter,public               :: iname_veg =1001                    ! named variable to denote a vegetation state variable
  integer(i4b),parameter,public               :: iname_soil=1002                    ! named variable to denote a soil layer
  integer(i4b),parameter,public               :: iname_snow=1003                    ! named variable to denote a snow layer
  integer(i4b),parameter,public               :: iname_glce=1004                     ! named variable to denote an ice layer
  integer(i4b),parameter,public               :: iname_lake=1005                    ! named variable to denote a lake layer
  integer(i4b),parameter,public               :: iname_aquifer=1006                 ! named variable to denote an aquifer state variable
  ! define named variables to describe the state variable type
  integer(i4b),parameter,public               :: iname_nrgCanair=2001               ! named variable defining the energy of the canopy air space
  integer(i4b),parameter,public               :: iname_nrgCanopy=2002               ! named variable defining the energy of the vegetation canopy
  integer(i4b),parameter,public               :: iname_watCanopy=2003               ! named variable defining the mass of total water on the vegetation canopy
  integer(i4b),parameter,public               :: iname_liqCanopy=2004               ! named variable defining the mass of liquid water on the vegetation canopy
  integer(i4b),parameter,public               :: iname_nrgLayer=3001                ! named variable defining the energy state variable for layers
  integer(i4b),parameter,public               :: iname_watLayer=3002                ! named variable defining the total water state variable for layers
  integer(i4b),parameter,public               :: iname_liqLayer=3003                ! named variable defining the liquid  water state variable for layers
  integer(i4b),parameter,public               :: iname_matLayer=3004                ! named variable defining the matric head state variable for soil layers
  integer(i4b),parameter,public               :: iname_lmpLayer=3005                ! named variable defining the liquid matric potential state variable for soil layers
  integer(i4b),parameter,public               :: iname_watAquifer=3006              ! named variable defining the water storage in the aquifer
  integer(i4b),parameter,public               :: iname_watSnow=3007                 ! named variable defining the water storage in the snowpack  
  integer(i4b),parameter,public               :: iname_watLake=3008                 ! named variable defining the water storage in the lake
  integer(i4b),parameter,public               :: iname_watGlce=3009                 ! named variable defining the water storage in the glacier ice
  integer(i4b),parameter,public               :: iname_watIce=3010                  ! named variable defining the water storage in the lake or glacier ice
  ! define named variables to describe the form and structure of the band-diagonal matrices used in the numerical solver
  ! NOTE: This indexing scheme provides the matrix structure expected by lapack and sundials. Specifically, they require kl extra rows for additional storage.
  !       Consequently, all indices are offset by kl and the total number of bands for storage is 2*kl+ku+1 instead of kl+ku+1.
  integer(i4b),parameter,public               :: nRHS=1                             ! number of unknown variables on the RHS of the linear system A.X=B
  integer(i4b),parameter,public               :: ku=3                               ! number of super-diagonal bands, ku>=3 to accommodate coupled layer above
  integer(i4b),parameter,public               :: kl=4                               ! number of sub-diagonal bands, kl>=4 to accommodate vegetation
  integer(i4b),parameter,public               :: ixDiag=kl+ku+1                     ! index for the diagonal band
  integer(i4b),parameter,public               :: nBands=2*kl+ku+1                   ! length of the leading dimension of the band diagonal matrix
  ! define named variables for the type of matrix used in the numerical solution.
  integer(i4b),parameter,public               :: ixFullMatrix=1001                  ! named variable for the full Jacobian matrix
  integer(i4b),parameter,public               :: ixBandMatrix=1002                  ! named variable for the band diagonal matrix
  ! define indices describing the first and last layers of the Jacobian to print (for debugging)
  integer(i4b),parameter,public               :: iJac1=1                            ! first layer of the Jacobian to print
  integer(i4b),parameter,public               :: iJac2=100                          ! last layer of the Jacobian to print 
  ! define limit checks
  real(rkind),parameter,public                :: iceResidWaterFrac=0.03_rkind       ! residual volumetric liquid water content in ice (-), 0.03 to keep from melting up to 0.05K below freezing, 0.03 for temperate glacier ice
  real(rkind),parameter,public                :: maxVolIceContent=0.7               ! snow maximum volumetric ice content to store water (-)
  real(rkind),parameter,public                :: verySmall=1.e-6_rkind              ! a small number
  real(rkind),parameter,public                :: verySmaller=1.e-12_rkind           ! a smaller number than verySmall
  real(rkind),parameter,public                :: veryBig=1.e+20_rkind               ! a very big number
  ! define summary information on all data structures
  integer(i4b),parameter                      :: nStruct=15                         ! number of data structures
  type(struct_info),parameter,public,dimension(nStruct) :: structInfo=(/&
                   struct_info('time',  'TIME' , maxvarTime ), &                    ! the time data structure
                   struct_info('forc',  'FORCE', maxvarForc ), &                    ! the forcing data structure
                   struct_info('attr',  'ATTR' , maxvarAttr ), &                    ! the attribute data structure
                   struct_info('type',  'TYPE' , maxvarType ), &                    ! the type data structure
                   struct_info('id'  ,  'ID'   , maxvarId   ), &                    ! the type data structure
                   struct_info('mpar',  'PARAM', maxvarMpar ), &                    ! the model parameter data structure
                   struct_info('bpar',  'BPAR' , maxvarBpar ), &                    ! the basin parameter data structure
                   struct_info('bvar',  'BVAR' , maxvarBvar ), &                    ! the basin variable data structure
                   struct_info('grid',  'GRID' , maxvarGrid ), &                    ! the glacier grid data structure
                   struct_info('indx',  'INDEX', maxvarIndx ), &                    ! the model index data structure
                   struct_info('prog',  'PROG',  maxvarProg ), &                    ! the prognostic (state) variable data structure
                   struct_info('diag',  'DIAG' , maxvarDiag ), &                    ! the diagnostic variable data structure
                   struct_info('flux',  'FLUX' , maxvarFlux ), &                    ! the flux data structure
                   struct_info('deriv', 'DERIV', maxvarDeriv), &                    ! the model derivative data structure
                   struct_info('lookup','LOOKUP',maxvarLookup)  /)                  ! the lookup table data structure
  ! fixed model decisions
  logical(lgt)          , parameter, public   :: overwriteRSMIN=.false.             ! flag to overwrite RSMIN

  ! ----------------------------------------------------------------------------------------------------------------
  ! * part 2: globally constant variables/structures that require initialization
  ! ----------------------------------------------------------------------------------------------------------------

 ! define Not-a-Number (NaN)
  real(rkind),save,public                     :: dNaN
  ! define default parameter values and parameter bounds
  type(par_info),save,public                  :: localParFallback(maxvarMpar) ! local column default parameters
  type(par_info),save,public                  :: basinParFallback(maxvarBpar) ! basin-average default parameters
  ! define vectors of metadata
  type(var_info),save,public                  :: time_meta(maxvarTime)        ! model time information
  type(var_info),save,public                  :: forc_meta(maxvarForc)        ! model forcing data
  type(var_info),save,public                  :: attr_meta(maxvarAttr)        ! local attributes
  type(var_info),save,public                  :: type_meta(maxvarType)        ! local classification of veg, soil, etc.
  type(var_info),save,public                  :: id_meta(maxvarId)            ! local classification of veg, soil, etc.
  type(var_info),save,public                  :: mpar_meta(maxvarMpar)        ! local model parameters for each HRU
  type(var_info),save,public                  :: indx_meta(maxvarIndx)        ! local model indices for each HRU
  type(var_info),save,public                  :: prog_meta(maxvarProg)        ! local state variables for each HRU
  type(var_info),save,public                  :: diag_meta(maxvarDiag)        ! local diagnostic variables for each HRU
  type(var_info),save,public                  :: flux_meta(maxvarFlux)        ! local model fluxes for each HRU
  type(var_info),save,public                  :: deriv_meta(maxvarDeriv)      ! local model derivatives for each HRU
  type(var_info),save,public                  :: lookup_meta(maxvarLookup)    ! local lookup tables for each HRU
  type(var_info),save,public                  :: bpar_meta(maxvarBpar)        ! basin parameters for aggregated processes
  type(var_info),save,public                  :: bvar_meta(maxvarBvar)        ! basin variables for aggregated processes
  type(var_info),save,public                  :: grid_meta(maxvarGrid)        ! glacier grid variables
  ! ancillary metadata structures
  type(flux2state),   save,public             :: flux2state_orig(maxvarFlux)  ! named variables for the states affected by each flux (original)
  type(flux2state),   save,public             :: flux2state_liq(maxvarFlux)   ! named variables for the states affected by each flux (liquid water)
  type(extended_info),save,public,allocatable :: averageFlux_meta(:)          ! timestep-average model fluxes
  ! mapping from original to child structures
  integer(i4b),save,public,allocatable        :: forcChild_map(:)             ! index of the child data structure: stats forc
  integer(i4b),save,public,allocatable        :: progChild_map(:)             ! index of the child data structure: stats prog
  integer(i4b),save,public,allocatable        :: diagChild_map(:)             ! index of the child data structure: stats diag
  integer(i4b),save,public,allocatable        :: fluxChild_map(:)             ! index of the child data structure: stats flux
  integer(i4b),save,public,allocatable        :: indxChild_map(:)             ! index of the child data structure: stats indx
  integer(i4b),save,public,allocatable        :: bvarChild_map(:)             ! index of the child data structure: stats bvar
  integer(i4b),save,public,allocatable        :: gridChild_map(:)             ! index of the child data structure: stats grid
  ! child metadata structures
  type(extended_info),save,public,allocatable :: statForc_meta(:)             ! child metadata for stats
  type(extended_info),save,public,allocatable :: statProg_meta(:)             ! child metadata for stats
  type(extended_info),save,public,allocatable :: statDiag_meta(:)             ! child metadata for stats
  type(extended_info),save,public,allocatable :: statFlux_meta(:)             ! child metadata for stats
  type(extended_info),save,public,allocatable :: statIndx_meta(:)             ! child metadata for stats
  type(extended_info),save,public,allocatable :: statBvar_meta(:)             ! child metadata for stats
  type(extended_info),save,public,allocatable :: statGrid_meta(:)             ! child metadata for stats

  ! ----------------------------------------------------------------------------------------------------------------
  ! * part 3: run time variables
  ! ----------------------------------------------------------------------------------------------------------------

  ! define the model decisions
  type(model_options),save,public                  :: model_decisions(maxvarDecisions)  ! the model decision structure
  ! define index variables describing the indices of the first and last HRUs in the forcing file
  integer(i4b),save,public                         :: ixHRUfile_min                     ! minimum index
  integer(i4b),save,public                         :: ixHRUfile_max                     ! maximum index
  ! define mapping structures
  type(gru2hru_map),allocatable,save,public        :: gru_struc(:)                      ! gru2hru map
  type(hru2gru_map),allocatable,save,public        :: index_map(:)                      ! hru2gru map
  ! define variables used for the vegetation phenology
  real(rkind),dimension(12),save,public            :: greenVegFrac_monthly              ! fraction of green vegetation in each month (0-1)
  real(rkind),save,public                          :: minExpLogHgtFac=0.02_rkind        ! factor for minimum height of transition from the exponential to the logarithmic wind profile
  ! define variable used to smooth the ice freezing curve
  real(rkind),save,public                          :: icefrz_mult=10._rkind             ! freezing curve scaling factor multipier of snow to ice, closer to a step function since ice does not hold water
  integer(i4b),save,public                         :: nLakeIceLayers_poss=1             ! number of ice layers in a lake that can accumulate 
  integer(i4b),save,public                         :: nMeltingIceLayers=1               ! number of glacier ice layers that can have a change in total water content 
  real(rkind),save,public                          :: thick4area=0.1                    ! an arbitrary small threshold for glacier thickness to be considered as glacier area (m)
  ! define variables used for horizontal domain type          
  integer(i4b),save,public                         :: upland=1                          ! upland domain
  integer(i4b),save,public                         :: glacCln1=2                        ! glacier clean first domain
  integer(i4b),save,public                         :: glacCln2=3                        ! glacier clean second domain
  integer(i4b),save,public                         :: glacDbr=4                         ! glacier debris domain
  integer(i4b),save,public                         :: wetland=5                         ! wetland domain
  ! define the model output file
  character(len=256),save,public                   :: fileout=''                        ! output filename
  character(len=256),save,public                   :: output_fileSuffix=''              ! suffix for the output file
  ! define controls on model output
  logical(lgt),dimension(maxvarFreq),save,public   :: finalizeStats=.false.             ! flags to finalize statistics
  logical(lgt),save,public                         :: allowRoutingOutput=.false.        ! flag to allow routing variable output (currently very large and slow to write, so turned off by default)
  integer(i4b),save,public                         :: maxLayers                         ! maximum number of layers
  integer(i4b),save,public                         :: maxSnowLayers                     ! maximum number of snow layers
  integer(i4b),save,public                         :: maxSoilLayers                     ! maximum number of soil layers
  integer(i4b),save,public                         :: maxGlaciers                       ! maximum number of glaciers in a GRU
  integer(i4b),save,public                         :: maxGrid                           ! maximum number of grids in a GRU
  integer(i4b),save,public                         :: maxGridX                          ! maximum number of grid cells in the x-direction
  integer(i4b),save,public                         :: maxGridY                          ! maximum number of grid cells in the y-direction
  integer(i4b),save,public                         :: maxGlceLayers                     ! maximum number of glacier ice layers on any glacier
  integer(i4b),save,public                         :: maxWetlands                       ! maximum number of wetlands in a GRU
  integer(i4b),save,public                         :: maxLakeLayers                     ! maximum number of lake layers
  ! define control variables
  integer(i4b),save,public                         :: startGRU                          ! index of the starting GRU for parallelization run
  integer(i4b),save,public                         :: checkHRU                          ! index of the HRU for a single HRU run
  integer(i4b),save,public                         :: iRunMode                          ! define the current running mode
  integer(i4b),save,public                         :: nThreads=1                        ! number of threads
  integer(i4b),save,public                         :: ixProgress=ixProgress_id          ! define frequency to write progress
  integer(i4b),save,public                         :: ixRestart=ixRestart_never         ! define frequency to write restart files
  integer(i4b),save,public                         :: newOutputFile=noNewFiles          ! define option for new output files
  ! define common variables
  integer(i4b),save,public                         :: numtim                            ! number of time steps
  integer(i4b),save,public                         :: maxDOM                            ! number of domains (every HRU may have multiple) in the run space
  integer(i4b),save,public                         :: nHRUrun                           ! number of HRUs in the run space
  integer(i4b),save,public                         :: nGRUrun                           ! number of GRUs in the run space
  real(rkind),save,public                          :: data_step                         ! length of the time_step
  real(rkind),save,public                          :: refJulDay                         ! reference time in fractional julian days
  real(rkind),save,public                          :: refJulDay_data                    ! reference time in fractional julian days (data files)
  real(rkind),save,public                          :: dJulianStart                      ! julian day of start time of simulation
  real(rkind),save,public                          :: dJulianFinsh                      ! julian day of end time of simulation
  integer(i4b),save,public                         :: nHRUfile                          ! number of HRUs in the file
  integer(i4b),save,public                         :: urbanVegCategory                  ! vegetation category for urban areas
  logical(lgt),save,public                         :: globalPrintFlag=.false.           ! flag to compute the Jacobian, residual, and step progress
  integer(i4b),save,public                         :: chunksize=1024                    ! chunk size for the netcdf read/write
  integer(i4b),save,public                         :: outputPrecision=nf90_double       ! variable type
  integer(i4b),save,public                         :: outputCompressionLevel=4          ! output netcdf file deflate level: 0-9. 0 is no compression.
  ! define data structures for the buffered read
  integer(i4b),save,public                         :: ixStartRead                       ! start index of the data read
  real(rkind),save,public,allocatable              :: fulltimeVec(:)                    ! full time vector in an input file (nRead)
  type(gru_hru_double),save,public,allocatable     :: fullforcingStruct(:)              ! x(:)%gru(:)%hru(:)%var(:) -- full model forcing data
  ! define data structures for the buffered write  
  type(gru_hru_dom_int),save,public,allocatable    :: fullIndxSave(:)                   ! x%gru(:)%hru(:)%dom(:)%var(:) -- saved output for indices
  type(gru_hru_double),save,public,allocatable     :: fullForcSave(:)                   ! x%gru(:)%hru(:)%var(:)        -- saved output for forcing
  type(gru_hru_dom_double),save,public,allocatable :: fullProgSave(:)                   ! x%gru(:)%hru(:)%dom(:)%var(:) -- saved output for prognostic variables
  type(gru_hru_dom_double),save,public,allocatable :: fullDiagSave(:)                   ! x%gru(:)%hru(:)%dom(:)%var(:) -- saved output for diagnostic variables
  type(gru_hru_dom_double),save,public,allocatable :: fullFluxSave(:)                   ! x%gru(:)%hru(:)%dom(:)%var(:) -- saved output for flux variables
  type(gru_double),save,public,allocatable         :: fullBvarSave(:)                   ! x%gru(:)%var(:)               -- saved output for basin variables
  ! define result from the time calls
  integer(i4b),dimension(8),save,public            :: startInit,endInit                 ! date/time for the start and end of the initialization
  integer(i4b),dimension(8),save,public            :: startSetup,endSetup               ! date/time for the start and end of the parameter setup
  integer(i4b),dimension(8),save,public            :: startRestart,endRestart           ! date/time for the start and end to read restart data
  integer(i4b),dimension(8),save,public            :: startRead,endRead                 ! date/time for the start and end of the data read
  integer(i4b),dimension(8),save,public            :: startWrite,endWrite               ! date/time for the start and end of the stats/write
  integer(i4b),dimension(8),save,public            :: startPhysics,endPhysics           ! date/time for the start and end of the physics
  ! define elapsed time
  real(rkind),save,public                          :: elapsedInit                       ! elapsed time for the initialization
  real(rkind),save,public                          :: elapsedSetup                      ! elapsed time for the parameter setup
  real(rkind),save,public                          :: elapsedRestart                    ! elapsed time to read restart data
  real(rkind),save,public                          :: elapsedRead                       ! elapsed time for the data read
  real(rkind),save,public                          :: elapsedWrite                      ! elapsed time for the stats/write
  real(rkind),save,public                          :: elapsedPhysics                    ! elapsed time for the physics
  real(rkind),save,public                          :: elapsedUpdateArea                 ! elapsed time for updating glacier and wetland area
  ! define ancillary data structures
  type(var_i),save,public                          :: startTime                         ! start time for the model simulation
  type(var_i),save,public                          :: finshTime                         ! end time for the model simulation
  type(var_i),save,public                          :: refTime                           ! reference time for the model simulation
  type(var_i),save,public                          :: oldTime                           ! time for the previous model time step
  ! output file information
  logical(lgt),dimension(maxvarFreq),save,public   :: outFreq                           ! true if the output frequency is desired
  integer(i4b),dimension(maxvarFreq),save,public   :: ncid                              ! netcdf output file id
  ! look-up values for the choice of the time zone information (formerly in modelDecisions module)
  integer(i4b),parameter,public                    :: ncTime=1                          ! time zone information from NetCDF file (timeOffset = longitude/15. - ncTimeOffset)
  integer(i4b),parameter,public                    :: utcTime=2                         ! all times in UTC (timeOffset = longitude/15. hours)
  integer(i4b),parameter,public                    :: localTime=3                       ! all times local (timeOffset = 0)
  ! define metadata for model forcing datafile non-Actors
  type(file_info),save,public,allocatable          :: forcFileInfo(:)                   ! file info for model forcing data
  ! define indices in the forcing data files non-Actors
  integer(i4b),save,public                         :: iFile=1                           ! index of current forcing file from forcing file list
  integer(i4b),save,public                         :: forcingStep=integerMissing        ! index of current time step in current forcing file
  integer(i4b),save,public                         :: forcNcid=integerMissing           ! netcdf id for current netcdf forcing file
  ! define controls on model output non-Actors
  integer(i4b),dimension(maxvarFreq),save,public   :: statCounter=0                     ! time counter for stats
  integer(i4b),dimension(maxvarFreq),save,public   :: outputTimeStep=0                  ! timestep in output files
  logical(lgt),dimension(maxvarFreq),save,public   :: resetStats=.true.                 ! flags to reset statistics
  ! define common variables non-Actors
  real(rkind),save,public                          :: fracJulDay                        ! fractional julian days since the start of year
  real(rkind),save,public                          :: tmZoneOffsetFracDay               ! time zone offset in fractional days
  integer(i4b),save,public                         :: yearLength                        ! number of days in the current year
  ! define fixed dimensions
  integer(i4b),parameter,public                    :: nSpecBand=2                       ! number of spectral bands
  integer(i4b),parameter,public                    :: nTimeDelay=2000                   ! number of time steps in the time delay histogram (default: ~1 season = 24*365/4)
  ! printing step frequency
  integer(i4b),parameter,public                    :: print_step_freq = 1000            ! frequency (in time steps) to print number of steps taken in solver
END MODULE globalData

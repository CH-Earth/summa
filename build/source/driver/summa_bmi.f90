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

module summabmi
  ! provides functions needed for summa driver routines adding BMI functions
  ! *****************************************************************************
  ! * use desired modules
  ! *****************************************************************************
  ! data types
  USE,intrinsic :: iso_c_binding, only: c_ptr, c_loc, c_f_pointer
  USE nr_type                                                 ! variable types, etc.
! NGEN_ACTIVE is to be set when running in the Nextgen framework
! https://github.com/NOAA-OWP/ngen
#ifdef NGEN_ACTIVE
  use bmif_2_0_iso                                            ! BMI libraries NextGen
#else
  use bmif_2_0                                                ! BMI libraries standard
#endif
  USE summa_type, only: summa1_type_dec                       ! master summa data type
  USE data_types, only: gru2hru_map                           ! mapping between the GRUs and HRUs
  USE data_types, only: hru2gru_map                           ! mapping between the GRUs and HRUs
  USE data_types, only: model_options                         ! the model decision structure
  USE data_types, only: var_info                              ! metadata for variables in each model structure
  USE data_types, only: extended_info                         ! extended metadata for variables in each model structure  
  USE data_types, only: file_info                             ! metadata for model forcing datafile
  USE data_types, only: var_i                                 ! vector of integers
  ! subroutines and functions: model setup
  USE summa_init, only: summa_initialize                      ! used to allocate/initialize summa data structures
  USE summa_setup, only: summa_paramSetup                     ! used to initialize parameter data structures (e.g. vegetation and soil parameters)
  USE summa_restart, only: summa_readRestart                  ! used to read restart data and reset the model state
  ! subroutines and functions: model simulation
  USE summa_forcing, only: summa_readForcing                  ! used to read forcing data
  USE summa_modelRun, only: summa_runPhysics                  ! used to run the summa physics for one time step
  USE summa_writeOutput, only: summa_writeOutputFiles         ! used to write the summa output files
  ! utility functions
  USE summa_util, only: stop_program                          ! used to stop the summa program (with errors)
  USE summa_util, only: handle_err                            ! used to process errors
  ! global data
  USE globalData, only: integerMissing                        ! missing integer
  USE globalData, only: print_step_freq                       ! frequency of printing model output
  USE globalData, only: dJulianStart                          ! julian day of start time of simulation
  USE globalData, only: dJulianFinsh                          ! julian day of end time of simulation
  USE globalData, only: data_step                             ! length of time steps for the outermost timeloop
  USE globalData, only: startGRU                              ! index of the starting GRU for run
  USE globalData, only: ixRestart_iy,ixRestart_im,ixRestart_id,ixRestart_end,ixRestart_never ! restart file frequency options
  USE globalData, only: gru_struc                             ! HRU information for given GRU
  USE globalData, only: index_map                             ! GRU information for given HRU
  USE globalData, only: model_decisions                       ! model decision structure 
  USE globalData, only: fileout, output_fileSuffix            ! output filename and suffix
  USE globalData, only: outFreq                               ! output frequency flags
  USE globalData, only: ncid                                  ! netcdf output file id
  USE globalData, only: maxLayers,maxSnowLayers,maxLakeLayers,maxGlceLayers,maxSoilLayers ! maximum number of layers
  USE globalData, only: maxGlaciers,maxWetlands,maxGrid,maxGridX,maxGridY ! maximum grids and grid dimensions
  USE globalData, only: ixProgress                            ! define frequency to write progress
  USE globalData, only: ixRestart                             ! define frequency to write restart files
  USE globalData, only: newOutputFile                         ! define option for new output files
  USE globalData, only: nHRUrun                               ! number of HRUs in the run domain
  USE globalData, only: nGRUrun                               ! number of GRUs in the run domain
  USE globalData, only: maxDOM                                ! maximum number of domains (every HRU may have multiple) in the run domain
  USE globalData, only: urbanVegCategory                      ! vegetation category for urban areas
#ifndef NGEN_FORCING_ACTIVE
  USE globalData, only: ixHRUfile_min, ixHRUfile_max          ! indices of the first and last HRUs in the forcing file
  USE globalData, only: nHRUfile                              ! number of HRUs in the forcing file
  USE globalData, only: forcFileInfo                          ! file info for model forcing data
  USE globalData, only: iFile                                 ! index of current forcing file from forcing file list
  USE globalData, only: forcingStep, forcNcid                 ! index of current time step in current forcing file and netcdf id
#endif
  USE globalData, only: statCounter, outputTimeStep           ! timestep in output files and time counter for stats
  USE globalData, only: resetStats, finalizeStats             ! flags to reset and finalize statistics
  USE globalData, only: oldTime                               ! time for the previous model time step
  USE globalData, only: elapsedInit                           ! elapsed time for the initialization
  USE globalData, only: elapsedSetup                          ! elapsed time for the parameter setup
  USE globalData, only: elapsedRestart                        ! elapsed time to read restart data
  USE globalData, only: elapsedRead                           ! elapsed time for the forcing data read
  USE globalData, only: elapsedWrite                          ! elapsed time for the stats/write
  USE globalData, only: elapsedPhysics                        ! elapsed time for the physics
  USE multiconst, only: secprday                              ! number of seconds in a day
  ! provide access to the named variables that describe elements of parent model structures
  USE var_lookup, only: iLookTIME                             ! named variables for time data structure
  USE var_lookup, only: iLookATTR                             ! named variables for real valued attribute data structure
  USE var_lookup, only: iLookFORCE                            ! named variables for forcing data structure
  USE var_lookup, only: iLookFLUX                             ! named variables for local flux variables
  USE var_lookup, only: iLookDIAG                             ! named variables for local diagnostic variables
  USE var_lookup, only: iLookPROG                             ! named variables for local prognostic variables
  USE var_lookup, only: iLookBVAR                             ! named variables for basin (GRU) variables
  USE var_lookup, only: maxvarDecisions                       ! maximum number of decisions
  USE var_lookup, only: maxvarFreq                            ! maximum number of output files

  implicit none

  ! Define the attributes of the model.
  type :: summa_model
     integer(i4b)                       :: timeStep                          ! index of model time step
     type(summa1_type_dec), allocatable :: summa1_struc(:)                   ! master summa data structure
     type(gru2hru_map), allocatable     :: gru_struc(:)                      ! HRU information for given GRU
     type(hru2gru_map), allocatable     :: index_map(:)                      ! GRU information for given HRU
     type(model_options)                :: model_decisions(maxvarDecisions)  ! the model decisions, could change if different decisions for different GRUs
     character(len=256)                 :: fileout, output_fileSuffix        ! output filename and suffix
     logical(lgt),dimension(maxvarFreq) :: outFreq                           ! true if the output frequency is desired
     integer(i4b),dimension(maxvarFreq) :: ncid                              ! netcdf output file id
     integer(i4b)                       :: maxLayers,maxSnowLayers           ! maximum number of layers and snow layers, could be different for different GRUs
     integer(i4b)                       :: maxLakeLayers,maxGlceLayers,maxSoilLayers ! maximum number of lake, soil, and glce layers
     integer(i4b)                       :: maxGlaciers,maxWetlands           ! maximum number of glaciers and wetlands, could be different for different GRUs
     integer(i4b)                       :: maxGrid,maxGridX,maxGridY         ! maximum grids and grid dimensions, could be different for different GRUs
     integer(i4b)                       :: ixProgress                        ! define frequency to write progress
     integer(i4b)                       :: ixRestart                         ! define frequency to write restart files
     integer(i4b)                       :: newOutputFile                     ! define option for new output files
     integer(i4b)                       :: nHRUrun                           ! number of HRUs in the run domain
     integer(i4b)                       :: nGRUrun                           ! number of GRUs in the run domain
     integer(i4b)                       :: maxDOM                            ! maximum number of domains (every HRU may have multiple) in the run domain
     integer(i4b)                       :: urbanVegCategory                  ! vegetation category for urban areas
#ifndef NGEN_FORCING_ACTIVE
     integer(i4b)                       :: ixHRUfile_min, ixHRUfile_max      ! indices of the first and last HRUs in the forcing file
     integer(i4b)                       :: nHRUfile                          ! number of HRUs in the forcing file
     type(file_info), allocatable       :: forcFileInfo(:)                   ! file info for model forcing data
     integer(i4b)                       :: iFile                             ! index of current forcing file from forcing file list
     integer(i4b)                       :: forcingStep, forcNcid             ! index of current time step in current forcing file and netcdf id
#endif
     integer(i4b),dimension(maxvarFreq) :: statCounter, outputTimeStep       ! time counter for stats and time step in output files
     logical(lgt),dimension(maxvarFreq) :: resetStats, finalizeStats         ! flags to reset and finalize statistics
     type(var_i)                        :: oldTime                           ! time for the previous model time step
     real(rkind)                        :: elapsedInit                       ! elapsed time for the initialization
     real(rkind)                        :: elapsedSetup                      ! elapsed time for the parameter setup
     real(rkind)                        :: elapsedRestart                    ! elapsed time to read restart data
     real(rkind)                        :: elapsedRead                       ! elapsed time for the forcing data read
     real(rkind)                        :: elapsedWrite                      ! elapsed time for the stats/write
     real(rkind)                        :: elapsedPhysics                    ! elapsed time for the physics
  end type summa_model

  type, extends (bmi) :: summa_bmi
     private
     type (summa_model) :: model
   contains
     procedure :: get_component_name => summa_component_name
     procedure :: get_input_item_count => summa_input_item_count
     procedure :: get_output_item_count => summa_output_item_count
     procedure :: get_input_var_names => summa_input_var_names
     procedure :: get_output_var_names => summa_output_var_names
     procedure :: initialize => summa_bmi_initialize
     procedure :: finalize => summa_finalize
     procedure :: get_start_time => summa_start_time
     procedure :: get_end_time => summa_end_time
     procedure :: get_current_time => summa_current_time
     procedure :: get_time_step => summa_time_step
     procedure :: get_time_units => summa_time_units
     procedure :: update => summa_update
     procedure :: update_until => summa_update_until
     procedure :: get_var_grid => summa_var_grid
     procedure :: get_grid_type => summa_grid_type
     procedure :: get_grid_rank => summa_grid_rank
     procedure :: get_grid_shape => summa_grid_shape
     procedure :: get_grid_size => summa_grid_size
     procedure :: get_grid_spacing => summa_grid_spacing
     procedure :: get_grid_origin => summa_grid_origin
     procedure :: get_grid_x => summa_grid_x
     procedure :: get_grid_y => summa_grid_y
     procedure :: get_grid_z => summa_grid_z
     procedure :: get_grid_node_count => summa_grid_node_count
     procedure :: get_grid_edge_count => summa_grid_edge_count
     procedure :: get_grid_face_count => summa_grid_face_count
     procedure :: get_grid_edge_nodes => summa_grid_edge_nodes
     procedure :: get_grid_face_edges => summa_grid_face_edges
     procedure :: get_grid_face_nodes => summa_grid_face_nodes
     procedure :: get_grid_nodes_per_face => summa_grid_nodes_per_face
     procedure :: get_var_type => summa_var_type
     procedure :: get_var_units => summa_var_units
     procedure :: get_var_itemsize => summa_var_itemsize
     procedure :: get_var_nbytes => summa_var_nbytes
     procedure :: get_var_location => summa_var_location
     procedure :: get_value_int => summa_get_int
     procedure :: get_value_float => summa_get_float
     procedure :: get_value_double => summa_get_double
     generic :: get_value => &
          get_value_int, &
          get_value_float, &
          get_value_double
     procedure :: get_value_ptr_int => summa_get_ptr_int
     procedure :: get_value_ptr_float => summa_get_ptr_float
     procedure :: get_value_ptr_double => summa_get_ptr_double
     generic :: get_value_ptr => &
          get_value_ptr_int, &
          get_value_ptr_float, &
          get_value_ptr_double
     procedure :: get_value_at_indices_int => summa_get_at_indices_int
     procedure :: get_value_at_indices_float => summa_get_at_indices_float
     procedure :: get_value_at_indices_double => summa_get_at_indices_double
     generic :: get_value_at_indices => &
          get_value_at_indices_int, &
          get_value_at_indices_float, &
          get_value_at_indices_double
     procedure :: set_value_int => summa_set_int
     procedure :: set_value_float => summa_set_float
     procedure :: set_value_double => summa_set_double
     generic :: set_value => &
          set_value_int, &
          set_value_float, &
          set_value_double
     procedure :: set_value_at_indices_int => summa_set_at_indices_int
     procedure :: set_value_at_indices_float => summa_set_at_indices_float
     procedure :: set_value_at_indices_double => summa_set_at_indices_double
     generic :: set_value_at_indices => &
          set_value_at_indices_int, &
          set_value_at_indices_float, &
          set_value_at_indices_double
  end type summa_bmi

  private
  public :: summa_bmi

  ! *****************************************************************************
  ! * variable definitions
  ! *****************************************************************************
  character (len=BMI_MAX_COMPONENT_NAME), target :: &
       component_name = "Structure for Unifying Multiple Modeling Alternatives: SUMMA"
  ! define parameters for the model simulation
  integer(i4b), parameter            :: n=1                        ! number of instantiations
  ! Exchange items
#ifdef NGEN_ACTIVE
  integer, parameter :: input_item_count = 8
#else
  integer, parameter :: input_item_count = 7
#endif
  integer, parameter :: output_item_count = 16
  character (len=BMI_MAX_VAR_NAME), target,dimension(input_item_count)  :: input_items
  character (len=BMI_MAX_VAR_NAME), target,dimension(output_item_count) :: output_items
  ! Buffers behind summa_get_ptr_int/float.  The BMI contract is that the returned pointer
  ! stays valid after the call, so these cannot be the callers' own automatic arrays.
  real,    target, allocatable :: ptr_buffer(:)
  integer, target              :: iptr_buffer
  ! ---------------------------------------------------------------------------------------

  contains

   ! *****************************************************************************
   ! * model setup/initialization
   ! *****************************************************************************
   function summa_bmi_initialize(this, config_file) result (bmi_status)
     class (summa_bmi), intent(out) :: this
     character (len=*), intent(in) :: config_file
     ! error control
     integer(i4b)                       :: err=0                      ! error code
     character(len=1024)                :: message=''                 ! error message
     character(len=1024)                :: file_manager
     character(len=16)                  :: restart_print_freq
     integer(i4b)                       :: attrib_file_HRU_order
     character(len=16)                  :: ixRestart_str
     integer  :: bmi_status,i,fu,rc
     ! namelist definition
     namelist /parameters/ file_manager, attrib_file_HRU_order, restart_print_freq

     ! initialize global variables
     this%model%timeStep = 0
     fileout = ''
     output_fileSuffix = ''
#ifndef NGEN_FORCING_ACTIVE
     iFile = 1
     forcingStep = integerMissing
     forcNcid = integerMissing
#endif
     statCounter = 0
     outputTimeStep = 0
     resetStats = .true.
     finalizeStats = .false.

     ! allocate space for the master summa structure
     allocate(this%model%summa1_struc(n), stat=err)
     if(err/=0) call stop_program(1, 'problem allocating master summa structure')

     ! if using the BMI interface, there is an argument pointing to the file manager file
     !  then make sure summaFileManagerFile is set before executing initialization
     if (len(config_file) > 0)then
#ifdef NGEN_ACTIVE
       ! with NGEN the argument gives the file manager file as an input parameter in a namelist
       open (action='read', file=config_file, iostat=rc, newunit=fu)
       read (nml=parameters, iostat=rc, unit=fu)
       this%model%summa1_struc(n)%summaFileManagerFile=trim(file_manager)
       startGRU = attrib_file_HRU_order
       ixRestart_str = trim(restart_print_freq)
       select case (ixRestart_str)
        case ('y' , 'year');  ixRestart = ixRestart_iy
        case ('m' , 'month'); ixRestart = ixRestart_im
        case ('d' , 'day');   ixRestart = ixRestart_id
        case ('e' , 'end');   ixRestart = ixRestart_end
        case ('n' , 'never'); ixRestart = ixRestart_never
        case default;  print*, 'unknown frequency to write restart files in NGEN parameters namelist'; err=1; return
       end select
       print *
       print *, "INFO: NGEN detected, using file manager file ", trim(file_manager), " and GRU index in files ", startGRU
#else
       ! without NGEN the argument gives the file manager file directly
       ! Note, if this is more than 80 characters the pre-built BMI libraries will fail
       this%model%summa1_struc(n)%summaFileManagerFile=trim(config_file)
#endif
     endif

     ! declare and allocate summa data structures and initialize model state to known values
     call summa_initialize(this%model%summa1_struc(n), err, message)
     call handle_err(err, message)

     ! initialize parameter data structures (e.g. vegetation and soil parameters)
     call summa_paramSetup(this%model%summa1_struc(n), err, message)
     call handle_err(err, message)

     ! read restart data and reset the model state
     call summa_readRestart(this%model%summa1_struc(n), err, message)
     call handle_err(err, message)

     ! get global variables that are constants throughout the model simulation
     this%model%gru_struc = gru_struc
     this%model%index_map = index_map
     this%model%model_decisions = model_decisions
     this%model%output_fileSuffix = output_fileSuffix
     this%model%maxLayers = maxLayers
     this%model%maxSnowLayers = maxSnowLayers
     this%model%maxLakeLayers = maxLakeLayers
     this%model%maxGlceLayers = maxGlceLayers
     this%model%maxSoilLayers = maxSoilLayers
     this%model%maxGlaciers = maxGlaciers
     this%model%maxWetlands = maxWetlands
     this%model%maxGrid = maxGrid
     this%model%maxGridX = maxGridX
     this%model%maxGridY = maxGridY
     this%model%urbanVegCategory = urbanVegCategory
     this%model%ixProgress = ixProgress
     this%model%ixRestart = ixRestart
     this%model%newOutputFile = newOutputFile
#ifndef NGEN_FORCING_ACTIVE
     this%model%ixHRUfile_min = ixHRUfile_min
     this%model%ixHRUfile_max = ixHRUfile_max
     this%model%forcFileInfo = forcFileInfo
     this%model%iFile = iFile
     this%model%forcingStep = forcingStep
     this%model%forcNcid = forcNcid
#endif
     ! update global variables in the model structure that change during the model simulation
     this%model%timeStep = 1
     this%model%oldTime = oldTime
     this%model%outFreq = outFreq
     this%model%ncid = ncid
     this%model%elapsedInit = elapsedInit
     this%model%elapsedSetup = elapsedSetup
     this%model%elapsedRestart = elapsedRestart
     this%model%elapsedRead = elapsedRead
     this%model%elapsedWrite = elapsedWrite
     this%model%elapsedPhysics = elapsedPhysics
     bmi_status = BMI_SUCCESS
   end function summa_bmi_initialize

   ! *****************************************************************************
   ! * advance model by one time step.
   ! *****************************************************************************
   function summa_update(this) result (bmi_status)
     class (summa_bmi), intent(inout) :: this
     ! error control
     integer(i4b)                       :: err=0                      ! error code
     character(len=1024)                :: message=''                 ! error message
     integer :: bmi_status

     ! get global variables that are constants throughout the model simulation
     gru_struc = this%model%gru_struc
     index_map = this%model%index_map
     model_decisions = this%model%model_decisions
     output_fileSuffix = this%model%output_fileSuffix
     newOutputFile = this%model%newOutputFile
     outFreq = this%model%outFreq
     maxLayers = this%model%maxLayers
     maxSnowLayers = this%model%maxSnowLayers
     maxLakeLayers = this%model%maxLakeLayers
     maxGlceLayers = this%model%maxGlceLayers
     maxSoilLayers = this%model%maxSoilLayers
     maxGlaciers = this%model%maxGlaciers
     maxWetlands = this%model%maxWetlands
     maxGrid = this%model%maxGrid
     maxGridX = this%model%maxGridX
     maxGridY = this%model%maxGridY
     nHRUrun = this%model%nHRUrun
     nGRUrun = this%model%nGRUrun
     maxDOM = this%model%maxDOM
     urbanVegCategory = this%model%urbanVegCategory
     ixProgress = this%model%ixProgress
     ixRestart = this%model%ixRestart
#ifndef NGEN_FORCING_ACTIVE
     ixHRUfile_min = this%model%ixHRUfile_min
     ixHRUfile_max = this%model%ixHRUfile_max
     forcFileInfo = this%model%forcFileInfo
     iFile = this%model%iFile
     forcingStep = this%model%forcingStep
     forcNcid = this%model%forcNcid
#endif
     ! initialize global variables that change during the model simulation
     fileout = this%model%fileout
     ncid = this%model%ncid
     statCounter = this%model%statCounter
     outputTimeStep = this%model%outputTimeStep
     resetStats = this%model%resetStats
     finalizeStats = this%model%finalizeStats
     oldTime = this%model%oldTime
     elapsedRead = this%model%elapsedRead
     elapsedWrite = this%model%elapsedWrite
     elapsedPhysics = this%model%elapsedPhysics

     ! read model forcing data
     call summa_readForcing(this%model%timeStep, this%model%summa1_struc(n), err, message)
     call handle_err(err, message)

#ifndef NGEN_ACTIVE
     if (mod(this%model%timeStep, print_step_freq) == 0)then
       print *, 'step ---> ', this%model%timeStep
     endif
#endif
     ! run the summa physics for one time step
     call summa_runPhysics(this%model%timeStep, this%model%summa1_struc(n), err, message)
     call handle_err(err, message)

     ! write the model output
#ifndef NGEN_OUTPUT_ACTIVE
     call summa_writeOutputFiles(this%model%timeStep, this%model%summa1_struc(n), err, message)
     call handle_err(err, message)
#endif

     ! update global variables that change during the model simulation
     this%model%timeStep = this%model%timeStep + 1
     this%model%fileout = fileout
     this%model%ncid = ncid
     this%model%nHRUrun = nHRUrun
     this%model%nGRUrun = nGRUrun
     this%model%maxDOM = maxDOM
#ifndef NGEN_FORCING_ACTIVE
     this%model%nHRUrun = nHRUfile
     this%model%iFile = iFile
     this%model%forcingStep = forcingStep
     this%model%forcNcid = forcNcid
     this%model%forcFileInfo = forcFileInfo
#endif
     this%model%statCounter = statCounter
     this%model%outputTimeStep = outputTimeStep
     this%model%resetStats = resetStats
     this%model%finalizeStats = finalizeStats
     this%model%oldTime = oldTime
     this%model%elapsedRead = elapsedRead
     this%model%elapsedWrite = elapsedWrite
     this%model%elapsedPhysics = elapsedPhysics
     bmi_status = BMI_SUCCESS
   end function summa_update

   ! ****************************************************************************
   ! * advance the model until the given time
   ! ****************************************************************************
   function summa_update_until(this, time) result (bmi_status)
     class (summa_bmi), intent(inout) :: this
     double precision, intent(in) :: time
     integer :: bmi_status, istat, n_steps, i
     double precision :: current

     istat = this%get_current_time(current) ! unit seconds
     if (time < current) then
       bmi_status = BMI_FAILURE
       return
     end if

     n_steps = nint( (time - current)/data_step ) + 1 ! model can only do a full data_step
     ! SUMMA runs the ending step (so start=end would still run a step)
     do i = 1, n_steps
       istat = this%update()
     end do
     bmi_status = BMI_SUCCESS
   end function summa_update_until

   ! *****************************************************************************
   ! * successful end (NOTE, NGEN does not call finalize)
   ! *****************************************************************************
   function summa_finalize(this) result (bmi_status)
     class (summa_bmi), intent(inout) :: this
     integer :: bmi_status

     ! get global variables that change during the model simulation
     ncid = this%model%ncid
     elapsedInit = this%model%elapsedInit
     elapsedSetup = this%model%elapsedSetup
     elapsedRestart = this%model%elapsedRestart
     elapsedRead = this%model%elapsedRead
     elapsedWrite = this%model%elapsedWrite
     elapsedPhysics = this%model%elapsedPhysics

     call stop_program(0, 'finished simulation successfully.')
     ! to prevent exiting before HDF5 has closed
     call sleep(2)
     bmi_status = BMI_SUCCESS
   end function summa_finalize

   ! *****************************************************************************
   ! * extra BMI functions
   ! *****************************************************************************

   ! Get the name of the model
   function summa_component_name(this, name) result (bmi_status)
     class (summa_bmi), intent(in) :: this
     character (len=*), pointer, intent(out) :: name
     integer :: bmi_status

     name => component_name
     bmi_status = BMI_SUCCESS
   end function summa_component_name

   ! Count the input variables
   function summa_input_item_count(this, count) result (bmi_status)
     class (summa_bmi), intent(in) :: this
     integer, intent(out) :: count
     integer :: bmi_status

     count = input_item_count
     bmi_status = BMI_SUCCESS
   end function summa_input_item_count

   ! Count the output variables
   function summa_output_item_count(this, count) result (bmi_status)
     class (summa_bmi), intent(in) :: this
     integer, intent(out) :: count
     integer :: bmi_status

     count = output_item_count
     bmi_status = BMI_SUCCESS
   end function summa_output_item_count

   ! List output variables standardized as "https://csdms.colorado.edu/wiki/CSDMS_Standard_Names"
   ! These are the inputs we will need if we do not want to call read_force inside summa_forcing.f90
   ! NGEN uses two component wind and a time vector that is not currently separable
   !   (compute wind speed from the two components and time from start time and hourly step assumption)
   function summa_input_var_names(this, names) result (bmi_status)
     class (summa_bmi), intent(in) :: this
     character (*), pointer, intent(out) :: names(:)
     integer :: bmi_status

     input_items(1) = 'atmosphere_water__precipitation_mass_flux'
     input_items(2) = 'land_surface_air__temperature'
     input_items(3) = 'atmosphere_air_water~vapor__relative_saturation'
#ifdef NGEN_ACTIVE
     input_items(4) = 'land_surface_wind__x_component_of_velocity'
     input_items(8) = 'land_surface_wind__y_component_of_velocity'
#else
     input_items(4) = 'land_surface_wind__speed'
#endif
     input_items(5) = 'land_surface_radiation~incoming~shortwave__energy_flux'
     input_items(6) = 'land_surface_radiation~incoming~longwave__energy_flux'
     input_items(7) = 'land_surface_air__pressure'

     names => input_items
     bmi_status = BMI_SUCCESS
   end function summa_input_var_names

   ! List output variables standardized as "https://csdms.colorado.edu/wiki/CSDMS_Standard_Names"
   function summa_output_var_names(this, names) result (bmi_status)
     class (summa_bmi), intent(in) :: this
     character (*), pointer, intent(out) :: names(:)
     integer :: bmi_status, i

     output_items(1) = 'land_surface_water__runoff_volume_flux'
     output_items(2) = 'land_surface_water__evaporation_mass_flux'
     output_items(3) = 'land_vegetation_water__evaporation_mass_flux'
     output_items(4) = 'land_vegetation_water__transpiration_mass_flux'
     output_items(5) = 'snowpack__sublimation_mass_flux'
     output_items(6) = 'land_vegetation_water__sublimation_mass_flux'
     output_items(7) = 'snowpack_mass'
     output_items(8) = 'soil_water__mass'
     output_items(9) = 'land_vegetation_water__mass'
     output_items(10)= 'land_surface_radiation~net~total__energy_flux'
     output_items(11)= 'land_atmosphere_heat~net~latent__energy_flux'   !(incoming to the *atmosphere*, since atmosphere is last)
     output_items(12)= 'land_atmosphere_heat~net~sensible__energy_flux' !(incoming to the *atmosphere*, since atmosphere is last)
     output_items(13)= 'atmosphere_energy~net~total__energy_flux'
     output_items(14)= 'land_vegetation_energy~net~total__energy_flux'
     output_items(15)= 'land_surface_energy~net~total__energy_flux'
     output_items(16)= 'land_surface_water__baseflow_volume_flux'
     names => output_items
     bmi_status = BMI_SUCCESS
   end function summa_output_var_names

   ! Model start time
   function summa_start_time(this, time) result (bmi_status)
     class (summa_bmi), intent(in) :: this
     double precision, intent(out) :: time
     integer :: bmi_status

     time = 0.0 ! unit seconds
     bmi_status = BMI_SUCCESS
   end function summa_start_time

   ! Model end time
   function summa_end_time(this, time) result (bmi_status)
     class (summa_bmi), intent(in) :: this
     double precision, intent(out) :: time
     integer :: bmi_status

     time = (dJulianFinsh - dJulianStart)*secprday ! unit seconds
     bmi_status = BMI_SUCCESS
   end function summa_end_time

   ! Model current time
   function summa_current_time(this, time) result (bmi_status)
     class (summa_bmi), intent(in) :: this
     double precision, intent(out) :: time
     integer :: bmi_status

     if(this%model%timeStep==0)then
       time = 0.0 ! unit seconds
     else
       time = (data_step*real(this%model%timeStep-1,dp)) ! unit seconds
     end if
     bmi_status = BMI_SUCCESS
   end function summa_current_time

   ! Model time step
   function summa_time_step(this, time_step) result (bmi_status)
     class (summa_bmi), intent(in) :: this
     double precision, intent(out) :: time_step
     integer :: bmi_status

     time_step = data_step  ! unit seconds
     bmi_status = BMI_SUCCESS
   end function summa_time_step

   ! Model time units
   function summa_time_units(this, units) result (bmi_status)
     class (summa_bmi), intent(in) :: this
     character (len=*), intent(out) :: units
     integer :: bmi_status

     units = "s"
     bmi_status = BMI_SUCCESS
   end function summa_time_units

   ! Get the grid id for a particular variable
   function summa_var_grid(this, name, grid) result (bmi_status)
     class (summa_bmi), intent(in) :: this
     character (len=*), intent(in) :: name
     integer, intent(out) :: grid
     integer :: bmi_status

     select case(name)
     case default
       grid = 0
       bmi_status = BMI_SUCCESS
     end select
   end function summa_var_grid

   ! The type of a variable's grid
   function summa_grid_type(this, grid, type) result (bmi_status)
     class (summa_bmi), intent(in) :: this
     integer, intent(in) :: grid
     character (len=*), intent(out) :: type
     integer :: bmi_status

     select case(grid)
     case(0)
       type = 'points'
       bmi_status = BMI_SUCCESS
     case default
       type = "-"
       bmi_status = BMI_FAILURE
     end select
   end function summa_grid_type

   ! The number of dimensions of a grid, latitude and longitude and elevation
   function summa_grid_rank(this, grid, rank) result (bmi_status)
     class (summa_bmi), intent(in) :: this
     integer, intent(in) :: grid
     integer, intent(out) :: rank
     integer :: bmi_status

     select case(grid)
     case(0)
       rank = 3
       bmi_status = BMI_SUCCESS
     case default
       rank = -1
       bmi_status = BMI_FAILURE
     end select
   end function summa_grid_rank

   ! The dimensions of a grid, not applicable to unstructured
   function summa_grid_shape(this, grid, shape) result (bmi_status)
     class (summa_bmi), intent(in) :: this
     integer, intent(in) :: grid
     integer, dimension(:), intent(out) :: shape
     integer :: bmi_status

     select case(grid)
     case default
       shape(:) = -1
       bmi_status = BMI_FAILURE
     end select
   end function summa_grid_shape

   ! The total number of elements in a grid
   function summa_grid_size(this, grid, size) result (bmi_status)
     class (summa_bmi), intent(in) :: this
     integer, intent(in) :: grid
     integer, intent(out) :: size
     integer :: bmi_status

     select case(grid)
     case(0)
       size = sum(gru_struc(:)%hruCount)
       bmi_status = BMI_SUCCESS
     case default
       size = -1
       bmi_status = BMI_FAILURE
     end select
   end function summa_grid_size

   ! The distance between nodes of a grid, not applicable to unstructured
   function summa_grid_spacing(this, grid, spacing) result (bmi_status)
     class (summa_bmi), intent(in) :: this
     integer, intent(in) :: grid
     double precision, dimension(:), intent(out) :: spacing
     integer :: bmi_status

     select case(grid)
     case default
       spacing(:) = -1.d0
       bmi_status = BMI_FAILURE
     end select
   end function summa_grid_spacing

   ! Coordinates of grid origin, not applicable to unstructured
   function summa_grid_origin(this, grid, origin) result (bmi_status)
     class (summa_bmi), intent(in) :: this
     integer, intent(in) :: grid
     double precision, dimension(:), intent(out) :: origin
     integer :: bmi_status

     select case(grid)
     case default
       origin(:) = -1.d0
       bmi_status = BMI_FAILURE
     end select
   end function summa_grid_origin

   ! X-coordinates of grid nodes, longitude (degrees east)
   function summa_grid_x(this, grid, x) result (bmi_status)
     class (summa_bmi), intent(in) :: this
     integer, intent(in) :: grid
     double precision, dimension(:), intent(out) :: x
     integer :: bmi_status, iGRU, jHRU

     summaVars: associate(attrStruct => this%model%summa1_struc(n)%attrStruct    & ! x%gru(:)%hru(:)%var(:)     -- local attributes for each HRU
      )
      select case(grid)
      case default
        do iGRU = 1, this%model%summa1_struc(n)%nGRU
          do jHRU = 1, gru_struc(iGRU)%hruCount
            x((iGRU-1) * gru_struc(iGRU)%hruCount + jHRU) = attrStruct%gru(iGRU)%hru(jHRU)%var(iLookATTR%longitude)
          end do
        end do
        bmi_status = BMI_SUCCESS
      end select
     end associate summaVars
   end function summa_grid_x

   ! Y-coordinates of grid nodes, latitude (degrees north)
   function summa_grid_y(this, grid, y) result (bmi_status)
     class (summa_bmi), intent(in) :: this
     integer, intent(in) :: grid
     double precision, dimension(:), intent(out) :: y
     integer :: bmi_status, iGRU, jHRU

     summaVars: associate(attrStruct => this%model%summa1_struc(n)%attrStruct    & ! x%gru(:)%hru(:)%var(:)     -- local attributes for each HRU
      )
      select case(grid)
      case default
        do iGRU = 1, this%model%summa1_struc(n)%nGRU
          do jHRU = 1, gru_struc(iGRU)%hruCount
            y((iGRU-1) * gru_struc(iGRU)%hruCount + jHRU) = attrStruct%gru(iGRU)%hru(jHRU)%var(iLookATTR%latitude)
          end do
        end do
        bmi_status = BMI_SUCCESS
      end select
     end associate summaVars
   end function summa_grid_y

   ! Z-coordinates of grid nodes, elevation (m)
   function summa_grid_z(this, grid, z) result (bmi_status)
     class (summa_bmi), intent(in) :: this
     integer, intent(in) :: grid
     double precision, dimension(:), intent(out) :: z
     integer :: bmi_status, iGRU, jHRU

     summaVars: associate(attrStruct => this%model%summa1_struc(n)%attrStruct    & ! x%gru(:)%hru(:)%var(:)     -- local attributes for each HRU
      )
      select case(grid)
      case default
        do iGRU = 1, this%model%summa1_struc(n)%nGRU
          do jHRU = 1, gru_struc(iGRU)%hruCount
            z((iGRU-1) * gru_struc(iGRU)%hruCount + jHRU) = attrStruct%gru(iGRU)%hru(jHRU)%var(iLookATTR%elevation)
          end do
        end do
        bmi_status = BMI_SUCCESS
      end select
     end associate summaVars
   end function summa_grid_z

   ! Get the number of nodes in an unstructured grid
   function summa_grid_node_count(this, grid, count) result(bmi_status)
     class(summa_bmi), intent(in) :: this
     integer, intent(in) :: grid
     integer, intent(out) :: count
     integer :: bmi_status

     select case(grid)
     case default
       count = sum(gru_struc(:)%hruCount)
       bmi_status = BMI_SUCCESS
     end select
   end function summa_grid_node_count

   ! Get the number of edges in an unstructured grid, points is 0 BUT USED FOR ROUTING FLOW BETWEEN GRUs
   function summa_grid_edge_count(this, grid, count) result(bmi_status)
     class(summa_bmi), intent(in) :: this
     integer, intent(in) :: grid
     integer, intent(out) :: count
     integer :: bmi_status

     select case(grid)
     case default
       count = 0
       bmi_status = BMI_SUCCESS
     end select
   end function summa_grid_edge_count

   ! Get the number of faces in an unstructured grid, points is 0 BUT COULD BE USED FOR GRUs
   function summa_grid_face_count(this, grid, count) result(bmi_status)
     class(summa_bmi), intent(in) :: this
     integer, intent(in) :: grid
     integer, intent(out) :: count
     integer :: bmi_status

     select case(grid)
     case default
       count = 0 ! could be this%model%summa1_struc(n)%nGRU
       bmi_status = BMI_SUCCESS
     end select
   end function summa_grid_face_count

   ! Get the edge-node connectivity, points is 0 BUT COULD BE USED FOR LATERAL FLOW BETWEEN HRUs
   function summa_grid_edge_nodes(this, grid, edge_nodes) result(bmi_status)
     class(summa_bmi), intent(in) :: this
     integer, intent(in) :: grid
     integer, dimension(:), intent(out) :: edge_nodes
     integer :: bmi_status

     select case(grid)
     case default
       edge_nodes(:) = -1
       bmi_status = BMI_FAILURE
     end select
   end function summa_grid_edge_nodes

  ! Get the face-edge connectivity, points is 0 BUT  COULD BE USED FOR LATERAL FLOW TO GRU
   function summa_grid_face_edges(this, grid, face_edges) result(bmi_status)
     class(summa_bmi), intent(in) :: this
     integer, intent(in) :: grid
     integer, dimension(:), intent(out) :: face_edges
     integer :: bmi_status

     select case(grid)
     case default
       face_edges(:) = -1
       bmi_status = BMI_FAILURE
     end select
   end function summa_grid_face_edges

   ! Get the face-node connectivity, points is 0 BUT COULD BE USED FOR HRU FLOW TO GRU
   function summa_grid_face_nodes(this, grid, face_nodes) result(bmi_status)
     class(summa_bmi), intent(in) :: this
     integer, intent(in) :: grid
     integer, dimension(:), intent(out) :: face_nodes
     integer :: bmi_status

     select case(grid)
     case default
       face_nodes(:) = -1 !gru_struc(iGRU)%hruInfo(iHRU)%hru_id for iHRU=1,hruCount
       bmi_status = BMI_FAILURE
     end select
   end function summa_grid_face_nodes

   ! Get the number of nodes for each face, points is 0 BUT COULD BE USED FOR GRUs
   function summa_grid_nodes_per_face(this, grid, nodes_per_face) result(bmi_status)
     class(summa_bmi), intent(in) :: this
     integer, intent(in) :: grid
     integer, dimension(:), intent(out) :: nodes_per_face
     integer :: bmi_status

     select case(grid)
     case default
       nodes_per_face(:) = -1 ! could be gru_struc(iGRU)%hruCount
       bmi_status = BMI_SUCCESS
     end select
   end function summa_grid_nodes_per_face

   ! The data type of the variable, as a string
   function summa_var_type(this, name, type) result (bmi_status)
     class (summa_bmi), intent(in) :: this
     character (len=*), intent(in) :: name
     character (len=*), intent(out) :: type
     integer :: bmi_status

     if(name(1:5)=='model')then ! not currently used, left in for future integer type needs
       type = "integer"
     else
       type = "real"
     endif
     bmi_status = BMI_SUCCESS
   end function summa_var_type

   ! The units of the given variable
   function summa_var_units(this, name, units) result (bmi_status)
     class (summa_bmi), intent(in) :: this
     character (len=*), intent(in) :: name
     character (len=*), intent(out) :: units
     integer :: bmi_status

     select case (name)
     ! input, note using the NGEN preferred unit definitions, equivalent to standard SUMMA definitions as noted
     case('atmosphere_water__precipitation_mass_flux')              ; units = 'mm s-1';   ; bmi_status = BMI_SUCCESS !equivalent kg m-2 s-1
     case('land_surface_air__temperature')                          ; units = 'K'         ; bmi_status = BMI_SUCCESS
     case('atmosphere_air_water~vapor__relative_saturation')        ; units = 'kg kg-1'   ; bmi_status = BMI_SUCCESS
#ifdef NGEN_ACTIVE
     case('land_surface_wind__x_component_of_velocity')             ; units = 'm s-1'     ; bmi_status = BMI_SUCCESS
     case('land_surface_wind__y_component_of_velocity')             ; units = 'm s-1'     ; bmi_status = BMI_SUCCESS
#else
     case('land_surface_wind__speed')                               ; units = 'm s-1'     ; bmi_status = BMI_SUCCESS
#endif
     case('land_surface_radiation~incoming~shortwave__energy_flux') ; units = 'W m-2'     ; bmi_status = BMI_SUCCESS
     case('land_surface_radiation~incoming~longwave__energy_flux')  ; units = 'W m-2'     ; bmi_status = BMI_SUCCESS
     case('land_surface_air__pressure')                             ; units = 'kg m-1 s-2'; bmi_status = BMI_SUCCESS

     ! output
     case('land_surface_water__runoff_volume_flux')        ; units = 'm s-1'     ; bmi_status = BMI_SUCCESS
     case('land_surface_water__evaporation_mass_flux')     ; units = 'mm s-1'    ; bmi_status = BMI_SUCCESS !equivalent kg m-2 s-1
     case('land_vegetation_water__evaporation_mass_flux')  ; units = 'mm s-1'    ; bmi_status = BMI_SUCCESS !equivalent kg m-2 s-1
     case('land_vegetation_water__transpiration_mass_flux'); units = 'mm s-1'    ; bmi_status = BMI_SUCCESS !equivalent kg m-2 s-1
     case('snowpack__sublimation_mass_flux')               ; units = 'mm s-1'    ; bmi_status = BMI_SUCCESS !equivalent kg m-2 s-1
     case('land_vegetation_water__sublimation_mass_flux')  ; units = 'mm s-1'    ; bmi_status = BMI_SUCCESS !equivalent kg m-2 s-1
     case('snowpack_mass')                                 ; units = 'kg m-2'    ; bmi_status = BMI_SUCCESS
     case('soil_water__mass')                              ; units = 'kg m-2'    ; bmi_status = BMI_SUCCESS
     case('land_vegetation_water__mass')                   ; units = 'kg m-2'    ; bmi_status = BMI_SUCCESS
     case('land_surface_radiation~net~total__energy_flux') ; units = 'W m-2'     ; bmi_status = BMI_SUCCESS
     case('land_atmosphere_heat~net~latent__energy_flux')  ; units = 'W m-2'     ; bmi_status = BMI_SUCCESS
     case('land_atmosphere_heat~net~sensible__energy_flux'); units = 'W m-2'     ; bmi_status = BMI_SUCCESS
     case('atmosphere_energy~net~total__energy_flux')      ; units = 'W m-2'     ; bmi_status = BMI_SUCCESS
     case('land_vegetation_energy~net~total__energy_flux') ; units = 'W m-2'     ; bmi_status = BMI_SUCCESS
     case('land_surface_energy~net~total__energy_flux')    ; units = 'W m-2'     ; bmi_status = BMI_SUCCESS
     case('land_surface_water__baseflow_volume_flux')      ; units = 'm s-1'     ; bmi_status = BMI_SUCCESS
     case default; units = "-"; bmi_status = BMI_FAILURE
     end select
   end function summa_var_units

   ! Memory use per array element
   function summa_var_itemsize(this, name, size) result (bmi_status)
     class (summa_bmi), intent(in) :: this
     character (len=*), intent(in) :: name
     integer, intent(out) :: size
     real, target    :: target_arr(sum(gru_struc(:)%hruCount))
     integer ,target :: itarget_arr
     integer :: bmi_status

     call get_basin_field(this, name, target_arr, itarget_arr) ! See near bottom of file
     ! use the real or integer target
     if(name(1:5)=='model')then ! not currently used, left in for future integer type needs
       size = sizeof(itarget_arr) ! 'sizeof' in gcc & ifort
     else
       size = sizeof(target_arr(1)) ! 'sizeof' in gcc & ifort
     endif
     bmi_status = BMI_SUCCESS
   end function summa_var_itemsize

   ! The size of the given variable
   function summa_var_nbytes(this, name, nbytes) result (bmi_status)
     class (summa_bmi), intent(in) :: this
     character (len=*), intent(in) :: name
     integer, intent(out) :: nbytes
     integer :: bmi_status
     integer :: s1, s2, s3, grid, grid_size, item_size

     s1 = this%get_var_grid(name, grid)
     s2 = this%get_grid_size(grid, grid_size)
     s3 = this%get_var_itemsize(name, item_size)
     if ((s1 == BMI_SUCCESS).and.(s2 == BMI_SUCCESS).and.(s3 == BMI_SUCCESS)) then
       nbytes = item_size * grid_size
       bmi_status = BMI_SUCCESS
     else
       nbytes = -1
       bmi_status = BMI_FAILURE
     end if
   end function summa_var_nbytes

   ! The location (node, face, edge) of the given variable
   function summa_var_location(this, name, location) result (bmi_status)
     class (summa_bmi), intent(in) :: this
     character (len=*), intent(in) :: name
     character (len=*), intent(out) :: location
     integer :: bmi_status

     select case(name)
     case default
        location = "node"
        bmi_status = BMI_SUCCESS
     end select
   end function summa_var_location

   ! Get a copy of a integer variable's values, flattened
   function summa_get_int(this, name, dest) result (bmi_status)
     class (summa_bmi), intent(in) :: this
     character (len=*), intent(in) :: name
     integer, intent(inout) :: dest(:)
     real, target    :: target_arr(sum(gru_struc(:)%hruCount))
     integer ,target :: itarget_arr
     integer :: bmi_status

     select case(name)
     case default
       call get_basin_field(this, name, target_arr, itarget_arr) ! See near bottom of file
       ! use the integer target
       dest = itarget_arr
       bmi_status = BMI_SUCCESS
     end select
   end function summa_get_int

   ! Get a copy of a real variable's values, flattened
   function summa_get_float(this, name, dest) result (bmi_status)
     class (summa_bmi), intent(in) :: this
     character (len=*), intent(in) :: name
     real, intent(inout) :: dest(:)
     real, target    :: target_arr(sum(gru_struc(:)%hruCount))
     integer ,target :: itarget_arr
     integer :: bmi_status

     select case(name)
     case default
       call get_basin_field(this, name, target_arr, itarget_arr) ! See near bottom of file
       ! use the real target
       dest = target_arr
       bmi_status = BMI_SUCCESS
     end select
   end function summa_get_float

   ! Get a copy of a double variable's values, flattened
   function summa_get_double(this, name, dest) result (bmi_status)
     class (summa_bmi), intent(in) :: this
     character (len=*), intent(in) :: name
     double precision, intent(inout) :: dest(:)
     integer :: bmi_status

     select case(name)
     case default
       dest(:) = -1.d0
       bmi_status = BMI_FAILURE
     end select
   end function summa_get_double

   ! Get a reference to an integer-valued variable, flattened
   function summa_get_ptr_int(this, name, dest_ptr) result (bmi_status)
     class (summa_bmi), intent(in) :: this
     character (len=*), intent(in) :: name
     integer, pointer, intent(inout) :: dest_ptr(:)
     integer :: bmi_status
     real    :: target_arr(sum(gru_struc(:)%hruCount))
     type (c_ptr) :: src

     select case(name)
     case default
       call get_basin_field(this, name, target_arr, iptr_buffer) ! See near bottom of file
       ! use the integer target, which is a scalar (the 'model*' fields are not per-HRU)
       src = c_loc(iptr_buffer)
       call c_f_pointer(src, dest_ptr, [1])
       bmi_status = BMI_SUCCESS
     end select
   end function summa_get_ptr_int

   ! Get a reference to a real-valued variable, flattened
   function summa_get_ptr_float(this, name, dest_ptr) result (bmi_status)
     class (summa_bmi), intent(in) :: this
     character (len=*), intent(in) :: name
     real, pointer, intent(inout) :: dest_ptr(:)
     integer :: bmi_status, n_elements
     integer :: itarget_arr
     type (c_ptr) :: src

     select case(name)
     case default
       n_elements = sum(gru_struc(:)%hruCount)
       if (.not. allocated(ptr_buffer)) allocate(ptr_buffer(n_elements))
       call get_basin_field(this, name, ptr_buffer, itarget_arr) ! See near bottom of file
       ! use the real target
       src = c_loc(ptr_buffer(1))
       call c_f_pointer(src, dest_ptr, [n_elements])
       bmi_status = BMI_SUCCESS
     end select
   end function summa_get_ptr_float

   ! Get a reference to an double-valued variable, flattened
   function summa_get_ptr_double(this, name, dest_ptr) result (bmi_status)
     class (summa_bmi), intent(in) :: this
     character (len=*), intent(in) :: name
     double precision, pointer, intent(inout) :: dest_ptr(:)
     integer :: bmi_status, n_elements
     type (c_ptr) :: src

     select case(name)
     case default
       bmi_status = BMI_FAILURE
     end select
   end function summa_get_ptr_double

   ! Get values of an integer variable at the given locations
   function summa_get_at_indices_int(this, name, dest, inds) result (bmi_status)
     class (summa_bmi), intent(in) :: this
     character (len=*), intent(in) :: name
     integer, intent(inout) :: dest(:)
     integer, intent(in) :: inds(:)
     integer :: bmi_status, i, n_elements
     real, target    :: target_arr(sum(gru_struc(:)%hruCount))
     integer ,target :: itarget_arr
     type (c_ptr) src
     integer, pointer :: src_flattened(:)

     select case(name)
     case default
       call get_basin_field(this, name, target_arr, itarget_arr) ! See near bottom of file
       ! use the integer target, which is a scalar (the 'model*' fields are not per-HRU)
       src = c_loc(itarget_arr)
       call c_f_pointer(src, src_flattened, [1])
       n_elements = size(inds)
       do i = 1, n_elements
          dest(i) = src_flattened(inds(i))
       end do
       bmi_status = BMI_SUCCESS
     end select
   end function summa_get_at_indices_int

   ! Get values of a real variable at the given locations
   function summa_get_at_indices_float(this, name, dest, inds) result (bmi_status)
     class (summa_bmi), intent(in) :: this
     character (len=*), intent(in) :: name
     real, intent(inout) :: dest(:)
     integer, intent(in) :: inds(:)
     integer :: bmi_status, i, n_elements
     real, target    :: target_arr(sum(gru_struc(:)%hruCount))
     integer ,target :: itarget_arr
     type (c_ptr) src
     real, pointer :: src_flattened(:)

     select case(name)
     case default
       call get_basin_field(this, name, target_arr, itarget_arr) ! See near bottom of file
       ! use the real target
       src = c_loc(target_arr(1))
       call c_f_pointer(src, src_flattened, [size(target_arr)])
       n_elements = size(inds)
       do i = 1, n_elements
          dest(i) = src_flattened(inds(i))
       end do
       bmi_status = BMI_SUCCESS
     end select
   end function summa_get_at_indices_float

   ! Get values of a double variable at the given locations
   function summa_get_at_indices_double(this, name, dest, inds) result (bmi_status)
     class (summa_bmi), intent(in) :: this
     character (len=*), intent(in) :: name
     double precision, intent(inout) :: dest(:)
     integer, intent(in) :: inds(:)
     integer :: bmi_status, i, n_elements
     type (c_ptr) src
     double precision, pointer :: src_flattened(:)

     select case(name)
     case default
       bmi_status = BMI_FAILURE
     end select
   end function summa_get_at_indices_double

   ! Set new integer values, ONLY FOR INPUT VARIABLES
   function summa_set_int(this, name, src) result (bmi_status)
     class (summa_bmi), intent(inout) :: this
     character (len=*), intent(in) :: name
     integer, intent(in) :: src(:)
     real :: rsrc(sum(gru_struc(:)%hruCount))
     integer :: bmi_status

     select case(name)
     case default
       rsrc = -999.0
       call assign_basin_field(this, name, rsrc, src(1)) ! See near bottom of file
       bmi_status = BMI_SUCCESS
     end select
   end function summa_set_int

   ! Set new real values, ONLY FOR INPUT VARIABLES
   function summa_set_float(this, name, src) result (bmi_status)
     class (summa_bmi), intent(inout) :: this
     character (len=*), intent(in) :: name
     real, intent(in) :: src(:)
     integer :: bmi_status, isrc

     select case(name)
     case default
       isrc = -999
       call assign_basin_field(this, name, src, isrc) ! See near bottom of file
       bmi_status = BMI_SUCCESS
     end select
   end function summa_set_float

   ! Set new double values
   function summa_set_double(this, name, src) result (bmi_status)
     class (summa_bmi), intent(inout) :: this
     character (len=*), intent(in) :: name
     double precision, intent(in) :: src(:)
     integer :: bmi_status

     select case(name)
     case default
       bmi_status = BMI_FAILURE
     end select
   end function summa_set_double

   ! Set integer values at particular locations
   function summa_set_at_indices_int(this, name, inds, src) result (bmi_status)
     class (summa_bmi), intent(inout) :: this
     character (len=*), intent(in) :: name
     integer, intent(in) :: inds(:)
     integer, intent(in) :: src(:)
     integer :: bmi_status
     type (c_ptr) dest
     integer, pointer :: dest_flattened(:)
     integer :: i

     select case(name)
     case default
       bmi_status = BMI_FAILURE
     end select
   end function summa_set_at_indices_int

   ! Set real values at particular locations, ONLY FOR INPUT VARIABLES
   function summa_set_at_indices_float(this, name, inds, src) result (bmi_status)
     class (summa_bmi), intent(inout) :: this
     character (len=*), intent(in) :: name
     integer, intent(in) :: inds(:)
     real, intent(in) :: src(:)
     integer :: bmi_status
     type (c_ptr) dest
     real, pointer :: dest_flattened(:)
     integer :: i

     select case(name)
     case default
       bmi_status = BMI_FAILURE
     end select
   end function summa_set_at_indices_float

   ! Set double values at particular locations
   function summa_set_at_indices_double(this, name, inds, src) result (bmi_status)
     class (summa_bmi), intent(inout) :: this
     character (len=*), intent(in) :: name
     integer, intent(in) :: inds(:)
     double precision, intent(in) :: src(:)
     integer :: bmi_status
     type (c_ptr) dest
     double precision, pointer :: dest_flattened(:)
     integer :: i

     select case(name)
     case default
       bmi_status = BMI_FAILURE
     end select
   end function summa_set_at_indices_double

#ifdef NGEN_ACTIVE
   function register_bmi(this) result(bmi_status) bind(C, name="register_bmi")
     use, intrinsic:: iso_c_binding, only: c_ptr, c_loc, c_int
     use iso_c_bmif_2_0
     implicit none
     type(c_ptr) :: this ! If not value, then from the C perspective `this` is a void**
     integer(kind=c_int) :: bmi_status
     !Create the model instance to use
     type(summa_bmi), pointer :: bmi_model
     !Create a simple pointer wrapper
     type(box), pointer :: bmi_box

     !allocate model
     allocate(summa_bmi::bmi_model)
     !allocate the pointer box
     allocate(bmi_box)

     !associate the wrapper pointer the created model instance
     bmi_box%ptr => bmi_model

     if( .not. associated( bmi_box ) .or. .not. associated( bmi_box%ptr ) ) then
       bmi_status = BMI_FAILURE
     else
       ! Return the pointer to box
       this = c_loc(bmi_box)
       bmi_status = BMI_SUCCESS
     endif
   end function register_bmi
#endif

   ! non-BMI helper function to assign input fields
   subroutine assign_basin_field(this, name, src_arr, isrc_arr)
     implicit none
     class (summa_bmi), intent(inout) :: this
     character (len=*), intent(in) :: name
     real, intent(in)    :: src_arr(sum(gru_struc(:)%hruCount))
     integer, intent(in) :: isrc_arr
     integer ::  iGRU, jHRU, i, iDOM

     summaVars: associate(&
      timeStruct           => this%model%summa1_struc(n)%timeStruct  , & ! x%var(:)                          -- model time data
      forcStruct           => this%model%summa1_struc(n)%forcStruct  , & ! x%gru(:)%hru(:)%var(:)            -- model forcing data
      diagStruct           => this%model%summa1_struc(n)%diagStruct    & ! x%gru(:)%hru(:)%dom(:)%var(:)%dat -- model diagnostic variables
      )

      if(name(1:5)=='model')then ! not currently used, left in for future integer type needs
        select case (name)
        ! input
        case('model__time_year')
          timeStruct%var(iLookTIME%iyyy) = isrc_arr
        end select
      else
        do iGRU = 1, this%model%summa1_struc(n)%nGRU
          do jHRU = 1, gru_struc(iGRU)%hruCount
            i = (iGRU-1) * gru_struc(iGRU)%hruCount + jHRU
            select case (name)
            ! input is same for all domains (for now)
            case('atmosphere_water__precipitation_mass_flux')
              forcStruct%gru(iGRU)%hru(jHRU)%var(iLookFORCE%pptrate) = src_arr(i)
            case('land_surface_air__temperature')
              forcStruct%gru(iGRU)%hru(jHRU)%var(iLookFORCE%airtemp) = src_arr(i)
            case('atmosphere_air_water~vapor__relative_saturation')
              forcStruct%gru(iGRU)%hru(jHRU)%var(iLookFORCE%spechum) = src_arr(i)
            case('land_surface_wind__x_component_of_velocity')
              do iDOM = 1, gru_struc(iGRU)%hruInfo(jHRU)%domCount ! same for all domains for now
                diagStruct%gru(iGRU)%hru(jHRU)%dom(iDOM)%var(iLookDIAG%windspd_x)%dat(1) = src_arr(i)
              end do
            case('land_surface_wind__y_component_of_velocity')
              do iDOM = 1, gru_struc(iGRU)%hruInfo(jHRU)%domCount ! same for all domains for now
                diagStruct%gru(iGRU)%hru(jHRU)%dom(iDOM)%var(iLookDIAG%windspd_y)%dat(1) = src_arr(i)
              end do
            case('land_surface_wind__speed')
              forcStruct%gru(iGRU)%hru(jHRU)%var(iLookFORCE%windspd) = src_arr(i)
            case('land_surface_radiation~incoming~shortwave__energy_flux')
              forcStruct%gru(iGRU)%hru(jHRU)%var(iLookFORCE%SWRadAtm) = src_arr(i)
            case('land_surface_radiation~incoming~longwave__energy_flux')
              forcStruct%gru(iGRU)%hru(jHRU)%var(iLookFORCE%LWRadAtm) = src_arr(i)
            case('land_surface_air__pressure')
              forcStruct%gru(iGRU)%hru(jHRU)%var(iLookFORCE%airpres) = src_arr(i)
            end select
          end do
        end do
      endif
     end associate summaVars
   end subroutine assign_basin_field

   ! non-BMI helper function to get fields, for every HRU in the run.
   ! NOTE: this used to take a do_nHRU argument and fill only the first do_nHRU entries,
   !       leaving the rest at -999.  Every caller but the itemsize probe wants them all,
   !       and all were passing 1, so callers were handed one value followed by -999s.
   subroutine get_basin_field(this, name, target_arr, itarget_arr)
     implicit none
     class (summa_bmi), intent(in) :: this
     character (len=*), intent(in) :: name
     real, target, intent(out)    :: target_arr(sum(gru_struc(:)%hruCount))
     integer, target, intent(out) :: itarget_arr
     integer ::  iGRU, jHRU, i, iDOM
     real :: fracDOM

     summaVars: associate(&
      timeStruct           => this%model%summa1_struc(n)%timeStruct  , & ! x%var(:)                          -- model time data
      forcStruct           => this%model%summa1_struc(n)%forcStruct  , & ! x%gru(:)%hru(:)%var(:)            -- model forcing data
      progStruct           => this%model%summa1_struc(n)%progStruct  , & ! x%gru(:)%hru(:)%dom(:)%var(:)%dat -- model prognostic (state) variables
      diagStruct           => this%model%summa1_struc(n)%diagStruct  , & ! x%gru(:)%hru(:)%dom(:)%var(:)%dat -- model diagnostic variables
      fluxStruct           => this%model%summa1_struc(n)%fluxStruct  , & ! x%gru(:)%hru(:)%dom(:)%var(:)%dat -- model fluxes
      bvarStruct           => this%model%summa1_struc(n)%bvarStruct    & ! x%gru(:)%var(:)%dat               -- model basin (GRU) variables
      )
      target_arr = -999.0
      itarget_arr = -999
      if(name(1:5)=='model')then ! not currently used, left in for future integer type needs
        select case (name)
        ! input
        case('model__time_year')
          itarget_arr = timeStruct%var(iLookTIME%iyyy)
        end select
      else
        do iGRU = 1, this%model%summa1_struc(n)%nGRU
          do jHRU = 1, gru_struc(iGRU)%hruCount
            i = (iGRU-1) * gru_struc(iGRU)%hruCount + jHRU
            target_arr(i) = 0._rkind
            do iDOM = 1, gru_struc(iGRU)%hruInfo(jHRU)%domCount
              fracDOM = progStruct%gru(iGRU)%hru(jHRU)%dom(iDOM)%var(iLookPROG%DOMarea)%dat(1)/ bvarStruct%gru(iGRU)%var(iLookBVAR%basin__totalArea)%dat(1)

              select case (name)
              ! input is same for all domains (for now)
              case('atmosphere_water__precipitation_mass_flux')
                target_arr(i) = forcStruct%gru(iGRU)%hru(jHRU)%var(iLookFORCE%pptrate); exit
              case('land_surface_air__temperature')
                target_arr(i) = forcStruct%gru(iGRU)%hru(jHRU)%var(iLookFORCE%airtemp); exit
              case('atmosphere_air_water~vapor__relative_saturation')
                target_arr(i) = forcStruct%gru(iGRU)%hru(jHRU)%var(iLookFORCE%spechum); exit
              case('land_surface_wind__x_component_of_velocity')
                target_arr(i) = diagStruct%gru(iGRU)%hru(jHRU)%dom(1)%var(iLookDIAG%windspd_x)%dat(1); exit ! same for all domains for now
              case('land_surface_wind__y_component_of_velocity')
                target_arr(i) = diagStruct%gru(iGRU)%hru(jHRU)%dom(1)%var(iLookDIAG%windspd_y)%dat(1); exit ! same for all domains for now
              case('land_surface_wind__speed')
                target_arr(i) = forcStruct%gru(iGRU)%hru(jHRU)%var(iLookFORCE%windspd); exit
              case('land_surface_radiation~incoming~shortwave__energy_flux')
                target_arr(i) = forcStruct%gru(iGRU)%hru(jHRU)%var(iLookFORCE%SWRadAtm); exit
              case('land_surface_radiation~incoming~longwave__energy_flux')
                target_arr(i) = forcStruct%gru(iGRU)%hru(jHRU)%var(iLookFORCE%LWRadAtm); exit
              case('land_surface_air__pressure')
                target_arr(i) = forcStruct%gru(iGRU)%hru(jHRU)%var(iLookFORCE%airpres); exit

              ! output is averaged over domains
              case('land_surface_water__runoff_volume_flux')
                target_arr(i) = bvarStruct%gru(iGRU)%var(iLookBVAR%averageRoutedRunoff)%dat(1); exit
              case('land_surface_water__evaporation_mass_flux')
                target_arr(i) = target_arr(i) + fluxStruct%gru(iGRU)%hru(jHRU)%dom(iDOM)%var(iLookFLUX%scalarGroundEvaporation)%dat(1) * fracDOM
              case('land_vegetation_water__evaporation_mass_flux')
                target_arr(i) = target_arr(i) + fluxStruct%gru(iGRU)%hru(jHRU)%dom(iDOM)%var(iLookFLUX%scalarCanopyEvaporation)%dat(1) * fracDOM
              case('land_vegetation_water__transpiration_mass_flux')
                target_arr(i) = target_arr(i) + fluxStruct%gru(iGRU)%hru(jHRU)%dom(iDOM)%var(iLookFLUX%scalarCanopyTranspiration)%dat(1) * fracDOM
              case('snowpack__sublimation_mass_flux')
                target_arr(i) = target_arr(i) + fluxStruct%gru(iGRU)%hru(jHRU)%dom(iDOM)%var(iLookFLUX%scalarGroundSublimation)%dat(1) * fracDOM
              case('land_vegetation_water__sublimation_mass_flux')
                target_arr(i) = target_arr(i) + fluxStruct%gru(iGRU)%hru(jHRU)%dom(iDOM)%var(iLookFLUX%scalarCanopySublimation)%dat(1) * fracDOM
              case('snowpack_mass')
                target_arr(i) = target_arr(i) + progStruct%gru(iGRU)%hru(jHRU)%dom(iDOM)%var(iLookPROG%scalarSWE)%dat(1) * fracDOM
              case('soil_water__mass')
                target_arr(i) = target_arr(i) + diagStruct%gru(iGRU)%hru(jHRU)%dom(iDOM)%var(iLookDIAG%scalarTotalSoilWat)%dat(1) * fracDOM
              case('land_vegetation_water__mass')
                target_arr(i) = target_arr(i) + progStruct%gru(iGRU)%hru(jHRU)%dom(iDOM)%var(iLookPROG%scalarCanopyWat)%dat(1) * fracDOM
              case('land_surface_radiation~net~total__energy_flux')
                target_arr(i) = target_arr(i) + fluxStruct%gru(iGRU)%hru(jHRU)%dom(iDOM)%var(iLookFLUX%scalarNetRadiation)%dat(1) * fracDOM
              case('land_atmosphere_heat~net~latent__energy_flux')
                target_arr(i) = target_arr(i) + fluxStruct%gru(iGRU)%hru(jHRU)%dom(iDOM)%var(iLookFLUX%scalarLatHeatTotal)%dat(1) * fracDOM
              case('land_atmosphere_heat~net~sensible__energy_flux')
                target_arr(i) = target_arr(i) + fluxStruct%gru(iGRU)%hru(jHRU)%dom(iDOM)%var(iLookFLUX%scalarSenHeatTotal)%dat(1) * fracDOM
              case('atmosphere_energy~net~total__energy_flux')
                target_arr(i) = target_arr(i) + fluxStruct%gru(iGRU)%hru(jHRU)%dom(iDOM)%var(iLookFLUX%scalarCanairNetNrgFlux)%dat(1) * fracDOM
              case('land_vegetation_energy~net~total__energy_flux')
                target_arr(i) = target_arr(i) + fluxStruct%gru(iGRU)%hru(jHRU)%dom(iDOM)%var(iLookFLUX%scalarCanopyNetNrgFlux)%dat(1) * fracDOM
              case('land_surface_energy~net~total__energy_flux')
                target_arr(i) = target_arr(i) + fluxStruct%gru(iGRU)%hru(jHRU)%dom(iDOM)%var(iLookFLUX%scalarGroundNetNrgFlux)%dat(1) * fracDOM
              case('land_surface_water__baseflow_volume_flux')
                target_arr(i) = target_arr(i) + fluxStruct%gru(iGRU)%hru(jHRU)%dom(iDOM)%var(iLookFLUX%scalarAquiferBaseflow)%dat(1) * fracDOM
              end select
            end do
          end do
        end do
      endif
     end associate summaVars
   end subroutine get_basin_field

end module summabmi

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
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License
! along with this program. If not, see <http://www.gnu.org/licenses/>.


! **************************************************************************************************
! SUMMA optimization driver
!
! Performs a common cold-start spinup and dynamically evaluates parameter
! samples across independent MPI model instances.
! **************************************************************************************************

program summa_driver_opt

  ! data types
  USE nr_type,    only: i4b,rkind,lgt
  USE summa_type, only: config_info
  USE summa_type, only: parallel_context_type

  ! logging
  USE globalData,      only: iulog
  USE iso_fortran_env, only: error_unit
  USE iso_fortran_env, only: output_unit

  ! parameter sampling
  USE summa_parameter_sampling, only: initialize_parameter_evaluation
  USE summa_parameter_sampling, only: dispatch_parameter_samples

  ! MPI
  USE mpi, only: MPI_Init,MPI_Finalize
  USE mpi, only: MPI_Bcast
  
  USE mpi, only: MPI_Comm_split_type
  USE mpi, only: MPI_Comm_split

  USE mpi, only: MPI_COMM_TYPE_SHARED
  USE mpi, only: MPI_INFO_NULL
  USE mpi, only: MPI_UNDEFINED

  USE mpi, only: MPI_INTEGER
  USE mpi, only: MPI_COMM_WORLD,MPI_COMM_SELF
  USE mpi, only: MPI_SUCCESS
  
  USE mpi_context, only: set_mpi_context
  USE error_utils, only: check_mpi,abort_mpi

  ! SUMMA subroutines/functions
  USE summa_config, only: read_manifest
  USE summa_util,   only: getCommandArguments
  USE summa_util,   only: stop_program

  implicit none

  type(config_info) :: config                              ! SUMMA configuration information

  type(parallel_context_type) :: world_parallel            ! all MPI ranks
  type(parallel_context_type) :: node_parallel             ! all MPI ranks on the same physical node
  type(parallel_context_type) :: leader_parallel           ! one leader rank from each physical node
  type(parallel_context_type) :: instance_parallel         ! ranks assigned to one calibration case
  type(parallel_context_type) :: domain_parallel           ! SUMMA domain-parallel context (one instance per core)

  integer(i4b) :: node_index                               ! Index of the physical node
  integer(i4b) :: nNodes                                   ! Number of physical nodes
  integer(i4b) :: case_group                               ! Case group on the local node
  integer(i4b) :: global_case_group                        ! Case group across all nodes
  integer(i4b) :: nCaseGroups                              ! Total number of case groups
  integer(i4b) :: ranks_per_case                           ! MPI ranks assigned to each case
  integer(i4b) :: leader_color                             ! Color used to construct node-leader communicator

  integer(i4b) :: iCase                                    ! Index of the current case
  integer(i4b) :: nCases                                   ! Total number of cases available for execution
  integer(i4b) :: first_case                               ! First case assigned to this case group
  integer(i4b) :: case_stride                              ! Interval between cases assigned to this case group

  integer(i4b),     parameter :: nSamples=5000             ! total number of parameter samples

  integer(i4b)        :: err=0                             ! SUMMA error code
  integer(i4b)        :: mpi_err=0                         ! MPI error code
  character(len=1024) :: message=''                        ! SUMMA error message
  character(len=256)  :: mpi_message=''                    ! MPI error message

  ! ---------------------------------------------------------------------------------------
  ! Initialize MPI
  ! ---------------------------------------------------------------------------------------

  call MPI_Init(mpi_err)
  call check_mpi(-1,mpi_err,'MPI_Init failed')

  ! establish the global MPI context
  world_parallel%comm=MPI_COMM_WORLD

  call set_mpi_context(world_parallel%comm,  & ! global MPI communicator
                       world_parallel%rank,  & ! rank in global communicator
                       world_parallel%size,  & ! number of global MPI ranks
                       mpi_err,mpi_message)    ! MPI error information

  if(mpi_err/=MPI_SUCCESS) &
    call abort_mpi(world_parallel%rank,trim(mpi_message))

  ! ---------------------------------------------------------------------------------------
  ! Read run configuration
  ! ---------------------------------------------------------------------------------------

  ! process command-line arguments once before configuring individual cases
  call getCommandArguments(config,err,message)
  if(err/=0) call abort_mpi(world_parallel%rank,trim(message))

  ! prevent command-line arguments from being reread during case initialization
  config%read_cli = .false.

  ! read the multi-case manifest when specified on the command line
  if(allocated(config%manifest_file))then
    call read_manifest(config%manifest_file,config,err,message)
    if(err/=0) call abort_mpi(world_parallel%rank,trim(message))
  endif

  ! ---------------------------------------------------------------------------------------
  ! Build the MPI communicator hierarchy
  ! ---------------------------------------------------------------------------------------

  ! ----- single-case run: use all available MPI ranks for one calibration -----

  if(.not.allocated(config%manifest_file))then

    instance_parallel%comm=MPI_COMM_WORLD

    call set_mpi_context(instance_parallel%comm,  &
                         instance_parallel%rank,  &
                         instance_parallel%size,  &
                         mpi_err,mpi_message)

    if(mpi_err/=MPI_SUCCESS) &
      call abort_mpi(world_parallel%rank,trim(mpi_message))

  ! ----- multi-case run: partition MPI ranks among independent calibrations -----

  else

    ! -------------------------------------------------------------------------
    ! Construct node-local communicator
    ! -------------------------------------------------------------------------

    ! construct a node-local communicator containing only ranks on the same physical node
    call MPI_Comm_split_type(world_parallel%comm,    & ! parent communicator
                             MPI_COMM_TYPE_SHARED,    & ! group ranks sharing physical memory
                             world_parallel%rank,     & ! rank ordering key
                             MPI_INFO_NULL,           & ! no additional MPI information
                             node_parallel%comm,      & ! node-local communicator
                             mpi_err)                   ! MPI error code

    call check_mpi(world_parallel%rank,mpi_err, &
                   'unable to create node-local communicator')

    ! establish rank and size within the physical node
    call set_mpi_context(node_parallel%comm,  & ! node-local MPI communicator
                         node_parallel%rank,  & ! rank within the physical node
                         node_parallel%size,  & ! number of MPI ranks on the physical node
                         mpi_err,mpi_message)   ! MPI error code and message

    if(mpi_err/=MPI_SUCCESS) &
      call abort_mpi(world_parallel%rank,trim(mpi_message))

    ! -------------------------------------------------------------------------
    ! Identify physical nodes
    ! -------------------------------------------------------------------------

    ! include only node-local rank zero in the node-leader communicator
    if(node_parallel%rank==0)then
      leader_color=0
    else
      leader_color=MPI_UNDEFINED
    endif

    call MPI_Comm_split(world_parallel%comm,     &
                        leader_color,            &
                        world_parallel%rank,      &
                        leader_parallel%comm,    &
                        mpi_err)

    call check_mpi(world_parallel%rank,mpi_err, &
                   'unable to create node-leader communicator')

    ! node leaders determine the node index and total number of nodes
    if(node_parallel%rank==0)then

      call set_mpi_context(leader_parallel%comm,  &
                           leader_parallel%rank,  &
                           leader_parallel%size,  &
                           mpi_err,mpi_message)

      if(mpi_err/=MPI_SUCCESS) &
        call abort_mpi(world_parallel%rank,trim(mpi_message))

      node_index=leader_parallel%rank
      nNodes=leader_parallel%size

    endif

    ! distribute node information to all ranks on the physical node
    call MPI_Bcast(node_index,1,MPI_INTEGER,0,node_parallel%comm,mpi_err)
    call check_mpi(world_parallel%rank,mpi_err, &
                   'unable to broadcast node index')

    call MPI_Bcast(nNodes,1,MPI_INTEGER,0,node_parallel%comm,mpi_err)
    call check_mpi(world_parallel%rank,mpi_err, &
                   'unable to broadcast number of nodes')

    ! -------------------------------------------------------------------------
    ! Partition ranks on each node among independent cases
    ! -------------------------------------------------------------------------

    ! require equal-sized calibration groups on each node
    if(mod(node_parallel%size,config%cases_per_node)/=0)then
      write(message,'(A,I0,A,I0,A)')                                  &
        'number of MPI ranks on node (',node_parallel%size,            &
        ') must be divisible by cases_per_node (',config%cases_per_node,')'
      call abort_mpi(world_parallel%rank,trim(message))
    endif

    ! determine the number of ranks assigned to each independent calibration
    ranks_per_case=node_parallel%size/config%cases_per_node

    ! assign contiguous node-local ranks to calibration groups
    case_group=node_parallel%rank/ranks_per_case

    ! identify this calibration group across all physical nodes
    global_case_group=node_index*config%cases_per_node + case_group

    ! determine the total number of calibration groups across all nodes
    nCaseGroups=nNodes*config%cases_per_node

    ! split the node communicator into independent case communicators
    call MPI_Comm_split(node_parallel%comm,      & ! parent node-local communicator
                        case_group,              & ! calibration group identifier
                        node_parallel%rank,      & ! preserve node-local rank ordering
                        instance_parallel%comm,  & ! communicator for one calibration
                        mpi_err)                   ! MPI error code

    call check_mpi(world_parallel%rank,mpi_err, &
                   'unable to create case communicator')

    ! establish rank and size within the assigned calibration group
    call set_mpi_context(instance_parallel%comm,  & ! case-specific MPI communicator
                         instance_parallel%rank,  & ! rank within the calibration group
                         instance_parallel%size,  & ! number of MPI ranks in the calibration group
                         mpi_err,mpi_message)       ! MPI error code and message

    if(mpi_err/=MPI_SUCCESS) &
      call abort_mpi(world_parallel%rank,trim(mpi_message))

  endif
 
  ! dynamic parameter evaluation requires one dispatcher and at least one worker
  if(instance_parallel%size < 2)then
    call abort_mpi(world_parallel%rank, &
                   'parallel parameter evaluation requires at least two MPI ranks per case')
  endif

  ! ---------------------------------------------------------------------------------------
  ! Define the domain-parallel communicator
  ! ---------------------------------------------------------------------------------------

  ! each MPI rank independently evaluates the complete spatial domain
  domain_parallel%comm=MPI_COMM_SELF

  call set_mpi_context(domain_parallel%comm,  & ! single-rank MPI communicator
                       domain_parallel%rank,  & ! rank within the domain communicator
                       domain_parallel%size,  & ! number of ranks in the domain communicator
                       mpi_err,mpi_message)      ! MPI error code and message

  if(mpi_err/=MPI_SUCCESS) &
    call abort_mpi(world_parallel%rank,trim(mpi_message))
 
  ! ---------------------------------------------------------------------------------------
  ! Log Processor layout 
  ! ---------------------------------------------------------------------------------------

  if(allocated(config%manifest_file))then

    ! report the MPI rank and calibration-group assignment
    write(*,'(A,I4,A,I3,A,I3,A,I3,A,I3,A,I3,A,I3)') &
      '  rank = '          , world_parallel%rank,    &
      ', node = '          , node_index,             &
      ', node_rank = '     , node_parallel%rank,     &
      ', case_group = '    , case_group,             &
      ', instance_rank = ' , instance_parallel%rank, &
      ', global_group = '  , global_case_group,      &
      ', nCaseGroups = '   , nCaseGroups
  
  endif

  ! ---------------------------------------------------------------------------------------
  ! Define case execution loop
  ! ---------------------------------------------------------------------------------------

  if(allocated(config%manifest_file))then

    ! multi-case run: assign each case group a subset of cases from the manifest
    nCases      = size(config%case_names)
    first_case  = global_case_group+1
    case_stride = nCaseGroups

  else

    ! single-case run: execute the case defined by the standard configuration
    nCases      = 1
    first_case  = 1
    case_stride = 1

  endif

  ! ---------------------------------------------------------------------------------------
  ! Run assigned SUMMA cases
  ! ---------------------------------------------------------------------------------------

  ! process the cases assigned to this case group sequentially
  do iCase=first_case,nCases,case_stride

    ! set the case name and config file from the manifest for multi-case runs
    if(allocated(config%manifest_file))then
      config%manifest_casename=trim(config%case_names(iCase))
      config%config_file=trim(config%template_path)//trim(config%template_file)
      if(instance_parallel%rank==0) write(output_unit,'(A)') 'Running case: '//trim(config%manifest_casename)
    endif
  
    ! initialize and execute parameter calibration for the current SUMMA case
    call run_case(config,                    & ! SUMMA configuration structure
                  domain_parallel,           & ! MPI context for domain parallelism
                  instance_parallel,         & ! MPI context for model-instance parallelism
                  nSamples,                  & ! total number of parameter samples
                  err,message)                 ! error code and message
    if(err/=0) call abort_mpi(instance_parallel%rank,trim(message))

  enddo

  ! ---------------------------------------------------------------------------------------
  ! Finalize MPI
  ! ---------------------------------------------------------------------------------------
  
  call MPI_Finalize(mpi_err)
  call check_mpi(world_parallel%rank,mpi_err,'MPI_Finalize failed')

  if(world_parallel%rank==0) &
    call stop_program(0,'finished parallel parameter evaluation successfully.')

contains

  ! --------------------------------------------------------------------------------------------------
  ! --------------------------------------------------------------------------------------------------
  ! --------------------------------------------------------------------------------------------------
  ! --------------------------------------------------------------------------------------------------
  ! --------------------------------------------------------------------------------------------------
  ! --------------------------------------------------------------------------------------------------
  ! --------------------------------------------------------------------------------------------------

  ! **************************************************************************************************
  ! internal subroutine run_case: initialize and execute parameter calibration for a single SUMMA case
  !
  ! This routine performs all case-specific initialization and parameter evaluation. It configures
  ! SUMMA from the specified case file, establishes rank-specific logging, performs a cold-start
  ! spinup, initializes parameter sampling, evaluates the requested parameter samples, and writes
  ! the calibration output. All MPI ranks in the instance-parallel communicator participate in the
  ! evaluation of the case.
  !
  ! The routine is intended to be called repeatedly to support sequential execution of multiple
  ! independent cases within a single MPI execution.
  !
  ! Notes:
  !   - Command-line arguments are assumed to have been processed before this routine is called.
  !   - Case-specific configuration is established by init_config using config%case_file.
  !   - Parameter-search state and calibration output are initialized independently for each case.
  !   - All ranks must call this routine collectively.
  !
  ! **************************************************************************************************
  subroutine run_case(config,domain_parallel,instance_parallel,nSamples,err,message)
  
    ! logging
    USE globalData,      only: iulog
    USE iso_fortran_env, only: error_unit
    USE iso_fortran_env, only: output_unit
  
    ! MPI
    USE mpi, only: MPI_Bcast
    USE mpi, only: MPI_CHARACTER
  
    USE error_utils, only: check_mpi,abort_mpi
  
    ! SUMMA globals
    USE globalData, only: restart_filename
  
    ! SUMMA paths/filenames
    USE summaFileManager, only: OUTPUT_PATH
    USE summaFileManager, only: MODEL_INITCOND
  
    ! SUMMA parameter information
    USE parameter_search, only: parameter_spec,parameter_search_info
  
    ! SUMMA subroutines/functions
    USE summa_init,   only: init_config
    USE summa_spinup, only: spinup_from_cold
  
    USE calibration_output_module, only: create_calibration_output
    USE calibration_output_module, only: close_calibration_output
  
    implicit none
  
    ! ---------------------------------------------------------------------------------------
    ! Dummy arguments
    ! ---------------------------------------------------------------------------------------
  
    type(config_info), intent(inout) :: config
    type(parallel_context_type), intent(in) :: domain_parallel
    type(parallel_context_type), intent(in) :: instance_parallel
  
    integer(i4b), intent(in) :: nSamples
    integer(i4b), intent(out) :: err
    character(len=*), intent(out) :: message
  
    ! ---------------------------------------------------------------------------------------
    ! Local variables
    ! ---------------------------------------------------------------------------------------
  
    character(len=4)   :: rankString
    character(len=256) :: log_file
  
    type(parameter_spec)        :: param_spec
    type(parameter_search_info) :: search
  
    character(len=64), allocatable :: param_name(:)
  
    real(rkind), allocatable :: x_best(:)
    real(rkind)              :: F_best
    integer(i4b)             :: sample_best
  
    integer(i4b)       :: ncid_calib
    character(len=256) :: calib_file
  
    integer(i4b) :: mpi_err
  
    ! ---------------------------------------------------------------------------------------
    ! Initialize error control
    ! ---------------------------------------------------------------------------------------
  
    err=0
    message='run_case/'
    mpi_err=0
  
    ! ---------------------------------------------------------------------------------------
    ! Configure SUMMA
    ! ---------------------------------------------------------------------------------------
   
    ! use stderr until the configuration and output paths are known
    iulog=error_unit
   
    ! read the configuration files to establish file paths, simulation settings, and calibration options
    call init_config(config,err,message)
    if(err/=0) call abort_mpi(instance_parallel%rank,trim(message))
    
    config%read_config = .false.
    
    ! configure rank-specific logging
    iulog=99
    config%iulog_summa=iulog
    
    write(rankString,'(I4.4)') instance_parallel%rank
    log_file=trim(OUTPUT_PATH)//'logs/'//trim(config%case_name)// &
             '_rank'//rankString//'.log'
    
    call execute_command_line('mkdir -p "'//trim(OUTPUT_PATH)//'logs"')
    open(unit=iulog,file=trim(log_file),status='replace',action='write')
  
    ! ---------------------------------------------------------------------------------------
    ! Spin up SUMMA from a cold state
    ! ---------------------------------------------------------------------------------------
  
    ! perform a one-year cold-start spinup to establish the initial model state
    call spinup_from_cold(config,             & ! SUMMA configuration structure
                          domain_parallel,    & ! MPI context for domain parallelism
                          instance_parallel,  & ! MPI context for model-instance parallelism
                          err,message)          ! error code and message
    if(err/=0) call abort_mpi(instance_parallel%rank,trim(message))
  
    ! broadcast the rank-0 restart filename to all model instances
    call MPI_Bcast(restart_filename,len(restart_filename),MPI_CHARACTER,0, &
                   instance_parallel%comm,mpi_err)
    call check_mpi(instance_parallel%rank,mpi_err, &
                   'unable to broadcast restart filename')
    
    ! use the common cold-start restart state as the initial conditions
    MODEL_INITCOND=trim(restart_filename)
  
    ! ---------------------------------------------------------------------------------------
    ! Initialize parameter sampling
    ! ---------------------------------------------------------------------------------------
  
    ! Initialize parameter-evaluation state that is invariant across samples and shared by all ranks
    call initialize_parameter_evaluation(config,            & ! SUMMA configuration
                                         instance_parallel, & ! MPI instance-parallel context
                                         param_spec,search, & ! parameter specification and search information
                                         param_name,        & ! complete SUMMA parameter-name vector
                                         x_best,F_best,     & ! current best solution and objective value
                                         sample_best,       & ! sample index associated with current best solution
                                         err,message)         ! error code and message
    if(err/=0) call abort_mpi(instance_parallel%rank,trim(message))
  
    ! ---------------------------------------------------------------------------------------
    ! Create calibration output file
    ! ---------------------------------------------------------------------------------------
    
    if(instance_parallel%rank == 0)then
      calib_file=trim(OUTPUT_PATH)//trim(config%case_name)//'_calibration.nc'
      call create_calibration_output(calib_file,param_spec,nSamples,instance_parallel%size-1,           &
                                     config%case_name,config%calib%metric,config%calib%obs_transform,   &
                                     ncid_calib,err,message)
      if(err/=0) call abort_mpi(instance_parallel%rank,trim(message))
    else
      ncid_calib=-1
    endif
  
    ! ---------------------------------------------------------------------------------------
    ! Evaluate parameter samples
    ! ---------------------------------------------------------------------------------------
  
    call dispatch_parameter_samples(config,                     & ! SUMMA configuration structure
                                    domain_parallel,            & ! MPI context for domain parallelism
                                    instance_parallel,          & ! MPI context for model-instance parallelism
                                    param_spec,search,          & ! parameter specification and search information
                                    param_name,ncid_calib,      & ! parameter names and calibration output
                                    x_best,F_best,sample_best,  & ! current best solution and objective value; sample index
                                    nSamples,err,message)         ! total number of parameter samples; error code and message
    if(err/=0) call abort_mpi(instance_parallel%rank,trim(message))
  
    ! ---------------------------------------------------------------------------------------
    ! Close calibration output
    ! ---------------------------------------------------------------------------------------
  
    if(instance_parallel%rank == 0)then
      call close_calibration_output(ncid_calib,err,message)
      if(err/=0) call abort_mpi(instance_parallel%rank,trim(message))
    endif
  
    ! close the rank-specific logging file
    close(iulog);  iulog=error_unit

  end subroutine run_case

end program summa_driver_opt

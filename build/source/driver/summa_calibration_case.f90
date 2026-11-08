module summa_calibration_case

  ! data types
  USE nr_type, only: i4b,rkind,lgt
  USE data_types, only: gru_double,gru_hru_dom_doublevec
  USE summa_type, only: config_info,parallel_context_type

  ! MPI
  USE mpi, only: MPI_COMM_SELF,MPI_SUCCESS
  USE mpi_context, only: set_mpi_context

  ! Parameter sampling and optimization
  USE summa_parameter_sampling, only: check_search_settings
  USE summa_parameter_sampling, only: initialize_parameter_evaluation
  USE summa_parameter_sampling, only: dispatch_parameter_samples
  USE summa_parameter_fields,   only: summa_parameter_map
  USE summa_parameter_fields,   only: summa_calibration_adapter
  USE summa_parameter_fields,   only: initialize_summa_calibration_adapter
  USE summa_parameter_fields,   only: apply_summa_calibration_adapter
  USE summa_parameter_export,   only: write_physical_trial_parameters
  USE summa_simulation,         only: n_calibration_targets

  implicit none
  private

  public :: initialize_summa_calibration_config
  public :: run_case

  public :: execute_summa_case

  ! State passed explicitly to the generic case dispatcher.
  type, public :: summa_case_context
    type(config_info), pointer  :: config => null()
    type(parallel_context_type) :: instance_parallel
    integer(i4b)                :: node_index = 0
    integer(i4b)                :: global_case_group = 0
  end type summa_case_context

contains

  ! **************************************************************************************************
  ! Execute a SUMMA calibration case assigned by the generic scheduler.
  !
  ! The context contains the configuration and MPI information required to
  ! execute one case. This module procedure does not capture host variables.
  ! **************************************************************************************************
  subroutine execute_summa_case(context,case_index,err,message)

    class(*),         intent(inout) :: context
    integer(i4b),     intent(in)    :: case_index
    integer(i4b),     intent(out)   :: err
    character(len=*), intent(out)   :: message

    select type (context)

    type is (summa_case_context)

      if(.not.associated(context%config))then
        message='execute_summa_case/configuration pointer is not associated'
        err=20; return
      endif

      call run_case(context%config,context%instance_parallel, &
                    context%node_index,context%global_case_group, &
                    case_index,err,message)

    class default
      err=20; message='execute_summa_case/invalid case context type'

    end select

  end subroutine execute_summa_case


  ! **************************************************************************************************
  ! initialize_summa_calibration_config: Read configuration for a SUMMA calibration execution.
  !
  ! Processes command-line arguments and, when specified, reads the multi-case manifest.
  ! Configuration is established once before MPI case grouping and individual case execution.
  ! **************************************************************************************************
  subroutine initialize_summa_calibration_config(config,err,message)

    USE summa_util,   only: getCommandArguments
    USE summa_config, only: read_manifest

    implicit none

    ! ---------------------------------------------------------------------------------------
    ! Dummy arguments
    ! ---------------------------------------------------------------------------------------
    type(config_info), intent(inout) :: config   ! SUMMA configuration structure
    integer(i4b),      intent(out)   :: err      ! Error code
    character(len=*),  intent(out)   :: message  ! Error message

    ! initialize error control
    err=0
    message='initialize_summa_calibration_config/'

    ! get the current working directory
    call get_environment_variable('PWD',config%cwd,status=err)
    if(err/=0)then
      message=trim(message)//'unable to determine current working directory'
      return
    endif

    ! process command-line arguments once before configuring individual cases
    call getCommandArguments(config,err,message)
    if(err/=0) return

    ! prevent command-line arguments from being reread during case initialization
    config%read_cli=.false.

    ! read the multi-case manifest when specified
    if(allocated(config%manifest_file))then
      call read_manifest(config,err,message)
      if(err/=0) return
    endif

  end subroutine initialize_summa_calibration_config

  ! **************************************************************************************************
  ! run_case: Initialize and execute parameter calibration for a single SUMMA case.
  !
  ! This routine manages the SUMMA-specific calibration workflow for one case. It initializes
  ! the case configuration, establishes rank-specific logging, performs a common cold-start
  ! spinup, and extracts the a-priori parameter fields. It then initializes the parameter
  ! search, constructs the mapping between calibration parameters and SUMMA data structures,
  ! dispatches parameter evaluations, and writes the calibration results.
  !
  ! Rank 0 coordinates parameter sampling and collects objective values, while the remaining
  ! ranks evaluate independent parameter trials. All ranks in the instance-parallel communicator
  ! participate in the calibration of the case.
  !
  ! The routine supports sequential calibration of multiple independent cases within a single
  ! MPI execution. Each case has its own configuration, parameter-search state, and output files.
  !
  ! Notes:
  !   - Command-line arguments must be processed before calling this routine.
  !   - Case-specific configuration is established by init_config using config%case_file.
  !   - A-priori parameter fields and their mappings are retained across parameter trials.
  !   - All ranks in the instance-parallel communicator must call this routine collectively.
  !
  ! **************************************************************************************************
  subroutine run_case(config,instance_parallel,node_index,global_case_group,case_index,err,message)
    
    ! logging
    USE globalData,      only: iulog
    USE iso_fortran_env, only: error_unit
    USE iso_fortran_env, only: output_unit
    ! MPI
    USE mpi, only: MPI_Bcast
    USE mpi, only: MPI_Barrier
    USE mpi, only: MPI_CHARACTER
    USE error_utils, only: check_mpi,abort_mpi
    ! SUMMA globals
    USE globalData, only: restart_filename
    USE globalData, only: output_fileSuffix
    ! SUMMA paths/filenames
    USE summaFileManager, only: OUTPUT_PATH
    USE summaFileManager, only: MODEL_INITCOND
    ! SUMMA subroutines/functions
    USE summa_init,   only: init_config
    USE summa_spinup, only: spinup_from_cold
    USE calibration_output_module, only: create_calibration_output
    USE calibration_output_module, only: close_calibration_output
    ! General parameter information
    USE parameter_search, only: parameter_spec,parameter_search_info
    USE parameter_search, only: search_state_type
    implicit none

    ! ---------------------------------------------------------------------------------------
    ! Dummy arguments
    ! ---------------------------------------------------------------------------------------
    type(config_info),           intent(inout) :: config             ! SUMMA configuration structure
    type(parallel_context_type), intent(in)    :: instance_parallel  ! MPI context for model-instance parallelism
    integer(i4b),                intent(in)    :: node_index         ! Index of the physical compute node
    integer(i4b),                intent(in)    :: global_case_group  ! Index of the case group across all nodes
    integer(i4b), optional,      intent(in)    :: case_index         ! Assigned case index; absent for single-case runs
    integer(i4b),                intent(out)   :: err                ! Error code
    character(len=*),            intent(out)   :: message            ! Error message

    ! ---------------------------------------------------------------------------------------
    ! Local variables
    ! ---------------------------------------------------------------------------------------

    ! Files and logging
    integer(i4b)                           :: manifest_case_index  ! Index of the basin in the manifest
    character(len=4)                       :: rankString    ! MPI rank formatted as a string
    character(len=4)                       :: runString     ! run index formatted as a string
    character(len=256)                     :: log_file      ! path to rank-specific log file
    integer(i4b)                           :: log_unit      ! unit number for rank-specific log file
    integer(i4b)                           :: istat         ! file I/O status
    integer(i4b)                           :: nSamples      ! number of parameter samples

    ! Initialized SUMMA parameters
    type(gru_hru_dom_doublevec)            :: mpar_init     ! initialized HRU parameters
    type(gru_double)                       :: bpar_init     ! initialized basin parameters
    logical(lgt)                           :: single_domain ! one GRU, one HRU and one domain

    ! Parameter specification and sampling
    type(parameter_spec)                   :: param_spec    ! complete SUMMA parameter specification
    type(summa_parameter_map), allocatable :: param_map(:)  ! Mapping from calibration parameters to SUMMA parameter locations
    type(summa_calibration_adapter)        :: adapter       ! SUMMA-specific calibration context
    type(parameter_search_info)            :: search        ! parameter sampling information
    character(len=64), allocatable         :: param_name(:) ! complete SUMMA parameter names

    ! Objective function and search state
    type(search_state_type)                :: search_state  ! state of the parameter search
    integer(i4b)                           :: sample_best   ! sample index associated with best objective

    ! Calibration output
    integer(i4b)                           :: ncid_calib    ! NetCDF identifier for calibration output
    character(len=1024)                    :: calib_file    ! calibration output file path
    character(len=1024)                    :: file_prefix   ! prefix for calib filename

    ! Reconstructed parameter fields
    type(gru_hru_dom_doublevec)            :: mpar_trial    ! Physical HRU parameters for export
    type(gru_double)                       :: bpar_trial    ! Physical basin parameters for export
    integer(i4b), allocatable              :: selected_ids(:)      ! Selected sample indices
    real(rkind), allocatable               :: selected_values(:,:) ! Native calibration values
    integer(i4b)                           :: iSelected     ! Selected sample index
    character(len=20)                      :: sampleString  ! Sample identifier for filename
    character(len=1024)                    :: trial_file    ! Trial parameter output file

    ! MPI
    type(parallel_context_type)            :: domain_parallel  ! MPI context for SUMMA domain parallelism
    integer(i4b)                           :: mpi_err       ! MPI error code
    character(len=256)                     :: mpi_message   ! MPI error message

    ! ---------------------------------------------------------------------------------------
    ! Initialize
    ! ---------------------------------------------------------------------------------------

    ! initialize error control
    err=0
    mpi_err=0
    message='run_case/'

    ! ---------------------------------------------------------------------------------------
    ! Initialize SUMMA domain parallelism
    ! ---------------------------------------------------------------------------------------

    ! each model instance independently evaluates the complete spatial domain
    domain_parallel%comm=MPI_COMM_SELF

    call set_mpi_context(domain_parallel%comm,  & ! Single-rank MPI communicator
                         domain_parallel%rank,  & ! Rank within the domain communicator
                         domain_parallel%size,  & ! Number of ranks in the domain communicator
                         mpi_err,mpi_message)     ! MPI error information

    if(mpi_err/=MPI_SUCCESS) call abort_mpi(instance_parallel%rank,trim(mpi_message))

    ! ---------------------------------------------------------------------------------------
    ! Provide the additional configuration needed for multi-case execution
    ! ---------------------------------------------------------------------------------------

    if(present(case_index))then

      ! map the assigned work item to a basin and independent optimization run
      manifest_case_index=(case_index-1)/config%n_runs+1
      config%run_index=mod(case_index-1,config%n_runs)+1

      ! configure the assigned SUMMA case
      config%manifest_casename=trim(config%case_names(manifest_case_index))
      config%config_file=trim(config%template_path)//trim(config%template_file)

      if(instance_parallel%rank==0)then
        write(output_unit,'(A,I0,A,A)') &
          'Case group ',global_case_group, &
          ' running case: ',trim(config%manifest_casename)
      endif

    endif

    ! ---------------------------------------------------------------------------------------
    ! Configure SUMMA
    ! ---------------------------------------------------------------------------------------

    ! use stderr until the configuration and output paths are known
    iulog=error_unit

    ! set case root (for printing)
    config%is_case_root = (instance_parallel%rank == 0)

    ! read the configuration files to establish file paths, simulation settings, and calibration options
    call init_config(config,err,message)
    if(err/=0) call abort_mpi(instance_parallel%rank,trim(message))
    config%read_config = .false.

    ! the sample budget is a configuration setting, so it is only known now
    nSamples = config%calib%n_samples

    ! refuse a search the configuration cannot run before the spinup is paid for
    call check_search_settings(config,err,message)
    if(err/=0) call abort_mpi(instance_parallel%rank,trim(message))
    if(instance_parallel%rank==0 .and. trim(config%calib%algorithm)=='dds' .and. n_calibration_targets(config)>1) &
      write(output_unit,'(A,I0,A)') 'WARNING: DDS searches the weighted sum of the ',n_calibration_targets(config), &
                                    ' calibration targets; nsga2 keeps them separate'

    ! get the run index as character string to avoid collision in filenames
    write(runString,'(I4.4)') config%run_index
    ! common prefix for calibration and exported parameter files.
    file_prefix=trim(OUTPUT_PATH)//trim(config%case_name)//'_run'//runString
    if(len_trim(output_fileSuffix)>0) file_prefix=trim(file_prefix)//'_'//trim(output_fileSuffix)
    ! rank-specific logging
    iulog=99
    log_unit=iulog
    config%iulog_summa=iulog
    write(rankString,'(I4.4)') instance_parallel%rank
    log_file=trim(OUTPUT_PATH)//'/logs/'//trim(config%case_name)//'_run'//runString// '_rank'//rankString//'.log'
    call execute_command_line('mkdir -p "'//trim(OUTPUT_PATH)//'/logs"')
    open(unit=iulog,file=trim(log_file),status='replace',action='write')

    ! ---------------------------------------------------------------------------------------
    ! Spin up SUMMA from a cold state
    ! ---------------------------------------------------------------------------------------
    ! perform a one-year cold-start spinup to establish the initial model state
    call spinup_from_cold(config,               & ! SUMMA configuration structure
                          domain_parallel,      & ! MPI context for domain parallelism
                          instance_parallel,    & ! MPI context for model-instance parallelism
                          mpar_init, bpar_init, & ! initialized HRU/basin parameters
                          err,message)            ! error code and message
    if(err/=0) call abort_mpi(instance_parallel%rank,trim(message))

    ! a single-domain run has one initialized value per parameter, which the search keeps for the
    ! parameters it does not sample; a larger domain falls back to the parameter-file defaults
    single_domain = (size(mpar_init%gru)==1 .and. size(bpar_init%gru)==1)
    if(single_domain) single_domain = (size(mpar_init%gru(1)%hru)==1)
    if(single_domain) single_domain = (size(mpar_init%gru(1)%hru(1)%dom)==1)

    ! broadcast the rank-0 restart filename to all model instances
    call MPI_Bcast(restart_filename,len(restart_filename),MPI_CHARACTER,0, instance_parallel%comm,mpi_err)
    call check_mpi(instance_parallel%rank,mpi_err, 'unable to broadcast restart filename')

    ! use the common cold-start restart state as the initial conditions
    MODEL_INITCOND=trim(restart_filename)

    ! ---------------------------------------------------------------------------------------
    ! Initialize parameter sampling
    ! ---------------------------------------------------------------------------------------
    ! Initialize parameter-evaluation state that is invariant across samples and shared by all ranks
    if(single_domain)then
      call initialize_parameter_evaluation(config,instance_parallel, &
                                           param_spec,param_map,search,param_name, &
                                           search_state,sample_best, &
                                           mpar_init%gru(1)%hru(1)%dom(1), &
                                           bpar_init%gru(1),err,message)
    else
      call initialize_parameter_evaluation(config,instance_parallel, &
                                           param_spec,param_map,search,param_name, &
                                           search_state,sample_best, &
                                           err=err,message=message)
    endif
    if(err/=0) call abort_mpi(instance_parallel%rank,trim(message))

    ! Rank 0 checks physical feasibility.
    ! Workers apply proposals already accepted by rank 0.
    call initialize_summa_calibration_adapter(param_spec,search,param_map, &
                                              mpar_init,bpar_init,adapter, &
                                              err,message)
    if(err/=0) call abort_mpi(instance_parallel%rank,trim(message))


    ! ---------------------------------------------------------------------------------------
    ! Export the initial physical parameter fields before calibration.
    ! ---------------------------------------------------------------------------------------
    if(instance_parallel%rank==0)then
    
      trial_file=trim(file_prefix)//'_sample0_trialParams.nc'
    
      call write_physical_trial_parameters(trial_file,config%trial_param_names, &
                                           param_spec,param_map,mpar_init,bpar_init, &
                                           0,err,message)
      if(err/=0) call abort_mpi(instance_parallel%rank,trim(message))
    
    endif

    ! ---------------------------------------------------------------------------------------
    ! Create calibration output file
    ! ---------------------------------------------------------------------------------------
    if(instance_parallel%rank == 0)then
      calib_file=trim(file_prefix)//'_calibration.nc'
      call create_calibration_output(calib_file,param_spec,nSamples,                                   &
                                     instance_parallel%size-1,node_index,global_case_group,            &
                                     config%case_name,config%calib%targets,                            &
                                     config%calib%algorithm,config%calib%nsga2%population_size,        &
                                     ncid_calib,err,message)
      if(err/=0) call abort_mpi(instance_parallel%rank,trim(message))
      ! A few concise provenance attributes; MPI layout is already written above.
      call write_minimal_calibration_provenance(ncid_calib,config,err,message)
      if(err/=0) call abort_mpi(instance_parallel%rank,trim(message))
    else
      ncid_calib=-1
    endif

    ! ---------------------------------------------------------------------------------------
    ! Evaluate parameter samples
    ! ---------------------------------------------------------------------------------------
    call dispatch_parameter_samples(config,                       & ! SUMMA configuration structure
                                    domain_parallel,              & ! MPI context for domain parallelism
                                    instance_parallel,            & ! MPI context for model-instance parallelism
                                    param_spec,search,adapter,    & ! parameter specification and search information
                                    param_name,ncid_calib,        & ! parameter names and calibration output
                                    search_state,                 & ! information from previously evaluated parameter sets
                                    sample_best,                  & ! sample index
                                    selected_ids,selected_values, & ! ids and values for selected samples
                                    nSamples,                     & ! total number of parameter samples
                                    err,message)                    ! error code and message
    if(err/=0) call abort_mpi(instance_parallel%rank,trim(message))

    ! ---------------------------------------------------------------------------------------
    ! Close calibration output
    ! ---------------------------------------------------------------------------------------
    if(instance_parallel%rank == 0)then
      call close_calibration_output(ncid_calib,err,message)
      if(err/=0) call abort_mpi(instance_parallel%rank,trim(message))
    endif

    ! ---------------------------------------------------------------------------------------
    ! Export the best DDS sample or the NSGA-II Pareto samples.
    ! ---------------------------------------------------------------------------------------

    if(instance_parallel%rank==0)then
      if(size(selected_ids)>0)then
   
        call initialize_summa_calibration_adapter(param_spec,search,param_map, &
                                                  mpar_init,bpar_init,adapter,err,message)
        if(err/=0) call abort_mpi(instance_parallel%rank,trim(message))
  
        do iSelected=1,size(selected_ids)
       
          adapter%calibration%values=selected_values(:,iSelected)
       
          mpar_trial=mpar_init
          bpar_trial=bpar_init
       
          ! Apply spatial multipliers and parameter ratios.
          call apply_summa_calibration_adapter(adapter,mpar_trial,bpar_trial,err,message)
          if(err/=0) call abort_mpi(instance_parallel%rank,trim(message))
       
          ! define the new parameter trial file
          write(sampleString,'(I0)') selected_ids(iSelected)
          trial_file=trim(file_prefix)//'_sample'//trim(sampleString)//'_trialParams.nc'
       
          ! Write the resulting physical parameter fields.
          call write_physical_trial_parameters(trial_file,config%trial_param_names,param_spec,param_map, &
                                               mpar_trial,bpar_trial,selected_ids(iSelected),err,message)
          if(err/=0) call abort_mpi(instance_parallel%rank,trim(message))
       
        enddo
   
      endif
    endif

    ! ---------------------------------------------------------------------------------------
    ! Finalize 
    ! ---------------------------------------------------------------------------------------

    ! close the rank-specific logging file
    ! NOTE: close the unit that was opened, not whatever iulog currently holds. A callee
    !       that redirects iulog would otherwise leave this file open and get stderr
    !       closed in its place.
    close(log_unit);  iulog=error_unit

    ! ensure all ranks have closed their files before case-level cleanup
    call MPI_Barrier(instance_parallel%comm,mpi_err)
    call check_mpi(instance_parallel%rank,mpi_err,'unable to synchronize case cleanup')

    ! ---------------------------------------------------------------------------------------
    ! Persist calibration results and remove temporary logs.
    ! ---------------------------------------------------------------------------------------
    if(instance_parallel%rank==0)then
    
      if(allocated(config%persistent_output))then
   
        ! Copy calibration output to persistent storage 
        call execute_command_line('mkdir -p "'//trim(config%persistent_output)//'"',exitstat=istat)
        if(istat/=0) call abort_mpi(instance_parallel%rank,'unable to create persistent output directory')
        
        call execute_command_line('cp "'//trim(calib_file)//'" "'//trim(config%persistent_output)//'/"',exitstat=istat)
        if(istat/=0) call abort_mpi(instance_parallel%rank,'unable to persist calibration output')

        ! Copy selected physical parameter files to persistent storage
        do iSelected=1,size(selected_ids)
          write(sampleString,'(I0)') selected_ids(iSelected)
          trial_file=trim(file_prefix)//'_sample'//trim(sampleString)//'_trialParams.nc'
          call execute_command_line('cp "'//trim(trial_file)//'" "'//trim(config%persistent_output)//'/"',exitstat=istat)
          if(istat/=0) call abort_mpi(instance_parallel%rank,'unable to persist trial parameters')
        enddo
    
      endif
   
      ! remove log files 
      call execute_command_line('rm -f "'//trim(OUTPUT_PATH)//'/logs/'//trim(config%case_name)// &
                                '_run'//runString//'_rank"*.log',exitstat=istat)
      if(istat/=0) call abort_mpi(instance_parallel%rank,'unable to remove calibration logs')
    
    endif

  end subroutine run_case


  ! **************************************************************************************************
  ! A minimal, human-readable provenance record for the calibration NetCDF file.
  ! Existing global attributes already include case_name, mpi_workers, node_index,
  ! global_case_group and algorithm. Do not embed entire input files or checksums.
  ! Called by the case-group leader after create_calibration_output.
  ! **************************************************************************************************
  subroutine write_minimal_calibration_provenance(ncid_calib,config,err,message)
    USE mpi, only: MPI_Get_processor_name,MPI_MAX_PROCESSOR_NAME,MPI_SUCCESS
    USE summaFileManager, only: SETTINGS_PATH,M_DECISIONS
    USE netcdf, only: NF90_GLOBAL,NF90_NOERR,NF90_REDEF,NF90_PUT_ATT,NF90_ENDDEF,NF90_STRERROR
    implicit none

    integer(i4b),     intent(in)  :: ncid_calib
    type(config_info), intent(in) :: config
    integer(i4b),     intent(out) :: err
    character(*),     intent(out) :: message

    character(len=MPI_MAX_PROCESSOR_NAME) :: hostname
    integer(i4b) :: hostname_length,mpi_err,nc_err

    ! Generated at build time; this identifies the compiled executable, not
    ! necessarily the current Git HEAD when the executable is run.
    include 'summaversion.inc'

    err=0
    message='write_minimal_calibration_provenance/'

    call MPI_Get_processor_name(hostname,hostname_length,mpi_err)
    if(mpi_err/=MPI_SUCCESS)then
      err=20; message=trim(message)//'MPI_Get_processor_name failed'; return
    endif

    nc_err=nf90_redef(ncid_calib)
    if(nc_err/=NF90_NOERR) goto 900

    ! Same resolved path components used to open the SUMMA decisions file.
    nc_err=nf90_put_att(ncid_calib,NF90_GLOBAL,'model_decisions_file', &
                       trim(SETTINGS_PATH)//trim(M_DECISIONS))
    if(nc_err/=NF90_NOERR) goto 900

    nc_err=nf90_put_att(ncid_calib,NF90_GLOBAL,'compute_hostname',hostname(:hostname_length))
    if(nc_err/=NF90_NOERR) goto 900

    nc_err=nf90_put_att(ncid_calib,NF90_GLOBAL,'git_hash',trim(gitHash))
    if(nc_err/=NF90_NOERR) goto 900

    nc_err=nf90_put_att(ncid_calib,NF90_GLOBAL,'run_index',config%run_index)
    if(nc_err/=NF90_NOERR) goto 900

    nc_err=nf90_put_att(ncid_calib,NF90_GLOBAL,'cases_per_node',config%cases_per_node)
    if(nc_err/=NF90_NOERR) goto 900

    nc_err=nf90_enddef(ncid_calib)
    if(nc_err/=NF90_NOERR) goto 900
    return

900 continue
    err=20
    message=trim(message)//trim(nf90_strerror(nc_err))
  end subroutine write_minimal_calibration_provenance

end module summa_calibration_case

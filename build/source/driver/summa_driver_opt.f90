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

  USE mpi, only: MPI_INTEGER,MPI_DOUBLE_PRECISION
  USE mpi, only: MPI_CHARACTER
  
  USE mpi, only: MPI_COMM_WORLD,MPI_COMM_SELF
  USE mpi, only: MPI_ANY_TAG,MPI_ANY_SOURCE,MPI_STATUS_SIZE
  USE mpi, only: MPI_SOURCE,MPI_TAG
  USE mpi, only: MPI_SUCCESS
  
  USE mpi_context, only: set_mpi_context
  USE error_utils, only: check_mpi,abort_mpi

  ! SUMMA subroutines/functions
  USE summa_util,   only: stop_program

  implicit none

  type(config_info) :: config                              ! SUMMA configuration information

  type(parallel_context_type) :: domain_parallel           ! MPI context for domain parallelism
  type(parallel_context_type) :: instance_parallel         ! MPI context for model-instance parallelism

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

  ! get the MPI context for the optimization ensemble
  instance_parallel%comm=MPI_COMM_WORLD

  call set_mpi_context(instance_parallel%comm,  &
                       instance_parallel%rank,  &
                       instance_parallel%size,  &
                       mpi_err,mpi_message)

  if(mpi_err/=MPI_SUCCESS) &
    call abort_mpi(instance_parallel%rank,trim(mpi_message))

  ! dynamic parameter evaluation requires one dispatcher and at least one worker
  if(instance_parallel%size < 2)then
    call abort_mpi(instance_parallel%rank, &
                   'parallel parameter evaluation requires at least two MPI ranks')
  endif

  ! each model instance runs independently without domain parallelism
  domain_parallel%comm=MPI_COMM_SELF
  domain_parallel%rank=0
  domain_parallel%size=1

  ! ---------------------------------------------------------------------------------------
  ! Run SUMMA case
  ! ---------------------------------------------------------------------------------------

  ! initialize and execute parameter calibration for the current SUMMA case
  call run_case(config,                    & ! SUMMA configuration structure
                domain_parallel,           & ! MPI context for domain parallelism
                instance_parallel,         & ! MPI context for model-instance parallelism
                nSamples,                  & ! total number of parameter samples
                err,message)                 ! error code and message
  if(err/=0) call abort_mpi( instance_parallel%rank, trim(message) )

  ! ---------------------------------------------------------------------------------------
  ! Finalize MPI
  ! ---------------------------------------------------------------------------------------
  
  call MPI_Finalize(mpi_err)
  call check_mpi(instance_parallel%rank,mpi_err,'MPI_Finalize failed')

  if(instance_parallel%rank==0) &
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
    USE globalData, only: initConfig
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
    message=''
    mpi_err=0
  
    ! ---------------------------------------------------------------------------------------
    ! Configure SUMMA
    ! ---------------------------------------------------------------------------------------
    
    ! use stderr until the configuration and output paths are known
    iulog=error_unit
   
    ! read the configuration files to establish file paths, simulation settings, and calibration options
    call init_config(config,err,message)
    if(err/=0) call abort_mpi(instance_parallel%rank,trim(message))
    
    initConfig=.false.
    
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

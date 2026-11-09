
module calibration_parallel

  USE nr_type, only: i4b
  USE parallel_context, only: parallel_context_type
  USE iso_fortran_env, only: output_unit

  ! MPI
  USE mpi, only: MPI_Comm_split_type,MPI_Comm_split,MPI_Bcast
  USE mpi, only: MPI_COMM_TYPE_SHARED,MPI_INFO_NULL,MPI_UNDEFINED
  USE mpi, only: MPI_INTEGER,MPI_SUCCESS
  USE mpi_context, only: set_mpi_context
  USE error_utils, only: check_mpi,abort_mpi

  implicit none
  private

  public :: init_calibration_parallel

contains

  ! **************************************************************************************************
  ! Initialize MPI communicators for independent calibration cases.
  !
  ! For a single case, all ranks participate in one calibration. For multiple cases, ranks are
  ! divided into node-local groups that independently execute calibration cases. Each group
  ! contains one dispatcher and one or more workers.
  !
  ! This routine must be called collectively by all ranks in the global communicator.
  ! **************************************************************************************************
  subroutine init_calibration_parallel(multi_case,cases_per_node,world_parallel, &
                                       node_parallel,instance_parallel, &
                                       node_index,case_group,global_case_group,nCaseGroups,err,message)

    ! ---------------------------------------------------------------------------------------
    ! Dummy arguments
    ! ---------------------------------------------------------------------------------------
    logical,                     intent(in)    :: multi_case         ! Enable multiple independent calibration cases
    integer(i4b),                intent(inout) :: cases_per_node     ! Number of concurrent cases per node
    type(parallel_context_type), intent(in)    :: world_parallel    ! Global MPI context
    type(parallel_context_type), intent(out)   :: node_parallel     ! Node-local MPI context
    type(parallel_context_type), intent(out)   :: instance_parallel ! Calibration instance MPI context
    integer(i4b),                intent(out)   :: node_index        ! Index of the physical node
    integer(i4b),                intent(out)   :: case_group        ! Case group on the local node
    integer(i4b),                intent(out)   :: global_case_group ! Case group across all nodes
    integer(i4b),                intent(out)   :: nCaseGroups       ! Total number of case groups
    integer(i4b),                intent(out)   :: err               ! Error code
    character(len=*),            intent(out)   :: message           ! Error message

    ! ---------------------------------------------------------------------------------------
    ! Local variables
    ! ---------------------------------------------------------------------------------------
    type(parallel_context_type) :: leader_parallel  ! MPI context for physical node leaders
    integer(i4b)                :: nNodes           ! Number of physical nodes
    integer(i4b)                :: ranks_per_case   ! MPI ranks assigned to each case
    integer(i4b)                :: leader_color     ! Color for constructing node-leader communicator
    integer(i4b)                :: requested_cases  ! Requested number of concurrent cases per node
    integer                     :: mpi_err          ! MPI error code
    character(len=256)          :: mpi_message      ! MPI error message

    err=0
    message=''

    node_index=0
    case_group=0
    global_case_group=0
    nCaseGroups=1

    ! ---------------------------------------------------------------------------------------
    ! Single-case run: use all available MPI ranks for one calibration
    ! ---------------------------------------------------------------------------------------
    if(.not.multi_case)then

      instance_parallel%comm=world_parallel%comm

      call set_mpi_context(instance_parallel%comm,  &
                           instance_parallel%rank,  &
                           instance_parallel%size,  &
                           mpi_err,mpi_message)

      if(mpi_err/=MPI_SUCCESS) call abort_mpi(world_parallel%rank,trim(mpi_message))

    ! ---------------------------------------------------------------------------------------
    ! Multi-case run: partition MPI ranks among independent calibrations
    ! ---------------------------------------------------------------------------------------
    else

      ! construct a node-local communicator
      call MPI_Comm_split_type(world_parallel%comm, MPI_COMM_TYPE_SHARED, &
                               world_parallel%rank,MPI_INFO_NULL,node_parallel%comm,mpi_err)

      call check_mpi(world_parallel%rank,mpi_err,'unable to create node-local communicator')

      call set_mpi_context(node_parallel%comm,  &
                           node_parallel%rank,  &
                           node_parallel%size,  &
                           mpi_err,mpi_message)

      if(mpi_err/=MPI_SUCCESS) call abort_mpi(world_parallel%rank,trim(mpi_message))

      ! construct the node-leader communicator
      if(node_parallel%rank==0)then
        leader_color=0
      else
        leader_color=MPI_UNDEFINED
      endif

      call MPI_Comm_split(world_parallel%comm,leader_color,world_parallel%rank, &
                          leader_parallel%comm,mpi_err)

      call check_mpi(world_parallel%rank,mpi_err,'unable to create node-leader communicator')

      ! determine physical node indices
      if(node_parallel%rank==0)then

        call set_mpi_context(leader_parallel%comm,  &
                             leader_parallel%rank,  &
                             leader_parallel%size,  &
                             mpi_err,mpi_message)

        if(mpi_err/=MPI_SUCCESS) call abort_mpi(world_parallel%rank,trim(mpi_message))

        node_index=leader_parallel%rank
        nNodes=leader_parallel%size

      endif

      ! distribute node information
      call MPI_Bcast(node_index,1,MPI_INTEGER,0,node_parallel%comm,mpi_err)
      call check_mpi(world_parallel%rank,mpi_err,'unable to broadcast node index')

      call MPI_Bcast(nNodes,1,MPI_INTEGER,0,node_parallel%comm,mpi_err)
      call check_mpi(world_parallel%rank,mpi_err,'unable to broadcast number of nodes')

      ! -------------------------------------------------------------------------------------
      ! Partition ranks on each node among independent calibration cases
      ! -------------------------------------------------------------------------------------
      requested_cases=cases_per_node

      if(node_parallel%size < 2)then
        write(message,'(A,I0,A)')                                                 &
          'parameter calibration needs at least 2 MPI ranks per node, but got ',  &
          node_parallel%size,'; one rank coordinates the search and the rest evaluate samples'
        call abort_mpi(world_parallel%rank,trim(message))
      endif

      ! limit concurrent cases to those supported by the available ranks
      if(cases_per_node > node_parallel%size/2) cases_per_node=node_parallel%size/2

      ! require equal-sized case groups
      do while(mod(node_parallel%size,cases_per_node)/=0)
        cases_per_node=cases_per_node-1
      enddo

      if(cases_per_node /= requested_cases .and. world_parallel%rank == 0)then
        write(output_unit,'(A,I0,A,I0,A)')                              &
          'WARNING: cases_per_node reduced from ',requested_cases,' to ', &
          cases_per_node,' to fit the ranks available'
        write(output_unit,'(A,I0,A)')                                   &
          '         with ',node_parallel%size,' ranks per node, cases_per_node must divide'
        write(output_unit,'(A)')                                        &
          '         that evenly and leave at least 2 ranks per case'
        write(output_unit,'(A,I0,A)')                                   &
          '         run with ',requested_cases*2,' ranks per node to get the cases you asked for'
      endif

      ! assign contiguous ranks to independent calibration cases
      ranks_per_case=node_parallel%size/cases_per_node
      case_group=node_parallel%rank/ranks_per_case

      ! identify case groups across all physical nodes
      global_case_group=node_index*cases_per_node + case_group
      nCaseGroups=nNodes*cases_per_node

      ! construct the instance-parallel communicator
      call MPI_Comm_split(node_parallel%comm,case_group,node_parallel%rank, &
                          instance_parallel%comm,mpi_err)

      call check_mpi(world_parallel%rank,mpi_err,'unable to create case communicator')

      call set_mpi_context(instance_parallel%comm,  &
                           instance_parallel%rank,  &
                           instance_parallel%size,  &
                           mpi_err,mpi_message)

      if(mpi_err/=MPI_SUCCESS) call abort_mpi(world_parallel%rank,trim(mpi_message))

    endif

    ! require one dispatcher and at least one worker per calibration
    if(instance_parallel%size < 2)then
      call abort_mpi(world_parallel%rank, &
                     'parallel parameter evaluation requires at least two MPI ranks per case')
    endif

  end subroutine init_calibration_parallel

end module calibration_parallel


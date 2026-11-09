program calibration_driver

  USE nr_type, only: rkind
  USE mpi
  USE test_model, only: initialize_model,evaluate_model

  implicit none

  integer :: ierr,rank
  real(rkind) :: theta(1)
  real(rkind), allocatable :: simulated(:)

  call MPI_Init(ierr)
  call MPI_Comm_rank(MPI_COMM_WORLD,rank,ierr)

  call initialize_model()

  theta(1)=2._rkind
  call evaluate_model(theta,simulated)

  if(rank==0) print *, 'Simulated:',simulated

  call MPI_Finalize(ierr)

end program calibration_driver

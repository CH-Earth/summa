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

program summa_driver_mpi
  
  ! **** Driver program for SUMMA simulations ****

  ! * module access *
  ! data types
  USE nr_type, only: i4b                                      ! variable types, etc.
  USE summa_type, only: summa1_type_dec                       ! master summa data type
  USE mpi, only : MPI_COMM_WORLD, MPI_SUCCESS                 ! MPI constants
  ! subroutines and functions: MPI 
  USE mpi, only : MPI_Init, MPI_Finalize                      ! MPI subroutine interfaces
  USE mpi_context, only : set_mpi_context                     ! define MPI communication context
  USE error_utils, only : check_mpi, abort_mpi                ! check MPI errors
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
  USE globalData, only: numtim                                ! number of model time steps
  USE globalData, only: isPrint                               ! flag to enable informational screen/log output
  USE globalData, only: print_step_freq

  implicit none

  ! * driver variables *

  ! master summa data structure
  type(summa1_type_dec), allocatable :: summa1_struc(:)

  ! timing information
  integer(i4b), parameter :: n = 1
  integer(i4b) :: modelTimeStep

  ! MPI
  integer(i4b) :: rank = 0
  integer(i4b) :: size = 1

  ! error codes
  integer(i4b) :: err = 0
  integer(i4b) :: mpi_err = 0

  ! error messages
  character(len=1024) :: message = ''
  character(len=256)  :: mpi_message = ''

  ! ---- Initialize MPI ----------------------------------------

  call MPI_Init(mpi_err)
  call check_mpi(-1, mpi_err, 'MPI_Init failed')

  call set_mpi_context(MPI_COMM_WORLD, rank, size, mpi_err, mpi_message)
  if (mpi_err /= MPI_SUCCESS) call abort_mpi(rank, trim(mpi_message)) 

  isPrint = (rank == 0)

  ! ---- Initialize SUMMA --------------------------------------
 
  call initialize_summa_driver

  ! ---- Update SUMMA ------------------------------------------
  
  call update_summa_driver

  ! ---- Finalize SUMMA ----------------------------------------
  
  call finalize_summa_driver

  ! ---- Finalize MPI ----------------------------------------
  
  call MPI_Finalize(mpi_err)
  if (mpi_err /= MPI_SUCCESS)then
    write(message,'(A,I0,A)') 'ERROR [rank ', rank, ']: MPI_Finalize failed'
    call handle_err(mpi_err, message)
  endif

  if (rank == 0) then
    call stop_program(0, 'finished simulation successfully.')
  end if

contains

  ! --------------------------------------------------------------------------------------------------
  ! --------------------------------------------------------------------------------------------------
  ! --------------------------------------------------------------------------------------------------
  ! --------------------------------------------------------------------------------------------------

  subroutine initialize_summa_driver
   ! *** Initial operations for SUMMA driver program ***

   allocate(summa1_struc(n), stat=err)
   if (err/=0) call handle_err(err, 'problem allocating master summa structure')

   summa1_struc(n)%parallel%comm = MPI_COMM_WORLD
   summa1_struc(n)%parallel%rank = rank
   summa1_struc(n)%parallel%size = size

   ! declare and allocate summa data structures and initialize model state to known values
   call summa_initialize(summa1_struc(n), err, message)
   call handle_err(err, message)

   ! initialize parameter data structures (e.g. vegetation and soil parameters)
   call summa_paramSetup(summa1_struc(n), err, message)
   call handle_err(err, message)

   ! read restart data and reset the model state
   call summa_readRestart(summa1_struc(n), err, message)
   call handle_err(err, message)

  end subroutine initialize_summa_driver

  ! --------------------------------------------------------------------------------------------------
  ! --------------------------------------------------------------------------------------------------

  subroutine update_summa_driver

    ! *** Update operations for SUMMA driver program ***

   ! loop through time
   do modelTimeStep=1,numtim
 
     ! read model forcing data
     call summa_readForcing(modelTimeStep, summa1_struc(n), err, message)
     call handle_err(err, message)
 
     if (mod(modelTimeStep, print_step_freq) == 0) then
       print *, 'step ---> ', modelTimeStep
     end if
 
     ! run the summa physics for one time step
     call summa_runPhysics(modelTimeStep, summa1_struc(n), err, message)
     call handle_err(err, message)
 
     ! write the model output
     call summa_writeOutputFiles(modelTimeStep, summa1_struc(n), err, message)
     call handle_err(err, message)
 
   end do  ! end looping through time
  end subroutine update_summa_driver

  ! --------------------------------------------------------------------------------------------------
  ! --------------------------------------------------------------------------------------------------

  subroutine finalize_summa_driver

  ! *** Final operations for SUMMA driver program ***

  ! Any SUMMA/BMI/interface cleanup goes here

  ! Allow output libraries to complete file closure
  call sleep(2)

  end subroutine finalize_summa_driver

end program summa_driver_mpi

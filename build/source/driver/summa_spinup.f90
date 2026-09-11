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

module summa_spinup

  ! data types
  USE nr_type,    only: i4b,rkind
  USE summa_type, only: config_info,parallel_context_type

  implicit none
  private

  public :: spinup_from_cold

contains


  ! **************************************************************************************************
  ! Spin up SUMMA from a cold state.
  !
  ! Performs a one-year cold-start spinup.
  ! Additional parameter-specific warmup will be required after this point.
  !
  ! This initial warmup initializes parameter data structures that can be accessed by other routines.
  ! **************************************************************************************************

  subroutine spinup_from_cold(config,             & ! SUMMA configuration structure
                              domain_parallel,    & ! MPI context for domain parallelism
                              instance_parallel,  & ! MPI context for model-instance parallelism
                              err,message)          ! error code and error message

    ! SUMMA global configuration
    USE globalData, only: ixRestart
    USE globalData, only: ixRestart_end,ixRestart_never
    USE globalData, only: output_fileSuffix
    
    ! filenames
    USE summaFileManager, only: SIM_START_TM,SIM_END_TM
    
    ! SUMMA simulation
    USE summa_simulation, only: run_simulation

    implicit none

    ! dummy variables
    type(config_info),           intent(inout) :: config             ! SUMMA configuration structure
    type(parallel_context_type), intent(in)    :: domain_parallel    ! MPI context for domain parallelism
    type(parallel_context_type), intent(in)    :: instance_parallel  ! MPI context for model-instance parallelism
    integer(i4b),                intent(out)   :: err                ! error code
    character(*),                intent(out)   :: message            ! error message

    ! strings
    character(len=4)   :: rankString
    character(len=256) :: cmessage

    character(len=:), allocatable :: outputFileSuffix_orig
    character(len=:), allocatable :: simStartOriginal,simEndOriginal,spinStart
    character(len=:), allocatable :: timeUnits,flowUnits

    ! baseline parameter vectors
    character(len=64), allocatable :: param_name(:)
    real(rkind),       allocatable :: param_value(:)

    ! simulated time and flow
    real(rkind), allocatable :: timeSim(:),flowSim(:)

    ! local variables
    integer(i4b) :: iyear

    err=0
    message='spinup_from_cold/'

    ! zero-length parameter vectors for baseline initialization
    allocate(param_name(0),param_value(0),stat=err)
    if(err/=0)then
      message=trim(message)//'unable to allocate baseline parameter vectors'
      return
    endif

    ! Modify dates to define a one-year cold-start spinup
    ! prior to the start of the simulation period

    simStartOriginal=trim(SIM_START_TM)
    simEndOriginal=trim(SIM_END_TM)

    spinStart=simStartOriginal

    read(spinStart(1:4),*) iyear
    write(spinStart(1:4),'(I4.4)') iyear-1

    SIM_START_TM=spinStart
    SIM_END_TM=simStartOriginal

    ! use rank-specific output filenames during the common spinup

    write(rankString,'(I4.4)') instance_parallel%rank

    outputFileSuffix_orig=trim(output_fileSuffix)
    output_fileSuffix=trim(outputFileSuffix_orig)//'_spinup_rank'//rankString

    ! write the common restart state on rank 0 only
    ixRestart=merge(ixRestart_end, ixRestart_never, instance_parallel%rank == 0)

    ! run SUMMA for one year following the cold start
    call run_simulation(config,                 & ! SUMMA configuration structure
                        domain_parallel,        & ! MPI context for domain parallelism
                        instance_parallel,      & ! MPI context for model-instance parallelism
                        timeSim,flowSim,        & ! simulated time and streamflow
                        timeUnits,flowUnits,    & ! time and streamflow units
                        param_name,param_value, & ! parameter names and values
                        err,cmessage)             ! error code and message
    
    ! restore original global simulation settings
    ! NOTE: Do this before processing the error code
    SIM_START_TM=simStartOriginal
    SIM_END_TM=simEndOriginal
    output_fileSuffix=outputFileSuffix_orig
    ixRestart=ixRestart_never

    if(err/=0)then; message=trim(message)//trim(cmessage); return; endif

  end subroutine spinup_from_cold

end module summa_spinup

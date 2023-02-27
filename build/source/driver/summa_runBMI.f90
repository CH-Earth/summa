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

program summa_runBMI
  ! driver program for summa simulations
  ! *****************************************************************************
  ! * use desired modules
  ! *****************************************************************************
  ! subroutines and functions: model simulation
  use summa_driver
  use, intrinsic :: iso_fortran_env, only : file_unit=>input_unit

  implicit none

  ! *****************************************************************************
  ! * variable definitions
  ! *****************************************************************************
  character (len=*), parameter :: output_file = "summa_bmi.out"
  character (len=*), parameter :: var_name = "land_surface_water__runoff_volume_flux"
  type (summa_bmi) :: model
  integer                      :: arg_count = 0
  character (len=80)           :: arg
  integer                      :: i, j, istat, grid_id, grid_size
  double precision             :: current_time, end_time
  real, allocatable            :: runoff(:)

  ! *****************************************************************************
  ! * model simulation
  ! *****************************************************************************
  do while (arg_count <= 1)
    call get_command_argument(arg_count, arg)
    arg_count = arg_count + 1
  end do

  if (len_trim(arg) == 0) then
     write(*,"(a)") "Usage: summa_bmi.exe CONFIGURATION_FILE<80char"
     write(*,"(a)")
     write(*,"(a)") "Run the summa model through its BMI with a configuration file."
     write(*,"(a)") "Output is written to the file `summa_bmi.out`."
     stop
  end if

  open(file_unit,file=output_file)

  write(file_unit,"(a)") "Initialize model."
  istat = model%initialize(arg)

  istat = model%get_current_time(current_time)
  istat = model%get_end_time(end_time)
! variable will need to be on a grid that is constant through simulation and set in initialization
  istat = model%get_var_grid(var_name, grid_id)
  istat = model%get_grid_size(grid_id, grid_size)

  allocate(runoff(grid_size))

  do while (current_time <= end_time)
    write(file_unit,"(a, f10.2)") "Model values at time = ", current_time
    istat = model%get_value(var_name, runoff)
    do j = 1, grid_size
      write (file_unit,"(e13.5)", advance="no") runoff(j)
    end do
    write (file_unit,*)
    istat = model%update()
    istat = model%get_current_time(current_time)
  end do

  deallocate(runoff)
  istat = model%finalize()
  write(file_unit,"(a)") "Finalize model."

  close(file_unit)

end program summa_runBMI

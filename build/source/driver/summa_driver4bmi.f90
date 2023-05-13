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

program summa_driver4bmi
  ! driver program for summa simulations
  ! *****************************************************************************
  ! * use desired modules
  ! *****************************************************************************
  USE nrtype                                                  ! variable types, etc.
  ! subroutines and functions: model simulation
  USE summa_bmi
  ! global data
  USE globalData,only:numtim                                  ! number of time steps

  implicit none

  ! *****************************************************************************
  ! * variable definitions
  ! *****************************************************************************
  type (summa_bmi) :: model
  integer(i4b)                       :: istat
  ! define timing information
  integer(i4b)                       :: modelTimeStep         ! index of model time step

  ! *****************************************************************************
  ! * model simulation
  ! *****************************************************************************
  ! give this a 0 length argument to use fileManager from summa standard command arguments
  istat = model%initialize('')

  ! loop through time where numtim has been already computed as
  ! numtim = nint( (dJulianFinsh - dJulianStart)*secprday/data_step ) + 1
  ! SUMMA runs the ending step (so start=end would still run a step)
  do modelTimeStep=1,numtim
    istat = model%update()
  end do  ! (looping through time)
  istat = model%finalize()

end program summa_driver4bmi

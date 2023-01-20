! SUMMA - Structure for Unifying Multiple Modeling Alternatives
! Copyright (C) 2014-2015 NCAR/RAL
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

program multi_driver
! used to evaluate different methods for simulating snow processes
! *****************************************************************************
! use desired modules
! *****************************************************************************
use,intrinsic :: iso_c_binding
USE,intrinsic :: ieee_arithmetic                            ! IEEE arithmetic (obviously)
use summa_bmi, only: initialize, update, finalize, modelTimeStep
USE globalData,only:numtim                                  ! number of time steps
implicit none

integer(kind=c_int) :: istat

!'main routine'
istat = initialize()

! loop through time
do modelTimeStep=1,numtim
    istat = update()
end do  ! (looping through time)
istat = finalize()

end program multi_driver

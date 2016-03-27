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

module paramCheck_module
! define numerical recipes data type
USE nrtype
! define look-up values for the choice of method to combine and sub-divide snow layers
USE mDecisions_module,only:&
 sameRulesAllLayers, & ! SNTHERM option: same combination/sub-dividion rules applied to all layers
 rulesDependLayerIndex ! CLM option: combination/sub-dividion rules depend on layer index
implicit none
private
public::paramCheck
contains


 ! ************************************************************************************************
 ! public subroutine paramCheck: check consistency of model parameters
 ! ************************************************************************************************
 subroutine paramCheck(mpar_data,err,message)
 ! model decisions
 USE globalData,only:model_decisions  ! model decision structure
 USE var_lookup,only:iLookDECISIONS   ! named variables for elements of the decision structure
 ! SUMMA look-up variables
 USE var_lookup,only:iLookPARAM       ! named variables for elements of the data structures
 implicit none
 ! define input
 real(dp),intent(in)            :: mpar_data(:)         ! parameter vector
 ! define output
 integer(i4b),intent(out)       :: err                  ! error code
 character(*),intent(out)       :: message              ! error message
 ! local variables
 integer(i4b)                   :: iLayer               ! index of model layers
 real(dp),dimension(5)          :: zminLayer            ! minimum layer depth in each layer (m)
 real(dp),dimension(4)          :: zmaxLayer_lower      ! lower value of maximum layer depth
 real(dp),dimension(4)          :: zmaxLayer_upper      ! upper value of maximum layer depth
 ! Start procedure here
 err=0; message="paramCheck/"

 ! *****
 ! * check that the snow layer bounds are OK...
 ! ********************************************

 ! select option for combination/sub-division of snow layers
 select case(model_decisions(iLookDECISIONS%snowLayers)%iDecision)
  ! SNTHERM option
  case(sameRulesAllLayers)
   if(mpar_data(iLookPARAM%zmax)/mpar_data(iLookPARAM%zmin) < 2.5_dp)then
    message=trim(message)//'zmax must be at least 2.5 times larger than zmin: this avoids merging layers that have just been divided'
    err=20; return
   endif
  ! CLM option
  case(rulesDependLayerIndex)
   ! (build vectors of min/max)
   zminLayer       = (/mpar_data(iLookPARAM%zminLayer1),&
                       mpar_data(iLookPARAM%zminLayer2),&
                       mpar_data(iLookPARAM%zminLayer3),&
                       mpar_data(iLookPARAM%zminLayer4),&
                       mpar_data(iLookPARAM%zminLayer5)/)
   zmaxLayer_lower = (/mpar_data(iLookPARAM%zmaxLayer1_lower),&
                       mpar_data(iLookPARAM%zmaxLayer2_lower),&
                       mpar_data(iLookPARAM%zmaxLayer3_lower),&
                       mpar_data(iLookPARAM%zmaxLayer4_lower)/)
   zmaxLayer_upper = (/mpar_data(iLookPARAM%zmaxLayer1_upper),&
                       mpar_data(iLookPARAM%zmaxLayer2_upper),&
                       mpar_data(iLookPARAM%zmaxLayer3_upper),&
                       mpar_data(iLookPARAM%zmaxLayer4_upper)/)
   ! (check consistency)
   do iLayer=1,4  ! NOTE: the lower layer does not have a maximum value
    ! ensure that we have higher maximum thresholds for sub-division when fewer number of layers
    if(zmaxLayer_lower(iLayer) < zmaxLayer_upper(iLayer))then
     write(message,'(a,2(i0,a))') trim(message)//'expect the maximum threshold for sub-division in the case where there is only ', &
                                  iLayer,' layer(s) is greater than the maximum threshold for sub-division in the case where there are > ',&
                                  iLayer,' layer(s)'
     err=20; return
    endif
    ! ensure that the maximum thickness is 3 times greater than the minimum thickness
    if(zmaxLayer_upper(iLayer)/zminLayer(iLayer) < 2.5_dp .or. zmaxLayer_upper(iLayer)/zminLayer(iLayer+1) < 2.5_dp)then
     write(*,'(a,1x,3(f20.10,1x))') 'zmaxLayer_upper(iLayer), zminLayer(iLayer), zminLayer(iLayer+1) = ', &
                                     zmaxLayer_upper(iLayer), zminLayer(iLayer), zminLayer(iLayer+1)
     write(message,'(a,3(i0,a))') trim(message)//'zmaxLayer_upper for layer ',iLayer,' must be 2.5 times larger than zminLayer for layers ',&
                                  iLayer,' and ',iLayer+1,': this avoids merging layers that have just been divided'
     err=20; return
    endif
   end do  ! loop through layers
  case default; err=20; message=trim(message)//'unable to identify option to combine/sub-divide snow layers'; return
 end select ! (option to combine/sub-divide snow layers)

 ! -------------------------------------------------------------------------------------------------------------------------------------------

 ! *****
 ! * check soil stress functionality...
 ! ************************************

 ! check that the maximum transpiration limit is within bounds
 if(mpar_data(iLookPARAM%critSoilTranspire)>mpar_data(iLookPARAM%theta_sat) .or. &
    mpar_data(iLookPARAM%critSoilTranspire)<mpar_data(iLookPARAM%theta_res))then
  message=trim(message)//'critSoilTranspire parameter is out of range '// &
                         '[NOTE: if overwriting Noah-MP soil table values in paramTrial, must overwrite all soil parameters]'
  err=20; return
 endif

 ! check that the soil wilting point is within bounds
 if(mpar_data(iLookPARAM%critSoilWilting)>mpar_data(iLookPARAM%theta_sat) .or. &
    mpar_data(iLookPARAM%critSoilWilting)<mpar_data(iLookPARAM%theta_res))then
  print*, 'mpar_data(iLookPARAM%theta_res) = ', mpar_data(iLookPARAM%theta_res)
  print*, 'mpar_data(iLookPARAM%theta_sat) = ', mpar_data(iLookPARAM%theta_sat)
  message=trim(message)//'critSoilWilting parameter is out of range '// &
                         '[NOTE: if overwriting Noah-MP soil table values in paramTrial, must overwrite all soil parameters]'
  err=20; return
 endif

 ! check that the field capacity is within bounds
 if(mpar_data(iLookPARAM%fieldCapacity)>mpar_data(iLookPARAM%theta_sat) .or. &
    mpar_data(iLookPARAM%fieldCapacity)<mpar_data(iLookPARAM%theta_res))then
  print*, 'mpar_data(iLookPARAM%fieldCapacity) = ', mpar_data(iLookPARAM%fieldCapacity)
  message=trim(message)//'fieldCapacity parameter is out of range '// &
                         '[NOTE: if overwriting Noah-MP soil table values in paramTrial, must overwrite all soil parameters]'
  err=20; return
 endif

 end subroutine paramCheck


end module paramCheck_module

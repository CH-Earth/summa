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
module summa_lakeCoupling
! lake forcing passed from SUMMA to mizuRoute
implicit none
private
public::set_lake_coupling
contains

 ! lake forcing of each water-class HRU for mizuRoute: the HRU's precipitation and Oudin potential evaporation (m s-1)
 subroutine set_lake_coupling(summa1_struc)
 USE nr_type
 USE summa_type, only:summa1_type_dec
 USE globalData,only:gru_struc                       ! gru-hru mapping structures
 USE var_lookup,only:iLookTYPE,iLookATTR,iLookFORCE,iLookPARAM ! named variables for HRU type, attributes, forcing and parameters
 USE module_sf_noahmplsm,only:isWater                ! water land cover type
 USE multiconst,only:iden_water                      ! intrinsic density of liquid water (kg m-3)
 USE derivforce_module,only:calcPotentialEvap_Oudin2005 ! potential evaporation (mm s-1)
 implicit none
 type(summa1_type_dec),intent(inout)   :: summa1_struc ! master summa data structure
 integer(i4b)                          :: iGRU,iHRU    ! GRU and HRU indices
 integer(i4b)                          :: iWater       ! index of the water-class HRU

 ! find the water-class HRUs once: attributes are not read yet when mizuRoute is initialized
 if(.not.allocated(summa1_struc%coupling_lake))then
  iWater = 0
  do iGRU=1,summa1_struc%nGRU_local
   iWater = iWater + count([(summa1_struc%typeStruct%gru(iGRU)%hru(iHRU)%var(iLookTYPE%vegTypeIndex)==isWater, iHRU=1,gru_struc(iGRU)%hruCount)])
  end do
  allocate(summa1_struc%coupling_lake(iWater))
  iWater = 0
  do iGRU=1,summa1_struc%nGRU_local
   do iHRU=1,gru_struc(iGRU)%hruCount
    if(summa1_struc%typeStruct%gru(iGRU)%hru(iHRU)%var(iLookTYPE%vegTypeIndex)/=isWater) cycle
    iWater = iWater + 1
    summa1_struc%coupling_lake(iWater)%id    = gru_struc(iGRU)%hruInfo(iHRU)%hru_id
    summa1_struc%coupling_lake(iWater)%ixGRU = iGRU
    summa1_struc%coupling_lake(iWater)%ixHRU = iHRU
    summa1_struc%coupling_lake(iWater)%area  = summa1_struc%attrStruct%gru(iGRU)%hru(iHRU)%var(iLookATTR%HRUarea)
   end do
  end do
 endif

 do iWater=1,size(summa1_struc%coupling_lake)
  associate(water => summa1_struc%coupling_lake(iWater))
  associate(forc => summa1_struc%forcStruct%gru(water%ixGRU)%hru(water%ixHRU)%var, &
            mpar => summa1_struc%mparStruct%gru(water%ixGRU)%hru(water%ixHRU)%dom(1)%var)
  water%precip = forc(iLookFORCE%pptrate)/iden_water
  water%evap   = calcPotentialEvap_Oudin2005(forc(iLookFORCE%SWRadAtm), forc(iLookFORCE%airtemp), &
                 mpar(iLookPARAM%oudinPETScaleK1)%dat(1), mpar(iLookPARAM%oudinPETTempThrK2)%dat(1))/iden_water
  end associate
  end associate
 end do

 end subroutine set_lake_coupling

end module summa_lakeCoupling

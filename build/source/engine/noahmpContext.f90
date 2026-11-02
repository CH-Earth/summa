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

module noahmpContext_module

! data types
USE nr_type

! global variables
USE globalData,only:&
                    model_decisions,      & ! model decision structure
                    greenVegFrac_monthly, & ! fraction of green vegetation in each month (0-1)
                    urbanVegCategory,     & ! vegetation category for urban areas
                    overwriteRSMIN          ! flag to overwrite RSMIN

! provide access to the derived types to define the data structures
USE data_types,only:&
                    var_i,            & ! data vector (i4b)
                    var_dlength         ! data vector with variable length dimension (rkind)

! named variables defining elements in the data structures
USE var_lookup,only:iLookTYPE,iLookPARAM  ! named variables for structure elements
USE var_lookup,only:iLookDECISIONS        ! named variables for elements of the decision structure

! look-up values for the choice of LAI method
USE mDecisions_module,only:specified      ! LAI/SAI computed from green vegetation fraction and winterSAI and summerLAI parameters

! Noah-MP parameters for a single HRU
USE noahmp_globals,only:noahmp_context

! privacy
implicit none
private
public::noahmpContext
contains


 ! ************************************************************************************************
 ! public subroutine noahmpContext: get the Noah-MP parameters for a single HRU
 ! ************************************************************************************************
 ! NOTE: copies the HRU's values out of the shared read-only tables, so concurrent HRUs do not overwrite each other
 subroutine noahmpContext(&
                          type_data,         & ! intent(in):  type of vegetation and soil
                          mpar_data,         & ! intent(in):  model parameters
                          context,           & ! intent(out): Noah-MP parameters for this HRU
                          err,message)         ! intent(out): error control
 ! -------------------------------------------------------------------------------------------------
 ! modules
 USE NOAHMP_ROUTINES,only:redprm             ! assign Noah-MP parameters from the tables
 implicit none
 ! -------------------------------------------------------------------------------------------------
 ! input/output
 type(var_i),intent(in)            :: type_data      ! type of vegetation and soil
 type(var_dlength),intent(in)      :: mpar_data      ! model parameters
 type(noahmp_context),intent(out)  :: context        ! Noah-MP parameters for this HRU
 integer(i4b),intent(out)          :: err            ! error code
 character(*),intent(out)          :: message        ! error message
 ! -------------------------------------------------------------------------------------------------
 ! initialize error control
 err=0; message='noahmpContext/'

 ! populate parameters from the Noah-MP tables
 ! Passing an unrealistically large number of soil layers in order to pass the check for NROOT, that is done to avoid making any changes to Noah-MP code.
 !  --> NROOT from Noah-MP veg tables (as read here) is not used in SUMMA
 call REDPRM(type_data%var(iLookTYPE%vegTypeIndex),      & ! vegetation type index
             type_data%var(iLookTYPE%soilTypeIndex),     & ! soil type
             type_data%var(iLookTYPE%slopeTypeIndex),    & ! slope type index
             10000_i4b,                                  & ! number of soil layers
             urbanVegCategory,                           & ! vegetation category for urban areas
             context)                                      ! Noah-MP parameters for this HRU

 ! overwrite the minimum resistance
 if(overwriteRSMIN) context%RSMIN = mpar_data%var(iLookPARAM%minStomatalResistance)%dat(1)

 ! overwrite the vegetation height
 context%HVT = mpar_data%var(iLookPARAM%heightCanopyTop)%dat(1)
 context%HVB = mpar_data%var(iLookPARAM%heightCanopyBottom)%dat(1)

 ! overwrite the tables for LAI and SAI
 if(model_decisions(iLookDECISIONS%LAI_method)%iDecision == specified)then
   context%SAIM(:) = mpar_data%var(iLookPARAM%winterSAI)%dat(1)
   context%LAIM(:) = mpar_data%var(iLookPARAM%summerLAI)%dat(1)*greenVegFrac_monthly
 end if

 end subroutine noahmpContext

end module noahmpContext_module

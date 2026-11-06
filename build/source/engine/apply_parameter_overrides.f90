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

module parameter_overrides_module

USE nr_type, only: i4b,rkind
USE globalData, only: integerMissing

USE data_types, only: gru_double
USE data_types, only: gru_hru_dom_doubleVec

implicit none
private

public :: apply_parameter_overrides
public :: split_domain_target

contains

 ! ************************************************************************************************
 ! public subroutine apply_overrides: apply user-specified parameter values
 ! ************************************************************************************************
 subroutine apply_parameter_overrides(nGRU_local,param_name,param_value, &
                                      mparStruct,bparStruct,err,message)
 USE get_ixname_module,only:get_ixParam,get_ixBpar           ! access function to find index of elements in structure
 USE globalData,only:gru_struc                               ! mapping from GRUs to HRUs
 implicit none

 integer(i4b)            , intent(in)    :: nGRU_local       ! number of GRUs assigned to this rank
 character(*)            , intent(in)    :: param_name(:)    ! parameter names specified on the command line
 real(rkind)             , intent(in)    :: param_value(:)   ! parameter values specified on the command line
 type(gru_hru_dom_doubleVec) , intent(inout) :: mparStruct   ! model parameters for each local HRU and domain
 type(gru_double)        , intent(inout) :: bparStruct       ! basin parameters for each local GRU
 integer(i4b)            , intent(out)   :: err              ! error code
 character(*)            , intent(out)   :: message          ! error message
 ! local variables
 integer(i4b)                            :: iParam           ! parameter index in command-line arrays
 integer(i4b)                            :: ixParam          ! index of parameter in data structure
 integer(i4b)                            :: iGRU             ! local GRU index
 integer(i4b)                            :: iHRU             ! local HRU index
 integer(i4b)                            :: iDOM             ! local domain index
 character(len=64)                       :: name             ! parameter name without its domain type
 integer(i4b)                            :: domType          ! domain type the value applies to, 0 for every domain

 ! Start procedure here
 err=0; message = "apply_parameter_overrides/"

 ! loop through the parameters specified on the command line
 do iParam=1,size(param_name)

   ! a name@domainType value goes only to the domains of that type
   call split_domain_target(param_name(iParam),name,domType)
   if(domType==integerMissing)then
     message=trim(message)//'unknown domain type in "'//trim(param_name(iParam))//'" (use upland, glacCln1, glacCln2, glacDbr, or wetland)'
     err=20; return
   endif

   ! check for a local HRU parameter
   ixParam=get_ixParam(trim(name))
   if(ixParam/=integerMissing)then

     ! apply parameter value to every domain of every local HRU, or to those of the named type
     do iGRU=1,nGRU_local
       do iHRU=1,gru_struc(iGRU)%hruCount
         do iDOM=1,gru_struc(iGRU)%hruInfo(iHRU)%domCount
           if(domType>0 .and. gru_struc(iGRU)%hruInfo(iHRU)%domInfo(iDOM)%dom_type/=domType) cycle
           mparStruct%gru(iGRU)%hru(iHRU)%dom(iDOM)%var(ixParam)%dat(:)=param_value(iParam)
         enddo
       enddo
     enddo
     cycle
   endif

   ! check for a basin parameter, which has no domains
   ixParam=get_ixBpar(trim(name))
   if(ixParam/=integerMissing)then
     if(domType>0)then
       message=trim(message)//'basin parameter "'//trim(name)//'" has no domain type to target'
       err=20; return
     endif

     ! apply parameter value to all local GRUs
     do iGRU=1,nGRU_local
       bparStruct%gru(iGRU)%var(ixParam)=param_value(iParam)
     enddo
     cycle
   endif

   ! parameter name was not recognized
   message=trim(message)//'unable to identify parameter "'//trim(param_name(iParam))//'"'
   err=20; return
 enddo

 end subroutine apply_parameter_overrides
 
 ! ************************************************************************************************
 ! public subroutine split_domain_target: split "name@domainType" into the parameter name and the domain type
 ! ************************************************************************************************
 ! domType is 0 with no suffix (every domain), and integerMissing for a suffix that names no domain type
 subroutine split_domain_target(target, name, domType)
 USE globalData,only:upland,glacCln1,glacCln2,glacDbr,wetland ! horizontal domain types
 implicit none
 character(*), intent(in)  :: target   ! parameter name, optionally followed by @ and a domain type
 character(*), intent(out) :: name     ! parameter name
 integer(i4b), intent(out) :: domType  ! domain type the value applies to, 0 for every domain
 integer(i4b)              :: iAt      ! position of the @
 iAt = index(target,'@')
 if(iAt==0)then
   name = trim(target); domType = 0; return
 endif
 name = target(1:iAt-1)
 select case(trim(target(iAt+1:)))
   case('upland');   domType = upland
   case('glacCln1'); domType = glacCln1
   case('glacCln2'); domType = glacCln2
   case('glacDbr');  domType = glacDbr
   case('wetland');  domType = wetland
   case default;     domType = integerMissing
 end select
 end subroutine split_domain_target

end module parameter_overrides_module

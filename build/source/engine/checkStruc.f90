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

module checkStruc_module
USE nrtype
implicit none
private
public::checkStruc
! define missing values
integer(i4b),parameter :: missingInteger=-9999
real(dp),parameter     :: missingDouble=-9999._dp
contains


 ! ************************************************************************************************
 ! public subroutine checkStruc: check data structures 
 ! ************************************************************************************************
 subroutine checkStruc(err,message)
 ! ascii utilities
 USE ascii_util_module,only:split_line
 ! number of variables in each data structure
 USE var_lookup,only:maxvarTime,maxvarForc,maxvarAttr,maxvarType    ! maximum number variables in each data structure
 USE var_lookup,only:maxvarState,maxvarDiag,maxvarFlux,maxvarDeriv  ! maximum number variables in each data structure
 USE var_lookup,only:maxvarMpar,maxvarMvar,maxvarIndx               ! maximum number variables in each data structure
 USE var_lookup,only:maxvarBpar,maxvarBvar                          ! maximum number variables in each data structure
 ! metadata structures
 USE data_struc,only:time_meta,forc_meta,attr_meta,type_meta        ! metadata structures
 USE data_struc,only:state_meta,diag_meta,flux_meta,deriv_meta      ! metadata structures
 USE data_struc,only:mpar_meta,mvar_meta,indx_meta                  ! metadata structures
 USE data_struc,only:bpar_meta,bvar_meta                            ! metadata structures
 ! named variables defining strructure elements
 USE var_lookup,only:iLookTIME,iLookFORCE,iLookATTR,iLookTYPE       ! named variables showing the elements of each data structure
 USE var_lookup,only:iLookSTATE,iLookDIAG,iLookFLUX,iLookDERIV      ! named variables showing the elements of each data structure
 USE var_lookup,only:iLookPARAM,iLookMVAR,iLookINDEX                ! named variables showing the elements of each data structure
 USE var_lookup,only:iLookBPAR,iLookBVAR                            ! named variables showing the elements of each data structure
 implicit none
 ! dummy variables
 integer(i4b),intent(out)             :: err         ! error code
 character(*),intent(out)             :: message     ! error message
 ! local variables
 integer(i4b)                         :: iStruct     ! index of data structure
 integer(i4b),parameter               :: nStruct=13  ! number of data structures
 character(len=8192)                  :: longString  ! string containing the indices defined in the structure constructor
 character(len=32),allocatable        :: words(:)    ! vector of words extracted from the long string
 integer(i4b)                         :: ix          ! index of the variable in the data structure
 integer(i4b)                         :: ixTest      ! test the structure constructor = (1,2,3,...,nVar)
 character(len=256)                   :: cmessage    ! error message of downwind routine
 ! -----------------------------------------------------------------------------------------------------------------------------------
 ! data structure information
 type info ! data structure information
  character(len=32)                   :: structName  ! name of the data structure
  character(len=32)                   :: lookName    ! name of the look-up variables
  integer(i4b)                        :: nVar        ! number of variables in each data structure
 end type info
 ! populate structure information
 type(info),parameter,dimension(nStruct) :: structInfo=(/&
                                             info('time',  'TIME' , maxvarTime ), & ! the time data structure   
                                             info('forc',  'FORCE', maxvarForc ), & ! the forcing data structure
                                             info('attr',  'ATTR' , maxvarAttr ), & ! the attribute data structure
                                             info('type',  'TYPE' , maxvarType ), & ! the type data structure
                                             info('mpar',  'PARAM', maxvarMpar ), & ! the model parameter data structure
                                             info('mvar',  'MVAR' , maxvarMvar ), & ! the model variable data structure
                                             info('bpar',  'BPAR' , maxvarBpar ), & ! the basin parameter data structure
                                             info('bvar',  'BVAR' , maxvarBvar ), & ! the basin variable data structure
                                             info('indx',  'INDEX', maxvarIndx ), & ! the model index data structure
                                             info('state', 'STATE', maxvarState), & ! the state variable data structure
                                             info('diag',  'DIAG' , maxvarDiag ), & ! the diagnostic variable data structure
                                             info('flux',  'FLUX' , maxvarFlux ), & ! the flux data structure
                                             info('deriv', 'DERIV', maxvarDeriv) /) ! the model derivative data structure
 ! -----------------------------------------------------------------------------------------------------------------------------------
 ! initialize errors
 err=0; message="checkStruc/"

 ! -----
 ! * check that the structure constructors are correct...
 ! ------------------------------------------------------

 ! loop through data structures 
 do iStruct=1,nStruct
  ! convert the lookup structures to a character string
  select case(trim(structInfo(iStruct)%structName))
   case('time');  write(longString,*) iLookTIME
   case('forc');  write(longString,*) iLookFORCE
   case('attr');  write(longString,*) iLookATTR
   case('type');  write(longString,*) iLookTYPE
   case('mpar');  write(longString,*) iLookPARAM
   case('mvar');  write(longString,*) iLookMVAR
   case('bpar');  write(longString,*) iLookBPAR
   case('bvar');  write(longString,*) iLookBVAR
   case('indx');  write(longString,*) iLookINDEX
   case('state'); write(longString,*) iLookSTATE
   case('diag');  write(longString,*) iLookDIAG
   case('flux');  write(longString,*) iLookFLUX
   case('deriv'); write(longString,*) iLookDERIV
   case default; err=20; message=trim(message)//'unable to identify lookup structure'; return
  end select
  ! convert the string to a character vector
  call split_line(longString,words,err,cmessage)
  if(err/=0)then; message=trim(message)//trim(cmessage); return; endif
  if(size(words)/=structInfo(iStruct)%nVar)then; err=20; message=trim(message)//'unexpected number of elements'; return; endif
  ! check that the integer value is appropriate
  do ix=1,structInfo(iStruct)%nVar
   read(words(ix),*) ixTest  ! convert character to integer; store in ixTest
   if(ixTest/=ix)then
    write(message,'(a,i0,a)')trim(message)//'problem with structure constructor iLook'//trim(structInfo(iStruct)%lookName)//' [element=',ix,']'
    err=20; return
   endif
  end do
 end do  ! looping through data structures

 ! -----
 ! * check that the metadata is fully populated...
 ! -----------------------------------------------

 ! loop through data structures
 do iStruct=1,nStruct
  ! check that the metadata is fully populated 
  select case(trim(structInfo(iStruct)%structName))
   case('time');  call checkPopulated(iStruct,time_meta,err,cmessage)
   case('forc');  call checkPopulated(iStruct,forc_meta,err,cmessage) 
   case('attr');  call checkPopulated(iStruct,attr_meta,err,cmessage) 
   case('type');  call checkPopulated(iStruct,type_meta,err,cmessage) 
   case('mpar');  call checkPopulated(iStruct,mpar_meta,err,cmessage) 
   case('mvar');  call checkPopulated(iStruct,mvar_meta,err,cmessage) 
   case('bpar');  call checkPopulated(iStruct,bpar_meta,err,cmessage) 
   case('bvar');  call checkPopulated(iStruct,bvar_meta,err,cmessage) 
   case('indx');  call checkPopulated(iStruct,indx_meta,err,cmessage) 
   case('state'); call checkPopulated(iStruct,state_meta,err,cmessage) 
   case('diag');  call checkPopulated(iStruct,diag_meta,err,cmessage) 
   case('flux');  call checkPopulated(iStruct,flux_meta,err,cmessage) 
   case('deriv'); call checkPopulated(iStruct,deriv_meta,err,cmessage) 
   case default; err=20; message=trim(message)//'unable to identify lookup structure'; return
  end select
  if(err/=0)then; message=trim(message)//trim(cmessage); return; endif  ! (check for errors) 
 end do  ! looping through data structures


 contains

  ! ************************************************************************************************
  ! internal subroutine checkPopulated: check that the metadata is fully populated...
  ! ************************************************************************************************
  subroutine checkPopulated(iStruct,metadata,err,message)
  ! access the data type for the metadata structures
  USE data_struc,only:var_info 
  ! get index from character string
  USE get_ixname_module,only: get_ixtime
  USE get_ixname_module,only: get_ixattr
  USE get_ixname_module,only: get_ixtype
  USE get_ixname_module,only: get_ixforce
  USE get_ixname_module,only: get_ixparam
  USE get_ixname_module,only: get_ixmvar
  USE get_ixname_module,only: get_ixindex
  USE get_ixname_module,only: get_ixbpar
  USE get_ixname_module,only: get_ixbvar
  USE get_ixname_module,only: get_ixstate
  USE get_ixname_module,only: get_ixdiag
  USE get_ixname_module,only: get_ixflux
  USE get_ixname_module,only: get_ixderiv
  implicit none
  ! dummy variables
  integer(i4b),intent(in)   :: iStruct     ! index of data structure
  type(var_info)            :: metadata(:) ! metadata structure 
  integer(i4b),intent(out)  :: err         ! error code
  character(*),intent(out)  :: message     ! error message
  ! local variables
  integer(i4b)              :: iVar        ! index of variable within a data structure
  integer(i4b)              :: jVar        ! index of variable within a data structure (returned from the variable name)
  integer(i4b)              :: jStruct     ! index of data structure
  ! initialize error control
  err=0; message='checkPopulated/'
 
  ! loop through variables
  do iVar=1,size(metadata)

   ! check that the variable is populated
   if(len(trim(metadata(iVar)%varname))==0)then
    write(message,'(a,i0,a)') trim(message)//trim(structInfo(iStruct)%structName)//'_meta structure is not populated for named variable # ',iVar, ' in structure iLook'//trim(structInfo(iStruct)%lookName)
    err=20; return
   endif

   ! check that the index-from-name lookup returns the correct variable
   do jStruct=1,nStruct

    ! (identify if the variable exists in a given structure)
    select case(trim(structInfo(jStruct)%structName))
     case('time');  jVar = get_ixtime(trim(metadata(iVar)%varname))
     case('forc');  jVar = get_ixforce(trim(metadata(iVar)%varname))  
     case('attr');  jVar = get_ixattr(trim(metadata(iVar)%varname)) 
     case('type');  jVar = get_ixtype(trim(metadata(iVar)%varname)) 
     case('mpar');  jVar = get_ixparam(trim(metadata(iVar)%varname))  
     case('mvar');  jVar = get_ixmvar(trim(metadata(iVar)%varname)) 
     case('bpar');  jVar = get_ixbpar(trim(metadata(iVar)%varname)) 
     case('bvar');  jVar = get_ixbvar(trim(metadata(iVar)%varname)) 
     case('indx');  jVar = get_ixindex(trim(metadata(iVar)%varname)) 
     case('state'); jVar = get_ixstate(trim(metadata(iVar)%varname)) 
     case('diag');  jVar = get_ixdiag(trim(metadata(iVar)%varname)) 
     case('flux');  jVar = get_ixflux(trim(metadata(iVar)%varname)) 
     case('deriv'); jVar = get_ixderiv(trim(metadata(iVar)%varname)) 
     case default; err=20; message=trim(message)//'unable to identify lookup structure'; return
    end select
    if(jVar>0)then   ! found the variable in structure jStruct

     ! --> check that the variable is in the correct structure
     if(jStruct/=iStruct)then
      if(trim(structInfo(iStruct)%structName)=='mvar' .or. trim(structInfo(jStruct)%structName)=='mvar') cycle
      message=trim(message)//'variable '//trim(metadata(iVar)%varname)//' from structure '//trim(structInfo(iStruct)%structName)//'_meta is in structure '//trim(structInfo(jStruct)%structName)//'_meta'
      err=20; return

     ! in the correct structure
     else

      ! --> check that the variable index is correct
      if(jVar/=iVar)then
       write(message,'(a,i0,a,i0,a)') trim(message)//'variable '//trim(metadata(iVar)%varname)//' has index ', iVar, ' (expect index ', jVar, '); problem possible in popMetadat, get_ix'//trim(structInfo(iStruct)%structName)//', or var_lookup'
       err=20; return
      endif

     endif  ! variable found in the correct structure

    ! variable does not exist in structure jStruct
    else

     ! --> check that we found the variable
     if(iStruct==jStruct)then
      message = trim(message)//'cannot find variable '//trim(metadata(iVar)%varname)//' in structure '//trim(structInfo(iStruct)%structName)//'_meta; you need to add variable to get_ix'//trim(structInfo(iStruct)%structName)
      err=20; return
     endif

    endif  ! if the variable exists in structure jStruct

   end do  ! looping through data structures

  end do  ! looping through variables in structure iStruct

  end subroutine checkPopulated

 end subroutine checkStruc

end module checkStruc_module

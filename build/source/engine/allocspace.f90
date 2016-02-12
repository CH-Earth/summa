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

module allocspace_module
USE nrtype
implicit none
private
public::init_metad
public::alloc_stim
public::alloc_time
public::alloc_forc
public::alloc_attr
public::alloc_type
public::alloc_mpar
public::alloc_mvar
public::alloc_indx
public::alloc_bpar
public::alloc_bvar
! define missing values
integer(i4b),parameter :: missingInteger=-9999
real(dp),parameter     :: missingDouble=-9999._dp
contains


 ! ************************************************************************************************
 ! public subroutine init_metad: initialize metadata structures
 ! ************************************************************************************************
 subroutine init_metad(err,message)
 ! ascii utilities
 USE ascii_util_module,only:split_line
 ! number of variables in each data structure
 USE var_lookup,only:maxvarTime,maxvarForc,maxvarAttr,maxvarType    ! maximum number variables in each data structure
 USE var_lookup,only:maxvarMpar,maxvarMvar,maxvarIndx               ! maximum number variables in each data structure
 USE var_lookup,only:maxvarDiag,maxvarFlux,maxvarDeriv              ! maximum number variables in each data structure
 USE var_lookup,only:maxvarBpar,maxvarBvar                          ! maximum number variables in each data structure
 ! metadata structures
 USE data_struc,only:time_meta,forc_meta,attr_meta,type_meta        ! metadata structures
 USE data_struc,only:mpar_meta,mvar_meta,indx_meta                  ! metadata structures
 USE data_struc,only:diag_meta,flux_meta,deriv_meta                 ! metadata structures
 USE data_struc,only:bpar_meta,bvar_meta                            ! metadata structures
 ! named variables defining strructure elements
 USE var_lookup,only:iLookTIME,iLookFORCE,iLookATTR,iLookTYPE       ! named variables showing the elements of each data structure
 USE var_lookup,only:iLookPARAM,iLookMVAR,iLookINDEX                ! named variables showing the elements of each data structure
 USE var_lookup,only:iLookDIAG,iLookFLUX,iLookDERIV                 ! named variables showing the elements of each data structure
 USE var_lookup,only:iLookBPAR,iLookBVAR                            ! named variables showing the elements of each data structure
 implicit none
 ! dummy variables
 integer(i4b),intent(out)             :: err         ! error code
 character(*),intent(out)             :: message     ! error message
 ! define look-up variables for data structures
 integer(i4b)                         :: iStruct     ! loop through data structures
 integer(i4b),parameter               :: nStruct=12  ! number of data structures
 integer(i4b),parameter               :: ixTime=1    ! the time data structure
 integer(i4b),parameter               :: ixForc=2    ! the time data structure
 integer(i4b),parameter               :: ixAttr=3    ! the time data structure
 integer(i4b),parameter               :: ixType=4    ! the time data structure
 integer(i4b),parameter               :: ixMpar=5    ! the time data structure
 integer(i4b),parameter               :: ixMvar=6    ! the time data structure
 integer(i4b),parameter               :: ixBpar=7    ! the time data structure
 integer(i4b),parameter               :: ixBvar=8    ! the time data structure
 integer(i4b),parameter               :: ixIndx=9    ! the time data structure
 integer(i4b),parameter               :: ixDiag=10   ! the time data structure
 integer(i4b),parameter               :: ixFlux=11   ! the time data structure
 integer(i4b),parameter               :: ixDeriv=12  ! the time data structure
 ! check that the structure constructors are correct
 character(len=8192)                  :: longString  ! string containing the indices defined in the structure constructor
 character(len=32),allocatable        :: words(:)    ! vector of words extracted from the long string
 character(len=32)                    :: cTry        ! character string of trial structure
 integer(i4b)                         :: ix          ! index of the variable in the data structure
 integer(i4b)                         :: ixTest      ! test the structure constructor = (1,2,3,...,nVar)
 integer(i4b)                         :: nVar        ! the number of variables in each data structure
 ! error control
 character(len=256)                   :: cmessage    ! error message of downwind routine
 ! initialize errors
 err=0; message="init_model/"

 ! ensure metadata structures are deallocated
 if (associated(time_meta))  deallocate(time_meta)   ! 1
 if (associated(forc_meta))  deallocate(forc_meta)   ! 2
 if (associated(attr_meta))  deallocate(attr_meta)   ! 3
 if (associated(type_meta))  deallocate(type_meta)   ! 4
 if (associated(mpar_meta))  deallocate(mpar_meta)   ! 5
 if (associated(mvar_meta))  deallocate(mvar_meta)   ! 6
 if (associated(bpar_meta))  deallocate(bpar_meta)   ! 7
 if (associated(bvar_meta))  deallocate(bvar_meta)   ! 8
 if (associated(indx_meta))  deallocate(indx_meta)   ! 9
 if (associated(diag_meta))  deallocate(diag_meta)   ! 10
 if (associated(flux_meta))  deallocate(flux_meta)   ! 11
 if (associated(deriv_meta)) deallocate(deriv_meta)  ! 12

 ! allocate metadata structures
 allocate(time_meta(maxvarTime),forc_meta(maxvarForc),attr_meta(maxvarAttr),type_meta(maxvarType),&
          mpar_meta(maxvarMpar),mvar_meta(maxvarMvar),indx_meta(maxvarIndx),&
          diag_meta(maxvarDiag),flux_meta(maxvarFlux),deriv_meta(maxvarDeriv),&
          bpar_meta(maxvarBpar),bvar_meta(maxvarBvar),stat=err)
 if(err/=0)then; err=20; message=trim(message)//"problemAllocateMetadata"; return; endif

 ! check that the structure constructors are correct
 do iStruct=1,nStruct
  ! convert the lookup structures to a character string
  select case(iStruct)
   case(ixTime);  write(longString,*) iLookTIME;  nVar=maxvarTime; cTry='iLookTIME'
   case(ixForc);  write(longString,*) iLookFORCE; nVar=maxvarForc; cTry='iLookFORCE'
   case(ixAttr);  write(longString,*) iLookATTR;  nVar=maxvarAttr; cTry='iLookATTR'
   case(ixType);  write(longString,*) iLookTYPE;  nVar=maxvarType; cTry='iLookTYPE'
   case(ixMpar);  write(longString,*) iLookPARAM; nVar=maxvarMpar; cTry='iLookMPAR'
   case(ixMvar);  write(longString,*) iLookMVAR;  nVar=maxvarMvar; cTry='iLookMVAR'
   case(ixBpar);  write(longString,*) iLookBPAR;  nVar=maxvarBpar; cTry='iLookBPAR'
   case(ixBvar);  write(longString,*) iLookBVAR;  nVar=maxvarBvar; cTry='iLookBVAR'
   case(ixIndx);  write(longString,*) iLookINDEX; nVar=maxvarIndx; cTry='iLookINDX'
   case(ixDiag);  write(longString,*) iLookDIAG;  nVar=maxvarDiag; cTry='iLookDIAG'
   case(ixFlux);  write(longString,*) iLookFLUX;  nVar=maxvarFlux; cTry='iLookFLUX'
   case(ixDeriv); write(longString,*) iLookDERIV; nVar=maxvarDeriv; cTry='iLookDERIV'
   case default; err=20; message=trim(message)//'unable to identify lookup structure'; return
  end select
  ! convert the string to a character vector
  call split_line(longString,words,err,cmessage)
  if(err/=0)then; message=trim(message)//trim(cmessage); return; endif
  if(size(words)/=nVar)then; err=20; message=trim(message)//'unexpected number of elements'; return; endif
  ! check that the integer value is appropriate
  do ix=1,nVar
   read(words(ix),*) ixTest  ! convert character to integer; store in ixTest
   if(ixTest/=ix)then
    write(message,'(a,i0,a)')trim(message)//'problem with structure constructor '//trim(cTry)//' [element=',ix,']'
    err=20; return
   endif
  end do
 end do  ! looping through data structures

 end subroutine init_metad

 ! ************************************************************************************************
 ! public subroutine alloc_stim: initialize data structures for scalar time structures
 ! ************************************************************************************************
 subroutine alloc_stim(datastr,err,message)
 ! used to initialize structure components for model variables
 USE data_struc,only:var_i,time_meta                 ! data structures
 implicit none
 ! dummy variables
 type(var_i),intent(out),pointer      :: datastr     ! data structure to allocate
 integer(i4b),intent(out)             :: err         ! error code
 character(*),intent(out)             :: message     ! error message
 ! initialize errors
 err=0; message="alloc_stim/"
 ! check that the metadata structure is allocated
 if(.not.associated(time_meta))then
  err=10; message=trim(message)//"metadataNotInitialized"; return
 endif
 ! initialize top-level data structure
 if(associated(datastr)) deallocate(datastr)
 allocate(datastr,stat=err)
 if(err/=0)then; err=20; message=trim(message)//"problemAllocateDataTopLevel"; return; endif
 ! initialize second level data structure
 allocate(datastr%var(size(time_meta)),stat=err)
 if(err/=0)then; err=20; message=trim(message)//"problemAllocateData2ndLevel"; return; endif
 ! set values to missing
 datastr%var(:) = missingInteger
 end subroutine alloc_stim

 ! ************************************************************************************************
 ! public subroutine alloc_time: initialize data structures for time structures
 ! ************************************************************************************************
 subroutine alloc_time(nHRU,err,message)
 ! used to initialize structure components for model variables
 USE data_struc,only:time_hru,time_meta              ! data structures
 implicit none
 ! dummy variables
 integer(i4b),intent(in)              :: nHRU        ! number of HRUs
 integer(i4b),intent(out)             :: err         ! error code
 character(*),intent(out)             :: message     ! error message
 ! local variables
 integer(i4b)                         :: iHRU        ! loop through HRUs
 integer(i4b)                         :: nVar        ! number of variables
 ! initialize errors
 err=0; message="alloc_time/"
 ! check that the metadata structure is allocated
 if(.not.associated(time_meta))then
  err=10; message=trim(message)//"metadataNotInitialized"; return
 endif
 ! initialize top-level data structure
 if(associated(time_hru)) deallocate(time_hru)
 allocate(time_hru(nHRU),stat=err)
 if(err/=0)then; err=20; message=trim(message)//"problemAllocateDataTopLevel"; return; endif
 ! initialize second level data structure
 nVar = size(time_meta)
 do iHRU=1,nHRU
  allocate(time_hru(iHRU)%var(nVar),stat=err)
  if(err/=0)then; err=20; message=trim(message)//"problemAllocateData2ndLevel"; return; endif
  ! set values to missing
  time_hru(iHRU)%var(:) = missingInteger
 end do
 end subroutine alloc_time


 ! ************************************************************************************************
 ! public subroutine alloc_forc: initialize data structures for model forcing data
 ! ************************************************************************************************
 subroutine alloc_forc(nHRU,err,message)
 ! used to initialize structure components for model variables
 USE data_struc,only:forc_hru,forc_meta             ! data structures
 implicit none
 ! dummy variables
 integer(i4b),intent(in)              :: nHRU        ! number of HRUs
 integer(i4b),intent(out)             :: err         ! error code
 character(*),intent(out)             :: message     ! error message
 ! local variables
 integer(i4b)                         :: iHRU        ! loop through HRUs
 integer(i4b)                         :: nVar        ! number of variables
 ! initialize errors
 err=0; message="alloc_forc/"
 ! check that the metadata structure is allocated
 if(.not.associated(forc_meta))then
  err=10; message=trim(message)//"metadataNotInitialized"; return
 endif
 ! initialize top-level data structure
 if(associated(forc_hru)) deallocate(forc_hru)
 allocate(forc_hru(nHRU),stat=err)
 if(err/=0)then; err=20; message=trim(message)//"problemAllocateDataTopLevel"; return; endif
 ! initialize second level data structure
 nVar = size(forc_meta)
 do iHRU=1,nHRU
  allocate(forc_hru(iHRU)%var(nVar),stat=err)
  if(err/=0)then; err=20; message=trim(message)//"problemAllocateData2ndLevel"; return; endif
  ! set values to missing
  forc_hru(iHRU)%var(:) = missingDouble
 end do
 end subroutine alloc_forc


 ! ************************************************************************************************
 ! public subroutine alloc_attr: initialize data structures for local attributes
 ! ************************************************************************************************
 subroutine alloc_attr(nHRU,err,message)
 ! used to initialize structure components for model variables
 USE data_struc,only:attr_meta,attr_hru   ! data structures
 implicit none
 ! dummy variables
 integer(i4b),intent(in)              :: nHRU        ! number of HRUs
 integer(i4b),intent(out)             :: err         ! error code
 character(*),intent(out)             :: message     ! error message
 ! local variables
 integer(i4b)                         :: iHRU        ! loop through HRUs
 integer(i4b)                         :: nVar        ! number of variables
 ! initialize errors
 err=0; message="alloc_attr/"
 ! check that the metadata structure is allocated
 if(.not.associated(attr_meta))then
  err=10; message=trim(message)//"metadataNotInitialized"; return
 endif
 ! initialize top-level data structure
 if(associated(attr_hru)) deallocate(attr_hru)
 allocate(attr_hru(nHRU),stat=err)
 if(err/=0)then; err=20; message=trim(message)//"problemAllocateDataTopLevel"; return; endif
 ! initialize second level data structure
 nVar = size(attr_meta)
 do iHRU=1,nHRU
  allocate(attr_hru(iHRU)%var(nVar),stat=err)
  if(err/=0)then; err=20; message=trim(message)//"problemAllocateData2ndLevel"; return; endif
  ! fill data with missing values
  attr_hru(iHRU)%var(:) = missingDouble
 end do
 end subroutine alloc_attr


 ! *************************************************************************************************
 ! public subroutine alloc_type: initialize data structures for local classification of veg, soil, etc.
 ! *************************************************************************************************
 subroutine alloc_type(nHRU,err,message)
 ! used to initialize structure components for model variables
 USE data_struc,only:type_hru,type_meta             ! data structures
 implicit none
 ! dummy variables
 integer(i4b),intent(in)              :: nHRU        ! number of HRUs
 integer(i4b),intent(out)             :: err         ! error code
 character(*),intent(out)             :: message     ! error message
 ! local variables
 integer(i4b)                         :: iHRU        ! loop through HRUs
 integer(i4b)                         :: nVar        ! number of variables
 ! initialize errors
 err=0; message="f-alloc_type/"
 ! check that the metadata structure is allocated
 if(.not.associated(type_meta))then
  err=10; message=trim(message)//"metadataNotInitialized"; return
 endif
 ! initialize top-level data structure
 if(associated(type_hru)) deallocate(type_hru)
 allocate(type_hru(nHRU),stat=err)
 if(err/=0)then; err=20; message=trim(message)//"problemAllocateDataTopLevel"; return; endif
 ! initialize second level data structure
 nVar = size(type_meta)
 do iHRU=1,nHRU
  allocate(type_hru(iHRU)%var(nVar),stat=err)
  if(err/=0)then; err=20; message=trim(message)//"problemAllocateData2ndLevel"; return; endif
  ! fill data with missing values
  type_hru(iHRU)%var(:) = missingInteger
 end do
 end subroutine alloc_type


 ! *************************************************************************************************
 ! public subroutine alloc_mpar: initialize data structures for model parameters
 ! *************************************************************************************************
 subroutine alloc_mpar(nHRU,err,message)
 ! used to initialize structure components for model variables
 USE data_struc,only:mpar_hru,mpar_meta             ! data structures
 implicit none
 ! dummy variables
 integer(i4b),intent(in)              :: nHRU        ! number of HRUs
 integer(i4b),intent(out)             :: err         ! error code
 character(*),intent(out)             :: message     ! error message
 ! local variables
 integer(i4b)                         :: iHRU        ! loop through HRUs
 integer(i4b)                         :: nPar        ! number of parameters
 ! initialize errors
 err=0; message="f-alloc_mpar/"
 ! check that the metadata structure is allocated
 if(.not.associated(mpar_meta))then
  err=10; message=trim(message)//"metadataNotInitialized"; return
 endif
 ! initialize top-level data structure
 if(associated(mpar_hru)) deallocate(mpar_hru)
 allocate(mpar_hru(nHRU),stat=err)
 if(err/=0)then; err=20; message=trim(message)//"problemAllocateDataTopLevel"; return; endif
 ! get the number of parameters
 nPar = size(mpar_meta)
 ! loop through HRUs
 do iHRU=1,nHRU
  ! initialize second level data structure
  allocate(mpar_hru(iHRU)%var(nPar),stat=err)
  if(err/=0)then; err=20; message=trim(message)//"problemAllocateData2ndLevel"; return; endif
  ! set values to missing
  mpar_hru(iHRU)%var(:) = missingDouble
 end do  ! looping through HRUs
 end subroutine alloc_mpar


 ! *************************************************************************************************
 ! public subroutine alloc_mvar: initialize data structures for model variables
 ! *************************************************************************************************
 subroutine alloc_mvar(nHRU,err,message)
 ! used to initialize structure components for model variables
 USE data_struc,only:mvar_hru,mvar_meta             ! data structures
 implicit none
 ! dummy variables
 integer(i4b),intent(in)              :: nHRU        ! number of HRUs
 integer(i4b),intent(out)             :: err         ! error code
 character(*),intent(out)             :: message     ! error message
 ! local variables
 integer(i4b)                         :: iHRU        ! loop through HRUs
 integer(i4b)                         :: nVar        ! number of variables
 ! initialize errors
 err=0; message="alloc_mvar/"
 ! check that the metadata structure is allocated
 if(.not.associated(mvar_meta))then
  err=10; message=trim(message)//"metadataNotInitialized"; return
 endif
 ! initialize top-level data structure
 if(associated(mvar_hru)) deallocate(mvar_hru)
 allocate(mvar_hru(nHRU),stat=err)
 if(err/=0)then; err=20; message=trim(message)//"problemAllocateDataTopLevel"; return; endif
 ! initialize second level data structure
 nVar = size(mvar_meta)
 do iHRU=1,nHRU
  allocate(mvar_hru(iHRU)%var(nVar),stat=err)
  if(err/=0)then; err=20; message=trim(message)//"problemAllocateData2ndLevel"; return; endif
 end do ! (looping through the HRUs)
 end subroutine alloc_mvar


 ! *************************************************************************************************
 ! public subroutine alloc_indx: initialize structure components for model indices
 ! *************************************************************************************************
 subroutine alloc_indx(nHRU,err,message)
 ! used to initialize structure components for model variables
 USE data_struc,only:indx_hru,indx_meta             ! data structures
 implicit none
 ! dummy variables
 integer(i4b),intent(in)              :: nHRU        ! number of HRUs
 integer(i4b),intent(out)             :: err         ! error code
 character(*),intent(out)             :: message     ! error message
 ! local variables
 integer(i4b)                         :: iHRU        ! loop through HRUs
 integer(i4b)                         :: nVar        ! number of variables
 ! initialize errors
 err=0; message="alloc_indx/"
 ! check that the metadata structure is allocated
 if(.not.associated(indx_meta))then
  err=10; message=trim(message)//"metadataNotInitialized"; return
 endif
 ! initialize top-level data structure
 if(associated(indx_hru)) deallocate(indx_hru)
 allocate(indx_hru(nHRU),stat=err)
 if(err/=0)then; err=20; message=trim(message)//"problemAllocateDataTopLevel"; return; endif
 ! initialize second level data structure
 nVar = size(indx_meta)
 do iHRU=1,nHRU
  allocate(indx_hru(iHRU)%var(nVar),stat=err)
  if(err/=0)then; err=20; message=trim(message)//"problemAllocateData2ndLevel"; return; endif
 end do ! (looping through HRUs in the data structure)
 end subroutine alloc_indx


 ! *************************************************************************************************
 ! public subroutine alloc_bpar: initialize data structures for basin-average model parameters
 ! *************************************************************************************************
 subroutine alloc_bpar(err,message)
 ! used to initialize structure components for model variables
 USE data_struc,only:bpar_data,bpar_meta             ! data structures
 implicit none
 ! dummy variables
 integer(i4b),intent(out)             :: err         ! error code
 character(*),intent(out)             :: message     ! error message
 ! local variables
 integer(i4b)                         :: nPar        ! number of parameters
 ! initialize errors
 err=0; message="alloc_bpar/"
 ! check that the metadata structure is allocated
 if(.not.associated(bpar_meta))then
  err=10; message=trim(message)//"metadataNotInitialized"; return
 endif
 ! initialize top-level data structure
 if(associated(bpar_data)) deallocate(bpar_data)
 allocate(bpar_data,stat=err)
 if(err/=0)then; err=20; message=trim(message)//"problemAllocateDataTopLevel"; return; endif
 ! get the number of parameters
 nPar = size(bpar_meta)
 ! initialize second level data structure
 allocate(bpar_data%var(nPar),stat=err)
 if(err/=0)then; err=20; message=trim(message)//"problemAllocateData2ndLevel"; return; endif
 ! set values to missing
 bpar_data%var(:) = missingDouble
 end subroutine alloc_bpar


 ! *************************************************************************************************
 ! public subroutine alloc_bvar: initialize data structures for basin-average model variables
 ! *************************************************************************************************
 subroutine alloc_bvar(err,message)
 ! used to initialize structure components for model variables
 USE data_struc,only:bvar_data,bvar_meta             ! data structures
 implicit none
 ! dummy variables
 integer(i4b),intent(out)             :: err         ! error code
 character(*),intent(out)             :: message     ! error message
 ! local variables
 integer(i4b)                         :: iVar        ! index of variables
 integer(i4b)                         :: nVar        ! number of variables
 integer(i4b),parameter               :: nTimeDelay=2000 ! number of elements in the time delay histogram
 ! initialize errors
 err=0; message="alloc_bvar/"
 ! check that the metadata structure is allocated
 if(.not.associated(bvar_meta))then
  err=10; message=trim(message)//"metadataNotInitialized"; return
 endif
 ! initialize top-level data structure
 if(associated(bvar_data)) deallocate(bvar_data)
 allocate(bvar_data,stat=err)
 if(err/=0)then; err=20; message=trim(message)//"problemAllocateDataTopLevel"; return; endif
 ! get the number of parameters
 nVar = size(bvar_meta)
 ! initialize second level data structure
 allocate(bvar_data%var(nVar),stat=err)
 if(err/=0)then; err=20; message=trim(message)//"problemAllocateData2ndLevel"; return; endif
 ! initialize third-level data structures
 do iVar=1,nVar
  select case(bvar_meta(ivar)%vartype)
   case('scalarv'); allocate(bvar_data%var(ivar)%dat(1),stat=err)
   case('routing'); allocate(bvar_data%var(ivar)%dat(nTimeDelay),stat=err)
   case default
    err=40; message=trim(message)//"unknownVariableType[name='"//trim(bvar_meta(ivar)%varname)//"'; &
                                   &type='"//trim(bvar_meta(ivar)%vartype)//"']"; return
  endselect
  bvar_data%var(ivar)%dat(:) = missingDouble
 end do ! (looping through model variables)
 end subroutine alloc_bvar


end module allocspace_module

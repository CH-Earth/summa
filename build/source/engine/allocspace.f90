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
! data types
USE nrtype
! number of variables in each data structure
USE var_lookup,only:maxvarTime,maxvarForc,maxvarAttr,maxvarType    ! maximum number variables in each data structure
USE var_lookup,only:maxvarState,maxvarDiag,maxvarFlux,maxvarDeriv  ! maximum number variables in each data structure
USE var_lookup,only:maxvarMpar,maxvarMvar,maxvarIndx               ! maximum number variables in each data structure
USE var_lookup,only:maxvarBpar,maxvarBvar                          ! maximum number variables in each data structure
USE var_lookup,only:maxvarDecisions                                ! maximum number of decisions
! metadata structures
USE data_struc,only:time_meta,forc_meta,attr_meta,type_meta        ! metadata structures
USE data_struc,only:state_meta,diag_meta,flux_meta,deriv_meta      ! metadata structures
USE data_struc,only:mpar_meta,mvar_meta,indx_meta                  ! metadata structures
USE data_struc,only:bpar_meta,bvar_meta                            ! metadata structures
implicit none
private
public::init_metad
public::initStruct
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
! define fixed dimensions
integer(i4b),parameter :: nBand=2         ! number of spectral bands
integer(i4b),parameter :: nTimeDelay=2000 ! number of elements in the time delay histogram
! -----------------------------------------------------------------------------------------------------------------------------------
! data structure information
integer(i4b),parameter               :: nStruct=13  ! number of data structures
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
contains


 ! ************************************************************************************************
 ! public subroutine init_metad: initialize metadata structures
 ! ************************************************************************************************
 subroutine init_metad(err,message)
 implicit none
 ! dummy variables
 integer(i4b),intent(out)             :: err         ! error code
 character(*),intent(out)             :: message     ! error message
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
 if (associated(state_meta)) deallocate(state_meta)  ! 10
 if (associated(diag_meta))  deallocate(diag_meta)   ! 11
 if (associated(flux_meta))  deallocate(flux_meta)   ! 12
 if (associated(deriv_meta)) deallocate(deriv_meta)  ! 13

 ! allocate metadata structures
 allocate(time_meta(maxvarTime),forc_meta(maxvarForc),attr_meta(maxvarAttr),type_meta(maxvarType),&
          state_meta(maxvarState),diag_meta(maxvarDiag),flux_meta(maxvarFlux),deriv_meta(maxvarDeriv),&
          mpar_meta(maxvarMpar),mvar_meta(maxvarMvar),indx_meta(maxvarIndx),&
          bpar_meta(maxvarBpar),bvar_meta(maxvarBvar),stat=err)
 if(err/=0)then; err=20; message=trim(message)//"problemAllocateMetadata"; return; endif

 end subroutine init_metad


 ! ************************************************************************************************
 ! public subroutine initStruct: initialize data structures
 ! ************************************************************************************************
 subroutine initStruct(&
                       ! input: model control
                       nHRU,       &    ! number of HRUs
                       nSnow,      &    ! number of snow layers for each HRU
                       nSoil,      &    ! number of soil layers for each HRU
                       ! input: data structures
                       timeStruct, &    ! model time data
                       forcStruct, &    ! model forcing data
                       attrStruct, &    ! local attributes for each HRU
                       typeStruct, &    ! local classification of soil veg etc. for each HRU
                       mparStruct, &    ! model parameters
                       mvarStruct, &    ! model variables
                       indxStruct, &    ! model indices
                       bparStruct, &    ! basin-average parameters
                       bvarStruct, &    ! basin-average variables
                       ! output: error control
                       err,message)
 ! provide access to the derived types to define the data structures
 USE data_struc,only:&
                     var_int,             & ! x%var(:)            (i4b)
                     var_double,          & ! x%var(:)            (dp)
                     var_intVec,          & ! x%var(:)%dat        (i4b)
                     var_doubleVec,       & ! x%var(:)%dat        (dp)
                     spatial_int,         & ! x%hru(:)%var(:)     (i4b)
                     spatial_double,      & ! x%hru(:)%var(:)     (dp)
                     spatial_intVec,      & ! x%hru(:)%var(:)%dat (i4b)
                     spatial_doubleVec      ! x%hru(:)%var(:)%dat (dp)
 implicit none
 ! input: model control
 integer(i4b),            intent(in)   :: nHRU          ! number of HRUs
 integer(i4b),            intent(in)   :: nSnow(:)      ! number of snow layers for each HRU
 integer(i4b),            intent(in)   :: nSoil(:)      ! number of soil layers for each HRU
 ! input: data structures
 type(var_int),           intent(out)  :: timeStruct    ! model time data
 type(spatial_double),    intent(out)  :: forcStruct    ! model forcing data
 type(spatial_double),    intent(out)  :: attrStruct    ! local attributes for each HRU
 type(spatial_int),       intent(out)  :: typeStruct    ! local classification of soil veg etc. for each HRU
 type(spatial_double),    intent(out)  :: mparStruct    ! model parameters
 type(spatial_doubleVec), intent(out)  :: mvarStruct    ! model variables
 type(spatial_intVec),    intent(out)  :: indxStruct    ! model indices
 type(var_double),        intent(out)  :: bparStruct    ! basin-average parameters
 type(var_doubleVec),     intent(out)  :: bvarStruct    ! basin-average variables
 ! output: error control
 integer(i4b),intent(out)              :: err           ! error code
 character(*),intent(out)              :: message       ! error message
 ! ---------------------------------------------------------------------------------------------------------------
 ! local variables
 integer(i4b)                          :: iHRU          ! index of HRUs 
 integer(i4b)                          :: iStruct       ! index of data structure 
 logical(lgt)                          :: check         ! flag to check if the structure is already allocated
 integer(i4b)                          :: nLayers       ! total number of layers in the snow+soil domian
 character(len=256)                    :: cMessage      ! error message of downwind routine
 ! initialize errors
 err=0; message="initStruct/"

 ! -----
 ! * allocate spatial dimensions...
 ! --------------------------------

 ! loop through data structures
 do iStruct=1,nStruct
  check=.false. ! initialize check
  ! allocate spatial dimension
  select case(trim(structInfo(iStruct)%structName))
   case('time'); ! do nothing for time: no spatial dimension 
   case('forc'); if(allocated(forcStruct%hru))then; check=.true.; else; allocate(forcStruct%hru(nHRU),stat=err); endif
   case('attr'); if(allocated(attrStruct%hru))then; check=.true.; else; allocate(attrStruct%hru(nHRU),stat=err); endif 
   case('type'); if(allocated(typeStruct%hru))then; check=.true.; else; allocate(typeStruct%hru(nHRU),stat=err); endif 
   case('mpar'); if(allocated(mparStruct%hru))then; check=.true.; else; allocate(mparStruct%hru(nHRU),stat=err); endif
   case('mvar'); if(allocated(mvarStruct%hru))then; check=.true.; else; allocate(mvarStruct%hru(nHRU),stat=err); endif
   case('indx'); if(allocated(indxStruct%hru))then; check=.true.; else; allocate(indxStruct%hru(nHRU),stat=err); endif
   case('bpar'); ! do nothing for basin-average parameters: no spatial dimension 
   case('bvar'); ! do nothing for basin-average variables: no spatial dimension
   case('state');! do nothing -- ADD LATER
   case('diag'); ! do nothing -- Likely will not define at the top level
   case('flux'); ! do nothing -- Likely will not define at the top level
   case('deriv');! do nothing -- Likely will not define at the top level
   case default; err=20; message=trim(message)//'unable to identify lookup structure'; return
  end select
  ! check errors
  if(check) then; err=20; message=trim(message)//'structure '//trim(structInfo(iStruct)%structName)//' was unexpectedly allocated already'; return; endif
  if(err/=0)then; err=20; message=trim(message)//'problem allocating the '//trim(structInfo(iStruct)%structName)//' structure'; return; endif
 end do  ! looping through data structures

 ! -----
 ! * allocate variable dimensions...
 ! ---------------------------------

 do iStruct=1,nStruct  ! loop through data structures
  do iHRU=1,nHRU  ! loop through HRUs
   check=.false. ! initialize check
   ! cycle for variables without a HRU dimension
   select case(trim(structInfo(iStruct)%structName))
    case('time','bpar','bvar'); if(iHRU>1) cycle
   end select
   ! allocate variable dimension
   select case(trim(structInfo(iStruct)%structName))
    case('time'); if(allocated(timeStruct%var)          )then; check=.true.; else; allocate(timeStruct%var(structInfo(iStruct)%nVar),stat=err); endif
    case('forc'); if(allocated(forcStruct%hru(iHRU)%var))then; check=.true.; else; allocate(forcStruct%hru(iHRU)%var(structInfo(iStruct)%nVar),stat=err); endif
    case('attr'); if(allocated(attrStruct%hru(iHRU)%var))then; check=.true.; else; allocate(attrStruct%hru(iHRU)%var(structInfo(iStruct)%nVar),stat=err); endif
    case('type'); if(allocated(typeStruct%hru(iHRU)%var))then; check=.true.; else; allocate(typeStruct%hru(iHRU)%var(structInfo(iStruct)%nVar),stat=err); endif
    case('mpar'); if(allocated(mparStruct%hru(iHRU)%var))then; check=.true.; else; allocate(mparStruct%hru(iHRU)%var(structInfo(iStruct)%nVar),stat=err); endif
    case('mvar'); if(allocated(mvarStruct%hru(iHRU)%var))then; check=.true.; else; allocate(mvarStruct%hru(iHRU)%var(structInfo(iStruct)%nVar),stat=err); endif
    case('indx'); if(allocated(indxStruct%hru(iHRU)%var))then; check=.true.; else; allocate(indxStruct%hru(iHRU)%var(structInfo(iStruct)%nVar),stat=err); endif
    case('bpar'); if(allocated(bparStruct%var)          )then; check=.true.; else; allocate(bparStruct%var(structInfo(iStruct)%nVar),stat=err); endif 
    case('bvar'); if(allocated(bvarStruct%var)          )then; check=.true.; else; allocate(bvarStruct%var(structInfo(iStruct)%nVar),stat=err); endif
    case('state');! do nothing -- ADD LATER
    case('diag'); ! do nothing -- Likely will not define at the top level
    case('flux'); ! do nothing -- Likely will not define at the top level
    case('deriv');! do nothing -- Likely will not define at the top level
    case default; err=20; message=trim(message)//'unable to identify lookup structure'; return
   end select
   ! check errors
   if(check) then; err=20;write(message,'(a,i0)') trim(message)//'structure '//trim(structInfo(iStruct)%structName)//' was unexpectedly allocated already for HRU ', iHRU; return; endif
   if(err/=0)then; err=20;write(message,'(a,i0)') trim(message)//'problem allocating the '//trim(structInfo(iStruct)%structName)//' structure for HRU ', iHRU; return; endif
  end do  ! looping through data structures
 end do  ! loop through HRUs

 ! -----
 ! * allocate specific variables...
 ! --------------------------------

 do iStruct=1,nStruct  ! loop through data structures
  do iHRU=1,nHRU  ! loop through HRUs
   ! get the number of layers in the snow+soil sub-domain
   nLayers = nSnow(iHRU) + nSoil(iHRU)
   ! cycle for variables without a HRU dimension
   select case(trim(structInfo(iStruct)%structName))
    case('time','bpar','bvar'); if(iHRU>1) cycle
   end select
   ! allocate data dimension
   select case(trim(structInfo(iStruct)%structName))
    case('time'); ! do nothing: no data dimension 
    case('forc'); ! do nothing: no data dimension
    case('attr'); ! do nothing: no data dimension
    case('type'); ! do nothing: no data dimension
    case('mpar'); ! do nothing: no data dimension
    case('mvar'); call allocateDat_dp( mvar_meta,nSnow(iHRU),nSoil(iHRU),nLayers,mvarStruct%hru(iHRU),err,cmessage)
    case('indx'); call allocateDat_int(indx_meta,nSnow(iHRU),nSoil(iHRU),nLayers,indxStruct%hru(iHRU),err,cmessage)
    case('bpar'); ! do nothing: no data dimension 
    case('bvar'); call allocateDat_dp( bvar_meta,nSnow(iHRU),nSoil(iHRU),nLayers,bvarStruct,err,cmessage)  ! NOTE: no HRU dimension
    case('state');! do nothing -- ADD LATER
    case('diag'); ! do nothing -- Likely will not define at the top level
    case('flux'); ! do nothing -- Likely will not define at the top level
    case('deriv');! do nothing -- Likely will not define at the top level
    case default; err=20; message=trim(message)//'unable to identify lookup structure'; return
   end select
   ! check errors
   if(err/=0)then
    write(message,'(a,i0)') trim(message)//trim(cmessage)//'; attempting to allocate the '//trim(structInfo(iStruct)%structName)//' structure for HRU ', iHRU
    err=20; return
   endif
  end do  ! looping through data structures
 end do  ! loop through HRUs

 end subroutine initStruct

 ! ************************************************************************************************
 ! private subroutine allocateDat_dp: initialize data dimension of the data structures
 ! ************************************************************************************************
 subroutine allocateDat_dp(metadata,nSnow,nSoil,nLayers, & ! input
                           varData,err,message)            ! output
 ! provide access to derived data types
 USE data_struc,only:var_info       ! metadata structure
 USE data_struc,only:var_doubleVec  ! x%var(:)%dat (dp)
 implicit none
 ! input variables
 type(var_info),intent(in)         :: metadata(:) ! metadata structure
 integer(i4b),intent(in)           :: nSnow       ! number of snow layers
 integer(i4b),intent(in)           :: nSoil       ! number of soil layers
 integer(i4b),intent(in)           :: nLayers     ! total number of soil layers in the snow+soil domian (nSnow+nSoil)
 ! output variables
 type(var_doubleVec),intent(inout) :: varData     ! model variables for a local HRU
 integer(i4b),intent(out)          :: err         ! error code
 character(*),intent(out)          :: message     ! error message
 ! local variables
 integer(i4b)                      :: iVar        ! variable index
 ! initialize error control
 err=0; message='allocateDat_dp/'

 ! loop through variables in the data structure
 do iVar=1,size(metadata)

  ! check allocated
  if(allocated(varData%var(iVar)%dat))then
   message=trim(message)//'variable '//trim(metadata(iVar)%varname)//' is unexpectedly allocated'
   err=20; return

  ! allocate dimension
  else
   select case(trim(metadata(iVar)%vartype))
    case('scalarv'); allocate(varData%var(iVar)%dat(1),stat=err)
    case('wLength'); allocate(varData%var(iVar)%dat(nBand),stat=err)
    case('midSnow'); allocate(varData%var(iVar)%dat(nSnow),stat=err)
    case('midSoil'); allocate(varData%var(iVar)%dat(nSoil),stat=err)
    case('midToto'); allocate(varData%var(iVar)%dat(nLayers),stat=err)
    case('ifcSnow'); allocate(varData%var(iVar)%dat(0:nSnow),stat=err)
    case('ifcSoil'); allocate(varData%var(iVar)%dat(0:nSoil),stat=err)
    case('ifcToto'); allocate(varData%var(iVar)%dat(0:nLayers),stat=err)
    case('routing'); allocate(varData%var(iVar)%dat(nTimeDelay),stat=err)
    case('unknown'); allocate(varData%var(iVar)%dat(0),stat=err)  ! unknown=initialize with zero-length vector
    case default; err=40; message=trim(message)//"unknownVariableType[name='"//trim(metadata(iVar)%varname)//"'; type='"//trim(metadata(iVar)%vartype)//"']"; return
   endselect
   if(err/=0)then; err=20; message=trim(message)//'problem allocating variable '//trim(metadata(iVar)%varname); return; endif
   varData%var(iVar)%dat(:) = missingDouble
  endif  ! if not allocated

 end do  ! looping through variables

 end subroutine allocateDat_dp

 ! ************************************************************************************************
 ! private subroutine allocateDat_int: initialize data dimension of the data structures
 ! ************************************************************************************************
 subroutine allocateDat_int(metadata,nSnow,nSoil,nLayers, & ! input
                            varData,err,message)            ! output
 ! provide access to derived data types
 USE data_struc,only:var_info       ! metadata structure
 USE data_struc,only:var_intVec     ! x%var(:)%dat (i4b)
 implicit none
 ! input variables
 type(var_info),intent(in)         :: metadata(:) ! metadata structure
 integer(i4b),intent(in)           :: nSnow       ! number of snow layers
 integer(i4b),intent(in)           :: nSoil       ! number of soil layers
 integer(i4b),intent(in)           :: nLayers     ! total number of soil layers in the snow+soil domian (nSnow+nSoil)
 ! output variables
 type(var_intVec),intent(inout)    :: varData     ! model variables for a local HRU
 integer(i4b),intent(out)          :: err         ! error code
 character(*),intent(out)          :: message     ! error message
 ! local variables
 integer(i4b)                      :: iVar        ! variable index
 ! initialize error control
 err=0; message='allocateDat_int/'

 ! loop through variables in the data structure
 do iVar=1,size(metadata)

  ! check allocated
  if(allocated(varData%var(iVar)%dat))then
   message=trim(message)//'variable '//trim(metadata(iVar)%varname)//' is unexpectedly allocated'
   err=20; return

  ! allocate dimension
  else
   select case(trim(metadata(iVar)%vartype))
    case('scalarv'); allocate(varData%var(iVar)%dat(1),stat=err)
    case('wLength'); allocate(varData%var(iVar)%dat(nBand),stat=err)
    case('midSnow'); allocate(varData%var(iVar)%dat(nSnow),stat=err)
    case('midSoil'); allocate(varData%var(iVar)%dat(nSoil),stat=err)
    case('midToto'); allocate(varData%var(iVar)%dat(nLayers),stat=err)
    case('ifcSnow'); allocate(varData%var(iVar)%dat(0:nSnow),stat=err)
    case('ifcSoil'); allocate(varData%var(iVar)%dat(0:nSoil),stat=err)
    case('ifcToto'); allocate(varData%var(iVar)%dat(0:nLayers),stat=err)
    case('routing'); allocate(varData%var(iVar)%dat(nTimeDelay),stat=err)
    case('unknown'); allocate(varData%var(iVar)%dat(0),stat=err)  ! unknown=initialize with zero-length vector
    case default; err=40; message=trim(message)//"unknownVariableType[name='"//trim(metadata(iVar)%varname)//"'; type='"//trim(metadata(iVar)%vartype)//"']"; return
   endselect
   if(err/=0)then; err=20; message=trim(message)//'problem allocating variable '//trim(metadata(iVar)%varname); return; endif
   varData%var(iVar)%dat(:) = missingInteger
  endif  ! if not allocated

 end do  ! looping through variables

 end subroutine allocateDat_int

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
